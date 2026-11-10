# Simulate from the two-component mixture used in all simulation studies:
#   x ~ p1 * N(0, 1) + (1 - p1) * N(delta, 1)
# Returns list(x, z) where z = 0 for component 1 and z = 1 for component 2.
sim_mixture <- function(n, delta, p1 = 0.5) {
    z <- stats::rbinom(n, 1, 1 - p1)
    list(x = stats::rnorm(n, mean = delta * z, sd = 1), z = z)
}

# Initialize parameters from hard cluster labels
init_from_labels <- function(x, cl, K) {
    n <- length(x)
    # Fallback for empty clusters to prevent NaNs
    counts <- tabulate(cl, K)
    pi <- pmax(counts, 1e-6) / n
    pi <- pi / sum(pi)
    mu <- vapply(seq_len(K), function(k) {
        # x[sample.int(...)] instead of sample(x, 1): sample() treats a
        # length-1 numeric x as 1:x
        if (counts[k] == 0) x[sample.int(n, 1)] else mean(x[cl == k])
    }, numeric(1))
    sigma2 <- vapply(seq_len(K), function(k) {
        if (counts[k] > 1) stats::var(x[cl == k]) else stats::var(x)
    }, numeric(1))
    list(pi = pi, mu = mu, sigma2 = sigma2)
}

# Initialize parameters using K-Means
init_kmeans <- function(x, K) {
    km <- stats::kmeans(x, centers = K, nstart = 10)
    init_from_labels(x, km$cluster, K)
}

# Initialize parameters randomly
# Means are drawn from the DISTINCT data values. Drawing positions with
# sample(x, K) can return equal values when the data have ties (Faithful
# waiting times are integers), and equal means + equal weights + equal
# variances is an exact saddle point that EM never leaves.
init_random <- function(x, K) {
    ux <- unique(x)
    if (length(ux) < K) {
        stop("init_random: fewer distinct values than components")
    }
    list(
        pi = rep(1 / K, K),
        mu = ux[sample.int(length(ux), K)],
        sigma2 = rep(stats::var(x), K)
    )
}

# Fit GMM using multiple initializations (best log-likelihood selection)
#
#  x Observations vector
#  K Number of components
#  n_random Number of random initializations in addition to k-means
#  ... Forwarded arguments to em_gmm (tol, max_iter, var_floor)
em_best <- function(x, K, n_random = 10, ...) {
    n <- length(x)
    # Closed-form MLE for single component (no EM needed)
    if (K == 1) {
        s2 <- mean((x - mean(x))^2) # MLE biased variance
        ll <- sum(stats::dnorm(x, mean = mean(x), sd = sqrt(s2), log = TRUE))
        return(list(
            pi = 1,
            mu = mean(x),
            sigma2 = s2,
            K = 1,
            iterations = 0,
            loglik = ll,
            loglik_hist = ll,
            status = "closed form",
            gamma = matrix(1, nrow = n, ncol = 1),
            monotone = TRUE
        ))
    }

    # Generate runs (1 K-means run + n_random random starts)
    fits <- c(
        list(em_gmm(x, init_kmeans(x, K), ...)),
        lapply(seq_len(n_random), function(r) em_gmm(x, init_random(x, K), ...))
    )
    # Filter out degenerated runs
    ok <- vapply(fits, function(f) f$status != "degenerate", logical(1))
    valid_fits <- fits[ok]
    # Fallback if all initializations degenerate. em_gmm rejects the degenerate
    # M-step, so these fits still hold the last valid parameters and the
    # log-likelihood matching them.
    if (length(valid_fits) == 0) {
        warning(
            "All initializations degenerated. Returning the highest-loglik fit ",
            "(last valid parameters before collapse)."
        )
        valid_fits <- fits
    }
    best_idx <- which.max(vapply(valid_fits, `[[`, numeric(1), "loglik"))
    valid_fits[[best_idx]]
}

# Evaluate continuous GMM density over a grid (for plotting)
#  grid Numeric vector of evaluation points
#  fit Model output from em_gmm or em_best
gmm_density <- function(grid, fit) {
    dens <- vapply(seq_len(fit$K), function(k) {
        fit$pi[k] *
            stats::dnorm(grid, mean = fit$mu[k], sd = sqrt(fit$sigma2[k]))
    }, numeric(length(grid)))
    # keep the matrix shape for a one-point grid
    rowSums(matrix(dens, nrow = length(grid)))
}

# ------------------------------------------------------------------
# Model selection: AIC and BIC
# ------------------------------------------------------------------

# Number of free parameters of a K-component univariate Gaussian mixture:
#   (K - 1) mixing weights (they sum to 1) + K means + K variances = 3K - 1
gmm_npar <- function(K) 3 * K - 1

# Akaike information criterion:  AIC = -2 * loglik + 2 * p   (smaller is better)
gmm_aic <- function(loglik, K) {
    -2 * loglik + 2 * gmm_npar(K)
}

# Bayesian information criterion: BIC = -2 * loglik + p * log(n)   (smaller is better)
gmm_bic <- function(loglik, K, n) {
    -2 * loglik + gmm_npar(K) * log(n)
}

# One-row summary of a fitted model (from em_gmm or em_best) with its criteria.
# The log-likelihood is recomputed from the parameters so that the function
# also works for fits produced by other code.
gmm_criteria <- function(x, fit) {
    n <- length(x)
    ll <- gmm_loglikelihood(x, fit$pi, fit$mu, fit$sigma2)
    data.frame(
        K = fit$K,
        npar = gmm_npar(fit$K),
        loglik = ll,
        AIC = gmm_aic(ll, fit$K),
        BIC = gmm_bic(ll, fit$K, n),
        status = fit$status,
        stringsAsFactors = FALSE
    )
}

# Fit K = K_range with em_best and compare by AIC / BIC.
#  x Observations vector
#  K_range Candidate numbers of components
#  n_random Random initializations per K (forwarded to em_best)
#  ... Forwarded to em_gmm (tol, max_iter, var_floor)
# Returns a list:
#  table one row per K with loglik, AIC, BIC, a "valid" flag and differences
#        to the best valid model (dAIC, dBIC)
#  fits the fitted models, named "K1", "K2", ...
#  K_aic, K_bic the K chosen by each criterion (degenerate fits are excluded)
select_K <- function(x, K_range = 1:6, n_random = 10, ...) {
    fits <- lapply(K_range, function(K) em_best(x, K, n_random = n_random, ...))
    names(fits) <- paste0("K", K_range)
    tab <- do.call(rbind, lapply(fits, function(f) gmm_criteria(x, f)))
    rownames(tab) <- NULL
    # A degenerate fit sits next to a collapsing component, so it should not
    # compete in model selection
    tab$valid <- tab$status != "degenerate"
    if (!any(tab$valid)) {
        warning("All fits are degenerate; selecting among them anyway.")
        tab$valid <- TRUE
    }
    tab$dAIC <- tab$AIC - min(tab$AIC[tab$valid])
    tab$dBIC <- tab$BIC - min(tab$BIC[tab$valid])
    v <- tab[tab$valid, ]
    list(
        table = tab,
        fits = fits,
        K_aic = v$K[which.min(v$AIC)],
        K_bic = v$K[which.min(v$BIC)]
    )
}

# ------------------------------------------------------------------
# Standard errors and confidence intervals: nonparametric bootstrap
# ------------------------------------------------------------------
# Idea: resample the data with replacement B times, refit EM each time and look
# at how much the estimates move. No formulas needed.
#   SE = standard deviation of the bootstrap estimates
#   CI = percentile interval (e.g. 2.5% and 97.5% quantiles for level = 0.95)
#
# Each resample is refitted starting from the full-data estimate (warm start).
# This is fast, but it keeps every refit near the same solution, so with
# strongly overlapping components it can understate the variability.
#
#  x Observations vector
#  fit Fitted model (from em_gmm or em_best); also used as warm start
#  B Number of bootstrap resamples
#  level Confidence level
#  ... Forwarded to em_gmm (tol, max_iter, var_floor)
# Returns one row per parameter: estimate, se, lower, upper,
#   n_failed (degenerate resamples, excluded) and
#   n_not_converged (resamples that hit max_iter, included)
gmm_ci <- function(x, fit, B = 200, level = 0.95, ...) {
    n <- length(x)
    K <- fit$K
    init <- list(pi = fit$pi, mu = fit$mu, sigma2 = fit$sigma2)
    par_names <- c(
        paste0("pi", seq_len(K)),
        paste0("mu", seq_len(K)),
        paste0("sigma", seq_len(K)) # standard deviations, not variances
    )

    est <- matrix(NA_real_, nrow = B, ncol = 3 * K, dimnames = list(NULL, par_names))
    status_b <- character(B)
    for (b in seq_len(B)) {
        xb <- sample(x, n, replace = TRUE)
        fb <- em_gmm(xb, init, ...)
        status_b[b] <- fb$status
        # skip degenerate bootstrap fits (stay NA)
        if (fb$status != "degenerate") {
            est[b, ] <- c(fb$pi, fb$mu, sqrt(fb$sigma2))
        }
    }

    a <- (1 - level) / 2
    data.frame(
        parameter = par_names,
        estimate = c(fit$pi, fit$mu, sqrt(fit$sigma2)),
        se = apply(est, 2, stats::sd, na.rm = TRUE),
        lower = apply(est, 2, stats::quantile, probs = a, na.rm = TRUE),
        upper = apply(est, 2, stats::quantile, probs = 1 - a, na.rm = TRUE),
        n_failed = sum(status_b == "degenerate"),
        n_not_converged = sum(status_b == "max_iter reached"),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

# ------------------------------------------------------------------
# Convergence rate of EM
# ------------------------------------------------------------------
# EM converges LINEARLY: near the optimum theta*, the parameter error shrinks
# by a roughly constant factor r per iteration,
#   theta_{t+1} - theta* ~ r * (theta_t - theta*),   0 < r < 1,
# where r is the largest eigenvalue of the Jacobian of the EM map (the
# largest "fraction of missing information", Dempster, Laird & Rubin 1977).
# The parameter steps ||theta_{t+1} - theta_t|| shrink by the same factor, so
# their successive ratios estimate r.
# r close to 0 = fast, r close to 1 = very slow.

# Run EM step by step and record the full history (loglik and parameters).
# Returns a matrix with one row per iteration (row 1 = starting values, so
# nrow - 1 = number of M-steps); components are sorted by mean.
# attr(out, "status") is "converged", "max_iter reached" or "degenerate".
#  x Observations vector
#  init List with starting values pi, mu, sigma2
#  max_iter Maximum number of M-steps
#  tol Relative log-likelihood tolerance (set very small to see the asymptotic regime)
#  var_floor Minimum variance (as in em_gmm)
em_trace <- function(
    x, init, max_iter = 2000, tol = 1e-14,
    var_floor = (0.01 * stats::sd(x))^2
) {
    K <- length(init$pi)
    par <- init
    rows <- vector("list", max_iter + 1)
    status <- "max_iter reached"
    ll_prev <- NA_real_
    n_done <- 0
    for (i in 0:max_iter) {
        E <- e_step(x, par$pi, par$mu, par$sigma2)
        o <- order(par$mu)
        rows[[i + 1]] <- c(i, E$loglik, par$pi[o], par$mu[o], par$sigma2[o])
        n_done <- i + 1
        if (i > 0 && abs(E$loglik - ll_prev) < tol * (1 + abs(ll_prev))) {
            status <- "converged"
            break
        }
        if (i == max_iter) {
            break
        }
        ll_prev <- E$loglik
        M <- m_step(x, E$gamma, var_floor)
        if (M$degenerate) {
            status <- "degenerate"
            break
        }
        par <- list(pi = M$pi, mu = M$mu, sigma2 = M$sigma2)
    }
    out <- do.call(rbind, rows[seq_len(n_done)])
    colnames(out) <- c(
        "iter", "loglik",
        paste0("pi", seq_len(K)),
        paste0("mu", seq_len(K)),
        paste0("sigma2_", seq_len(K))
    )
    attr(out, "status") <- status
    out
}

# Estimate the linear convergence rate from an em_trace matrix.
#  floor Ignore parameter steps smaller than this (numerical noise)
#  n_last Rate = median of the last n_last valid ratios (the asymptotic regime)
# Returns list(rate, ratio):
#  rate    estimate of r from successive parameter steps
#  ratio   data frame (iter, ratio) of the parameter-step ratios
#          ||theta_{t+1} - theta_t|| / ||theta_t - theta_{t-1}||
# rate is NA when fewer than 3 usable ratios exist (very fast convergence).
em_rate <- function(tr, floor = 1e-9, n_last = 20) {
    if (nrow(tr) < 5) {
        return(list(rate = NA_real_, ratio = data.frame(iter = integer(0), ratio = numeric(0))))
    }
    P <- tr[, -(1:2), drop = FALSE]
    step <- sqrt(rowSums(diff(P)^2)) # step[t] = ||theta_t - theta_{t-1}||
    ratio <- step[-1] / step[-length(step)]
    ok <- which(step[-length(step)] > floor & step[-1] > floor)
    ratio_df <- data.frame(iter = ok, ratio = ratio[ok])
    rate <- if (nrow(ratio_df) >= 3) {
        stats::median(utils::tail(ratio_df$ratio, n_last))
    } else {
        NA_real_
    }
    list(rate = rate, ratio = ratio_df)
}
