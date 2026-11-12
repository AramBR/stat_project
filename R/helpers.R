# Initialize parameters from hard cluster labels
init_from_labels <- function(x, cl, K) {
    n <- length(x)
    # Fallback for empty clusters to prevent NaNs
    counts <- tabulate(cl, K)
    pi <- pmax(counts, 1e-6) / n
    pi <- pi / sum(pi)
    mu <- vapply(seq_len(K), function(k) {
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

gmm_density <- function(grid, fit) {
    dens <- vapply(seq_len(fit$K), function(k) {
        fit$pi[k] *
            stats::dnorm(grid, mean = fit$mu[k], sd = sqrt(fit$sigma2[k]))
    }, numeric(length(grid)))
    # keep the matrix shape for a one-point grid
    rowSums(matrix(dens, nrow = length(grid)))
}

gmm_npar <- function(K) 3 * K - 1

# Akaike information criterion:  AIC = -2 * loglik + 2 * p   (smaller is better)
gmm_aic <- function(loglik, K) {
    -2 * loglik + 2 * gmm_npar(K)
}

# Bayesian information criterion: BIC = -2 * loglik + p * log(n)   (smaller is better)
gmm_bic <- function(loglik, K, n) {
    -2 * loglik + gmm_npar(K) * log(n)
}

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

select_K <- function(x, K_range = 1:6, n_random = 10, ...) {
    fits <- lapply(K_range, function(K) em_best(x, K, n_random = n_random, ...))
    names(fits) <- paste0("K", K_range)
    tab <- do.call(rbind, lapply(fits, function(f) gmm_criteria(x, f)))
    rownames(tab) <- NULL
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

em_rate <- function(tr, floor = 1e-9, n_last = 20) {
    if (nrow(tr) < 5) {
        return(NA_real_)
    }
    P <- tr[, -(1:2), drop = FALSE]
    step <- sqrt(rowSums(diff(P)^2)) # step[t] = ||theta_t - theta_{t-1}||
    ratio <- step[-1] / step[-length(step)]
    ok <- step[-length(step)] > floor & step[-1] > floor
    if (sum(ok) < 3) {
        return(NA_real_)
    }
    stats::median(utils::tail(ratio[ok], n_last))
}
