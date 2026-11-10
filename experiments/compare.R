# ==================================================================
# Comparison of our scratch EM with a published EM implementation and with a
# general-purpose optimizer
# Data: Old Faithful waiting times, K = 2
#
#   EM (ours)  - the algorithm of this project
#   mixtools   - mixtools::normalmixEM, the standard EM for normal mixtures
#                in R (install.packages("mixtools"))
#   BFGS       - a quasi-Newton optimizer (stats::optim); R estimates the
#                gradient numerically from extra likelihood evaluations
#
# Note on mixtools: with its default settings normalmixEM does NOT run plain
# EM but ECM (Expectation / Conditional Maximization, Meng & Rubin 1993).
# Its code switches ECM on whenever the variances are allowed to differ
# (the line `ECM <- ECM || ... || any(sdcat != 1)` is always TRUE then).
# One ECM iteration:
#   1. E-step, then update the weights and the means
#   2. a second E-step with the new means, then update the weights again
#      and the variances
# So each iteration costs two E-steps instead of one. ECM still never lowers
# the log-likelihood and reaches the same maximum, but it follows a
# different path. (Verified: coding these steps by hand reproduces
# normalmixEM's output to machine precision.)
# ==================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))

if (!requireNamespace("mixtools", quietly = TRUE)) {
    stop("Please install mixtools first: install.packages(\"mixtools\")")
}

method_cols <- c(
    "EM (ours)" = "#1b9e77",
    "mixtools (ECM)" = "#e6ab02",
    "BFGS" = "#7570b3"
)
methods <- names(method_cols)

# ------------------------------------------------------------------
# 0. Utilities
# ------------------------------------------------------------------

# General optimizers work on unconstrained numbers, but pi must lie in (0, 1)
# and sigma2 must be positive. So they optimize over
#   theta = (log(pi_2 / pi_1), mu_1, mu_2, log sigma_1, log sigma_2)
# and we convert back and forth.
pack_par <- function(pi, mu, sigma2) {
    c(log(pi[-1] / pi[1]), mu, 0.5 * log(sigma2))
}

unpack_par <- function(theta, K) {
    a <- c(0, theta[seq_len(K - 1)])
    pi <- exp(a - max(a))
    pi <- pi / sum(pi)
    list(
        pi = pi,
        mu = theta[(K - 1) + seq_len(K)],
        sigma2 = exp(2 * theta[(2 * K - 1) + seq_len(K)])
    )
}

# Order components by mean (same convention as em_gmm)
sort_par <- function(par) {
    o <- order(par$mu)
    list(pi = par$pi[o], mu = par$mu[o], sigma2 = par$sigma2[o])
}

# Negative log-likelihood with bookkeeping.
# "cost" counts likelihood evaluations (each costs about as much as one
# E-step), so EM iterations and optimizer calls share the same x-axis.
make_traced <- function(x, K, var_floor) {
    env <- new.env()
    env$cost <- 0
    env$cost_v <- numeric(0)
    env$ll_v <- numeric(0)
    env$theta_v <- list()

    f <- function(theta) {
        env$cost <- env$cost + 1
        if (any(!is.finite(theta))) {
            return(1e10)
        }
        p <- unpack_par(theta, K)
        ll <- gmm_loglikelihood(x, p$pi, p$mu, pmax(p$sigma2, var_floor))
        if (!is.finite(ll)) {
            ll <- -1e10
        }
        env$cost_v <- c(env$cost_v, env$cost)
        env$ll_v <- c(env$ll_v, ll)
        env$theta_v[[length(env$theta_v) + 1]] <- theta
        -ll
    }
    list(f = f, env = env)
}

# EM path in (mu_1, mu_2) space (K = 2), re-using e_step / m_step
em_path <- function(x, init, n_iter = 100, var_floor = (0.01 * stats::sd(x))^2) {
    par <- init
    out <- vector("list", n_iter + 1)
    for (i in 0:n_iter) {
        m <- sort(par$mu)
        out[[i + 1]] <- data.frame(iter = i, mu1 = m[1], mu2 = m[2])
        if (i == n_iter) {
            break
        }
        E <- e_step(x, par$pi, par$mu, par$sigma2)
        M <- m_step(x, E$gamma, var_floor)
        par <- list(pi = M$pi, mu = M$mu, sigma2 = M$sigma2)
    }
    do.call(rbind, out)
}

failed_run <- function(method) {
    list(
        method = method,
        par = list(pi = rep(NA_real_, 2), mu = rep(NA_real_, 2), sigma2 = rep(NA_real_, 2)),
        loglik = NA_real_, cost = NA_real_, iterations = NA_real_,
        converged = FALSE, restarted = NA,
        trace = data.frame(cost = numeric(0), loglik = numeric(0)),
        path = data.frame(mu1 = numeric(0), mu2 = numeric(0))
    )
}

# Run mixtools::normalmixEM quietly (it prints the iteration count itself)
run_mixtools <- function(x, init, K, tol, maxit) {
    invisible(utils::capture.output(
        fit <- mixtools::normalmixEM(
            x, lambda = init$pi, mu = init$mu, sigma = sqrt(init$sigma2),
            k = K, epsilon = tol, maxit = maxit, verb = FALSE
        )
    ))
    fit
}

# Path of mixtools in (mu_1, mu_2) space. normalmixEM does not store its
# parameters per iteration, so it is re-run with maxit = 1, 2, ... from the
# same start (it is deterministic, so this reproduces the same path).
mixtools_path <- function(x, init, n_iter, K = 2, tol = 1e-10) {
    out <- vector("list", n_iter + 1)
    out[[1]] <- data.frame(mu1 = min(init$mu), mu2 = max(init$mu))
    for (i in seq_len(n_iter)) {
        m <- sort(run_mixtools(x, init, K, tol, maxit = i)$mu)
        out[[i + 1]] <- data.frame(mu1 = m[1], mu2 = m[2])
    }
    do.call(rbind, out)
}

# Run one method from given starting values.
# Returns the final parameters, the log-likelihood recomputed with our own
# function (so every method is judged by the same yardstick), the cost, the
# number of iterations (EM methods only) and a log-likelihood trace (best
# value so far vs cost).
# "restarted" is TRUE when mixtools abandoned our starting values because a
# variance collapsed and restarted from its own random values; such a run is
# not a same-start comparison.
run_method <- function(
    method, x, init, K = 2, maxit = 5000, tol = 1e-10,
    var_floor = (0.01 * stats::sd(x))^2, record_path = FALSE
) {
    theta0 <- pack_par(init$pi, init$mu, init$sigma2)
    tryCatch(
        {
            path <- data.frame(mu1 = numeric(0), mu2 = numeric(0))
            iterations <- NA_real_
            restarted <- FALSE
            if (method == "EM (ours)") {
                fit <- em_gmm(x, init, tol = tol, max_iter = maxit, var_floor = var_floor)
                par <- list(pi = fit$pi, mu = fit$mu, sigma2 = fit$sigma2)
                # one E-step (likelihood evaluation) per entry of loglik_hist,
                # including the one at the starting values
                trace <- data.frame(cost = seq_along(fit$loglik_hist), loglik = fit$loglik_hist)
                cost <- length(fit$loglik_hist)
                iterations <- fit$iterations
                converged <- fit$status == "converged"
                if (record_path) {
                    path <- em_path(x, init, n_iter = min(fit$iterations, 200), var_floor = var_floor)
                    path <- path[, c("mu1", "mu2")]
                }
            } else if (method == "mixtools (ECM)") {
                fit <- run_mixtools(x, init, K, tol, maxit)
                par <- list(pi = fit$lambda, mu = fit$mu, sigma2 = fit$sigma^2)
                iterations <- length(fit$all.loglik) - 1
                # one E-step at the start + two E-steps per ECM iteration
                trace <- data.frame(
                    cost = 1 + 2 * (seq_along(fit$all.loglik) - 1),
                    loglik = fit$all.loglik
                )
                cost <- 1 + 2 * iterations
                converged <- iterations < maxit
                restarted <- fit$restarts > 0
                if (record_path) {
                    path <- mixtools_path(x, init, n_iter = min(iterations, 200), K = K, tol = tol)
                }
            } else {
                tr <- make_traced(x, K, var_floor)
                res <- switch(method,
                    "BFGS" = stats::optim(
                        theta0, tr$f, method = "BFGS",
                        control = list(maxit = maxit, reltol = tol)
                    ),
                    stop("unknown method: ", method)
                )
                converged <- res$convergence == 0
                p <- unpack_par(res$par, K)
                par <- list(pi = p$pi, mu = p$mu, sigma2 = pmax(p$sigma2, var_floor))
                trace <- data.frame(cost = tr$env$cost_v, loglik = cummax(tr$env$ll_v))
                cost <- tr$env$cost
                if (record_path && length(tr$env$ll_v) > 0) {
                    # keep only the evaluations that improved on the best so far
                    idx <- which(diff(c(-Inf, cummax(tr$env$ll_v))) > 0)
                    path <- do.call(rbind, lapply(tr$env$theta_v[idx], function(th) {
                        m <- sort(unpack_par(th, K)$mu)
                        data.frame(mu1 = m[1], mu2 = m[2])
                    }))
                }
            }
            par <- sort_par(par)
            list(
                method = method, par = par,
                loglik = gmm_loglikelihood(x, par$pi, par$mu, par$sigma2),
                cost = cost, iterations = iterations,
                converged = converged, restarted = restarted,
                trace = trace, path = path
            )
        },
        error = function(e) failed_run(method)
    )
}

summarise_run <- function(r) {
    p <- r$par
    data.frame(
        method = r$method,
        loglik = r$loglik,
        pi1 = p$pi[1], mu1 = p$mu[1], mu2 = p$mu[2],
        sigma1 = sqrt(p$sigma2[1]), sigma2 = sqrt(p$sigma2[2]),
        cost = r$cost, iterations = r$iterations, converged = r$converged,
        restarted = r$restarted,
        row.names = NULL
    )
}

# ------------------------------------------------------------------
# 1. One common starting point: results table
# ------------------------------------------------------------------
data(faithful)
x <- faithful$waiting
set.seed(42)
init_km <- init_kmeans(x, 2)

runs <- lapply(methods, function(m) run_method(m, x, init_km))
names(runs) <- methods

res_tab <- bind_rows(lapply(runs, summarise_run))
ll_star <- max(res_tab$loglik, na.rm = TRUE)
res_tab$`gap to best` <- ll_star - res_tab$loglik

kable(
    res_tab %>% select(-restarted),
    digits = 4,
    caption = paste(
        "Table 7. Final fit of each method from the same k-means start (Faithful waiting times).",
        "Cost = number of likelihood evaluations / E-steps (EM: 1 per iteration, mixtools ECM: 2 per iteration, plus 1 at the start)."
    )
)

# ------------------------------------------------------------------
# 2. Convergence race: log-likelihood gap vs computational cost
# ------------------------------------------------------------------
trace_df <- bind_rows(lapply(runs, function(r) {
    data.frame(method = r$method, r$trace)
})) %>%
    mutate(gap = pmax(ll_star - loglik, 1e-12))

ggplot(trace_df, aes(cost, gap, colour = method)) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() +
    scale_y_log10(
        breaks = 10^seq(-12, 2, by = 2),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_colour_manual(values = method_cols) +
    labs(
        x = "Computational cost (likelihood evaluations, log scale)",
        y = "Log-likelihood gap to best (log scale)",
        colour = NULL,
        title = "Convergence race from the same starting point",
        subtitle = "Lower-left is better. Both EM variants improve at every step;\nmixtools (ECM) pays 2 E-steps per iteration."
    )

# ------------------------------------------------------------------
# 3. Do all methods agree? Estimates and cost by method
# ------------------------------------------------------------------
est_long <- res_tab %>%
    select(method, pi1, mu1, mu2, sigma1, sigma2) %>%
    pivot_longer(-method, names_to = "parameter", values_to = "estimate") %>%
    mutate(parameter = factor(parameter, levels = c("pi1", "mu1", "mu2", "sigma1", "sigma2")))

ggplot(est_long, aes(estimate, reorder(method, desc(method)), colour = method)) +
    geom_point(size = 3) +
    facet_wrap(~parameter, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = method_cols, guide = "none") +
    labs(
        x = "Estimate", y = NULL,
        title = "Parameter estimates by method",
        subtitle = "If all points line up vertically, every method found the same maximum."
    )

ggplot(res_tab, aes(reorder(method, cost), cost, fill = method)) +
    geom_col() +
    geom_text(aes(label = cost), hjust = -0.2) +
    coord_flip() +
    expand_limits(y = max(res_tab$cost, na.rm = TRUE) * 1.15) +
    scale_fill_manual(values = method_cols, guide = "none") +
    labs(
        x = NULL, y = "Likelihood evaluations",
        title = "Cost to reach the optimum (one run, k-means start)"
    )

# ------------------------------------------------------------------
# 4. Optimizer paths on the likelihood surface (slice in the (mu1, mu2) plane)
#    pi and sigma are held at the MLE, so the contour is a 2D slice of the
#    full likelihood; paths are projected onto it.
# ------------------------------------------------------------------
fit_star <- em_gmm(x, init_km)
ll_fit <- gmm_loglikelihood(x, fit_star$pi, fit_star$mu, fit_star$sigma2)

mu_grid <- seq(min(x), max(x), length.out = 90)
surf <- expand.grid(mu1 = mu_grid, mu2 = mu_grid)
surf$ll <- mapply(
    function(a, b) gmm_loglikelihood(x, fit_star$pi, c(a, b), fit_star$sigma2),
    surf$mu1, surf$mu2
)
surf$gap <- pmax(ll_fit - surf$ll, 0)

# A deliberately poor start so that the paths are long enough to see
init_poor <- list(pi = c(0.5, 0.5), mu = c(50, 60), sigma2 = c(8, 8)^2)

path_df <- bind_rows(lapply(methods, function(m) {
    r <- run_method(m, x, init_poor, record_path = TRUE)
    if (nrow(r$path) == 0) {
        return(NULL)
    }
    data.frame(method = m, step = seq_len(nrow(r$path)), r$path)
})) %>%
    mutate(method = factor(method, levels = methods))

ggplot() +
    geom_contour_filled(
        data = surf, aes(mu1, mu2, z = gap),
        breaks = c(0, 2, 5, 10, 25, 50, 100, 250, 500, Inf)
    ) +
    geom_path(
        data = path_df, aes(mu1, mu2),
        colour = "white", linewidth = 1.6
    ) +
    geom_path(
        data = path_df, aes(mu1, mu2, colour = method),
        linewidth = 0.8
    ) +
    geom_point(
        data = path_df %>% group_by(method) %>% slice(1),
        aes(mu1, mu2), shape = 8, size = 3, colour = "white"
    ) +
    geom_point(
        data = path_df %>% group_by(method) %>% slice_tail(n = 1),
        aes(mu1, mu2), shape = 21, size = 3, fill = "white", colour = "black"
    ) +
    facet_wrap(~method) +
    scale_fill_viridis_d(direction = -1, name = "Log-lik gap") +
    scale_colour_manual(values = method_cols, guide = "none") +
    coord_cartesian(xlim = range(x), ylim = range(x), expand = FALSE) +
    labs(
        x = expression(mu[1]), y = expression(mu[2]),
        title = "Paths of each optimizer on the likelihood surface",
        subtitle = "Slice with pi and sigma fixed at the MLE. Star = start, circle = end.\nOnly improving steps are drawn."
    )

# ------------------------------------------------------------------
# 5. Robustness: 100 random starts
# ------------------------------------------------------------------
R_starts <- 100
set.seed(2026)
inits <- lapply(seq_len(R_starts), function(i) init_random(x, 2))

rob <- bind_rows(lapply(seq_len(R_starts), function(i) {
    bind_rows(lapply(methods, function(m) {
        s <- summarise_run(run_method(m, x, inits[[i]]))
        s$start <- i
        s
    }))
}))
set.seed(42)

best_ll <- max(rob$loglik, na.rm = TRUE)
rob <- rob %>%
    mutate(
        gap = pmax(best_ll - loglik, 1e-10),
        success = !is.na(loglik) & (best_ll - loglik) < 1e-3
    )

rob_sum <- rob %>%
    group_by(method) %>%
    summarise(
        `success %` = 100 * mean(success),
        `converged %` = 100 * mean(converged),
        `median cost` = median(cost, na.rm = TRUE),
        `median iterations` = median(iterations, na.rm = TRUE),
        .groups = "drop"
    )
n_restarted <- sum(rob$restarted, na.rm = TRUE)
kable(
    rob_sum,
    digits = 1,
    caption = sprintf(
        paste(
            "Table 8. Robustness over %d random starts (success = final log-likelihood within 1e-3 of the best found by any method).",
            "mixtools abandoned the given start (collapsing variance) and restarted from its own random values in %d runs."
        ),
        R_starts, n_restarted
    )
)

ggplot(rob, aes(method, gap, colour = method)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.5) +
    geom_boxplot(fill = NA, outlier.shape = NA, colour = "grey30", width = 0.5) +
    scale_y_log10(
        breaks = 10^seq(-10, 3, by = 3),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_colour_manual(values = method_cols, guide = "none") +
    coord_flip() +
    labs(
        x = NULL, y = "Final log-likelihood gap to best (log scale)",
        title = "Where does each method end up from 100 random starts?",
        subtitle = "Points at the bottom reached the global maximum; points above are local optima or premature stops."
    )

ggplot(rob_sum, aes(reorder(method, `success %`), `success %`, fill = method)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.0f%%", `success %`)), hjust = -0.1) +
    coord_flip(ylim = c(0, 110)) +
    scale_fill_manual(values = method_cols, guide = "none") +
    labs(x = NULL, y = "Runs reaching the global maximum (%)", title = "Success rate from random starts")

ggplot(rob, aes(reorder(method, cost, median), cost, fill = method)) +
    geom_boxplot(outlier.alpha = 0.4) +
    scale_y_log10() +
    scale_fill_manual(values = method_cols, guide = "none") +
    coord_flip() +
    labs(x = NULL, y = "Cost (likelihood evaluations, log scale)", title = "Cost from random starts")
