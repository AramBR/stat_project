library(dplyr)
library(ggplot2)
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

make_traced <- function(x, K, var_floor) {
    env <- new.env()
    env$cost <- 0
    env$cost_v <- numeric(0)
    env$ll_v <- numeric(0)

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
        -ll
    }
    list(f = f, env = env)
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

run_method <- function(
    method, x, init, K = 2, maxit = 5000, tol = 1e-10,
    var_floor = (0.01 * stats::sd(x))^2
) {
    if (method == "EM (ours)") {
        fit <- em_gmm(x, init, tol = tol, max_iter = maxit, var_floor = var_floor)
        par <- list(pi = fit$pi, mu = fit$mu, sigma2 = fit$sigma2)
        # one E-step per entry of loglik_hist, including the one at the start
        trace <- data.frame(cost = seq_along(fit$loglik_hist), loglik = fit$loglik_hist)
    } else if (method == "mixtools (ECM)") {
        fit <- run_mixtools(x, init, K, tol, maxit)
        par <- list(pi = fit$lambda, mu = fit$mu, sigma2 = fit$sigma^2)
        # one E-step at the start + two E-steps per ECM iteration
        trace <- data.frame(
            cost = 1 + 2 * (seq_along(fit$all.loglik) - 1),
            loglik = fit$all.loglik
        )
    } else if (method == "BFGS") {
        tr <- make_traced(x, K, var_floor)
        res <- stats::optim(
            pack_par(init$pi, init$mu, init$sigma2), tr$f,
            method = "BFGS", control = list(maxit = maxit, reltol = tol)
        )
        p <- unpack_par(res$par, K)
        par <- list(pi = p$pi, mu = p$mu, sigma2 = pmax(p$sigma2, var_floor))
        # best log-likelihood found so far, after each evaluation
        trace <- data.frame(cost = tr$env$cost_v, loglik = cummax(tr$env$ll_v))
    } else {
        stop("unknown method: ", method)
    }
    list(
        method = method,
        loglik = gmm_loglikelihood(x, par$pi, par$mu, par$sigma2),
        trace = trace
    )
}

data(faithful)
x <- faithful$waiting
set.seed(42)
init_km <- init_kmeans(x, 2)

runs <- lapply(methods, function(m) run_method(m, x, init_km))
ll_star <- max(vapply(runs, `[[`, numeric(1), "loglik"))

# Log-likelihood gap to the best value vs computational cost
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
