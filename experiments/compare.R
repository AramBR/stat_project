# ==================================================================
# Comparison of our scratch EM with general-purpose optimizers
# and with library EM implementations (mixtools, mclust)
# Data: Old Faithful waiting times, K = 2
#
# Optional packages (sections are skipped if missing):
#   install.packages(c("mixtools", "mclust"))
# ==================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))

has_mixtools <- requireNamespace("mixtools", quietly = TRUE)
has_mclust <- requireNamespace("mclust", quietly = TRUE)

method_cols <- c(
    "EM (ours)" = "#1b9e77",
    "Nelder-Mead" = "#d95f02",
    "BFGS (num. grad)" = "#7570b3",
    "BFGS (analytic grad)" = "#e7298a",
    "L-BFGS-B" = "#66a61e",
    "nlm" = "#e6ab02",
    "mixtools::normalmixEM" = "#a6761d",
    "mclust::Mclust" = "#666666"
)

core_methods <- c(
    "EM (ours)", "Nelder-Mead", "BFGS (num. grad)",
    "BFGS (analytic grad)", "L-BFGS-B", "nlm"
)
lib_methods <- c(
    if (has_mixtools) "mixtools::normalmixEM",
    if (has_mclust) "mclust::Mclust"
)
all_methods <- c(core_methods, lib_methods)
# methods that accept our starting values and report a trace
traced_methods <- c(core_methods, if (has_mixtools) "mixtools::normalmixEM")

# ------------------------------------------------------------------
# 0. Utilities
# ------------------------------------------------------------------

# Unconstrained parametrisation for the general optimizers:
# theta = (logits of pi_2..pi_K, mu_1..mu_K, log sigma_1..log sigma_K)
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

# Negative log-likelihood (and its gradient) with bookkeeping.
# "cost" counts one unit per likelihood evaluation or gradient evaluation
# (each is roughly one E-step), so that EM iterations and optimizer calls
# can be put on the same horizontal axis.
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

    # Analytic gradient of the NEGATIVE log-likelihood in theta.
    # Uses the E-step responsibilities (Fisher's identity).
    gr <- function(theta) {
        env$cost <- env$cost + 1
        p <- unpack_par(theta, K)
        s2 <- pmax(p$sigma2, var_floor)
        G <- e_step(x, p$pi, p$mu, s2)$gamma
        g_a <- colSums(G)[-1] - length(x) * p$pi[-1]
        g_mu <- sapply(seq_len(K), function(k) {
            sum(G[, k] * (x - p$mu[k])) / s2[k]
        })
        g_ls <- sapply(seq_len(K), function(k) {
            sum(G[, k] * ((x - p$mu[k])^2 / s2[k] - 1))
        })
        -c(g_a, g_mu, g_ls)
    }
    list(f = f, gr = gr, env = env)
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
        loglik = NA_real_, cost = NA_real_, time_ms = NA_real_,
        converged = FALSE,
        trace = data.frame(cost = numeric(0), loglik = numeric(0)),
        path = data.frame(mu1 = numeric(0), mu2 = numeric(0))
    )
}

# Run one method from given starting values.
# Returns the final parameters, the TRUE log-likelihood recomputed with our own
# function (so every method is judged by the same yardstick), the cost, the
# elapsed time and a log-likelihood trace (best so far vs cost).
run_method <- function(
    method, x, init, K = 2, maxit = 5000, tol = 1e-10,
    var_floor = (0.01 * stats::sd(x))^2, record_path = FALSE
) {
    theta0 <- pack_par(init$pi, init$mu, init$sigma2)
    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch(
        {
            path <- data.frame(mu1 = numeric(0), mu2 = numeric(0))
            if (method == "EM (ours)") {
                fit <- em_gmm(x, init, tol = tol, max_iter = maxit, var_floor = var_floor)
                elapsed <- proc.time()[["elapsed"]] - t0
                par <- list(pi = fit$pi, mu = fit$mu, sigma2 = fit$sigma2)
                trace <- data.frame(cost = seq_along(fit$loglik_hist), loglik = fit$loglik_hist)
                cost <- fit$iterations
                converged <- fit$status == "converged"
                if (record_path) {
                    path <- em_path(x, init, n_iter = min(fit$iterations, 200), var_floor = var_floor)
                    path <- path[, c("mu1", "mu2")]
                }
            } else if (method == "mixtools::normalmixEM") {
                fit <- mixtools::normalmixEM(
                    x, lambda = init$pi, mu = init$mu, sigma = sqrt(init$sigma2),
                    k = K, epsilon = tol, maxit = maxit, verb = FALSE
                )
                elapsed <- proc.time()[["elapsed"]] - t0
                par <- list(pi = fit$lambda, mu = fit$mu, sigma2 = fit$sigma^2)
                trace <- data.frame(cost = seq_along(fit$all.loglik), loglik = fit$all.loglik)
                cost <- length(fit$all.loglik) - 1
                converged <- cost < maxit
            } else if (method == "mclust::Mclust") {
                # mclust chooses its own (hierarchical) initialisation
                fit <- mclust::Mclust(x, G = K, modelNames = "V", verbose = FALSE)
                elapsed <- proc.time()[["elapsed"]] - t0
                par <- list(
                    pi = fit$parameters$pro,
                    mu = fit$parameters$mean,
                    sigma2 = fit$parameters$variance$sigmasq
                )
                trace <- data.frame(cost = numeric(0), loglik = numeric(0))
                cost <- NA_real_
                converged <- TRUE
            } else {
                tr <- make_traced(x, K, var_floor)
                if (method == "nlm") {
                    res <- stats::nlm(tr$f, theta0, iterlim = maxit)
                    theta <- res$estimate
                    converged <- res$code %in% c(1, 2)
                } else {
                    res <- switch(method,
                        "Nelder-Mead" = stats::optim(
                            theta0, tr$f, method = "Nelder-Mead",
                            control = list(maxit = maxit, reltol = tol)
                        ),
                        "BFGS (num. grad)" = stats::optim(
                            theta0, tr$f, method = "BFGS",
                            control = list(maxit = maxit, reltol = tol)
                        ),
                        "BFGS (analytic grad)" = stats::optim(
                            theta0, tr$f, tr$gr, method = "BFGS",
                            control = list(maxit = maxit, reltol = tol)
                        ),
                        "L-BFGS-B" = stats::optim(
                            theta0, tr$f, method = "L-BFGS-B",
                            control = list(maxit = maxit, factr = tol / .Machine$double.eps)
                        ),
                        stop("unknown method: ", method)
                    )
                    theta <- res$par
                    converged <- res$convergence == 0
                }
                elapsed <- proc.time()[["elapsed"]] - t0
                p <- unpack_par(theta, K)
                par <- list(pi = p$pi, mu = p$mu, sigma2 = pmax(p$sigma2, var_floor))
                trace <- data.frame(cost = tr$env$cost_v, loglik = cummax(tr$env$ll_v))
                cost <- tr$env$cost
                if (record_path && length(tr$env$ll_v) > 0) {
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
                cost = cost, time_ms = 1000 * elapsed, converged = converged,
                trace = trace, path = path
            )
        },
        error = function(e) failed_run(method)
    )
    out
}

summarise_run <- function(r) {
    p <- r$par
    data.frame(
        method = r$method,
        loglik = r$loglik,
        pi1 = p$pi[1], mu1 = p$mu[1], mu2 = p$mu[2],
        sigma1 = sqrt(p$sigma2[1]), sigma2 = sqrt(p$sigma2[2]),
        cost = r$cost, time_ms = r$time_ms, converged = r$converged
    )
}

# ------------------------------------------------------------------
# 1. One common starting point: results table
# ------------------------------------------------------------------
data(faithful)
x <- faithful$waiting
set.seed(42)
init_km <- init_kmeans(x, 2)

runs <- lapply(all_methods, function(m) run_method(m, x, init_km))
names(runs) <- all_methods

res_tab <- bind_rows(lapply(runs, summarise_run))
ll_star <- max(res_tab$loglik, na.rm = TRUE)
res_tab$`gap to best` <- ll_star - res_tab$loglik

kable(
    res_tab,
    digits = 4,
    caption = "Table 7. Final fit of each method from the same k-means start (Faithful waiting times). Cost = number of E-step-equivalent evaluations."
)

# ------------------------------------------------------------------
# 2. Convergence race: log-likelihood gap vs computational cost
# ------------------------------------------------------------------
trace_df <- bind_rows(lapply(runs[traced_methods], function(r) {
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
        x = "Computational cost (likelihood / gradient evaluations, log scale)",
        y = "Log-likelihood gap to best (log scale)",
        colour = NULL,
        title = "Convergence race from the same starting point",
        subtitle = "Lower-left is better. EM is monotone; general optimizers can jump around."
    ) +
    guides(colour = guide_legend(nrow = 2))

# ------------------------------------------------------------------
# 3. Do all methods agree? Estimates by method
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

# Cost and time per method
cost_df <- res_tab %>%
    select(method, cost, time_ms) %>%
    pivot_longer(c(cost, time_ms), names_to = "measure", values_to = "value") %>%
    mutate(measure = recode(measure, cost = "Evaluations (E-step equivalents)", time_ms = "Elapsed time (ms)"))

ggplot(cost_df, aes(reorder(method, value), value, fill = method)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~measure, scales = "free_x") +
    scale_fill_manual(values = method_cols, guide = "none") +
    labs(x = NULL, y = NULL, title = "Cost to reach the optimum (one run)")

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

path_df <- bind_rows(lapply(core_methods, function(m) {
    r <- run_method(m, x, init_poor, record_path = TRUE)
    if (nrow(r$path) == 0) {
        return(NULL)
    }
    data.frame(method = m, step = seq_len(nrow(r$path)), r$path)
})) %>%
    mutate(method = factor(method, levels = core_methods))

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
        subtitle = "Slice with pi and sigma fixed at the MLE. Star = start, circle = end. Only improving steps are drawn."
    )

# ------------------------------------------------------------------
# 5. Robustness: 100 random starts
# ------------------------------------------------------------------
R_starts <- 100
set.seed(2026)
inits <- lapply(seq_len(R_starts), function(i) init_random(x, 2))

rob <- bind_rows(lapply(seq_len(R_starts), function(i) {
    bind_rows(lapply(traced_methods, function(m) {
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
        `median time (ms)` = median(time_ms, na.rm = TRUE),
        .groups = "drop"
    )
kable(
    rob_sum,
    digits = 1,
    caption = "Table 8. Robustness over 100 random starts (success = final log-likelihood within 1e-3 of the best found by any method)."
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
    labs(x = NULL, y = "Cost (evaluations, log scale)", title = "Cost from random starts")

# ------------------------------------------------------------------
# 6. Timing benchmark (30 repeated fits from the same start)
# ------------------------------------------------------------------
n_reps <- 30
timing <- bind_rows(lapply(seq_len(n_reps), function(i) {
    bind_rows(lapply(all_methods, function(m) {
        r <- run_method(m, x, init_km)
        data.frame(method = m, rep = i, time_ms = r$time_ms)
    }))
}))

ggplot(timing, aes(reorder(method, time_ms, median), time_ms, fill = method)) +
    geom_boxplot(outlier.alpha = 0.4) +
    scale_y_log10() +
    scale_fill_manual(values = method_cols, guide = "none") +
    coord_flip() +
    labs(
        x = NULL, y = "Elapsed time per fit (ms, log scale)",
        title = "Speed benchmark (30 repetitions, Faithful data)"
    )

# ------------------------------------------------------------------
# 7. Our EM vs library EM implementations on simulated data
#    (same starting values for mixtools; mclust uses its own start)
# ------------------------------------------------------------------
if (length(lib_methods) > 0) {
    lib_cmp <- bind_rows(lapply(c(0.5, 1.5, 3, 5), function(d) {
        bind_rows(lapply(1:30, function(r) {
            set.seed(42 + r)
            z <- rbinom(500, 1, 0.5)
            xs <- rnorm(500, mean = d * z)
            init <- init_kmeans(xs, 2)
            ours <- run_method("EM (ours)", xs, init)
            bind_rows(lapply(lib_methods, function(m) {
                theirs <- run_method(m, xs, init)
                data.frame(
                    delta = d, rep = r, library = m,
                    ll_diff = ours$loglik - theirs$loglik,
                    par_diff = max(abs(c(
                        ours$par$pi - theirs$par$pi,
                        ours$par$mu - theirs$par$mu,
                        sqrt(ours$par$sigma2) - sqrt(theirs$par$sigma2)
                    ))),
                    iter_ours = ours$cost,
                    iter_lib = theirs$cost
                )
            }))
        }))
    }))
    set.seed(42)

    print(
        ggplot(lib_cmp, aes(factor(delta), ll_diff, fill = library)) +
            geom_hline(yintercept = 0, linetype = 2) +
            geom_boxplot(outlier.alpha = 0.4) +
            scale_fill_manual(values = method_cols) +
            labs(
                x = expression(Delta * mu), y = "loglik(ours) - loglik(library)", fill = NULL,
                title = "Our EM vs library implementations: log-likelihood difference",
                subtitle = "Positive = ours found a higher likelihood. Zero = same optimum."
            )
    )

    print(
        ggplot(lib_cmp, aes(factor(delta), par_diff + 1e-16, fill = library)) +
            geom_boxplot(outlier.alpha = 0.4) +
            scale_y_log10() +
            scale_fill_manual(values = method_cols) +
            labs(
                x = expression(Delta * mu), y = "Max abs. parameter difference (log scale)", fill = NULL,
                title = "Our EM vs library implementations: parameter agreement"
            )
    )

    if (has_mixtools) {
        print(
            ggplot(
                lib_cmp %>% filter(library == "mixtools::normalmixEM"),
                aes(iter_lib, iter_ours, colour = factor(delta))
            ) +
                geom_abline(slope = 1, intercept = 0, linetype = 2) +
                geom_point(size = 2, alpha = 0.7) +
                scale_x_log10() +
                scale_y_log10() +
                scale_colour_brewer(palette = "Dark2") +
                labs(
                    x = "mixtools iterations", y = "Our EM iterations", colour = expression(Delta * mu),
                    title = "Iterations to converge: ours vs mixtools",
                    subtitle = "Stopping rules differ (relative vs absolute log-likelihood change), so exact equality is not expected."
                )
        )
    }
}
