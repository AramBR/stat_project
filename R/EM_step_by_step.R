# ==================================================================
# How EM works, step by step (Old Faithful waiting times, K = 2)
#
# Each EM iteration has two steps:
#   E-step: for every point, compute the probability that it belongs to
#           each component (its "responsibility"), given the current fit.
#   M-step: refit each component as a weighted Gaussian, using those
#           probabilities as weights (weighted mean, weighted variance,
#           share of total weight = mixing proportion).
# Repeating the two steps never lowers the log-likelihood.
#
# The figures show the fit at the start, at the points where a quarter, half,
# three quarters and 95% of the total improvement has been made, and at the
# end, starting from a deliberately poor guess so the changes are visible.
# ==================================================================
library(dplyr)
library(ggplot2)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))
col_c1 <- "#1b9e77" # component 1 (short waits)
col_c2 <- "#d95f02" # component 2 (long waits)

data(faithful)
x <- faithful$waiting

# Poor start: two almost identical components in the middle of the data
start <- list(pi = c(0.5, 0.5), mu = c(70, 72), sigma2 = c(13, 13)^2)

# em_trace keeps the parameters after every iteration (row 1 = start)
tr <- em_trace(x, start, tol = 1e-8)
n_iter <- nrow(tr) - 1

# Which iterations to draw: the first iteration that has made 25%, 50%, 75%
# and 95% of the total log-likelihood gain (plus the start and the end).
# From this start EM barely moves at first and then changes quickly, so fixed
# iteration numbers would show several identical panels.
ll <- tr[, "loglik"]
progress <- (ll - ll[1]) / (ll[n_iter + 1] - ll[1])
shown <- unique(c(
    0,
    vapply(c(0.25, 0.5, 0.75, 0.95), function(p) tr[which(progress >= p)[1], "iter"], numeric(1)),
    n_iter
))

# Parameters of the fit after iteration t, as a list
fit_at <- function(t) {
    p <- tr[t + 1, ]
    list(
        pi = unname(p[c("pi1", "pi2")]),
        mu = unname(p[c("mu1", "mu2")]),
        sigma2 = unname(p[c("sigma2_1", "sigma2_2")]),
        loglik = unname(p["loglik"])
    )
}

panel_label <- function(t) {
    f <- fit_at(t)
    lab <- if (t == n_iter) sprintf("Iteration %d (converged)", t) else sprintf("Iteration %d", t)
    sprintf("%s\nlog-likelihood = %.1f", lab, f$loglik)
}
panel_levels <- vapply(shown, panel_label, character(1))

# Fitted curves for every shown iteration
grid <- seq(35, 105, length.out = 400)
curves <- bind_rows(lapply(shown, function(t) {
    f <- fit_at(t)
    c1 <- f$pi[1] * dnorm(grid, f$mu[1], sqrt(f$sigma2[1]))
    c2 <- f$pi[2] * dnorm(grid, f$mu[2], sqrt(f$sigma2[2]))
    bind_rows(
        data.frame(x = grid, y = c1, curve = "Component 1"),
        data.frame(x = grid, y = c2, curve = "Component 2"),
        data.frame(x = grid, y = c1 + c2, curve = "Mixture")
    ) %>% mutate(panel = panel_label(t))
})) %>%
    mutate(panel = factor(panel, levels = panel_levels))

# Responsibilities (E-step) for every shown iteration:
# probability that each point belongs to component 2
resp <- bind_rows(lapply(shown, function(t) {
    f <- fit_at(t)
    g <- e_step(x, f$pi, f$mu, f$sigma2)$gamma
    data.frame(x = x, p2 = g[, 2], panel = panel_label(t))
})) %>%
    mutate(panel = factor(panel, levels = panel_levels))

# ------------------------------------------------------------------
# Figure 1. The fit after each shown iteration
# ------------------------------------------------------------------
ggplot() +
    geom_histogram(
        data = faithful, aes(waiting, after_stat(density)),
        binwidth = 2, fill = "grey88", colour = "white"
    ) +
    geom_line(
        data = curves %>% filter(curve != "Mixture"),
        aes(x, y, group = curve), colour = "grey35", linetype = 2, linewidth = 0.6
    ) +
    geom_line(
        data = curves %>% filter(curve == "Mixture"),
        aes(x, y), colour = "black", linewidth = 1
    ) +
    geom_point(
        data = resp, aes(x, y = -0.005, colour = p2),
        position = position_jitter(width = 0, height = 0.002, seed = 1),
        size = 1.3, alpha = 0.8
    ) +
    facet_wrap(~panel, ncol = 3) +
    scale_colour_gradient2(
        low = col_c1, mid = "grey75", high = col_c2, midpoint = 0.5,
        limits = c(0, 1), name = "P(component 2)"
    ) +
    labs(
        x = "Waiting time (minutes)", y = "Density",
        title = "EM, step by step",
        subtitle = paste(
            "Dashed = the two components, solid = their sum. Dots = data, coloured by the E-step:",
            "green = surely short wait, orange = surely long wait, grey = unsure.",
            sep = "\n"
        )
    )

# ------------------------------------------------------------------
# Figure 2. The log-likelihood rises at every step
# ------------------------------------------------------------------
ll_df <- data.frame(iter = tr[, "iter"], loglik = tr[, "loglik"])

ggplot(ll_df, aes(iter, loglik)) +
    geom_line(linewidth = 0.9, colour = "grey40") +
    geom_point(size = 1.2, colour = "grey40") +
    geom_point(
        data = ll_df %>% filter(iter %in% shown),
        size = 3.5, shape = 21, fill = col_c2, colour = "black"
    ) +
    labs(
        x = "Iteration", y = "Log-likelihood",
        title = "The log-likelihood never goes down",
        subtitle = "Orange points = the iterations drawn in the previous figure."
    )
