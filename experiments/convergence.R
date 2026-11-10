# ==================================================================
# Convergence rate analysis of EM
#
# Key idea: EM converges linearly. Near the optimum, every iteration shrinks
# the remaining error in the parameters by (roughly) the same factor r:
#     error_{t+1} ~ r * error_t        (0 < r < 1)
# so on a log scale the error is a straight line with slope log(r).
#   r = 0.5  -> ~3.3 iterations per extra correct digit
#   r = 0.99 -> ~230 iterations per extra correct digit
# Overlapping components make the class labels uncertain, which slows EM down.
# ==================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))

# ------------------------------------------------------------------
# 1. Faithful: error curves from several starting points
# ------------------------------------------------------------------
data(faithful)
x <- faithful$waiting
set.seed(42)

starts <- list(
    "good guess" = list(pi = c(.35, .65), mu = c(55, 80), sigma2 = c(6, 6)^2),
    "k-means" = init_kmeans(x, 2),
    "wide & far" = list(pi = c(.5, .5), mu = c(30, 110), sigma2 = c(20, 20)^2),
    "close means" = list(pi = c(.5, .5), mu = c(70, 72), sigma2 = c(13, 13)^2),
    "unbalanced" = list(pi = c(.98, .02), mu = c(78, 96), sigma2 = c(14, 1)^2)
)

traces <- lapply(starts, function(s) em_trace(x, s, max_iter = 3000))
rates <- lapply(traces, em_rate)

# Common reference optimum theta* = final parameters of the run with the
# highest final log-likelihood
final_ll <- vapply(traces, function(tr) unname(tr[nrow(tr), "loglik"]), numeric(1))
best_tr <- traces[[which.max(final_ll)]]
theta_star <- best_tr[nrow(best_tr), -(1:2)]

# Parameter error at every iteration: distance between theta_t and theta*
trace_df <- bind_rows(lapply(names(traces), function(nm) {
    tr <- traces[[nm]]
    P <- tr[, -(1:2), drop = FALSE]
    err <- sqrt(rowSums(sweep(P, 2, theta_star)^2))
    data.frame(start = nm, iter = tr[, "iter"], error = err)
}))

# Summary table: iterations, estimated rate, iterations per extra digit
rate_tab <- bind_rows(lapply(names(traces), function(nm) {
    r <- rates[[nm]]$rate
    data.frame(
        Start = nm,
        Iterations = nrow(traces[[nm]]) - 1,
        Status = attr(traces[[nm]], "status"),
        `Rate r` = r,
        `Iterations per digit` = -1 / log10(r),
        `Final loglik` = unname(final_ll[nm]),
        check.names = FALSE,
        row.names = NULL
    )
}))
kable(
    rate_tab,
    digits = 4,
    caption = "Table 11. Estimated linear convergence rate of EM on the Faithful data from five starting points."
)

# Straight-line asymptotes implied by the estimated rates
asym_df <- bind_rows(lapply(names(traces), function(nm) {
    r <- rates[[nm]]$rate
    d <- trace_df %>% filter(start == nm, error > 1e-6)
    if (is.na(r) || nrow(d) < 5) {
        return(NULL)
    }
    anchor <- d[round(nrow(d) * 0.7), ]
    data.frame(
        start = nm,
        iter = d$iter,
        error = anchor$error * r^(d$iter - anchor$iter)
    )
}))

# 1a. Error vs iteration (log y): straight lines = linear convergence
# (runs that stopped immediately, like the degenerate "unbalanced" start,
# have no line to draw and are left out)
plot_df <- trace_df %>%
    group_by(start) %>%
    filter(n() > 1, error > 1e-7) %>%
    ungroup()
ggplot(plot_df, aes(iter, error, colour = start)) +
    geom_line(linewidth = 1) +
    geom_line(
        data = asym_df, aes(iter, error, group = start),
        colour = "black", linetype = 3, linewidth = 0.6
    ) +
    scale_y_log10(
        breaks = 10^seq(-6, 3, by = 3),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
        x = "Iteration", y = "Distance to the optimal parameters (log scale)", colour = "Start",
        title = "EM converges linearly",
        subtitle = "Straight lines on a log scale. Dotted = slope implied by the estimated rate r."
    )

# 1b. Ratio of successive parameter steps settles at the rate
ratio_df <- bind_rows(lapply(names(rates), function(nm) {
    if (nrow(rates[[nm]]$ratio) == 0) {
        return(NULL)
    }
    data.frame(start = nm, rates[[nm]]$ratio)
}))
rate_lines <- rate_tab %>%
    filter(!is.na(`Rate r`)) %>%
    transmute(start = Start, rate = `Rate r`)

ggplot(ratio_df, aes(iter, ratio, colour = start)) +
    geom_hline(
        data = rate_lines, aes(yintercept = rate, colour = start),
        linetype = 2, linewidth = 0.5
    ) +
    geom_line(linewidth = 0.9) +
    coord_cartesian(ylim = c(0, 1.05)) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
        x = "Iteration",
        y = "Step size / previous step size",
        colour = "Start",
        title = "Ratio of successive parameter steps",
        subtitle = "It levels off at the convergence rate r (dashed). The rate belongs to the optimum, not to the start."
    )

# ------------------------------------------------------------------
# 2. Simulation: how does the rate depend on the separation?
# ------------------------------------------------------------------
deltas <- c(0.5, 1, 1.5, 2, 3, 5)
n_rate <- 500
R_rate <- 30 # replicates per delta  -> lower while testing
max_iter_rate <- 5000
true_pi <- 0.5

run_rate_rep <- function(d, r) {
    set.seed(42 + r)
    z <- rbinom(n_rate, 1, 1 - true_pi)
    xs <- rnorm(n_rate, mean = d * z)
    tr <- em_trace(xs, init_kmeans(xs, 2), max_iter = max_iter_rate)
    rt <- em_rate(tr)
    last <- tr[nrow(tr), ]
    gam <- e_step(
        xs,
        pi = unname(last[c("pi1", "pi2")]),
        mu = unname(last[c("mu1", "mu2")]),
        sigma2 = unname(last[c("sigma2_1", "sigma2_2")])
    )$gamma
    data.frame(
        delta = d, rep = r,
        rate = rt$rate,
        iterations = nrow(tr) - 1,
        converged = attr(tr, "status") == "converged",
        # average uncertainty of the class labels: 0 = labels clear, 0.5 = coin flips
        uncert = mean(1 - apply(gam, 1, max))
    )
}

rate_sim <- bind_rows(lapply(deltas, function(d) {
    bind_rows(lapply(seq_len(R_rate), function(r) run_rate_rep(d, r)))
}))
set.seed(42)

# The rate is NA when EM converges too fast to leave 3 usable ratios.
# Those runs are dropped from the medians, so the share is reported.
rate_sum <- rate_sim %>%
    group_by(`Delta mu` = delta) %>%
    summarise(
        `median rate` = median(rate, na.rm = TRUE),
        `iter. per digit` = -1 / log10(median(rate, na.rm = TRUE)),
        `median iterations` = median(iterations),
        `label uncertainty` = mean(uncert),
        `% not converged` = 100 * mean(!converged),
        `% rate NA` = 100 * mean(is.na(rate)),
        .groups = "drop"
    )
kable(
    rate_sum,
    digits = 3,
    caption = sprintf(
        "Table 12. EM convergence rate by component separation (n = %d, R = %d replicates per separation).",
        n_rate, R_rate
    )
)

# 2a. Rate vs separation
ggplot(rate_sim %>% filter(!is.na(rate)), aes(factor(delta), rate)) +
    geom_boxplot(fill = "grey85", outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, colour = "#1b9e77") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
        x = expression(Delta * mu ~ "(separation)"), y = "Convergence rate r (closer to 1 = slower)",
        title = "EM slows down when components overlap",
        subtitle = sprintf("n = %d, one point per replicate", n_rate)
    )

# 2b. Why: rate vs uncertainty about the class labels
ggplot(rate_sim %>% filter(!is.na(rate)), aes(uncert, rate, colour = factor(delta))) +
    geom_point(size = 2, alpha = 0.7) +
    scale_colour_brewer(palette = "Dark2") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
        x = "Average uncertainty of the class labels (1 - max responsibility)",
        y = "Convergence rate r", colour = expression(Delta * mu),
        title = "More label uncertainty = slower EM",
        subtitle = "When the E-step cannot tell which component a point belongs to, each step learns less."
    )

# 2c. What the rate means in practice: iterations to converge
ggplot(rate_sim, aes(factor(delta), iterations)) +
    geom_boxplot(fill = "grey85", outlier.alpha = 0.4) +
    scale_y_log10() +
    labs(
        x = expression(Delta * mu), y = "Iterations to converge (log scale)",
        title = sprintf("Iterations needed (tolerance 1e-14, n = %d)", n_rate)
    )
