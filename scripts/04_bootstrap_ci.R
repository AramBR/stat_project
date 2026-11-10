# ==================================================================
# Simulation 3: bootstrap standard errors, confidence intervals, coverage
# Truth: x ~ p1 * N(0, 1) + (1 - p1) * N(delta, 1)
# ==================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))

# ------------------------------------------------------------------
# 0. Quick demo on the real data (Old Faithful)
# ------------------------------------------------------------------
set.seed(42)
x <- faithful$waiting
fit_f <- em_gmm(x, init_kmeans(x, 2))
ci_f <- gmm_ci(x, fit_f, B = 500)
kable(
    ci_f,
    digits = 3,
    caption = "Table 9. Bootstrap standard errors and 95% percentile CIs for the Faithful fit (B = 500)."
)

ggplot(ci_f %>% filter(parameter != "pi2"), aes(estimate, parameter)) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
    geom_point(size = 3, colour = "#1b9e77") +
    facet_wrap(~parameter, scales = "free", ncol = 1) +
    labs(
        x = "Estimate with 95% bootstrap CI", y = NULL,
        title = "Faithful: parameter estimates with bootstrap intervals"
    ) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())

# ------------------------------------------------------------------
# 1. Coverage simulation
#    For each dataset: fit EM, build bootstrap CIs, check if truth is inside.
# ------------------------------------------------------------------
deltas <- c(0.5, 1.5, 3.0, 5.0)
R_ci <- 50 # replicates per delta (raise for smoother coverage)
B_ci <- 100 # bootstrap resamples per replicate
n_ci <- 500
true_pi <- 0.5
level <- 0.95

# total EM fits = length(deltas) * R_ci * B_ci  -> lower R_ci / B_ci while testing

run_ci_rep <- function(d, r) {
    set.seed(42 + r)
    z <- rbinom(n_ci, 1, 1 - true_pi)
    xs <- rnorm(n_ci, mean = d * z)
    f <- em_gmm(xs, init_kmeans(xs, 2), max_iter = 5000)
    # A degenerate fit is not an estimate; skip the replicate (counted below)
    if (f$status == "degenerate") {
        return(NULL)
    }
    ci <- gmm_ci(xs, f, B = B_ci, level = level, max_iter = 1000)
    truth <- c(pi1 = true_pi, pi2 = 1 - true_pi, mu1 = 0, mu2 = d, sigma1 = 1, sigma2 = 1)
    ci$truth <- unname(truth[ci$parameter])
    ci$delta <- d
    ci$rep <- r
    ci
}

ci_sim <- bind_rows(lapply(deltas, function(d) {
    bind_rows(lapply(seq_len(R_ci), function(r) run_ci_rep(d, r)))
})) %>%
    filter(parameter != "pi2") %>% # redundant with pi1
    mutate(
        covered = lower <= truth & truth <= upper,
        parameter = factor(parameter, levels = c("pi1", "mu1", "mu2", "sigma1", "sigma2"))
    )
set.seed(42)

# ------------------------------------------------------------------
# 2. Summary table
#    coverage  : share of intervals containing the truth (target: 95%)
#    SE ratio  : average bootstrap SE / empirical SD of the estimates
#                (empirical SD over replicates is the "true" SE; ratio ~ 1 is good)
# ------------------------------------------------------------------
ci_tab <- ci_sim %>%
    group_by(`Delta mu` = delta, parameter) %>%
    summarise(
        R = n(), # replicates used (degenerate full-data fits are skipped)
        `boot failed %` = 100 * mean(n_failed) / B_ci,
        `boot not conv. %` = 100 * mean(n_not_converged) / B_ci,
        `coverage %` = 100 * mean(covered),
        `mean CI width` = mean(upper - lower),
        `mean boot SE` = mean(se),
        `empirical SD` = sd(estimate),
        `SE ratio` = mean(se) / sd(estimate),
        .groups = "drop"
    )
kable(
    ci_tab,
    digits = 3,
    caption = sprintf(
        "Table 10. Bootstrap CI performance (up to %d replicates, B = %d, n = %d, nominal level %.0f%%).",
        R_ci, B_ci, n_ci, 100 * level
    )
)

# ------------------------------------------------------------------
# 3. Plots
# ------------------------------------------------------------------

# 3a. Coverage vs separation, with the Monte Carlo noise band around 95%
mc_band <- 1.96 * sqrt(level * (1 - level) / R_ci)
ggplot(ci_tab, aes(`Delta mu`, `coverage %` / 100, colour = parameter)) +
    annotate(
        "rect", xmin = -Inf, xmax = Inf,
        ymin = level - mc_band, ymax = min(level + mc_band, 1),
        fill = "grey80", alpha = 0.5
    ) +
    geom_hline(yintercept = level, linetype = 2) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
        x = expression(Delta * mu), y = "Coverage of 95% bootstrap CI", colour = NULL,
        title = "Do the bootstrap intervals cover the truth?",
        subtitle = "Dashed line = nominal 95%. Grey band = Monte Carlo noise from a finite number of replicates."
    )

# 3b. Are the bootstrap SEs the right size?
ggplot(ci_tab, aes(`Delta mu`, `SE ratio`, colour = parameter)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
        x = expression(Delta * mu), y = "Bootstrap SE / empirical SD", colour = NULL,
        title = "Size of bootstrap standard errors",
        subtitle = "Above 1: SEs too large (conservative). Below 1: SEs too small (intervals too narrow)."
    )

# 3c. Caterpillar plot: every interval for one separation, miss = red
delta_show <- 1.5
cat_df <- ci_sim %>%
    filter(delta == delta_show) %>%
    group_by(parameter) %>%
    arrange(estimate, .by_group = TRUE) %>%
    mutate(order = row_number()) %>%
    ungroup()

ggplot(cat_df, aes(y = order)) +
    geom_vline(aes(xintercept = truth), linetype = 2) +
    geom_errorbar(
        aes(xmin = lower, xmax = upper, colour = covered),
        width = 0, linewidth = 0.6, orientation = "y"
    ) +
    geom_point(aes(x = estimate, colour = covered), size = 1.2) +
    facet_wrap(~parameter, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = c(`TRUE` = "grey55", `FALSE` = "#d7191c"), name = "Covers truth") +
    labs(
        x = "Estimate and 95% bootstrap CI", y = "Replicate (sorted by estimate)",
        title = sprintf("Every interval for Delta mu = %.1f", delta_show),
        subtitle = "Dashed line = true value. Red intervals miss it."
    ) +
    theme(axis.text.y = element_blank())
