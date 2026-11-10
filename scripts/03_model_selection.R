# ==================================================================
# How many components? Choosing K with AIC and BIC
#
# Adding components always raises the log-likelihood (a bigger model can
# always fit at least as well), so the log-likelihood alone would always
# pick the largest K. AIC and BIC add a penalty for every extra parameter:
#   AIC = -2 * loglik + 2 * p
#   BIC = -2 * loglik + log(n) * p        (p = 3K - 1 parameters)
# Smaller is better. BIC's penalty is larger (log(500) = 6.2 > 2), so BIC
# prefers smaller models than AIC.
# ==================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))
crit_cols <- c(AIC = "#7570b3", BIC = "#d95f02")

# ------------------------------------------------------------------
# 1. Old Faithful: which K do AIC and BIC choose?
# ------------------------------------------------------------------
data(faithful)
x <- faithful$waiting
set.seed(42)

sel <- select_K(x, K_range = 1:6, n_random = 10)

kable(
    sel$table %>% select(K, npar, loglik, AIC, BIC, dAIC, dBIC, status),
    digits = 2,
    caption = sprintf(
        "Table 13. Model selection for the Faithful waiting times. AIC chooses K = %d, BIC chooses K = %d.",
        sel$K_aic, sel$K_bic
    )
)

# 1a. Criteria against K (the minimum is the chosen model)
crit_df <- sel$table %>%
    select(K, AIC, BIC) %>%
    pivot_longer(c(AIC, BIC), names_to = "criterion", values_to = "value")
best_df <- crit_df %>%
    group_by(criterion) %>%
    slice_min(value, n = 1)

ggplot(crit_df, aes(K, value, colour = criterion)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_point(data = best_df, size = 6, shape = 21, fill = NA, stroke = 1.3) +
    scale_x_continuous(breaks = 1:6) +
    scale_colour_manual(values = crit_cols) +
    labs(
        x = "Number of components K", y = "Criterion value (smaller is better)", colour = NULL,
        title = "Choosing K for the Faithful data",
        subtitle = "Circled = minimum = the K each criterion chooses."
    )

# 1b. Log-likelihood gain vs penalty: why the curve turns back up
gain_df <- sel$table %>%
    transmute(
        K,
        `-2 x log-likelihood` = -2 * loglik,
        `AIC penalty (2p)` = 2 * npar,
        `BIC penalty (log(n) p)` = log(length(x)) * npar
    ) %>%
    pivot_longer(-K, names_to = "part", values_to = "value") %>%
    mutate(panel = factor(
        ifelse(part == "-2 x log-likelihood", "Fit: keeps improving", "Penalty: keeps growing"),
        levels = c("Fit: keeps improving", "Penalty: keeps growing")
    ))

ggplot(gain_df, aes(K, value, colour = part)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~panel, scales = "free_y") +
    scale_x_continuous(breaks = 1:6) +
    scale_colour_manual(values = c("black", crit_cols[["AIC"]], crit_cols[["BIC"]])) +
    labs(
        x = "Number of components K", y = NULL, colour = NULL,
        title = "The two parts of AIC and BIC",
        subtitle = "After K = 2 the fit barely improves while the penalty keeps growing."
    ) +
    guides(colour = guide_legend(nrow = 1))

# 1c. What the fitted models look like
grid <- seq(35, 105, length.out = 400)
dens_df <- bind_rows(lapply(1:4, function(K) {
    data.frame(
        K = sprintf("K = %d", K),
        x = grid,
        y = gmm_density(grid, sel$fits[[K]])
    )
}))

ggplot() +
    geom_histogram(
        data = faithful, aes(waiting, after_stat(density)),
        binwidth = 2, fill = "grey85", colour = "white"
    ) +
    geom_line(data = dens_df, aes(x, y), linewidth = 1) +
    facet_wrap(~K, nrow = 1) +
    labs(
        x = "Waiting time (minutes)", y = "Density",
        title = "Fitted mixtures with 1 to 4 components",
        subtitle = "K = 1 misses the two modes; K = 3 and 4 add little visible improvement."
    )

# ------------------------------------------------------------------
# 2. Simulation: how often is the true K = 2 recovered?
#    Truth: 0.5 * N(0, 1) + 0.5 * N(delta, 1), candidates K = 1..4
# ------------------------------------------------------------------
deltas <- c(0.5, 1, 1.5, 2, 3)
R_sel <- 50 # replicates per delta (about 10 minutes in total; lower while testing)
n_sel <- 500
K_cand <- 1:4

run_sel_rep <- function(d, r) {
    set.seed(42 + r)
    z <- rbinom(n_sel, 1, 0.5)
    xs <- rnorm(n_sel, mean = d * z)
    # A looser tolerance (log-likelihood accurate to about 0.001) is plenty for
    # comparing AIC/BIC values, and much faster when components overlap
    s <- select_K(xs, K_range = K_cand, n_random = 3, tol = 1e-6, max_iter = 1000)
    data.frame(
        delta = d, rep = r,
        criterion = c("AIC", "BIC"),
        K_chosen = c(s$K_aic, s$K_bic)
    )
}

sel_sim <- bind_rows(lapply(deltas, function(d) {
    bind_rows(lapply(seq_len(R_sel), function(r) run_sel_rep(d, r)))
}))
set.seed(42)

sel_tab <- sel_sim %>%
    group_by(`Delta mu` = delta, criterion) %>%
    summarise(
        `% K = 1` = 100 * mean(K_chosen == 1),
        `% K = 2 (correct)` = 100 * mean(K_chosen == 2),
        `% K >= 3` = 100 * mean(K_chosen >= 3),
        .groups = "drop"
    )
kable(
    sel_tab,
    digits = 1,
    caption = sprintf(
        "Table 14. Share of replicates in which AIC / BIC choose each K (true K = 2, n = %d, R = %d).",
        n_sel, R_sel
    )
)

# 2a. Which K is chosen, by separation
ggplot(sel_sim, aes(factor(delta), fill = factor(K_chosen, levels = K_cand))) +
    geom_bar(position = "fill") +
    facet_wrap(~criterion) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(
        values = c(`1` = "#bdbdbd", `2` = "#1b9e77", `3` = "#e6ab02", `4` = "#d95f02"),
        drop = FALSE
    ) +
    labs(
        x = expression(Delta * mu ~ "(separation)"), y = "Share of replicates",
        fill = "K chosen",
        title = "Can AIC and BIC find the two components?",
        subtitle = "True K = 2 (green). Overlapping components look like a single Gaussian,\nso K = 1 wins when the separation is small."
    )

# 2b. Success rate side by side
ggplot(sel_tab, aes(`Delta mu`, `% K = 2 (correct)` / 100, colour = criterion)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_colour_manual(values = crit_cols) +
    labs(
        x = expression(Delta * mu), y = "Correct K chosen", colour = NULL,
        title = "BIC is cautious, AIC is generous",
        subtitle = "BIC needs more separation to detect the 2nd component;\nAIC finds it earlier but sometimes picks too many components."
    )
