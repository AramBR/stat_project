library(dplyr)
library(tidyr)
library(ggplot2)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))
crit_cols <- c(AIC = "#7570b3", BIC = "#d95f02")

data(faithful)
x <- faithful$waiting
set.seed(42)

sel <- select_K(x, K_range = 1:6, n_random = 10)

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
