library(dplyr)
library(ggplot2)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))

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

# Full EM path from every start, and the estimated rate r of each
traces <- lapply(starts, function(s) em_trace(x, s, max_iter = 3000))
rates <- vapply(traces, em_rate, numeric(1))

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

# Straight-line asymptotes implied by the estimated rates
asym_df <- bind_rows(lapply(names(traces), function(nm) {
    r <- rates[[nm]]
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
