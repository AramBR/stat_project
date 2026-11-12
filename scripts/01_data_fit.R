library(dplyr)
library(ggplot2)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))
set.seed(42) # init_kmeans uses random k-means starts

# Load the data and pull out the waiting times as a vector
data(faithful)
x <- faithful$waiting

# Histogram + kernel density + rug of the raw waiting times (shows the two modes)
ggplot(faithful, aes(waiting)) +
    geom_histogram(
        aes(y = after_stat(density)),
        binwidth = 2,
        fill = "grey80",
        colour = "white"
    ) +
    geom_density(colour = "black", linewidth = 1) +
    geom_rug(alpha = 0.3) +
    labs(
        x = "Waiting time to next eruption (minutes)",
        y = "Density",
        title = "Old Faithful waiting times"
    )

# Six different starting configurations (good, k-means, far apart, close means, unbalanced, saddle point)
init_km <- init_kmeans(x, 2)
starts <- list(
    "S1 good guess" = list(
        pi = c(.35, .65),
        mu = c(55, 80),
        sigma2 = c(6, 6)^2
    ),
    "S2 k-means" = init_km,
    "S3 wide & far" = list(
        pi = c(.5, .5),
        mu = c(30, 110),
        sigma2 = c(20, 20)^2
    ),
    "S4 close means" = list(
        pi = c(.5, .5),
        mu = c(70, 72),
        sigma2 = c(13, 13)^2
    ),
    "S5 unbalanced" = list(
        pi = c(.98, .02),
        mu = c(78, 96),
        sigma2 = c(14, 1)^2
    ),
    "S6 saddle" = list(
        pi = c(.5, .5),
        mu = rep(mean(x), 2),
        sigma2 = rep(var(x), 2)
    )
)

# Run EM from every starting point (up to 2000 iterations each)
runs1 <- lapply(starts, function(s) em_gmm(x, s, max_iter = 2000))

# Log-likelihood trajectories for each start (log x-axis) to compare convergence speed
traj <- bind_rows(lapply(names(runs1), function(nm) {
    data.frame(
        Start = nm,
        iteration = seq_along(runs1[[nm]]$loglik_hist),
        loglik = runs1[[nm]]$loglik_hist
    )
}))
ggplot(traj, aes(iteration, loglik, colour = Start)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 0.8) +
    scale_x_log10() +
    coord_cartesian(
        ylim = c(max(-1200, min(traj$loglik)), max(traj$loglik) + 5)
    ) +
    labs(
        x = "Iteration + 1 (log scale; 1 = starting values)",
        y = "Observed log-likelihood",
        colour = "Starting values"
    ) +
    scale_colour_brewer(palette = "Dark2")
