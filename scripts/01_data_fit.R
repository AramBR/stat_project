library(dplyr)
library(knitr)
library(ggplot2)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))
set.seed(42) # init_kmeans uses random k-means starts

col2 <- c("Component 1" = "#1b9e77", "Component 2" = "#d95f02")

# Load the data and pull out the waiting times as a vector
data(faithful)
x <- faithful$waiting
n <- length(x)

# Descriptive statistics table (n, mean, sd, quartiles, skewness)
desc <- faithful %>%
    summarise(
        n = n(),
        mean = mean(waiting),
        sd = sd(waiting),
        min = min(waiting),
        Q1 = quantile(waiting, .25),
        median = median(waiting),
        Q3 = quantile(waiting, .75),
        max = max(waiting),
        # moment skewness: third central moment over (second central moment)^1.5,
        # both with divisor n
        skewness = mean((waiting - mean(waiting))^3) /
            mean((waiting - mean(waiting))^2)^1.5
    )
kable(
    desc,
    digits = 2,
    caption = "Table 1. Descriptive statistics of the waiting time (minutes)."
)

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

# Initialise with k-means (2 clusters) and run the scratch EM algorithm
init_km <- init_kmeans(x, 2)
fit2 <- em_gmm(x, init_km)

# Table of fitted mixture weights, means, and standard deviations
est <- data.frame(
    Component = c("1 (short waits)", "2 (long waits)"),
    pi = fit2$pi,
    mu = fit2$mu,
    sigma = sqrt(fit2$sigma2)
)
kable(
    est,
    digits = 3,
    caption = sprintf(
        "Table 2. Scratch EM estimates (log-likelihood = %.3f, %d iterations, status: %s).",
        fit2$loglik,
        fit2$iterations,
        fit2$status
    )
)

# Sanity check: EM guarantees the log-likelihood never decreases
cat("Log-likelihood never decreased:", fit2$monotone, "\n")

# Build weighted component densities on a grid and hard-assign each point to a component
grid <- seq(35, 105, length.out = 500)
comp_df <- bind_rows(lapply(1:2, function(k) {
    data.frame(
        x = grid,
        y = fit2$pi[k] * dnorm(grid, fit2$mu[k], sqrt(fit2$sigma2[k])),
        component = paste("Component", k)
    )
}))
cls <- data.frame(
    x = x,
    component = paste("Component", max.col(fit2$gamma, ties.method = "first"))
)

# Plot histogram, each component (dashed), the total mixture (solid), and coloured class assignments
ggplot() +
    geom_histogram(
        data = faithful,
        aes(waiting, after_stat(density)),
        binwidth = 2,
        fill = "grey85",
        colour = "white"
    ) +
    geom_line(
        data = comp_df,
        aes(x, y, colour = component),
        linewidth = 1,
        linetype = 2
    ) +
    geom_line(
        data = data.frame(x = grid, y = gmm_density(grid, fit2)),
        aes(x, y),
        linewidth = 1.1
    ) +
    geom_point(
        data = cls,
        aes(x, y = -0.002, colour = component),
        position = position_jitter(height = 0.0012, width = 0),
        alpha = 0.6,
        size = 1.2
    ) +
    scale_colour_manual(values = col2) +
    labs(x = "Waiting time (minutes)", y = "Density", colour = NULL)

# Six different starting configurations (good, k-means, far apart, close means, unbalanced, saddle point)
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

# Table comparing final log-likelihood, iterations, convergence status, and estimates per start
tab1 <- bind_rows(lapply(names(runs1), function(nm) {
    f <- runs1[[nm]]
    data.frame(
        Start = nm,
        `Final loglik` = f$loglik,
        Iterations = f$iterations,
        Status = f$status,
        Monotone = f$monotone,
        pi1 = f$pi[1],
        mu1 = f$mu[1],
        mu2 = f$mu[2],
        sigma1 = sqrt(f$sigma2[1]),
        sigma2 = sqrt(f$sigma2[2]),
        check.names = FALSE
    )
}))
kable(
    tab1,
    digits = 3,
    caption = "Table 3. EM results from six starting configurations (faithful waiting times)."
)

# Log-likelihood trajectories for each start (log x-axis) to compare convergence speed
traj <- bind_rows(lapply(names(runs1), function(nm) {
    data.frame(
        Start = nm,
        # loglik_hist[1] is at the starting values, so this axis is
        # (number of M-steps + 1); the + 1 keeps the log scale defined
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
