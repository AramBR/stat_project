# ==================================================================
# Simulation: does EM recover the true parameters?
#
# We generate data where the truth is KNOWN,
#   x ~ 0.5 * N(0, 1) + 0.5 * N(delta, 1),
# fit EM, and compare the estimates with the truth. delta is the distance
# between the two means: small delta = the groups overlap, large delta =
# the groups are clearly separate.
# ==================================================================
library(dplyr)
library(tidyr)
library(knitr)
library(ggplot2)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))

# ---- One dataset ---------------------------------------------------

# Simulate one dataset, fit EM from a k-means start, return one row of results
run_rep <- function(d, r, n, p1) {
    set.seed(42 + r) # dataset r uses the same random draws for every delta
    z <- rbinom(n, 1, 1 - p1) # true group of each point: 0 = group 1, 1 = group 2
    xs <- rnorm(n, mean = d * z, sd = 1)

    f <- em_gmm(xs, init_kmeans(xs, 2), max_iter = 5000)

    # Assign every point to the component it most probably belongs to
    pred <- max.col(f$gamma, ties.method = "first") - 1
    mis <- mean(pred != z)

    # Log-likelihood at the TRUE parameters, for comparison with EM's fit
    ll_true <- gmm_loglikelihood(xs, c(p1, 1 - p1), c(0, d), c(1, 1))

    data.frame(
        delta = d,
        rep = r,
        status = f$status,
        iterations = f$iterations,
        pi1 = f$pi[1],
        mu1 = f$mu[1],
        mu2 = f$mu[2],
        s1 = sqrt(f$sigma2[1]),
        s2 = sqrt(f$sigma2[2]),
        # EM names its components by their means, so with heavy overlap it can
        # call group 1 "component 2" (all labels reversed). min(mis, 1 - mis)
        # counts the errors without caring which name each group got.
        misclass = min(mis, 1 - mis),
        # how much better EM's fit explains this data than the true parameters
        ll_gap = f$loglik - ll_true,
        monotone = f$monotone
    )
}

# ---- Run -----------------------------------------------------------

deltas <- c(0.5, 1.5, 3.0, 5.0)
R2 <- 100 # datasets per delta
n2 <- 500 # points per dataset
true_pi <- 0.5

grid <- expand.grid(d = deltas, r = seq_len(R2))
sim2 <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    run_rep(grid$d[i], grid$r[i], n2, true_pi)
}))
set.seed(42) # restore the global seed

# ---- Tables --------------------------------------------------------

# RMSE ("typical error") = sqrt(mean((estimate - truth)^2)):
# how far the estimate is from the truth on average, in the units of the parameter.
rmse <- function(est, truth) sqrt(mean((est - truth)^2))

# Table 4: how accurate are the estimates?
tab_est <- sim2 %>%
    group_by(`Delta mu` = delta) %>%
    summarise(
        `error pi1` = rmse(pi1, true_pi),
        `error mu1` = rmse(mu1, 0),
        `error mu2` = rmse(mu2, first(delta)),
        `error sigma1` = rmse(s1, 1),
        `error sigma2` = rmse(s2, 1),
        .groups = "drop"
    )
kable(
    tab_est,
    digits = 3,
    caption = sprintf(
        "Table 4. Typical error (RMSE) of the EM estimates by separation (R = %d datasets, n = %d). Truth: pi1 = 0.5, mu1 = 0, mu2 = Delta mu, sigma1 = sigma2 = 1.",
        R2, n2
    )
)

# Table 5: classification, fit and convergence
tab_fit <- sim2 %>%
    group_by(`Delta mu` = delta) %>%
    summarise(
        `misclassified %` = 100 * mean(misclass),
        `log-lik gap` = mean(ll_gap),
        `median iterations` = median(iterations),
        `% not converged` = 100 * mean(status != "converged"),
        `% loglik decreased` = 100 * mean(!monotone),
        .groups = "drop"
    )
kable(
    tab_fit,
    digits = 2,
    caption = paste(
        "Table 5. Share of points put in the wrong group, log-likelihood of EM's fit minus that of the true parameters,",
        "iterations, and convergence. '% loglik decreased' should always be 0."
    )
)

# ---- Plots ---------------------------------------------------------

# Estimation errors (estimate - truth) per parameter and separation
err <- sim2 %>%
    mutate(
        delta_f = factor(delta),
        `pi1` = pi1 - true_pi,
        `mu1` = mu1 - 0,
        `mu2` = mu2 - delta,
        `sigma1` = s1 - 1,
        `sigma2` = s2 - 1
    ) %>%
    select(delta_f, rep, pi1, mu1, mu2, sigma1, sigma2) %>%
    pivot_longer(-c(delta_f, rep), names_to = "parameter", values_to = "error")

ggplot(err, aes(delta_f, error)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_boxplot(fill = "grey85", outlier.alpha = 0.4) +
    facet_wrap(~parameter, scales = "free_y") +
    labs(
        x = "\u0394\u03bc (true separation)",
        y = "Estimate - truth",
        title = "EM estimation error by component separation",
        subtitle = "Dashed line = exactly right. Narrow boxes around 0 = accurate estimates."
    )

# Misclassification rate by separation
ggplot(tab_fit, aes(`Delta mu`, `misclassified %`)) +
    geom_line(linewidth = 1, colour = "#d95f02") +
    geom_point(size = 2.5, colour = "#d95f02") +
    scale_y_continuous(limits = c(0, 50)) +
    labs(
        x = "\u0394\u03bc",
        y = "Points in the wrong group (%)",
        title = "EM classification error by component separation",
        subtitle = "50% would be a coin flip."
    )

# Iterations needed to converge (slower when components overlap)
ggplot(sim2, aes(factor(delta), iterations)) +
    geom_boxplot(fill = "grey85", outlier.alpha = 0.4) +
    scale_y_log10() +
    labs(
        x = "\u0394\u03bc",
        y = "EM iterations (log scale)",
        title = "Convergence speed by component separation"
    )

# Outcome of each run (converged / degenerate / max_iter reached)
ggplot(sim2, aes(factor(delta), fill = status)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(palette = "Set2") +
    labs(
        x = "\u0394\u03bc",
        y = "Share of datasets",
        fill = "Status",
        title = "Convergence status by component separation"
    )
