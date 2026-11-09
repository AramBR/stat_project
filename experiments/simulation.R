library(dplyr)
library(tidyr)
library(knitr)
library(ggplot2)
source("R/core.R")
source("R/helpers.R")

# ------------------------------------------------------------------
# Simulation 2: parameter recovery of EM by component separation
# Truth: x ~ p1 * N(0, 1) + (1 - p1) * N(delta, 1)
# ------------------------------------------------------------------

# ---- Helpers -------------------------------------------------------

# Bayes-optimal misclassification rate for the 2-component mixture
# (best any classifier can do knowing the true parameters)
bayes_error <- function(delta, p1, s1 = 1, s2 = 1) {
    f <- function(x) {
        pmin(p1 * dnorm(x, 0, s1), (1 - p1) * dnorm(x, delta, s2))
    }
    integrate(f, lower = -10 * s1, upper = delta + 10 * s2, subdivisions = 1000)$value
}

# Adjusted Rand Index between two label vectors
adj_rand <- function(a, b) {
    tab <- table(a, b)
    choose2 <- function(m) m * (m - 1) / 2
    sum_ij <- sum(choose2(tab))
    sum_a <- sum(choose2(rowSums(tab)))
    sum_b <- sum(choose2(colSums(tab)))
    expected <- sum_a * sum_b / choose2(sum(tab))
    max_idx <- (sum_a + sum_b) / 2
    if (isTRUE(all.equal(max_idx, expected))) {
        return(1)
    }
    (sum_ij - expected) / (max_idx - expected)
}

# One replicate: simulate data, fit EM, return one row of results
run_rep <- function(d, r, n, p1, init = c("kmeans", "multistart")) {
    init <- match.arg(init)
    set.seed(42 + r) # reproducible dataset r (same z, noise across deltas)
    z <- rbinom(n, 1, 1 - p1) # z = 1 -> second component
    xs <- rnorm(n, mean = d * z, sd = 1)

    f <- tryCatch(
        if (init == "kmeans") {
            em_gmm(xs, init_kmeans(xs, 2), max_iter = 5000)
        } else {
            em_best(xs, 2, max_iter = 5000)
        },
        error = function(e) NULL
    )
    if (is.null(f)) {
        return(data.frame(
            delta = d, rep = r, n_obs = n, p1 = p1, status = "error",
            iterations = NA_integer_, pi1 = NA_real_, mu1 = NA_real_,
            mu2 = NA_real_, s1 = NA_real_, s2 = NA_real_,
            misclass = NA_real_, misclass_perm = NA_real_, ari = NA_real_,
            ll_gap = NA_real_, monotone = NA
        ))
    }

    pred <- max.col(f$gamma) - 1 # 0 = component 1, 1 = component 2
    mis <- mean(pred != z)
    # log-likelihood gap: fitted maximum minus likelihood at the true parameters
    ll_true <- gmm_loglikelihood(xs, c(p1, 1 - p1), c(0, d), c(1, 1))

    data.frame(
        delta = d,
        rep = r,
        n_obs = n,
        p1 = p1,
        status = f$status,
        iterations = f$iterations,
        pi1 = f$pi[1],
        mu1 = f$mu[1],
        mu2 = f$mu[2],
        s1 = sqrt(f$sigma2[1]),
        s2 = sqrt(f$sigma2[2]),
        misclass = mis,
        misclass_perm = min(mis, 1 - mis), # invariant to label swap
        ari = adj_rand(pred, z),
        ll_gap = f$loglik - ll_true,
        monotone = f$monotone
    )
}

# Full simulation over a grid of deltas (and optionally n, p1)
run_sim <- function(deltas, R, n, p1, init = "kmeans") {
    grid <- expand.grid(d = deltas, n = n, p1 = p1, r = seq_len(R))
    out <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
        run_rep(grid$d[i], grid$r[i], grid$n[i], grid$p1[i], init)
    }))
    set.seed(42) # restore the global seed
    out
}

# Summary table: bias, RMSE, convergence, classification + Monte Carlo SEs
summarise_sim <- function(sim) {
    sim %>%
        group_by(n_obs, p1, delta) %>%
        summarise(
            R = n(),
            `bias pi` = mean(pi1 - p1, na.rm = TRUE),
            `RMSE pi` = sqrt(mean((pi1 - p1)^2, na.rm = TRUE)),
            `bias mu1` = mean(mu1 - 0, na.rm = TRUE),
            `RMSE mu1` = sqrt(mean((mu1 - 0)^2, na.rm = TRUE)),
            `bias mu2` = mean(mu2 - delta, na.rm = TRUE),
            `RMSE mu2` = sqrt(mean((mu2 - delta)^2, na.rm = TRUE)),
            `bias sigma1` = mean(s1 - 1, na.rm = TRUE),
            `bias sigma2` = mean(s2 - 1, na.rm = TRUE),
            `MCSE bias pi` = sd(pi1, na.rm = TRUE) / sqrt(sum(!is.na(pi1))),
            `MCSE bias mu1` = sd(mu1, na.rm = TRUE) / sqrt(sum(!is.na(mu1))),
            `MCSE bias mu2` = sd(mu2, na.rm = TRUE) / sqrt(sum(!is.na(mu2))),
            `misclass. %` = 100 * mean(misclass_perm, na.rm = TRUE),
            `Bayes err. %` = 100 * bayes_error(first(delta), first(p1)),
            `mean ARI` = mean(ari, na.rm = TRUE),
            `% label swaps` = 100 * mean(misclass > 0.5, na.rm = TRUE),
            `mean ll gap` = mean(ll_gap, na.rm = TRUE),
            `median iter` = median(iterations, na.rm = TRUE),
            `mean iter` = mean(iterations, na.rm = TRUE),
            `% not conv./degen.` = 100 * mean(status != "converged"),
            `% non-monotone` = 100 * mean(!monotone, na.rm = TRUE),
            .groups = "drop"
        )
}

# ---- Run -----------------------------------------------------------

deltas <- c(0.5, 1.5, 3.0, 5.0)
R2 <- 100
n2 <- 500
true_pi <- 0.5

# init = "kmeans" reproduces the original study; "multistart" uses em_best()
sim2 <- run_sim(deltas, R2, n2, true_pi, init = "kmeans")

tab2 <- summarise_sim(sim2)

# Table 4: parameter recovery and convergence (same layout as before)
kable(
    tab2 %>%
        select(
            `Delta mu` = delta, `bias pi`, `RMSE pi`, `bias mu1`, `RMSE mu1`,
            `bias mu2`, `RMSE mu2`, `bias sigma1`, `bias sigma2`,
            `misclass. %`, `median iter`, `mean iter`, `% not conv./degen.`
        ),
    digits = 3,
    caption = "Table 4. Parameter recovery and convergence of EM by component separation (R = 100 replicates, n = 500)."
)

# Table 5: classification quality against the Bayes-optimal benchmark
kable(
    tab2 %>%
        transmute(
            `Delta mu` = delta,
            `misclass. %`, `Bayes err. %`,
            `excess %` = `misclass. %` - `Bayes err. %`,
            `mean ARI`, `% label swaps`, `mean ll gap`
        ),
    digits = 3,
    caption = "Table 5. Classification accuracy of EM vs the Bayes-optimal error, label swaps and log-likelihood gap to the truth."
)

# Table 6: Monte Carlo standard errors (how precise are the bias estimates?)
kable(
    tab2 %>%
        select(
            `Delta mu` = delta, R,
            `bias pi`, `MCSE bias pi`,
            `bias mu1`, `MCSE bias mu1`,
            `bias mu2`, `MCSE bias mu2`,
            `% non-monotone`
        ),
    digits = 4,
    caption = "Table 6. Monte Carlo standard errors of the bias estimates and monotonicity check."
)

# ---- Plots ---------------------------------------------------------

# Estimation errors (estimate - truth) per parameter and separation
err <- sim2 %>%
    mutate(
        delta_f = factor(delta),
        `pi` = pi1 - p1,
        `mu1` = mu1 - 0,
        `mu2` = mu2 - delta,
        `sigma1` = s1 - 1,
        `sigma2` = s2 - 1
    ) %>%
    select(delta_f, rep, `pi`, mu1, mu2, sigma1, sigma2) %>%
    pivot_longer(-c(delta_f, rep), names_to = "parameter", values_to = "error")

ggplot(err, aes(delta_f, error)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_boxplot(fill = "grey85", outlier.alpha = 0.4) +
    facet_wrap(~parameter, scales = "free_y") +
    labs(
        x = expression(Delta * mu ~ "(true separation)"),
        y = "Estimate - truth",
        title = "EM estimation error by component separation"
    )

# Misclassification vs Bayes-optimal error
mis_df <- tab2 %>%
    select(delta, EM = `misclass. %`, `Bayes-optimal` = `Bayes err. %`) %>%
    pivot_longer(-delta, names_to = "classifier", values_to = "error")

ggplot(mis_df, aes(delta, error, colour = classifier)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
        x = expression(Delta * mu),
        y = "Misclassification rate (%)",
        colour = NULL,
        title = "EM classification error vs Bayes-optimal error"
    )

# Iterations needed to converge (slower when components overlap)
ggplot(sim2, aes(factor(delta), iterations)) +
    geom_boxplot(fill = "grey85", outlier.alpha = 0.4) +
    scale_y_log10() +
    labs(
        x = expression(Delta * mu),
        y = "EM iterations (log scale)",
        title = "Convergence speed by component separation"
    )

# Outcome of each run (converged / degenerate / max_iter / error)
ggplot(sim2, aes(factor(delta), fill = status)) +
    geom_bar(position = "fill") +
    scale_fill_brewer(palette = "Set2") +
    labs(
        x = expression(Delta * mu),
        y = "Share of replicates",
        fill = "Status",
        title = "Convergence status by component separation"
    )

# ---- Optional extensions (uncomment to run) -------------------------
# Effect of sample size at fixed separation:
# sim_n <- run_sim(deltas = c(1.5, 3), R = 100, n = c(50, 200, 500, 2000), p1 = 0.5)
# kable(summarise_sim(sim_n), digits = 3)
#
# Unbalanced mixing proportions:
# sim_p <- run_sim(deltas = c(1.5, 3), R = 100, n = 500, p1 = c(0.1, 0.3, 0.5))
# kable(summarise_sim(sim_p), digits = 3)
#
# Single k-means start vs multi-start:
# sim_ms <- run_sim(deltas, R2, n2, true_pi, init = "multistart")
# kable(summarise_sim(sim_ms), digits = 3)
