# Joint log-density matrix
# Computes the n x K matrix of log-weights + log-densities:
# L[i, k] = log(pi_k) + log(phi(x_i; mu_k, sigma_k^2))
# Always returns an n x K matrix (also for n = 1 or K = 1).
log_dens_matrix <- function(x, pi, mu, sigma2) {
    n <- length(x)
    K <- length(pi)
    L <- matrix(NA_real_, nrow = n, ncol = K)
    for (k in seq_len(K)) {
        L[, k] <- log(pi[k]) +
            stats::dnorm(x, mean = mu[k], sd = sqrt(sigma2[k]), log = TRUE)
    }
    L
}

# Row-wise Log-Sum-Exp
# Computes log(sum(exp(L[i, ]))) stably using the max-trick.
# Rows that are entirely -Inf give -Inf (instead of NaN).
log_sum_exp <- function(L) {
    L <- as.matrix(L)
    # Row maximum for scaling
    m <- apply(L, 1, max)
    # A non-finite maximum (all -Inf, or +Inf) would give Inf - Inf = NaN
    m[!is.finite(m)] <- 0
    # sweep ensures explicit row-wise subtraction across matrix columns
    m + log(rowSums(exp(sweep(L, 1, m, "-"))))
}

# Observed-data log-likelihood
gmm_loglikelihood <- function(x, pi, mu, sigma2) {
    sum(log_sum_exp(log_dens_matrix(x, pi, mu, sigma2)))
}

# E-Step (Expectation)
# Computes posterior responsibilities gamma_ik and the log-likelihood
# at current parameter values.
e_step <- function(x, pi, mu, sigma2) {
    L <- log_dens_matrix(x, pi, mu, sigma2)
    lse <- log_sum_exp(L)
    # Normalize in log-space, then exponentiate
    gamma <- exp(sweep(L, 1, lse, "-"))
    list(
        gamma = gamma,
        loglik = sum(lse)
    )
}

# M-Step (Maximization)
# Updates mixture weights, means, and variances.
# Flags cluster collapse: a component with (almost) no weight, a non-finite
# estimate, or a variance below the floor.
m_step <- function(x, gamma, var_floor) {
    n <- length(x)
    nk <- colSums(gamma)
    pi <- nk / n
    mu <- colSums(gamma * x) / nk
    sigma2 <- vapply(seq_along(mu), function(k) {
        sum(gamma[, k] * (x - mu[k])^2) / nk[k]
    }, numeric(1))
    bad <- !is.finite(mu) | !is.finite(sigma2)
    hit_floor <- any(sigma2[!bad] < var_floor)
    sigma2 <- pmax(sigma2, var_floor) # Clamp variance
    list(
        pi = pi,
        mu = mu,
        sigma2 = sigma2,
        degenerate = any(bad) || hit_floor || any(nk < 1)
    )
}

# EM Algorithm for 1D Gaussian Mixture Models
# x Vector of observations
# init List containing starting values: pi, mu, sigma2
# tol Relative convergence tolerance
# max_iter Maximum number of EM iterations (M-steps)
# var_floor Minimum allowable variance (prevents singularity)
#
# Conventions of the returned object:
#   iterations  number of parameter updates (M-steps) actually accepted
#   loglik_hist log-likelihood at the start values (first entry) and after
#               every accepted M-step, so length = iterations + 1
#   loglik      log-likelihood at exactly the returned parameters
#   status      "converged", "max_iter reached" or "degenerate"; a degenerate
#               M-step is NOT accepted, so the returned parameters are the
#               last valid ones
em_gmm <- function(
    x,
    init,
    tol = 1e-8,
    max_iter = 2000,
    var_floor = (0.01 * stats::sd(x))^2
) {
    pi <- init$pi
    mu <- init$mu
    sigma2 <- init$sigma2
    # Pre-allocate to avoid slow vector growth in loops
    loglik_hist <- numeric(max_iter + 1)
    status <- "max_iter reached"
    n_updates <- 0

    # E-step at the starting values
    E <- e_step(x, pi, mu, sigma2)
    loglik_hist[1] <- E$loglik

    for (iter in seq_len(max_iter)) {
        # 1. M-Step (rejected if it collapses a component)
        M <- m_step(x, E$gamma, var_floor)
        if (M$degenerate) {
            status <- "degenerate"
            break
        }
        pi <- M$pi
        mu <- M$mu
        sigma2 <- M$sigma2
        n_updates <- iter

        # 2. E-Step at the new parameters
        E <- e_step(x, pi, mu, sigma2)
        loglik_hist[iter + 1] <- E$loglik

        # 3. Check convergence (relative change of the log-likelihood)
        change <- loglik_hist[iter + 1] - loglik_hist[iter]
        if (abs(change) < tol * (1 + abs(loglik_hist[iter]))) {
            status <- "converged"
            break
        }
    }
    # Trim history to completed iterations
    loglik_hist <- loglik_hist[seq_len(n_updates + 1)]
    # Sort components by mean (identifiability / label-switching resolution)
    o <- order(mu)
    prev <- utils::head(loglik_hist, -1)
    list(
        pi = pi[o],
        mu = mu[o],
        sigma2 = sigma2[o],
        K = length(mu),
        # E was computed at exactly these parameters, so loglik and gamma match them
        loglik = loglik_hist[n_updates + 1],
        loglik_hist = loglik_hist,
        iterations = n_updates,
        status = status,
        gamma = E$gamma[, o, drop = FALSE],
        # relative tolerance: allows for floating-point noise only
        monotone = all(diff(loglik_hist) > -1e-10 * (1 + abs(prev)))
    )
}
