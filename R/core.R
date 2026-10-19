# Joint log-density matrix
# Computes the n x K matrix of log-weights + log-densities:
# L[i, k] = log(pi_k) + log(phi(x_i; mu_k, sigma_k^2))
log_dens_matrix <- function(x, pi, mu, sigma2) {
  K <- length(pi)
  sapply(seq_len(K), function(k) {
    log(pi[k]) + stats::dnorm(x, mean = mu[k], sd = sqrt(sigma2[k]), log = TRUE)
  })
}

# Row-wise Log-Sum-Exp
# Computes log(sum(exp(L[i, ]))) stably using the max-trick.
log_sum_exp <- function(L) {
  L <- as.matrix(L)
  # Row maximum for scaling
  m <- apply(L, 1, max)
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
  L   <- log_dens_matrix(x, pi, mu, sigma2)
  lse <- log_sum_exp(L)
  
  # Normalize in log-space, then exponentiate
  gamma <- exp(sweep(L, 1, lse, "-"))
  
  list(
    gamma  = gamma,
    loglik = sum(lse)
  )
}

# M-Step (Maximization)
# Updates mixture weights, means, and variances.
# Includes safeguards for cluster collapse.
m_step <- function(x, gamma, var_floor) {
  n  <- length(x)
  nk <- colSums(gamma)
  
  pi <- nk / n
  mu <- colSums(gamma * x) / nk
  
  sigma2 <- sapply(seq_along(mu), function(k) {
    sum(gamma[, k] * (x - mu[k])^2) / nk[k]
  })
  
  hit_floor <- any(sigma2 < var_floor)
  sigma2    <- pmax(sigma2, var_floor) # Clamp variance
  
  list(
    pi         = pi,
    mu         = mu,
    sigma2     = sigma2,
    degenerate = hit_floor || any(nk < 1)
  )
}

# EM Algorithm for 1D Gaussian Mixture Models
# x Vector of observations
# init List containing starting values: pi, mu, sigma2
# tol Relative convergence tolerance
# max_iter Maximum number of EM iterations
# var_floor Minimum allowable variance (prevents singularity)
em_gmm <- function(x, init, tol = 1e-8, max_iter = 2000,
                   var_floor = (0.01 * stats::sd(x))^2) {
  pi     <- init$pi
  mu     <- init$mu
  sigma2 <- init$sigma2
  
  # Pre-allocate to avoid slow vector growth in loops
  loglik_hist <- numeric(max_iter)
  status      <- "max_iter reached"
  actual_iter <- 0
  
  for (iter in seq_len(max_iter)) {
    actual_iter <- iter
    
    # 1. E-Step
    E <- e_step(x, pi, mu, sigma2)
    loglik_hist[iter] <- E$loglik
    
    # 2. Check convergence
    if (iter > 1) {
      change <- loglik_hist[iter] - loglik_hist[iter - 1]
      if (abs(change) < tol * (1 + abs(loglik_hist[iter - 1]))) {
        status <- "converged"
        break
      }
    }
    
    # 3. M-Step
    M <- m_step(x, E$gamma, var_floor)
    pi     <- M$pi
    mu     <- M$mu
    sigma2 <- M$sigma2
    
    if (M$degenerate) {
      status <- "degenerate"
      break
    }
  }
  
  # Trim history to completed iterations
  loglik_hist <- loglik_hist[seq_len(actual_iter)]
  
  # Sort components by mean (identifiability / label-switching resolution)
  o <- order(mu)
  pi     <- pi[o]
  mu     <- mu[o]
  sigma2 <- sigma2[o]
  
  # Compute final responsibilities with sorted parameters
  final_gamma <- e_step(x, pi, mu, sigma2)$gamma
  
  list(
    pi          = pi,
    mu          = mu,
    sigma2      = sigma2,
    K           = length(mu),
    loglik      = utils::tail(loglik_hist, 1),
    loglik_hist = loglik_hist,
    iterations  = actual_iter,
    status      = status,
    gamma       = final_gamma,
    monotone    = all(diff(loglik_hist) > -1e-8)
  )
}
