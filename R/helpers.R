# Initialize parameters from hard cluster labels
init_from_labels <- function(x, cl, K) {
  n <- length(x)
  
  # Fallback for empty clusters to prevent NaNs
  counts <- tabulate(cl, K)
  pi     <- pmax(counts, 1e-6) / n
  pi     <- pi / sum(pi)
  
  mu <- sapply(seq_len(K), function(k) {
    if (counts[k] == 0) sample(x, 1) else mean(x[cl == k])
  })
  
  sigma2 <- sapply(seq_len(K), function(k) {
    if (counts[k] > 1) stats::var(x[cl == k]) else stats::var(x)
  })
  
  list(pi = pi, mu = mu, sigma2 = sigma2)
}

# Initialize parameters using K-Means
init_kmeans <- function(x, K) {
  km <- stats::kmeans(x, centers = K, nstart = 10)
  init_from_labels(x, km$cluster, K)
}

# Initialize parameters randomly
init_random <- function(x, K) {
  list(
    pi     = rep(1 / K, K),
    mu     = sample(x, K, replace = FALSE),
    sigma2 = rep(stats::var(x), K)
  )
}

# Fit GMM using multiple initializations (best log-likelihood selection)
# 
#  x Observations vector
#  K Number of components
#  n_random Number of random initializations in addition to k-means
#  ... Forwarded arguments to em_gmm (tol, max_iter, var_floor)
em_best <- function(x, K, n_random = 10, ...) {
  n <- length(x)
  
  # Closed-form MLE for single component (no EM needed)
  if (K == 1) {
    s2 <- mean((x - mean(x))^2) # MLE biased variance
    ll <- sum(stats::dnorm(x, mean = mean(x), sd = sqrt(s2), log = TRUE))
    
    return(list(
      pi          = 1,
      mu          = mean(x),
      sigma2      = s2,
      K           = 1,
      iterations  = 0,
      loglik      = ll,
      loglik_hist = ll,
      status      = "closed form",
      gamma       = matrix(1, nrow = n, ncol = 1),
      monotone    = TRUE
    ))
  }
  
  # Generate runs (1 K-means run + n_random random starts)
  fits <- c(
    list(em_gmm(x, init_kmeans(x, K), ...)),
    lapply(seq_len(n_random), function(r) em_gmm(x, init_random(x, K), ...))
  )
  
  # Filter out degenerated runs
  ok <- sapply(fits, function(f) f$status != "degenerate")
  valid_fits <- fits[ok]
  
  # Fallback if all initializations degenerate
  if (length(valid_fits) == 0) {
    warning("All initializations degenerated. Returning highest loglik degenerate fit.")
    valid_fits <- fits
  }
  
  best_idx <- which.max(sapply(valid_fits, `[[`, "loglik"))
  valid_fits[[best_idx]]
}

# Evaluate continuous GMM density over a grid (for plotting)
#  grid Numeric vector of evaluation points
#  fit Model output from em_gmm or em_best
gmm_density <- function(grid, fit) {
  rowSums(sapply(seq_len(fit$K), function(k) {
    fit$pi[k] * stats::dnorm(grid, mean = fit$mu[k], sd = sqrt(fit$sigma2[k]))
  }))
}
