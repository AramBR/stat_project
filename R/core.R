# log_dens_matrix(): n x K matrix with entries log(pi_k) + log phi(x_i; mu_k, sigma_k^2)
# 1-dim case
log_dens_matrix <- function(x, pi, mu, sigma2) {
  K <- length(pi)
  sapply(seq_len(K), function(k)
    log(pi[k]) + dnorm(x, mean = mu[k], sd = sqrt(sigma2[k]), log = TRUE))
}

# log_sum_exp(): row-wise log(sum(exp(.))) computed stably.
# Subtracting the row maximum before exponentiating prevents underflow.
log_sum_exp <- function(L) {
  L <- as.matrix(L)
  m <- apply(L, 1, max)
  m + log(rowSums(exp(L - m)))
}

# EM-LOOP

# 1) Observed-data log-likelihood  l(theta) = sum_i log sum_k pi_k phi(x_i)
gmm_loglikelihood <- function(x, pi, mu, sigma2) {
  sum(log_sum_exp(log_dens_matrix(x, pi, mu, sigma2)))
}

# 2) E-step: posterior responsibilities gamma_ik.
#    Returns the n x K matrix gamma and the log-likelihood at the current params
e_step <- function(x, pi, mu, sigma2) {
  L   <- log_dens_matrix(x, pi, mu, sigma2)   
  lse <- log_sum_exp(L)                      
  list(gamma  = exp(L - lse),               
       loglik = sum(lse))
}

# 3) M-step:
#    var_floor : lower bound for every variance(for infinities)
#    Returns new parameters and a flag if infinities are reached
m_step <- function(x, gamma, var_floor) {
  n  <- length(x)
  nk <- colSums(gamma)                        
  pi <- nk / n
  mu <- colSums(gamma * x) / nk
  sigma2 <- sapply(seq_along(mu), function(k)
    sum(gamma[, k] * (x - mu[k])^2) / nk[k])
  hit_floor <- any(sigma2 < var_floor)
  sigma2 <- pmax(sigma2, var_floor)            # variance safeguard
  list(pi = pi, mu = mu, sigma2 = sigma2,
       degenerate = hit_floor || any(nk < 1))
}

# 4) Full EM loop
#    init      : list(pi, mu, sigma2) of starting values
#    tol       : stop when |l(t) - l(t-1)| < tol * (1 + |l(t-1)|)  (relative)
#    max_iter  : hard limit on the number of iterations
#    Returns final estimates (sorted so that mu_1 < mu_2 < ...), the
#    log-likelihood history, number of iterations and convergence status.
# ---------------------------------------------------------------------------
em_gmm <- function(x, init, tol = 1e-8, max_iter = 2000,
                   var_floor = (0.01 * sd(x))^2) {
  pi <- init$pi; mu <- init$mu; sigma2 <- init$sigma2
  loglik_hist <- numeric(0)
  status <- "max_iter reached"
  for (iter in 1:max_iter) {
    E <- e_step(x, pi, mu, sigma2)
    loglik_hist <- c(loglik_hist, E$loglik)
    if (iter > 1) {
      change <- loglik_hist[iter] - loglik_hist[iter - 1]
      if (abs(change) < tol * (1 + abs(loglik_hist[iter - 1]))) {
        status <- "converged"; break
      }
    }
    M <- m_step(x, E$gamma, var_floor)
    pi <- M$pi; mu <- M$mu; sigma2 <- M$sigma2
    if (M$degenerate) {
      status <- "degenerate"
      loglik_hist <- c(loglik_hist, gmm_loglikelihood(x, pi, mu, sigma2))
      break
    }
  }
  o <- order(mu)                   
  list(pi = pi[o], mu = mu[o], sigma2 = sigma2[o], K = length(mu),
       loglik = tail(loglik_hist, 1), loglik_hist = loglik_hist,
       iterations = length(loglik_hist) - 1, status = status,
       gamma = e_step(x, pi[o], mu[o], sigma2[o])$gamma,
       monotone = all(diff(loglik_hist) > -1e-8))
}

# Starting-value helpers
# moments of a hard partition (cluster labels cl = 1..K)
init_from_labels <- function(x, cl, K) {
  list(pi     = tabulate(cl, K) / length(x),
       mu     = sapply(1:K, function(k) mean(x[cl == k])),
       sigma2 = sapply(1:K, function(k) if (sum(cl == k) > 1) var(x[cl == k]) else var(x)))
}
init_kmeans <- function(x, K) init_from_labels(x, kmeans(x, K, nstart = 10)$cluster, K)
init_random <- function(x, K) list(pi = rep(1 / K, K), mu = sample(x, K),
                                   sigma2 = rep(var(x), K))

# best of several starts: k-means + n_random random starts, keep the highest
# log-likelihood among non-degenerate runs
em_best <- function(x, K, n_random = 10, ...) {
  if (K == 1) {                                 # closed-form MLE, no EM needed
    s2 <- mean((x - mean(x))^2)
    return(list(pi = 1, mu = mean(x), sigma2 = s2, K = 1, iterations = 0,
                loglik = sum(dnorm(x, mean(x), sqrt(s2), log = TRUE)),
                status = "closed form"))
  }
  fits <- c(list(em_gmm(x, init_kmeans(x, K), ...)),
            lapply(1:n_random, function(r) em_gmm(x, init_random(x, K), ...)))
  ok <- sapply(fits, function(f) f$status != "degenerate")
  fits[ok][[which.max(sapply(fits[ok], `[[`, "loglik"))]]
}

# mixture density on a grid (for plots)
gmm_density <- function(grid, fit)
  rowSums(sapply(seq_len(fit$K), function(k)
    fit$pi[k] * dnorm(grid, fit$mu[k], sqrt(fit$sigma2[k]))))
