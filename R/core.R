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
