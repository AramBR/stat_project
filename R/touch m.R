# log_dens_matrix(): n x K matrix with entries log(pi_k) + log phi(x_i; mu_k, sigma_k^2)
# 1-dim case
log_dens_matrix <- function(x, pi, mu, sigma2) {
  K <- length(pi)
  sapply(seq_len(K), function(k)
    log(pi[k]) + dnorm(x, mean = mu[k], sd = sqrt(sigma2[k]), log = TRUE))
}
