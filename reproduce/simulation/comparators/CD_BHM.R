# ============================================================
# CD-BHM: equal randomization, no PDO data, non-adaptive borrowing across subgroups by a
# Bayesian hierarchical model. The trial flow is simulate_trial_noPDO() in trial_noPDO.R.
#
# Efficacy model, separately for each treatment A_j (control included), over the
# subgroups in which A_j is evaluated:
#   Y_R,i | theta_jk ~ Bernoulli(pi_R,jk),  pi_R,jk = Phi(theta_jk)        (probit)
#   theta_jk | mu_j, tau_j ~ N(mu_j, tau_j^2)                             (exchangeable)
#   mu_j ~ N(0, 10^2),  tau_j^2 ~ IG(0.01, 0.01)   (same IG prior as the variances of BOAT)
# Posterior by Gibbs sampling with Albert-Chib latent variables. Every full conditional
# is standard (tau_j^2 is drawn from its inverse-gamma full conditional, as in BOAT).
# ============================================================


# ============================================================
# BHM sampler
# - returns post_pE_mean, post_prob_eff, post_tau_mean
# - post_prob_eff[k, j+1] = P(pi_R,jk > pi_R,0k + delta_E | data)
# ============================================================
mcmc_bhm_efficacy <- function(
    df_pat,
    K = 3,
    J = 3,
    Sk_list,
    n_iter = 10000,
    n_burn = 5000,
    n_thin = 1,

    delta_E = 0.0,

    prior_mu_mean = 0, prior_mu_sd = 1e1,
    prior_tau2_a = 0.01, prior_tau2_b = 0.01,

    seed = 1e3
){
  set.seed(seed)

  N <- nrow(df_pat)
  stopifnot(N >= 1)
  k_pat   <- df_pat$k
  col_pat <- df_pat$j + 1L
  lower <- ifelse(df_pat$YE == 1, 0, -Inf)
  upper <- ifelse(df_pat$YE == 1, Inf, 0)

  # subgroups in which each treatment is evaluated
  Kj_list <- vector("list", J + 1)
  for (trt in 0:J) Kj_list[[trt + 1]] <- which(sapply(Sk_list, function(x) trt %in% x))

  pat_idx <- vector("list", K)
  for (k in 1:K) {
    pat_idx[[k]] <- vector("list", J + 1)
    for (col in 1:(J + 1)) pat_idx[[k]][[col]] <- which(k_pat == k & col_pat == col)
  }

  # initial values: smoothed empirical probit means
  theta <- matrix(NA_real_, nrow = K, ncol = J + 1)
  for (col in 1:(J + 1)) {
    for (k in Kj_list[[col]]) {
      idx <- pat_idx[[k]][[col]]
      p_tilde <- (sum(df_pat$YE[idx]) + 0.5) / (length(idx) + 1)
      theta[k, col] <- qnorm(p_tilde)
    }
  }
  mu  <- sapply(1:(J + 1), function(col) if (length(Kj_list[[col]]) > 0) mean(theta[Kj_list[[col]], col]) else NA_real_)
  tau <- rep(1, J + 1)

  keep_seq <- seq(n_burn + 1, n_iter, n_thin)
  n_keep <- length(keep_seq)
  pE_sum  <- matrix(0, nrow = K, ncol = J + 1)
  eff_sum <- matrix(0, nrow = K, ncol = J + 1)
  tau_sum <- numeric(J + 1)

  for (iter in 1:n_iter) {

    # 1) latent variables
    Z <- rtnorm(N, theta[cbind(k_pat, col_pat)], 1, lower = lower, upper = upper)

    for (col in 1:(J + 1)) {
      subs <- Kj_list[[col]]
      if (length(subs) == 0) next

      # 2) subgroup-specific probit means
      for (k in subs) {
        idx <- pat_idx[[k]][[col]]
        v <- 1 / (length(idx) + 1 / tau[col]^2)
        m <- v * (sum(Z[idx]) + mu[col] / tau[col]^2)
        theta[k, col] <- rnorm(1, m, sqrt(v))
      }
      th <- theta[subs, col]

      # 3) common mean
      v <- 1 / (length(subs) / tau[col]^2 + 1 / prior_mu_sd^2)
      m <- v * (sum(th) / tau[col]^2 + prior_mu_mean / prior_mu_sd^2)
      mu[col] <- rnorm(1, m, sqrt(v))

      # 4) between-subgroup variance: tau^2 | . ~ IG(a + m/2, b + sum((theta - mu)^2)/2)
      tau[col] <- sqrt(1 / rgamma(1, prior_tau2_a + length(th) / 2,
                                  prior_tau2_b + sum((th - mu[col])^2) / 2))
    }

    if (iter >= n_burn + 1 && ((iter - (n_burn + 1)) %% n_thin == 0)) {
      pE <- pnorm(theta)
      pE_sum  <- pE_sum + ifelse(is.na(pE), 0, pE)
      tau_sum <- tau_sum + tau
      for (k in 1:K) {
        p0 <- pE[k, 1]
        if (is.na(p0)) next
        for (col in 1:(J + 1)) {
          if (is.na(pE[k, col])) next
          eff_sum[k, col] <- eff_sum[k, col] + as.integer(pE[k, col] > p0 + delta_E)
        }
      }
    }
  }

  post_pE_mean <- pE_sum / n_keep
  post_pE_mean[is.na(theta)] <- NA
  post_prob_eff <- eff_sum / n_keep
  post_prob_eff[is.na(theta)] <- NA
  dimnames(post_prob_eff) <- list(paste0("k", 1:K), paste0("T", 0:J))

  list(
    post_pE_mean  = post_pE_mean,
    post_prob_eff = post_prob_eff,
    post_tau_mean = tau_sum / n_keep
  )
}


# ============================================================
# simulate_trial_CD_BHM(): one simulated CD-BHM trial
# Each analysis uses seed + 1000 * analysis_id for its MCMC, as in BOAT.
# ============================================================
simulate_trial_CD_BHM <- function(
    pE_mat,
    piT_vec,
    max_n_k,
    interim_fracs = c(0.5, 0.75),
    delta_E = 0.0,
    delta_T = 0.0,
    lambda = 0.95,
    mcmc_iter = 10000,
    mcmc_burn = 5000,
    mcmc_thin = 1,
    prior_tau2_a = 0.01,
    prior_tau2_b = 0.01,
    prior_piT_a = 1,
    prior_piT_b = 1,
    seed = 123,
    verbose = FALSE
) {
  eff_posterior <- function(df_pat, K, J, Sk_list, delta_E, analysis_id) {
    mcmc_bhm_efficacy(df_pat, K = K, J = J, Sk_list = Sk_list,
                      n_iter = mcmc_iter, n_burn = mcmc_burn, n_thin = mcmc_thin,
                      delta_E = delta_E, prior_tau2_a = prior_tau2_a, prior_tau2_b = prior_tau2_b,
                      seed = seed + 1000L * analysis_id)
  }
  simulate_trial_noPDO(pE_mat, piT_vec, max_n_k, eff_posterior = eff_posterior,
                       interim_fracs = interim_fracs, delta_E = delta_E, delta_T = delta_T,
                       lambda = lambda, prior_piT_a = prior_piT_a, prior_piT_b = prior_piT_b,
                       seed = seed, verbose = verbose)
}
