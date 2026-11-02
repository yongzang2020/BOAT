# ============================================================
# CD-SEP: conventional umbrella design. Each subgroup is analyzed separately with beta posteriors of the
# response probabilities (no borrowing across subgroups), and no PDO data are used. The trial flow
# is simulate_trial_noPDO() in trial_noPDO.R.
# ============================================================


# ============================================================
# Efficacy posterior (by subgroup k, independent across k)
#   p_{k,j} | data ~ Beta(a_{k,j}, b_{k,j})
#   post_prob_eff(k,j) = P(p_{k,j} > p_{k,0} + delta_E | data)
#     = ∫ f_{k,0}(t) * [1 - F_{k,j}(t + delta_E)] dt
#   df_pat may be empty
# ============================================================
eff_posterior_beta_by_k <- function(df_pat, K, J, Sk_list,
                                    prior_pE_a = 1, prior_pE_b = 1,
                                    delta_E = 0.0) {
  stopifnot(K >= 1, J >= 1)
  stopifnot(is.data.frame(df_pat))
  if (nrow(df_pat) > 0) stopifnot(all(c("k", "j", "YE") %in% names(df_pat)))
  
  a_post <- matrix(NA_real_, nrow = K, ncol = J + 1)
  b_post <- matrix(NA_real_, nrow = K, ncol = J + 1)
  
  for (k in 1:K) {
    for (j in Sk_list[[k]]) {
      if (nrow(df_pat) > 0) {
        idx <- which(df_pat$k == k & df_pat$j == j)
        y <- if (length(idx) > 0) sum(df_pat$YE[idx]) else 0
        n <- length(idx)
      } else {
        y <- 0
        n <- 0
      }
      a_post[k, j + 1] <- prior_pE_a + y
      b_post[k, j + 1] <- prior_pE_b + n - y
    }
  }
  
  post_pE_mean <- a_post / (a_post + b_post)
  
  post_prob_eff <- matrix(NA_real_, nrow = K, ncol = J + 1,
                          dimnames = list(paste0("k", 1:K), paste0("T", 0:J)))
  
  for (k in 1:K) {
    if (!(0 %in% Sk_list[[k]])) next
    
    a0 <- a_post[k, 1]
    b0 <- b_post[k, 1]
    if (!is.finite(a0) || !is.finite(b0)) next
    
    for (j in setdiff(Sk_list[[k]], 0)) {
      aj <- a_post[k, j + 1]
      bj <- b_post[k, j + 1]
      if (!is.finite(aj) || !is.finite(bj)) next
      
      if (delta_E >= 1) {
        post_prob_eff[k, j + 1] <- 0.0
        next
      }
      
      integrand <- function(t) {
        x <- t + delta_E
        Fj <- ifelse(x <= 0, 0,
                     ifelse(x >= 1, 1, pbeta(x, aj, bj)))
        dbeta(t, a0, b0) * (1 - Fj)
      }
      
      post_prob_eff[k, j + 1] <- tryCatch(
        integrate(integrand, lower = 0, upper = 1,
                  rel.tol = 1e-6, abs.tol = 0)$value,
        error = function(e) NA_real_
      )
    }
  }
  
  list(
    post_pE_mean  = post_pE_mean,
    post_prob_eff = post_prob_eff,
    a_post = a_post,
    b_post = b_post
  )
}


# ============================================================
# simulate_trial_CD_SEP(): one simulated CD-SEP trial
# ============================================================
simulate_trial_CD_SEP <- function(
    pE_mat,
    piT_vec,
    max_n_k,
    interim_fracs = c(0.5, 0.75),
    delta_E = 0.0,
    delta_T = 0.0,
    lambda = 0.95,
    prior_pE_a = 1,
    prior_pE_b = 1,
    prior_piT_a = 1,
    prior_piT_b = 1,
    seed = 123,
    verbose = FALSE
) {
  eff_posterior <- function(df_pat, K, J, Sk_list, delta_E, analysis_id) {
    eff_posterior_beta_by_k(df_pat = df_pat, K = K, J = J, Sk_list = Sk_list,
                            prior_pE_a = prior_pE_a, prior_pE_b = prior_pE_b, delta_E = delta_E)
  }
  simulate_trial_noPDO(pE_mat, piT_vec, max_n_k, eff_posterior = eff_posterior,
                       interim_fracs = interim_fracs, delta_E = delta_E, delta_T = delta_T,
                       lambda = lambda, prior_piT_a = prior_piT_a, prior_piT_b = prior_piT_b,
                       seed = seed, verbose = verbose)
}
