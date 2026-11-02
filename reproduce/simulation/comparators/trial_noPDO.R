# ============================================================
# Trial flow shared by the comparators CD-SEP, CD-BHM and CD-AB, which use no PDO data: equal
# randomization within subgroup, the interim analyses, drop rule and final rule of Section 3,
# and a beta toxicity model pooled across subgroups. The designs differ only in the efficacy
# analysis (CD_SEP.R, CD_BHM.R, CD_AB.R).
# A design supplies eff_posterior(df_pat, K, J, Sk_list, delta_E, analysis_id), which must
# return a list with post_prob_eff, a K x (J+1) matrix of P(pi_R,jk > pi_R,0k + delta_E | data).
# analysis_id is 1, ..., L for the interim analyses and L + 1 for the final analysis.
# ============================================================
library(msm)   # rtnorm for truncated normal draws


# ============================================================
# Toxicity posterior (pooled across subgroup k)
#   pi_j | data ~ Beta(a_j, b_j)
#   post_prob_safe(j): P(pi_j < pi_0 + delta_T | data)
#     = ∫ f0(t) * F_j(t + delta_T) dt
# ============================================================
tox_posterior_beta <- function(df_pat, J,
                               prior_piT_a = 1, prior_piT_b = 1,
                               delta_T = 0.0) {
  stopifnot(J >= 1)
  stopifnot(is.data.frame(df_pat))
  if (nrow(df_pat) > 0) stopifnot(all(c("j", "YT") %in% names(df_pat)))

  a_post <- numeric(J + 1)
  b_post <- numeric(J + 1)

  for (j in 0:J) {
    if (nrow(df_pat) > 0) {
      idx <- which(df_pat$j == j)
      y <- if (length(idx) > 0) sum(df_pat$YT[idx]) else 0
      n <- length(idx)
    } else {
      y <- 0
      n <- 0
    }
    a_post[j + 1] <- prior_piT_a + y
    b_post[j + 1] <- prior_piT_b + n - y
  }

  post_piT_mean <- a_post / (a_post + b_post)

  post_prob_safe <- rep(NA_real_, J + 1)
  a0 <- a_post[1]
  b0 <- b_post[1]

  for (j in 0:J) {
    if (j == 0) {
      post_prob_safe[j + 1] <- NA_real_
      next
    }
    if (delta_T >= 1) {
      post_prob_safe[j + 1] <- 1.0
      next
    }

    aj <- a_post[j + 1]
    bj <- b_post[j + 1]

    integrand <- function(t) {
      x <- t + delta_T
      Fj <- ifelse(x <= 0, 0,
                   ifelse(x >= 1, 1, pbeta(x, aj, bj)))
      dbeta(t, a0, b0) * Fj
    }

    post_prob_safe[j + 1] <- tryCatch(
      integrate(integrand, lower = 0, upper = 1,
                rel.tol = 1e-6, abs.tol = 0)$value,
      error = function(e) NA_real_
    )
  }

  list(
    post_piT_mean  = post_piT_mean,
    post_prob_safe = post_prob_safe,
    a_post = a_post,
    b_post = b_post
  )
}


# ============================================================
# simulate_trial_noPDO(): one simulated trial of a design without PDO data
#   - maximum sample size N_k of each subgroup: max_n_k
#   - interim analyses when every active subgroup reaches its interim sample size
#   - accrual balances n_k / N_k across the active subgroups
#   - equal randomization among the open treatments of the subgroup
#   - treatments are dropped for futility or toxicity only
# ============================================================
simulate_trial_noPDO <- function(
    pE_mat,
    piT_vec,
    max_n_k,
    eff_posterior,

    interim_fracs = c(0.5, 0.75),

    delta_E = 0.0,
    delta_T = 0.0,
    lambda = 0.95,

    prior_piT_a = 1,
    prior_piT_b = 1,

    seed = 123,
    verbose = FALSE
) {
  set.seed(seed)

  K <- nrow(pE_mat)
  J <- ncol(pE_mat) - 1
  stopifnot(length(piT_vec) == (J + 1))

  if (length(max_n_k) == 1) {
    max_n_k <- rep(as.integer(max_n_k), K)
  } else {
    stopifnot(length(max_n_k) == K)
    max_n_k <- as.integer(max_n_k)
  }
  stopifnot(all(is.finite(max_n_k)) && all(max_n_k >= 1))

  Sk_list <- vector("list", K)
  for (k in 1:K) Sk_list[[k]] <- which(!is.na(pE_mat[k, ])) - 1

  P_L_k <- function(nk, Nk) lambda * sqrt(nk / Nk)

  df_pat <- data.frame(
    id = integer(0), k = integer(0), j = integer(0),
    treat = character(0),
    YT = integer(0), YE = integer(0)
  )

  active_list <- lapply(Sk_list, function(x) x)
  names(active_list) <- as.character(1:K)
  active_groups <- as.integer(1:K)
  stopped_groups <- integer(0)

  n_by_k <- rep(0L, K)
  enrolled_count <- 0L

  interim_fracs <- sort(unique(interim_fracs))
  if (length(interim_fracs) > 0) {
    interim_targets <- sapply(interim_fracs, function(fr) pmin(as.integer(ceiling(fr * max_n_k)), max_n_k))
    if (is.vector(interim_targets)) interim_targets <- matrix(interim_targets, ncol = 1)
    colnames(interim_targets) <- paste0("I", seq_len(ncol(interim_targets)))
    rownames(interim_targets) <- paste0("k", 1:K)
  } else {
    interim_targets <- matrix(0L, nrow = K, ncol = 0)
  }
  next_interim_id <- 1L

  while (TRUE) {

    repeat {
      if (next_interim_id <= length(interim_fracs) && length(active_groups) > 0) {
        targets_k <- interim_targets[, next_interim_id]
        cond_vec <- sapply(active_groups, function(k) n_by_k[k] >= targets_k[k])

        if (all(cond_vec)) {

          if (nrow(df_pat) > 0) {

            postE <- eff_posterior(df_pat = df_pat, K = K, J = J, Sk_list = Sk_list,
                                   delta_E = delta_E, analysis_id = next_interim_id)

            tox <- tox_posterior_beta(
              df_pat = df_pat,
              J = J,
              prior_piT_a = prior_piT_a,
              prior_piT_b = prior_piT_b,
              delta_T = delta_T
            )

            for (k in as.integer(active_groups)) {
              k_chr <- as.character(k)
              drop_list <- c()

              nk_now <- n_by_k[k]
              PLk <- P_L_k(nk_now, max_n_k[k])

              for (j in active_list[[k_chr]]) {
                if (j == 0) next

                post_eff  <- postE$post_prob_eff[k, j + 1]
                post_safe <- tox$post_prob_safe[j + 1]
                if (is.na(post_eff) || is.na(post_safe)) next

                if (post_eff < PLk || post_safe < PLk) {
                  drop_list <- c(drop_list, j)
                }
              }

              if (length(drop_list) > 0) {
                active_list[[k_chr]] <- setdiff(active_list[[k_chr]], drop_list)
              }

              if (!(0 %in% active_list[[k_chr]])) {
                active_list[[k_chr]] <- union(active_list[[k_chr]], 0)
              }

              remaining_exp <- setdiff(active_list[[k_chr]], 0)
              if (length(remaining_exp) == 0) {
                active_groups  <- as.integer(setdiff(active_groups, k))
                stopped_groups <- union(stopped_groups, k)
              }
            }
          }

          if (length(active_groups) == 0) break
          next_interim_id <- next_interim_id + 1L
          next
        }
      }
      break
    }

    if (length(active_groups) == 0) break

    accrual_eligible <- active_groups[n_by_k[active_groups] < max_n_k[active_groups]]
    if (length(accrual_eligible) == 0) break

    frac_vec <- n_by_k[accrual_eligible] / max_n_k[accrual_eligible]
    min_frac <- min(frac_vec)
    cand <- accrual_eligible[which(frac_vec == min_frac)]
    k_i <- if (length(cand) == 1) cand else sample(cand, 1)

    enrolled_count <- enrolled_count + 1L
    new_id <- enrolled_count

    choices <- active_list[[as.character(k_i)]]
    j_choice <- sample(choices, 1)

    YT <- rbinom(1, 1, piT_vec[j_choice + 1])
    p_true <- pE_mat[k_i, j_choice + 1]
    YE <- rbinom(1, 1, p_true)

    df_pat <- rbind(
      df_pat,
      data.frame(
        id = new_id, k = k_i, j = j_choice,
        treat = paste0("T", j_choice),
        YT = YT, YE = YE
      )
    )

    n_by_k[k_i] <- n_by_k[k_i] + 1L

    if (verbose && enrolled_count %% 50 == 0) {
      cat("Enrolled:", enrolled_count,
          " Active groups:", paste(active_groups, collapse=","), "\n")
    }
  }

  final_selection <- vector("list", K)
  names(final_selection) <- as.character(1:K)
  for (k in 1:K) final_selection[[as.character(k)]] <- integer(0)

  if (nrow(df_pat) > 0) {
    not_stopped <- setdiff(1:K, stopped_groups)

    if (length(not_stopped) > 0) {
      postF <- eff_posterior(df_pat = df_pat, K = K, J = J, Sk_list = Sk_list,
                             delta_E = delta_E, analysis_id = length(interim_fracs) + 1L)

      toxF <- tox_posterior_beta(
        df_pat = df_pat,
        J = J,
        prior_piT_a = prior_piT_a,
        prior_piT_b = prior_piT_b,
        delta_T = delta_T
      )

      for (k in not_stopped) {
        k_chr <- as.character(k)
        open_exp <- setdiff(active_list[[k_chr]], 0)
        if (length(open_exp) == 0) next

        nk_now <- n_by_k[k]
        PLk <- P_L_k(nk_now, max_n_k[k])

        for (j in open_exp) {
          pe <- postF$post_prob_eff[k, j + 1]
          ps <- toxF$post_prob_safe[j + 1]
          if (is.na(pe) || is.na(ps)) next

          if (pe >= PLk && ps >= PLk) {
            final_selection[[k_chr]] <- union(final_selection[[k_chr]], j)
          }
        }
      }
    }
  }

  for (k in 1:K) {
    k_chr <- as.character(k)
    if (length(final_selection[[k_chr]]) > 0) {
      final_selection[[k_chr]] <- sort(unique(final_selection[[k_chr]]))
    } else {
      final_selection[[k_chr]] <- integer(0)
    }
  }

  n_by_trt <- setNames(rep(0L, J + 1), paste0("T", 0:J))
  if (nrow(df_pat) > 0) {
    tt <- table(df_pat$j)
    for (j in 0:J) {
      key <- as.character(j)
      if (key %in% names(tt)) n_by_trt[paste0("T", j)] <- as.integer(tt[[key]])
    }
  }

  n_by_trt_by_k <- matrix(0L, nrow = K, ncol = J + 1,
                          dimnames = list(paste0("k", 1:K), paste0("T", 0:J)))
  if (nrow(df_pat) > 0) {
    tab_kj <- table(df_pat$k, df_pat$j)
    for (k in 1:K) {
      for (j in 0:J) {
        if (as.character(k) %in% rownames(tab_kj) &&
            as.character(j) %in% colnames(tab_kj)) {
          n_by_trt_by_k[k, j + 1] <- as.integer(tab_kj[as.character(k), as.character(j)])
        }
      }
    }
  }

  return(list(
    final_selection = final_selection,
    n_by_trt = n_by_trt,
    n_by_trt_by_k = n_by_trt_by_k,
    n_by_k = n_by_k,
    max_n_k = max_n_k,
    stopped_groups = stopped_groups
  ))
}
