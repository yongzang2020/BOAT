# ============================================================
# CD-AB: equal randomization, no PDO data, adaptive borrowing across subgroups by the
# treatment-specific subgroup clustering of BOAT. The trial flow is simulate_trial_noPDO() in
# trial_noPDO.R.
#
# The sampler is the clustering sampler of BOAT with
# the PDO parts removed: no muO, rho, sigma_eps or sigma_mu, and no patient random
# effect b, which the binary outcome alone cannot identify. The model is then
#   Y_R,i = I(Z_i > 0),  Z_i ~ N(mu_R,jk, 1),  pi_R,jk = Phi(mu_R,jk),
# with mu_R,jk equal within a cluster, mu ~ N(0, 10^2) per cluster, and a uniform
# prior over the cluster configurations of each treatment. The RJ moves, proposals
# and acceptance ratios are those of BOAT without the PDO link term.
# ============================================================


# ============================================================
# RJMCMC sampler without PDO (Uniform prior on cluster configurations)
# - returns post_pE_mean, post_prob_eff, post_z_mean, post_part_prob
# ============================================================
mcmc_cluster_efficacy_noPDO <- function(
    df_pat,
    K=3,
    J=3,
    Sk_list,
    n_iter = 10000,
    n_burn = 5000,
    n_thin = 1,

    # decision margin for efficacy
    delta_E = 0.0,

    # Hyperparameters
    prior_muE_mean = 0, prior_muE_sd = 1e1,

    seed = 1e3
){

  set.seed(seed)

  # uniform prior => cancels in acceptance ratio
  get_prior_prob <- function(trt_idx, z_vec, subgroups) 1.0

  # ============================================================
  # Minimal helpers
  # ============================================================
  canonize_col <- function(z_mat, col, subs){
    zcol <- z_mat[, col]
    zs   <- zcol[subs]
    if(all(is.na(zs))) return(z_mat)

    labs <- unique(zs[!is.na(zs)])          # first appearance order
    map  <- setNames(seq_along(labs), labs) # old -> 1..C

    idx <- !is.na(zs)
    zs[idx] <- unname(map[as.character(zs[idx])])
    zcol[subs] <- zs

    z_mat[, col] <- zcol
    z_mat
  }

  omega_from_zcol <- function(zcol){
    out <- ifelse(is.na(zcol), NA_integer_, ifelse(duplicated(zcol), 0L, 1L))
    as.integer(out)
  }

  key_col <- function(z_mat, col, subs){
    ztmp <- canonize_col(z_mat, col, subs)
    paste0(ztmp[subs, col], collapse = "")
  }

  # ============================================================
  # Precompute indices
  # ============================================================
  N <- nrow(df_pat)
  stopifnot(N >= 1)

  k_pat   <- df_pat$k
  col_pat <- df_pat$j + 1L

  pat_idx <- vector("list", K)
  for(k in 1:K){
    pat_idx[[k]] <- vector("list", J+1)
    for(col in 1:(J+1)){
      pat_idx[[k]][[col]] <- which(k_pat == k & col_pat == col)
    }
  }

  idx_pat_subs_col <- function(subs, col){
    unlist(lapply(subs, function(kk) pat_idx[[kk]][[col]]), use.names = FALSE)
  }

  loglik_YE_subs_col <- function(muE_mat, subs, col, YE_star){
    idx <- idx_pat_subs_col(subs, col)
    if(length(idx) == 0) return(0.0)
    mu_vec <- muE_mat[cbind(k_pat[idx], rep.int(col, length(idx)))]
    sum(dnorm(YE_star[idx] - mu_vec, 0, 1, log = TRUE))
  }

  # ============================================================
  # Initialize parameters
  # ============================================================
  muE <- matrix(NA_real_, nrow = K, ncol = J+1)
  muE_est <- matrix(0, nrow = K, ncol = J+1)

  for(k in 1:K){
    for(j in Sk_list[[k]]){
      col <- j + 1L
      muE[k, col] <- 0

      idx_kj <- pat_idx[[k]][[col]]
      if(length(idx_kj) > 0){
        n_kj <- length(idx_kj)
        y_kj <- sum(df_pat$YE[idx_kj] == 1)

        # Jeffreys smoothing: p_tilde = (y + 0.5) / (n + 1)
        p_tilde <- (y_kj + 0.5) / (n_kj + 1)
        p_tilde <- pmin(pmax(p_tilde, 1e-8), 1 - 1e-8)

        muE_est[k, col] <- qnorm(p_tilde)
      }
    }
  }

  # for each treatment column, which subgroups have it
  Kj_list <- vector("list", J+1)
  for(trt in 0:J) {
    Kj_list[[trt+1]] <- which(sapply(Sk_list, function(x) trt %in% x))
  }

  # init z: all separate within each column
  z <- matrix(NA_integer_, nrow = K, ncol = J+1)
  for(col in 1:(J+1)) {
    subs <- Kj_list[[col]]
    if(length(subs) > 0) z[subs, col] <- seq_along(subs)
    z <- canonize_col(z, col, subs)
  }

  omega <- matrix(NA_integer_, nrow = K, ncol = J+1)
  for(col in 1:(J+1)) omega[, col] <- omega_from_zcol(z[, col])

  YE_star <- ifelse(df_pat$YE==0, -0.1, 0.1)

  # proposal sd in split moves
  sigma_new <- 1

  # ============================================================
  # RJ update helper: length-2 column (2 partitions: 11 vs 12)
  # ============================================================
  rj_update_col_len2 <- function(muE, z, omega, col, subs, trt_idx){

    z <- canonize_col(z, col, subs)
    omega[, col] <- omega_from_zcol(z[, col])

    zc <- z[subs, col]
    ncl <- length(unique(zc[!is.na(zc)]))

    # defensive
    if(ncl <= 0) return(list(muE=muE, z=z, omega=omega))
    if(ncl > 2)  return(list(muE=muE, z=z, omega=omega))

    # enforce current consistency by direct assignment
    if(ncl == 1) {
      muE[subs, col] <- muE[subs[1], col]
    }

    muE_new <- muE
    z_new <- z
    omega_new <- omega
    log_extra <- 0

    if(ncl == 2){
      # merge 2 -> 1
      z_new[subs, col] <- 1L
      z_new <- canonize_col(z_new, col, subs)
      omega_new[, col] <- omega_from_zcol(z_new[, col])

      removed_param <- muE[subs[2], col]
      lp_param_old <- dnorm(removed_param, prior_muE_mean, prior_muE_sd, log=TRUE)
      lprop_param  <- dnorm(removed_param, muE_est[subs[2], col], sigma_new, log=TRUE)
      log_extra <- - lp_param_old + lprop_param

      # direct assignment after merge
      muE_new[subs[2], col] <- muE_new[subs[1], col]

    } else if(ncl == 1){
      # split 1 -> 2
      z_new[subs[1], col] <- 1L
      z_new[subs[2], col] <- 2L
      z_new <- canonize_col(z_new, col, subs)
      omega_new[, col] <- omega_from_zcol(z_new[, col])

      mu_prop <- rnorm(1, muE_est[subs[2], col], sigma_new)

      # direct assignment after split
      muE_new[subs[1], col] <- muE[subs[1], col]
      muE_new[subs[2], col] <- mu_prop

      lp_param_new <- dnorm(mu_prop, prior_muE_mean, prior_muE_sd, log=TRUE)
      lprop_param  <- dnorm(mu_prop, muE_est[subs[2], col], sigma_new, log=TRUE)
      log_extra <- lp_param_new - lprop_param
    }

    if(!identical(z_new[, col], z[, col])){
      lp_old <- log(get_prior_prob(trt_idx, z[, col], subs))
      lp_new <- log(get_prior_prob(trt_idx, z_new[, col], subs))

      ll_old <- loglik_YE_subs_col(muE,     subs, col, YE_star)
      ll_new <- loglik_YE_subs_col(muE_new, subs, col, YE_star)

      accept_log <- (ll_new - ll_old) + (lp_new - lp_old) + log_extra

      if(log(runif(1)) < accept_log){
        muE <- muE_new
        z <- z_new
        omega <- omega_new
      }
    }

    list(muE=muE, z=z, omega=omega)
  }

  # ============================================================
  # RJ update helper: length-3 column (5 partitions)
  # ============================================================
  rj_update_col_len3 <- function(muE, z, omega, col, subs, trt_idx){

    z <- canonize_col(z, col, subs)
    omega[, col] <- omega_from_zcol(z[, col])

    zc  <- z[subs, col]
    ncl <- length(unique(zc[!is.na(zc)]))

    # defensive
    if(ncl <= 0) return(list(muE=muE, z=z, omega=omega))
    if(ncl > 3)  return(list(muE=muE, z=z, omega=omega))

    # enforce current consistency by direct assignment
    if(ncl >= 1){
      labs_now <- unique(zc[!is.na(zc)])
      for(lb in labs_now){
        mem <- subs[which(zc == lb)]
        muE[mem, col] <- muE[mem[1], col]
      }
    }

    muE_new <- muE
    z_new <- z
    omega_new <- omega
    log_extra <- 0

    if(ncl == 3){
      # merge one pair: 3 -> 2
      pairs <- list(c(1,2), c(1,3), c(2,3))
      pr <- pairs[[sample(1:3, 1)]]
      a <- subs[pr[1]]
      bsub <- subs[pr[2]]
      csub <- setdiff(subs, c(a, bsub))

      z_new[subs, col] <- NA_integer_
      z_new[a, col]    <- 1L
      z_new[bsub, col] <- 1L
      z_new[csub, col] <- 2L

      z_new <- canonize_col(z_new, col, subs)
      omega_new[, col] <- omega_from_zcol(z_new[, col])

      removed_param <- muE[bsub, col]
      lp_param_old <- dnorm(removed_param, prior_muE_mean, prior_muE_sd, log=TRUE)
      lprop_param  <- dnorm(removed_param, muE_est[bsub, col], sigma_new, log=TRUE)
      log_extra <- - lp_param_old + lprop_param

      # direct assignment
      muE_new[bsub, col] <- muE_new[a, col]
      muE_new[csub, col] <- muE[csub, col]

    } else if(ncl == 1){
      # split one singleton: 1 -> 2
      singleton <- sample(subs, 1)
      others <- setdiff(subs, singleton)

      z_new[subs, col] <- NA_integer_
      z_new[others, col] <- 1L
      z_new[singleton, col] <- 2L

      z_new <- canonize_col(z_new, col, subs)
      omega_new[, col] <- omega_from_zcol(z_new[, col])

      mu_prop <- rnorm(1, muE_est[singleton, col], sigma_new)
      theta_all <- muE[subs[1], col]

      # direct assignment
      muE_new[others, col] <- theta_all
      muE_new[singleton, col] <- mu_prop

      lp_param_new <- dnorm(mu_prop, prior_muE_mean, prior_muE_sd, log=TRUE)
      lprop_param  <- dnorm(mu_prop, muE_est[singleton, col], sigma_new, log=TRUE)
      log_extra <- lp_param_new - lprop_param

    } else if(ncl == 2){
      # move among 2-cluster states: 2->1, 2->3, or 2->2(swap)
      members1 <- subs[which(zc == 1)]
      members2 <- subs[which(zc == 2)]
      if(length(members1) == 1){
        singleton <- members1
        pair <- members2
      } else {
        singleton <- members2
        pair <- members1
      }

      theta_single <- muE[singleton, col]
      theta_pair   <- muE[pair[1], col]

      move_type <- sample(c("merge_to_1", "split_to_3", "swap_2to2"), 1)

      if(move_type == "merge_to_1"){
        rep_pair <- min(pair)

        removed_param <- muE[singleton, col]
        lp_param_old <- dnorm(removed_param, prior_muE_mean, prior_muE_sd, log=TRUE)
        lprop_param  <- dnorm(removed_param, muE_est[singleton, col], sigma_new, log=TRUE)
        log_extra <- - lp_param_old + lprop_param

        z_new[subs, col] <- 1L
        z_new <- canonize_col(z_new, col, subs)
        omega_new[, col] <- omega_from_zcol(z_new[, col])

        # direct assignment: all merge to pair value
        muE_new[subs, col] <- muE_new[rep_pair, col]

      } else if(move_type == "split_to_3"){
        pair <- sort(pair)
        new_one  <- max(pair)
        keep_one <- min(pair)

        z_new[subs, col] <- NA_integer_
        z_new[keep_one, col] <- 1L
        z_new[new_one, col]  <- 2L
        z_new[singleton, col] <- 3L

        z_new <- canonize_col(z_new, col, subs)
        omega_new[, col] <- omega_from_zcol(z_new[, col])

        mu_prop <- rnorm(1, muE_est[new_one, col], sigma_new)

        # direct assignment
        muE_new[keep_one, col] <- theta_pair
        muE_new[singleton, col] <- theta_single
        muE_new[new_one, col] <- mu_prop

        lp_param_new <- dnorm(mu_prop, prior_muE_mean, prior_muE_sd, log=TRUE)
        lprop_param  <- dnorm(mu_prop, muE_est[new_one, col], sigma_new, log=TRUE)
        log_extra <- lp_param_new - lprop_param

      } else {
        # swap member in size-2 cluster: 2->2
        moved <- sample(pair, 1)
        remain <- setdiff(pair, moved)

        z_new[subs, col] <- NA_integer_
        z_new[singleton, col] <- 1L
        z_new[moved, col]     <- 1L
        z_new[remain, col]    <- 2L

        z_new <- canonize_col(z_new, col, subs)
        omega_new[, col] <- omega_from_zcol(z_new[, col])

        # direct assignment with explicit mapping
        muE_new[singleton, col] <- theta_single
        muE_new[moved, col] <- theta_single
        muE_new[remain, col] <- theta_pair

        log_extra <- 0
      }
    }

    if(!identical(z_new[, col], z[, col])){

      lp_old <- log(get_prior_prob(trt_idx, z[, col], subs))
      lp_new <- log(get_prior_prob(trt_idx, z_new[, col], subs))

      ll_old <- loglik_YE_subs_col(muE,     subs, col, YE_star)
      ll_new <- loglik_YE_subs_col(muE_new, subs, col, YE_star)

      accept_log <- (ll_new - ll_old) + (lp_new - lp_old) + log_extra

      if(log(runif(1)) < accept_log){
        muE <- muE_new
        z <- z_new
        omega <- omega_new
      }
    }

    list(muE=muE, z=z, omega=omega)
  }

  # ============================================================
  # Storage: only keep post-burn & thin
  # ============================================================
  keep_seq <- seq(n_burn + 1, n_iter, n_thin)
  n_keep <- length(keep_seq)
  keep_id <- 0L

  pE_sum  <- matrix(0, nrow = K, ncol = J+1)
  eff_sum <- matrix(0, nrow = K, ncol = J+1)
  z_sum   <- matrix(0, nrow = K, ncol = J+1)

  part_key <- matrix(NA_character_, nrow = n_keep, ncol = J+1)
  colnames(part_key) <- paste0("T", 0:J)

  # ============================================================
  # MCMC loop
  # ============================================================
  for(iter in 1:n_iter){

    # 1) Update muE cluster-wise per column
    for(col in 1:(J+1)){
      subs <- Kj_list[[col]]
      if(length(subs) == 0) next

      z <- canonize_col(z, col, subs)
      omega[, col] <- omega_from_zcol(z[, col])

      zc <- z[subs, col]
      labs <- unique(zc[!is.na(zc)])

      for(lb in labs){
        mem <- subs[which(zc == lb)]
        repk <- mem[1]

        idx_all <- idx_pat_subs_col(mem, col)
        n_total <- length(idx_all)
        sum_resid <- if(n_total > 0) sum(YE_star[idx_all]) else 0

        var.muE <- 1 / (n_total + 1/(prior_muE_sd^2))
        mean.muE <- var.muE * (sum_resid + prior_muE_mean/(prior_muE_sd^2))

        muE[repk, col] <- rnorm(1, mean.muE, sqrt(var.muE))
        if(length(mem) > 1) muE[mem, col] <- muE[repk, col]
      }
    }

    # 2) Update YE_star
    mean_ystar <- muE[cbind(k_pat, col_pat)]
    YE_star <- rtnorm(N, mean_ystar, 1,
                      lower = ifelse(df_pat$YE==1, 0, -Inf),
                      upper = ifelse(df_pat$YE==1, Inf, 0))

    # ========================================================
    # 3) RJMCMC clustering moves (supports len=2 or len=3 for any treatment)
    # ========================================================
    for(trt in 0:J){
      col <- trt + 1L
      subs <- Kj_list[[col]]
      if(length(subs) <= 1) next

      if(length(subs) == 2){
        out <- rj_update_col_len2(muE, z, omega, col = col, subs = subs, trt_idx = trt)
        muE <- out$muE
        z <- out$z
        omega <- out$omega
      } else if(length(subs) == 3){
        out <- rj_update_col_len3(muE, z, omega, col = col, subs = subs, trt_idx = trt)
        muE <- out$muE
        z <- out$z
        omega <- out$omega
      }
    }

    # ========================================================
    # Store samples
    # ========================================================
    if(iter >= n_burn + 1 && ((iter - (n_burn + 1)) %% n_thin == 0)){
      keep_id <- keep_id + 1L

      for(col in 1:(J+1)){
        subs <- Kj_list[[col]]
        if(length(subs) > 0) z <- canonize_col(z, col, subs)
      }

      pE <- matrix(NA_real_, nrow = K, ncol = J+1)
      pE[!is.na(muE)] <- pnorm(muE[!is.na(muE)])

      pE_sum <- pE_sum + ifelse(is.na(pE), 0, pE)
      z_sum  <- z_sum  + ifelse(is.na(z), 0, z)

      # efficacy indicator: compare vs control within subgroup
      for(k in 1:K){
        p0 <- pE[k, 1]
        if(is.na(p0)) next
        for(col in 1:(J+1)){
          if(is.na(pE[k, col])) next
          eff_sum[k, col] <- eff_sum[k, col] + as.integer(pE[k, col] > p0 + delta_E)
        }
      }

      for(col in 1:(J+1)){
        subs <- Kj_list[[col]]
        part_key[keep_id, col] <- if(length(subs) > 0) key_col(z, col, subs) else NA_character_
      }
    }
  }

  # posterior partition probabilities per treatment
  post_part_prob <- vector("list", J+1)
  names(post_part_prob) <- paste0("T", 0:J)

  for(col in 1:(J+1)){
    keys <- part_key[, col]
    keys <- keys[!is.na(keys)]
    if(length(keys) == 0){
      post_part_prob[[col]] <- data.frame(pattern=character(0), prob=numeric(0))
    } else {
      tt <- sort(table(keys), decreasing = TRUE)
      post_part_prob[[col]] <- data.frame(
        pattern = names(tt),
        prob = as.numeric(tt) / sum(tt),
        row.names = NULL
      )
    }
  }

  post_prob_eff <- eff_sum / n_keep
  post_prob_eff[is.na(muE)] <- NA

  list(
    post_pE_mean   = pE_sum / n_keep,
    post_prob_eff  = post_prob_eff,
    post_z_mean    = z_sum  / n_keep,
    post_part_prob = post_part_prob
  )
}


# ============================================================
# simulate_trial_CD_AB(): one simulated CD-AB trial
# Each analysis uses seed + 1000 * analysis_id for its MCMC, as in BOAT.
# ============================================================
simulate_trial_CD_AB <- function(
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
    prior_piT_a = 1,
    prior_piT_b = 1,
    seed = 123,
    verbose = FALSE
) {
  eff_posterior <- function(df_pat, K, J, Sk_list, delta_E, analysis_id) {
    mcmc_cluster_efficacy_noPDO(df_pat, K = K, J = J, Sk_list = Sk_list,
                                n_iter = mcmc_iter, n_burn = mcmc_burn, n_thin = mcmc_thin,
                                delta_E = delta_E,
                                seed = seed + 1000L * analysis_id)
  }
  simulate_trial_noPDO(pE_mat, piT_vec, max_n_k, eff_posterior = eff_posterior,
                       interim_fracs = interim_fracs, delta_E = delta_E, delta_T = delta_T,
                       lambda = lambda, prior_piT_a = prior_piT_a, prior_piT_b = prior_piT_b,
                       seed = seed, verbose = verbose)
}
