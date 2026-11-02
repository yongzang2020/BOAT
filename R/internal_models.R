# Internal computational routines of the BOAT design: the data-generating response means, the
# toxicity posterior, the RJMCMC sampler of the joint response model and the probit model used
# for adaptive randomization.

# Probit-scale response means of the data-generating model
.boat_pE_to_muE <- function(pE_mat, Sk_list, sigma_b) {
  K <- nrow(pE_mat)
  muE_kj <- vector("list", K)
  
  total_sd <- sqrt(1 + sigma_b^2)
  pE_clip <- pmin(pmax(pE_mat, 1e-6), 1 - 1e-6)
  muE_calc <- total_sd * qnorm(pE_clip)
  
  for (k in 1:K) {
    arms <- Sk_list[[k]]
    vals <- muE_calc[k, arms + 1]
    names(vals) <- paste0("T", arms)
    muE_kj[[k]] <- vals
  }
  return(muE_kj)
}

# Toxicity: independent Beta posteriors and Pr(pi_T,j < pi_T,0 + delta_T | data)
.boat_tox_posterior <- function(df_pat, J,
                               prior_piT_a = 1, prior_piT_b = 1,
                               delta_T = 0.0) {
  stopifnot(all(c("j", "YT") %in% names(df_pat)))
  stopifnot(J >= 1)
  
  a_post <- numeric(J + 1)
  b_post <- numeric(J + 1)
  
  for (j in 0:J) {
    idx <- which(df_pat$j == j)
    y <- if (length(idx) > 0) sum(df_pat$YT[idx]) else 0
    n <- length(idx)
    a_post[j + 1] <- prior_piT_a + y
    b_post[j + 1] <- prior_piT_b + n - y
  }
  
  post_piT_mean <- a_post / (a_post + b_post)
  
  post_prob_safe <- rep(NA_real_, J + 1)
  a0 <- a_post[1]; b0 <- b_post[1]
  
  for (j in 0:J) {
    if (j == 0) {
      post_prob_safe[j + 1] <- NA_real_
      next
    }
    if (delta_T >= 1) {
      post_prob_safe[j + 1] <- 1.0
      next
    }
    
    aj <- a_post[j + 1]; bj <- b_post[j + 1]
    
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

# Response: joint probit model of patient responses and PDO responses with RJMCMC clustering
# of each treatment's response mean across subgroups (uniform prior on the partitions)
.boat_mcmc_efficacy <- function(
    df_pat,
    df_pdo,
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
    prior_sigma_b_alpha = 0.01, prior_sigma_b_beta = 0.01,
    prior_sigma_eps_alpha = 0.01, prior_sigma_eps_beta = 0.01,
    prior_sigma_mu_alpha = 0.01, prior_sigma_mu_beta = 0.01,
    
    # rho priors
    prior_rho_mean = 0,
    prior_rho_sd = 1e1,
    
    seed = 1e3,
    verbose = TRUE
){
  
  set.seed(seed)
  if(verbose) cat("Using uniform prior on cluster configurations (no empirical prior).\n")
  
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
  
  id_pat  <- df_pat$id
  k_pat   <- df_pat$k
  col_pat <- df_pat$j + 1L
  
  pat_idx <- vector("list", K)
  for(k in 1:K){
    pat_idx[[k]] <- vector("list", J+1)
    for(col in 1:(J+1)){
      pat_idx[[k]][[col]] <- which(k_pat == k & col_pat == col)
    }
  }
  
  df_pdo$col <- df_pdo$j + 1L
  k_pdo   <- df_pdo$k
  col_pdo <- df_pdo$col
  id_pdo  <- df_pdo$id
  
  pdo_idx_by_kcol <- vector("list", K)
  for(k in 1:K){
    pdo_idx_by_kcol[[k]] <- vector("list", J+1)
    for(col in 1:(J+1)){
      pdo_idx_by_kcol[[k]][[col]] <- which(k_pdo == k & col_pdo == col)
    }
  }
  
  idx_pat_subs_col <- function(subs, col){
    unlist(lapply(subs, function(kk) pat_idx[[kk]][[col]]), use.names = FALSE)
  }
  
  loglik_YE_subs_col <- function(muE_mat, subs, col, YE_star, b){
    idx <- idx_pat_subs_col(subs, col)
    if(length(idx) == 0) return(0.0)
    mu_vec <- muE_mat[cbind(k_pat[idx], rep.int(col, length(idx)))]
    sum(dnorm(YE_star[idx] - b[idx] - mu_vec, 0, 1, log = TRUE))
  }
  
  log_org_link_col <- function(muE_mat, muO, col_idx, sigma_mu) {
    sum(dnorm(muO[,col_idx] - muE_mat[,col_idx], 0, sigma_mu, log=TRUE), na.rm=TRUE)
  }
  
  # ============================================================
  # Initialize parameters
  # ============================================================
  muE <- matrix(NA_real_, nrow = K, ncol = J+1)
  muE_est <- matrix(0, nrow = K, ncol = J+1)
  muO <- matrix(NA_real_, nrow = K, ncol = J+1)
  
  for(k in 1:K){
    for(j in Sk_list[[k]]){
      col <- j + 1L
      muE[k, col] <- 0
      muO[k, col] <- 0
      
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
  
  b <- rep(0, N)
  sigma_b <- 1
  rho <- 0
  sigma_eps <- 1
  sigma_mu <- 1
  YE_star <- ifelse(df_pat$YE==0, -0.1, 0.1)
  
  # PDO indices per patient id (robust to id not equal to 1..N)
  pdo_indices_list <- split(seq_len(nrow(df_pdo)), df_pdo$id)
  
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
      
      ll_old <- loglik_YE_subs_col(muE,     subs, col, YE_star, b)
      ll_new <- loglik_YE_subs_col(muE_new, subs, col, YE_star, b)
      
      l_org_old <- log_org_link_col(muE,     muO, col, sigma_mu)
      l_org_new <- log_org_link_col(muE_new, muO, col, sigma_mu)
      
      accept_log <- (ll_new - ll_old) + (l_org_new - l_org_old) + (lp_new - lp_old) + log_extra
      
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
        
        # direct assignment with explicit mapping:
        # moved joins singleton's cluster -> takes singleton parameter
        # the remaining member keeps the pair's parameter
        muE_new[singleton, col] <- theta_single
        muE_new[moved, col] <- theta_single
        muE_new[remain, col] <- theta_pair
        
        log_extra <- 0
      }
    }
    
    if(!identical(z_new[, col], z[, col])){
      
      lp_old <- log(get_prior_prob(trt_idx, z[, col], subs))
      lp_new <- log(get_prior_prob(trt_idx, z_new[, col], subs))
      
      ll_old <- loglik_YE_subs_col(muE,     subs, col, YE_star, b)
      ll_new <- loglik_YE_subs_col(muE_new, subs, col, YE_star, b)
      
      l_org_old <- log_org_link_col(muE,     muO, col, sigma_mu)
      l_org_new <- log_org_link_col(muE_new, muO, col, sigma_mu)
      
      accept_log <- (ll_new - ll_old) + (l_org_new - l_org_old) + (lp_new - lp_old) + log_extra
      
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
  
  if(verbose) cat("Starting MCMC sampling...\n")
  
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
        sum_resid <- if(n_total > 0) sum(YE_star[idx_all] - b[idx_all]) else 0
        
        m_cluster <- length(mem)
        sum_muO <- sum(muO[mem, col], na.rm = TRUE)
        
        var.muE <- 1 / (n_total + 1/(prior_muE_sd^2) + m_cluster/(sigma_mu^2))
        mean.muE <- var.muE * (sum_resid + prior_muE_mean/(prior_muE_sd^2) + sum_muO/(sigma_mu^2))
        
        muE[repk, col] <- rnorm(1, mean.muE, sqrt(var.muE))
        if(m_cluster > 1) muE[mem, col] <- muE[repk, col]
      }
    }
    
    # 2) Update muO
    for(k in 1:K){
      for(j in Sk_list[[k]]){
        col <- j + 1L
        idx_kj <- pdo_idx_by_kcol[[k]][[col]]
        n_kj <- length(idx_kj)
        
        var.muO <- 1 / (n_kj/(sigma_eps^2) + 1/(sigma_mu^2))
        
        sum_pdo_resid <- if(n_kj > 0){
          sum(df_pdo$YO[idx_kj] - rho*b[id_pdo[idx_kj]])
        } else 0
        
        mean.muO <- var.muO * (sum_pdo_resid/(sigma_eps^2) + muE[k,col]/(sigma_mu^2))
        muO[k,col] <- rnorm(1, mean.muO, sqrt(var.muO))
      }
    }
    
    # 3) Update b
    for(i in 1:N){
      pid <- id_pat[i]
      k_i <- k_pat[i]
      col <- col_pat[i]
      
      pdo_idx <- pdo_indices_list[[as.character(pid)]]
      M_i <- if(is.null(pdo_idx)) 0L else length(pdo_idx)
      
      var.b <- 1 / (1 + 1/(sigma_b^2) + M_i*rho^2/(sigma_eps^2))
      term1 <- YE_star[i] - muE[k_i, col]
      
      term2 <- if(M_i > 0){
        muO_vals <- muO[cbind(k_pdo[pdo_idx], col_pdo[pdo_idx])]
        rho * sum(df_pdo$YO[pdo_idx] - muO_vals)
      } else 0
      
      mean.b <- var.b * (term1 + term2/(sigma_eps^2))
      b[i] <- rnorm(1, mean.b, sqrt(var.b))
    }
    
    # 4) Update variances
    sigma_b <- sqrt(1/rgamma(1, prior_sigma_b_alpha + N/2,
                             prior_sigma_b_beta + sum(b^2)/2))
    
    resid_pdo <- df_pdo$YO - muO[cbind(k_pdo, col_pdo)] - rho*b[id_pdo]
    sigma_eps <- sqrt(1/rgamma(1,
                               prior_sigma_eps_alpha + nrow(df_pdo)/2,
                               prior_sigma_eps_beta + sum(resid_pdo^2)/2))
    
    diff_mu <- (muE - muO)[!is.na(muE)]
    sigma_mu <- sqrt(1/rgamma(1,
                              prior_sigma_mu_alpha + length(diff_mu)/2,
                              prior_sigma_mu_beta + sum(diff_mu^2)/2))
    
    # 5) Update rho > 0
    denom_rho <- sum(b[id_pdo]^2)
    var.rho <- 1 / (denom_rho/(sigma_eps^2) + 1/(prior_rho_sd^2))
    lik_term <- sum((df_pdo$YO - muO[cbind(k_pdo, col_pdo)]) * b[id_pdo])/(sigma_eps^2)
    prior_term <- prior_rho_mean/(prior_rho_sd^2)
    mean.rho <- var.rho*(lik_term + prior_term)
    rho <- rtnorm(1, mean.rho, sqrt(var.rho), lower=0, upper=Inf)
    
    # 6) Update YE_star
    mean_ystar <- muE[cbind(k_pat, col_pat)] + b
    YE_star <- rtnorm(N, mean_ystar, 1,
                      lower = ifelse(df_pat$YE==1, 0, -Inf),
                      upper = ifelse(df_pat$YE==1, Inf, 0))
    
    # ========================================================
    # 7) RJMCMC clustering moves (supports len=2 or len=3 for any treatment)
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
      } else {
        # current implementation handles columns appearing in 2 or 3 subgroups
        # (no RJ move for other sizes)
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
      
      sd_total <- sqrt(1 + sigma_b^2)
      pE <- matrix(NA_real_, nrow = K, ncol = J+1)
      pE[!is.na(muE)] <- pnorm(muE[!is.na(muE)] / sd_total)
      
      pE_sum <- pE_sum + pE
      z_sum  <- z_sum  + z
      
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
    
    if(verbose && iter %% 500 == 0) cat("Iteration:", iter, "\n")
  }
  
  if(verbose) cat("MCMC Complete.\n")
  
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
  
  list(
    post_pE_mean   = pE_sum / n_keep,
    post_prob_eff  = eff_sum / n_keep,
    post_z_mean    = z_sum  / n_keep,
    post_part_prob = post_part_prob
  )
}

# Albert-Chib probit regression of patient response on the patient's mean PDO response, fitted
# separately for every open (subgroup, treatment) pair; used for response-adaptive randomization
.boat_ac_probit <- function(df_pat, df_pdo, K,
                                            active_groups,
                                            active_list,          # list by k: active arms
                                            n_iter_pp = 2000,     # total iterations for AC-probit
                                            n_burn_pp = 1000,     # burn-in iterations
                                            n_thin_pp = 1,        # thinning interval
                                            prior_mean = c(0, 0),
                                            prior_sd   = c(10, 10),
                                            standardize_x = TRUE,
                                            verbose = FALSE) {
  
  # organoid mean per (id,k,j) across replicates
  pdo_mean <- aggregate(YO ~ id + k + j, data = df_pdo, FUN = mean)
  
  gibbs_one <- function(y, x, n_iter_pp, n_burn_pp, n_thin_pp, m0, V0){
    n <- length(y)
    X <- cbind(1, x)
    p <- 2
    
    V0inv <- solve(V0)
    XtX   <- crossprod(X)
    
    alpha <- m0
    ystar <- ifelse(y == 1, 0.5, -0.5)
    
    total_iter <- n_iter_pp
    
    # keep from (burn+1), every thin steps
    keep_seq <- if (total_iter >= n_burn_pp + 1) {
      seq(n_burn_pp + 1, total_iter, by = n_thin_pp)
    } else {
      integer(0)
    }
    n_keep <- length(keep_seq)
    
    a0 <- numeric(n_keep)
    a1 <- numeric(n_keep)
    keep_id <- 0L
    
    for(it in 1:total_iter){
      mu <- as.vector(X %*% alpha)
      
      for(i in 1:n){
        if(y[i] == 1){
          ystar[i] <- msm::rtnorm(1, mean = mu[i], sd = 1, lower = 0, upper = Inf)
        } else {
          ystar[i] <- msm::rtnorm(1, mean = mu[i], sd = 1, lower = -Inf, upper = 0)
        }
      }
      
      Vn <- solve(XtX + V0inv)
      mn <- Vn %*% (crossprod(X, ystar) + V0inv %*% m0)
      alpha <- as.vector(mn + t(chol(Vn)) %*% rnorm(p))
      
      if(it >= n_burn_pp + 1 && ((it - (n_burn_pp + 1)) %% n_thin_pp == 0)){
        keep_id <- keep_id + 1L
        a0[keep_id] <- alpha[1]
        a1[keep_id] <- alpha[2]
      }
    }
    list(a0 = a0, a1 = a1)
  }
  
  m0 <- as.numeric(prior_mean)
  V0 <- diag(prior_sd^2, 2, 2)
  
  # output structure: list by k, inside list by j
  alpha_draws <- vector("list", K)
  names(alpha_draws) <- as.character(1:K)
  for(k in 1:K) alpha_draws[[as.character(k)]] <- list()
  
  for(k in active_groups){
    k_chr <- as.character(k)
    arms_k <- active_list[[k_chr]]
    if(length(arms_k) == 0) next
    
    for(j in arms_k){
      j_chr <- as.character(j)
      
      dat_pat_kj <- df_pat[df_pat$k == k & df_pat$j == j, c("id","k","j","YE")]
      if(nrow(dat_pat_kj) < 3) next
      if(length(unique(dat_pat_kj$YE)) < 2) next
      
      # merge barYO for (id,k,j)
      dat <- merge(dat_pat_kj, pdo_mean, by = c("id","k","j"), all.x = TRUE)
      dat <- dat[!is.na(dat$YO), ]
      if(nrow(dat) < 3) next
      if(length(unique(dat$YE)) < 2) next
      
      y <- as.integer(dat$YE)
      x_raw <- as.numeric(dat$YO)
      
      if(standardize_x){
        x_mean <- mean(x_raw)
        x_sd   <- sd(x_raw)
        if(!is.finite(x_sd) || x_sd < 1e-8) next
        x <- (x_raw - x_mean) / x_sd
      } else {
        x_mean <- 0
        x_sd   <- 1
        x <- x_raw
      }
      
      if(verbose) cat(sprintf("AC probit by (k,j): k=%d, j=%d, n=%d\n", k, j, length(y)))
      
      out <- gibbs_one(y = y, x = x,
                       n_iter_pp = n_iter_pp,
                       n_burn_pp = n_burn_pp,
                       n_thin_pp = n_thin_pp,
                       m0 = m0, V0 = V0)
      
      alpha_draws[[k_chr]][[j_chr]] <- list(
        a0 = out$a0,
        a1 = out$a1,
        x_mean = x_mean,
        x_sd   = x_sd
      )
    }
  }
  
  alpha_draws
}
