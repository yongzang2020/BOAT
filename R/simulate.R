#' Simulate one BOAT trial
#'
#' Simulates a complete trial under a scenario: patients enter the subgroups in a balanced way,
#' each patient's PDOs are tested with every open treatment of the subgroup, the patient is
#' randomized ([boat_randomize()]), interim analyses ([boat_fit()], [boat_decide()],
#' [boat_rar_model()]) are done when every active subgroup reaches its interim sample size, and
#' the final analysis gives the recommended treatments.
#'
#' Data are generated from the BOAT model. The patient random effect is
#' \eqn{b_i \sim N(0, \sigma_E^2)}, the clinical response is \eqn{Y_{R,i} = I(\mu_{R,j,k} + b_i + e_i > 0)}
#' with \eqn{e_i \sim N(0, 1)}, and the PDO replicates are
#' \eqn{Y_{O} \sim N(\mu_{O,j,k} + \rho b_i, \sigma_O^2)} with \eqn{\mu_{O,j,k} \sim N(\mu_{R,j,k}, \sigma_{OR}^2)}.
#' The response means \eqn{\mu_{R,j,k}} are set so that the response probabilities equal those of
#' the scenario.
#'
#' @param design A `boat_design` object.
#' @param scenario A `boat_scenario` object whose evaluated treatments match `design$arms`.
#' @param dgm Data-generating parameters from [boat_dgm()].
#' @param seed Random seed of the trial.
#' @param keep_data Keep the simulated patient and PDO data in the result.
#' @return An object of class `boat_trial` with elements `final_selection` (recommended
#'   treatments by subgroup), `n_by_trt` (patients per treatment), `n_by_trt_by_k` (patients per
#'   subgroup and treatment), `n_by_k`, `stopped_groups`, `interims` (decision tables of the
#'   interim analyses) and, if `keep_data = TRUE`, `pat` and `pdo`.
#' @examples
#' \donttest{
#' fast <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)), N = 30, lambda = 0.70,
#'                     mcmc = list(iter = 400, burn = 200), rar_mcmc = list(iter = 200, burn = 100))
#' truth <- boat_scenario(pR = rbind(c(0.50, 0.30, NA, 0.70),
#'                                   c(0.50, 0.70, 0.50, 0.30),
#'                                   c(0.30,   NA, 0.50, 0.30)),
#'                        piT = c(0.25, 0.35, 0.20, 0.15))
#' boat_simulate_trial(fast, truth, seed = 1001)
#' }
#' @export
boat_simulate_trial <- function(design, scenario, dgm = boat_dgm(), seed = 1001, keep_data = FALSE) {
  .check_design(design)
  .check_scenario(scenario, design)
  if (!inherits(dgm, "boat_dgm")) stop("dgm must be created by boat_dgm()")
  adaptive <- design$randomization == "adaptive"
  pE_mat <- scenario$pR; piT_vec <- scenario$piT
  sigma_b <- dgm$sigma_E; rho <- dgm$rho; sigma_eps <- dgm$sigma_O; sigma_mu <- dgm$sigma_OR; M <- dgm$M
  lambda <- design$lambda; interim_fracs <- design$interim

  set.seed(seed)
  K <- nrow(pE_mat)
  J <- ncol(pE_mat) - 1
  max_n_k <- design$N
  Sk_list <- vector("list", K)
  for (k in 1:K) Sk_list[[k]] <- which(!is.na(pE_mat[k, ])) - 1

  muE_kj <- .boat_pE_to_muE(pE_mat, Sk_list, sigma_b)
  muE_mat <- matrix(NA, nrow = K, ncol = J + 1)
  for (k in 1:K) {
    for (name in names(muE_kj[[k]])) {
      j_idx <- as.integer(sub("T", "", name))
      muE_mat[k, j_idx + 1] <- muE_kj[[k]][name]
    }
  }
  muO_mat <- muE_mat
  for (k in 1:K) for (j in Sk_list[[k]]) muO_mat[k, j + 1] <- rnorm(1, muE_mat[k, j + 1], sigma_mu)

  # data stored as growing vectors
  pat <- list(id = integer(0), k = integer(0), j = numeric(0), YT = integer(0), YE = numeric(0))
  pdo <- list(id = integer(0), k = integer(0), j = numeric(0), YO = numeric(0))
  df_pat <- function() data.frame(pat)
  df_pdo <- function() data.frame(pdo)

  active_list <- lapply(Sk_list, function(x) x)
  names(active_list) <- as.character(1:K)
  active_groups <- as.integer(1:K)
  stopped_groups <- integer(0)
  n_by_k <- rep(0L, K)
  interim_targets <- .interim_targets(design)
  next_interim_id <- 1L
  alpha_draws_current <- NULL
  enrolled_count <- 0L
  interims <- list()

  while (TRUE) {
    # interim analyses: whenever every active subgroup has reached the next interim size
    repeat {
      if (next_interim_id <= length(interim_fracs) && length(active_groups) > 0) {
        targets_k <- interim_targets[, next_interim_id]
        cond_vec <- sapply(active_groups, function(k) n_by_k[k] >= targets_k[k])
        if (all(cond_vec)) {
          if (length(pat$id) > 0) {
            dp <- df_pat(); dd <- df_pdo()
            postE <- .fit_efficacy(design, dp, dd, seed = seed + 1000L * next_interim_id)
            tox <- .fit_toxicity(design, dp)
            tab <- .decision_table(design, postE$post_prob_eff, tox$post_prob_safe, active_list, active_groups, n_by_k)
            res <- .interim_decide(postE$post_prob_eff, tox$post_prob_safe, active_list, active_groups,
                                   n_by_k, max_n_k, lambda, stopped_groups)
            tab$decision <- ifelse(is.na(tab$pr_eff) | is.na(tab$pr_safe), "continue",
                                   ifelse(tab$pr_eff < tab$cutoff | tab$pr_safe < tab$cutoff, "drop", "continue"))
            interims[[next_interim_id]] <- tab
            active_list <- res$active_list; active_groups <- res$active_groups; stopped_groups <- res$stopped_groups
          }
          if (length(active_groups) == 0) break
          if (adaptive) {
            alpha_draws_current <- if (length(pat$id) > 0)
              .fit_rar(design, df_pat(), df_pdo(), active_groups, active_list) else NULL
          }
          next_interim_id <- next_interim_id + 1L
          next
        }
      }
      break
    }
    if (length(active_groups) == 0) break

    # next patient: the eligible subgroup with the smallest accrued fraction (ties at random)
    accrual_eligible <- active_groups[n_by_k[active_groups] < max_n_k[active_groups]]
    if (length(accrual_eligible) == 0) break
    frac_vec <- n_by_k[accrual_eligible] / max_n_k[accrual_eligible]
    min_frac <- min(frac_vec)
    cand <- accrual_eligible[which(frac_vec == min_frac)]
    k_i <- if (length(cand) == 1) cand else sample(cand, 1)
    enrolled_count <- enrolled_count + 1L
    new_id <- enrolled_count
    b_i <- rnorm(1, 0, sigma_b)

    # PDOs tested with every open treatment of the subgroup
    choices <- active_list[[as.character(k_i)]]
    barYO_map <- setNames(rep(NA_real_, length(choices)), as.character(choices))
    for (j in choices) {
      YO_vals <- rnorm(M, mean = muO_mat[k_i, j + 1] + rho * b_i, sd = sigma_eps)
      barYO_map[as.character(j)] <- mean(YO_vals)
      pdo$id <- c(pdo$id, rep(new_id, M)); pdo$k <- c(pdo$k, rep(k_i, M))
      pdo$j <- c(pdo$j, rep(j, M)); pdo$YO <- c(pdo$YO, YO_vals)
    }

    # randomization
    if (adaptive) {
      probs <- .alloc_probs(design, alpha_draws_current, k_i, barYO_map, choices)
      j_choice <- .assign_adaptive(probs, choices)
    } else {
      j_choice <- .assign_equal(choices)
    }

    # outcomes
    YT <- rbinom(1, 1, piT_vec[j_choice + 1])
    latent_Y <- muE_mat[k_i, j_choice + 1] + b_i + rnorm(1, 0, 1)
    YE <- ifelse(latent_Y > 0, 1, 0)
    pat$id <- c(pat$id, new_id); pat$k <- c(pat$k, k_i); pat$j <- c(pat$j, j_choice)
    pat$YT <- c(pat$YT, YT); pat$YE <- c(pat$YE, YE)
    n_by_k[k_i] <- n_by_k[k_i] + 1L
  }

  # final analysis
  final_selection <- setNames(rep(list(integer(0)), K), as.character(1:K))
  final_table <- NULL
  if (length(pat$id) > 0) {
    not_stopped <- setdiff(1:K, stopped_groups)
    if (length(not_stopped) > 0) {
      dp <- df_pat(); dd <- df_pdo()
      postF <- .fit_efficacy(design, dp, dd, seed = seed + 1000L * (length(interim_fracs) + 1))
      toxF <- .fit_toxicity(design, dp)
      final_selection <- .final_select(postF$post_prob_eff, toxF$post_prob_safe, active_list, stopped_groups,
                                       n_by_k, max_n_k, lambda, K)
      final_table <- .decision_table(design, postF$post_prob_eff, toxF$post_prob_safe, active_list,
                                     not_stopped, n_by_k)
    }
  }

  n_by_trt <- setNames(rep(0L, J + 1), paste0("T", 0:J))
  n_by_trt_by_k <- matrix(0L, nrow = K, ncol = J + 1, dimnames = list(paste0("k", 1:K), paste0("T", 0:J)))
  for (i in seq_along(pat$id)) {
    n_by_trt[pat$j[i] + 1] <- n_by_trt[pat$j[i] + 1] + 1L
    n_by_trt_by_k[pat$k[i], pat$j[i] + 1] <- n_by_trt_by_k[pat$k[i], pat$j[i] + 1] + 1L
  }
  out <- list(final_selection = final_selection, n_by_trt = n_by_trt, n_by_trt_by_k = n_by_trt_by_k,
              n_by_k = n_by_k, stopped_groups = stopped_groups, interims = interims, final = final_table,
              seed = seed, randomization = design$randomization, arms = design$arms)
  if (keep_data) {
    out$pat <- data.frame(id = pat$id, k = pat$k, j = pat$j, YT = pat$YT, YR = pat$YE)
    out$pdo <- data.frame(id = pdo$id, k = pdo$k, j = pdo$j, YO = pdo$YO)
  }
  structure(out, class = "boat_trial")
}

#' Generate trial data
#'
#' Generates data of `n` patients per subgroup under a scenario, with equal randomization among
#' the treatments of each subgroup and the data-generating model of [boat_simulate_trial()].
#' Useful for trying the functions that conduct a trial.
#'
#' @inheritParams boat_simulate_trial
#' @param n Number of patients per subgroup.
#' @return A list with the patient data `pat` (columns `id`, `k`, `j`, `YT`, `YR`) and the PDO
#'   data `pdo` (columns `id`, `k`, `j`, `YO`) in the format of [boat_fit()].
#' @examples
#' design <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)))
#' truth <- boat_scenario(pR = rbind(c(0.50, 0.30, NA, 0.70),
#'                                   c(0.50, 0.70, 0.50, 0.30),
#'                                   c(0.30,   NA, 0.50, 0.30)),
#'                        piT = c(0.25, 0.35, 0.20, 0.15))
#' dat <- boat_generate_data(design, truth, n = 10, seed = 1)
#' head(dat$pat)
#' @export
boat_generate_data <- function(design, scenario, n, dgm = boat_dgm(), seed = NULL) {
  .check_design(design)
  .check_scenario(scenario, design)
  if (!is.null(seed)) set.seed(seed)
  K <- design$K
  muE <- matrix(NA_real_, K, design$J + 1)
  muE[!is.na(scenario$pR)] <- sqrt(1 + dgm$sigma_E^2) * qnorm(pmin(pmax(scenario$pR[!is.na(scenario$pR)], 1e-6), 1 - 1e-6))
  muO <- muE
  muO[!is.na(muE)] <- rnorm(sum(!is.na(muE)), muE[!is.na(muE)], dgm$sigma_OR)
  k <- rep(seq_len(K), each = n)
  pat <- data.frame(id = seq_along(k), k = k, j = NA_integer_, YT = NA_integer_, YR = NA_integer_)
  pdo <- vector("list", length(k))
  for (i in seq_along(k)) {
    arms <- design$arms[[k[i]]]
    b <- rnorm(1, 0, dgm$sigma_E)
    pdo[[i]] <- do.call(rbind, lapply(arms, function(j)
      data.frame(id = i, k = k[i], j = j, YO = rnorm(dgm$M, muO[k[i], j + 1] + dgm$rho * b, dgm$sigma_O))))
    j <- arms[sample.int(length(arms), 1)]
    pat$j[i] <- j
    pat$YT[i] <- rbinom(1, 1, scenario$piT[j + 1])
    pat$YR[i] <- as.integer(muE[k[i], j + 1] + b + rnorm(1) > 0)
  }
  list(pat = pat, pdo = do.call(rbind, pdo))
}

#' Operating characteristics
#'
#' Simulates `nsim` trials with [boat_simulate_trial()] and summarizes them by subgroup and
#' treatment.
#'
#' For each subgroup and treatment the summary gives the percentage of trials in which the
#' treatment is recommended (`final_rec`; for the control, the percentage of trials with no
#' experimental treatment recommended), the mean number of patients (`n_patients`) and their
#' percentage (`pct_patients`). For a subgroup with admissible treatments it gives the power
#' `power1` (\eqn{P_{k,1}}: some treatment is recommended and all recommended treatments are
#' admissible) and `power2` (\eqn{P_{k,2}}: at least one admissible treatment is recommended). For a
#' subgroup without admissible treatment it gives the family-wise error rate `fwer`.
#'
#' @inheritParams boat_simulate_trial
#' @param nsim Number of simulated trials.
#' @param seeds Seeds of the trials.
#' @param cores Number of cores (forked processes; 1 on Windows).
#' @return An object of class `boat_oc` with the summary `table` (long format) and the
#'   simulated trials `trials`.
#' @examples
#' \donttest{
#' fast <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)), N = 30, lambda = 0.70,
#'                     mcmc = list(iter = 400, burn = 200), rar_mcmc = list(iter = 200, burn = 100))
#' truth <- boat_scenario(pR = rbind(c(0.50, 0.30, NA, 0.70),
#'                                   c(0.50, 0.70, 0.50, 0.30),
#'                                   c(0.30,   NA, 0.50, 0.30)),
#'                        piT = c(0.25, 0.35, 0.20, 0.15))
#' boat_get_oc(fast, truth, nsim = 2)
#' }
#' @export
boat_get_oc <- function(design, scenario, nsim = 1000, seeds = 1000 + seq_len(nsim), dgm = boat_dgm(),
                        cores = 1) {
  .check_design(design)
  .check_scenario(scenario, design)
  trials <- .run_parallel(seeds, function(s) boat_simulate_trial(design, scenario, dgm, seed = s), cores)
  structure(list(table = .summarize_trials(trials, scenario, design), trials = trials, design = design,
                 scenario = scenario, dgm = dgm, seeds = seeds),
            class = "boat_oc")
}

#' Calibrate the tuning parameter lambda
#'
#' For each value of `lambda`, simulates the scenarios and computes the family-wise error rate
#' \eqn{\mathrm{FWER}_k} of every subgroup without admissible treatment. Larger `lambda` gives
#' stricter cutoffs and a lower FWER. The chosen \eqn{\lambda} is the value whose largest FWER is
#' closest to `alpha`. Ties go to the smallest value with FWER at most `alpha`.
#'
#' @inheritParams boat_get_oc
#' @param scenarios List of `boat_scenario` objects with at least one null subgroup.
#' @param lambda Values of \eqn{\lambda} to evaluate.
#' @param alpha Nominal FWER.
#' @return An object of class `boat_calibration` with the FWER of every null subgroup for every
#'   \eqn{\lambda} (`table`), the largest FWER per \eqn{\lambda} (`max_fwer`) and the chosen
#'   value (`lambda`).
#' @examples
#' \dontrun{
#' design <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)))
#' # subgroup 2 has no admissible treatment
#' truth <- boat_scenario(pR = rbind(c(0.50, 0.30, NA, 0.70),
#'                                   c(0.50, 0.70, 0.50, 0.30),
#'                                   c(0.30,   NA, 0.50, 0.30)),
#'                        piT = c(0.25, 0.35, 0.20, 0.15))
#' boat_calibrate(design, list(truth), lambda = seq(0.64, 0.76, by = 0.02), nsim = 1000, cores = 8)
#' }
#' @export
boat_calibrate <- function(design, scenarios, lambda = seq(0.60, 0.80, by = 0.02), nsim = 1000,
                           seeds = 1000 + seq_len(nsim), alpha = 0.10, dgm = boat_dgm(), cores = 1) {
  .check_design(design)
  if (inherits(scenarios, "boat_scenario")) scenarios <- list(scenarios)
  rows <- list()
  for (lam in lambda) {
    d <- design; d$lambda <- lam
    for (s in seq_along(scenarios)) {
      tab <- boat_get_oc(d, scenarios[[s]], seeds = seeds, dgm = dgm, cores = cores)$table
      f <- tab[tab$measure == "fwer", ]
      if (nrow(f)) rows[[length(rows) + 1]] <- data.frame(lambda = lam, scenario = s, subgroup = f$subgroup,
                                                          fwer = f$value / 100)
    }
  }
  if (!length(rows)) stop("none of the scenarios has a subgroup without admissible treatment")
  tab <- do.call(rbind, rows)
  mx <- aggregate(fwer ~ lambda, data = tab, FUN = max)
  structure(list(table = tab, max_fwer = mx, lambda = .choose_lambda(mx, alpha), alpha = alpha),
            class = "boat_calibration")
}

# the lambda whose largest FWER is closest to alpha (mx: data frame with columns lambda, fwer)
.choose_lambda <- function(mx, alpha) {
  dist <- abs(mx$fwer - alpha)
  best <- mx[abs(dist - min(dist)) < 1e-12, ]
  if (any(best$fwer <= alpha)) min(best$lambda[best$fwer <= alpha]) else min(best$lambda)
}

#' @export
print.boat_trial <- function(x, ...) {
  cat(sprintf("Simulated %s trial (seed %d): %d patients\n",
              if (x$randomization == "adaptive") "BOAT" else "BOAT (equal randomization)", x$seed, sum(x$n_by_k)))
  for (k in seq_along(x$final_selection)) {
    rec <- x$final_selection[[k]]
    arms <- x$arms[[k]]
    cat(sprintf("  B%d: %s; patients per treatment %s\n", k,
                if (length(rec)) paste("recommend", paste0("A", rec, collapse = ", ")) else "no treatment recommended",
                paste(sprintf("A%d %d", arms, x$n_by_trt_by_k[k, arms + 1]), collapse = ", ")))
  }
  if (length(x$stopped_groups)) cat("  accrual stopped early in", paste0("B", sort(x$stopped_groups), collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.boat_oc <- function(x, ...) {
  cat(sprintf("Operating characteristics of %s from %d simulated trials\n",
              if (x$design$randomization == "adaptive") "BOAT" else "BOAT (equal randomization)", length(x$trials)))
  .print_oc_table(x$table, x$design)
  cat("% Final Rec: trials recommending the treatment (A0: no treatment recommended)\n")
  cat("# Patients, % Patients: mean number and percentage of patients\n")
  cat("(P_k1, P_k2): power of a subgroup with admissible treatments; FWER_k: of a subgroup without\n")
  invisible(x)
}

#' @export
print.boat_calibration <- function(x, ...) {
  cat(sprintf("Largest FWER over the null subgroups (nominal %.2f):\n", x$alpha))
  print(data.frame(lambda = x$max_fwer$lambda, max_fwer = round(x$max_fwer$fwer, 3)), row.names = FALSE)
  cat(sprintf("Chosen lambda (largest FWER closest to %.2f): %s\n", x$alpha, format(x$lambda)))
  invisible(x)
}

# ---- internal ----
.run_parallel <- function(x, f, cores) {
  if (cores > 1 && .Platform$OS.type != "windows") {
    res <- parallel::mclapply(x, f, mc.cores = cores, mc.preschedule = FALSE)
    bad <- vapply(res, inherits, logical(1), "try-error")
    if (any(bad)) stop("simulation failed: ", as.character(res[[which(bad)[1]]]))
    res
  } else lapply(x, f)
}

# summary by subgroup and treatment
.summarize_trials <- function(trials, scenario, design) {
  R <- length(trials); K <- design$K; J <- design$J
  rec <- matrix(0, K, J + 1); n_sum <- matrix(0, K, J + 1); hit1 <- hit2 <- numeric(K)
  for (tr in trials) {
    for (k in 1:K) {
      ch <- tr$final_selection[[as.character(k)]]; ch <- sort(unique(ch[ch %in% 1:J]))
      if (length(ch) == 0) rec[k, 1] <- rec[k, 1] + 1 else rec[k, ch + 1] <- rec[k, ch + 1] + 1
      A <- .admissible(scenario, k, design$delta_R, design$delta_T)
      if (length(ch) > 0 && length(A) > 0) {
        hit1[k] <- hit1[k] + all(ch %in% A)
        hit2[k] <- hit2[k] + any(ch %in% A)
      }
    }
    n_sum <- n_sum + tr$n_by_trt_by_k
  }
  mean_n <- n_sum / R; pct_n <- 100 * mean_n / rowSums(mean_n)
  out <- list()
  for (k in 1:K) {
    for (a in 0:J) {
      if (is.na(scenario$pR[k, a + 1])) next
      out[[length(out) + 1]] <- data.frame(subgroup = k, treatment = a,
                                           measure = c("final_rec", "n_patients", "pct_patients"),
                                           value = round(c(100 * rec[k, a + 1] / R, mean_n[k, a + 1], pct_n[k, a + 1]), 1))
    }
    if (length(.admissible(scenario, k, design$delta_R, design$delta_T)) > 0) {
      out[[length(out) + 1]] <- data.frame(subgroup = k, treatment = NA, measure = c("power1", "power2"),
                                           value = round(100 * c(hit1[k], hit2[k]) / R, 1))
    } else {
      out[[length(out) + 1]] <- data.frame(subgroup = k, treatment = NA, measure = "fwer",
                                           value = round(100 * (1 - rec[k, 1] / R), 1))
    }
  }
  do.call(rbind, out)
}

# one block per subgroup: rows = measures, columns = treatments evaluated in the subgroup
.print_oc_table <- function(tab, design) {
  meas <- c(final_rec = "% Final Rec", n_patients = "# Patients", pct_patients = "% Patients")
  for (k in seq_len(design$K)) {
    t <- tab[tab$subgroup == k, ]
    arms <- sort(unique(t$treatment[!is.na(t$treatment)]))
    m <- sapply(arms, function(a) sapply(names(meas), function(ms) sprintf("%.1f", t$value[t$treatment %in% a & t$measure == ms])))
    m <- matrix(m, nrow = length(meas), dimnames = list(unname(meas), paste0("A", arms)))
    p <- t[is.na(t$treatment), ]
    extra <- if (all(c("power1", "power2") %in% p$measure))
      sprintf("(P_k1, P_k2) = (%.1f, %.1f)", p$value[p$measure == "power1"], p$value[p$measure == "power2"])
    else sprintf("FWER_k = %.1f", p$value[p$measure == "fwer"])
    cat(sprintf("\nB%d   %s\n", k, extra))
    print(noquote(m), right = TRUE)
  }
  cat("\n")
}
