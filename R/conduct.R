#' Fit the BOAT model to the trial data
#'
#' Fits the joint model of the patients' clinical responses and their PDO responses, with
#' adaptive clustering of each treatment's response mean across subgroups by reversible jump
#' MCMC, and the beta-binomial toxicity model. Used at every interim and at the final analysis.
#'
#' @param design A `boat_design` object.
#' @param pat Patient data: a data frame with one row per patient and columns `id`, `k`
#'   (subgroup, `1:K`), `j` (assigned treatment, `0` = control), `YT` (toxicity, 0/1) and `YR`
#'   (clinical response, 0/1).
#' @param pdo PDO data: a data frame with one row per PDO replicate and columns `id` (patient),
#'   `k`, `j` (treatment tested) and `YO` (PDO response). Each patient's PDOs are tested with the
#'   treatments open in the patient's subgroup at enrollment.
#' @param seed Seed set before the MCMC, or `NULL` to use the current random number stream.
#' @return An object of class `boat_fit` with elements
#'   \describe{
#'     \item{`summary`}{data frame with one row per subgroup and treatment: posterior mean
#'       response probability `pR`, \eqn{\Pr(\pi_{R,j,k} > \pi_{R,0,k} + \delta_R)} (`pr_eff`),
#'       posterior mean toxicity probability `pT` and \eqn{\Pr(\pi_{T,j} < \pi_{T,0} + \delta_T)}
#'       (`pr_safe`).}
#'     \item{`pr_eff`, `pR`}{the same quantities as `K` x `(J + 1)` matrices.}
#'     \item{`pr_safe`, `pT`}{vectors over the treatments.}
#'     \item{`partition`}{posterior probabilities of the partitions of the subgroups induced by
#'       the clustering of each treatment's response mean (a pattern such as `"112"` puts
#'       \eqn{B_1} and \eqn{B_2} in one cluster).}
#'   }
#' @examples
#' \donttest{
#' design <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)), N = 150, lambda = 0.70,
#'                       mcmc = list(iter = 1000, burn = 500))
#' truth <- boat_scenario(pR = rbind(c(0.50, 0.30, NA, 0.70),
#'                                   c(0.50, 0.70, 0.50, 0.30),
#'                                   c(0.30,   NA, 0.50, 0.30)),
#'                        piT = c(0.25, 0.35, 0.20, 0.15))
#' dat <- boat_generate_data(design, truth, n = 20, seed = 1)
#' fit <- boat_fit(design, dat$pat, dat$pdo, seed = 1)
#' fit
#' }
#' @export
boat_fit <- function(design, pat, pdo, seed = NULL) {
  .check_design(design)
  d <- .internal_data(design, pat, pdo)
  postE <- .fit_efficacy(design, d$pat, d$pdo, seed)
  tox <- .fit_toxicity(design, d$pat)
  K <- design$K; J <- design$J
  rows <- do.call(rbind, lapply(seq_len(K), function(k) {
    j <- design$arms[[k]]
    data.frame(k = k, j = j, n = sapply(j, function(jj) sum(d$pat$k == k & d$pat$j == jj)),
               pR = postE$post_pE_mean[k, j + 1], pr_eff = ifelse(j == 0, NA, postE$post_prob_eff[k, j + 1]),
               pT = tox$post_piT_mean[j + 1], pr_safe = tox$post_prob_safe[j + 1])
  }))
  dn <- list(paste0("k", seq_len(K)), paste0("T", 0:J))
  pr_eff <- postE$post_prob_eff; pr_eff[, 1] <- NA; dimnames(pr_eff) <- dn
  for (k in seq_len(K)) pr_eff[k, setdiff(0:J, design$arms[[k]]) + 1] <- NA   # treatments not evaluated
  pR <- postE$post_pE_mean; dimnames(pR) <- dn
  structure(list(summary = rows, pr_eff = pr_eff, pR = pR,
                 pr_safe = setNames(tox$post_prob_safe, paste0("T", 0:J)),
                 pT = setNames(tox$post_piT_mean, paste0("T", 0:J)),
                 partition = postE$post_part_prob, n = nrow(d$pat)),
            class = "boat_fit")
}

#' Interim and final decisions
#'
#' Applies the decision rule of [boat_boundary()] to a model fit. At an interim analysis a
#' treatment is dropped from a subgroup if its posterior probability of being effective or of
#' being safe is below the cutoff, and a subgroup with no experimental treatment left stops
#' accrual. At the final analysis the treatments still open are recommended if both
#' probabilities reach the cutoff.
#'
#' @param design A `boat_design` object.
#' @param fit A `boat_fit` object from [boat_fit()].
#' @param n Current number of patients in each subgroup (vector of length `K`).
#' @param active List of the treatments still open in each subgroup (control included).
#'   `NULL` means all treatments of the design.
#' @param final `TRUE` at the final analysis.
#' @return An object of class `boat_decision` with elements `table` (one row per open
#'   experimental treatment with the two probabilities, the cutoff and the decision), `active`
#'   (open treatments after the analysis), `stopped` (subgroups that stopped accrual) and, at the
#'   final analysis, `recommended` (list of recommended treatments per subgroup).
#' @examples
#' \donttest{
#' design <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)), N = 150, lambda = 0.70,
#'                       mcmc = list(iter = 1000, burn = 500))
#' truth <- boat_scenario(pR = rbind(c(0.50, 0.30, NA, 0.70),
#'                                   c(0.50, 0.70, 0.50, 0.30),
#'                                   c(0.30,   NA, 0.50, 0.30)),
#'                        piT = c(0.25, 0.35, 0.20, 0.15))
#' dat <- boat_generate_data(design, truth, n = 75, seed = 1)
#' fit <- boat_fit(design, dat$pat, dat$pdo, seed = 1)
#' boat_decide(design, fit, n = c(75, 75, 75))
#' }
#' @export
boat_decide <- function(design, fit, n, active = NULL, final = FALSE) {
  .check_design(design)
  if (!inherits(fit, "boat_fit")) stop("fit must be created by boat_fit()")
  K <- design$K
  if (length(n) != K) stop("n must give the number of patients of each subgroup")
  if (is.null(active)) active <- design$arms
  active_list <- setNames(lapply(active, function(a) a), as.character(seq_len(K)))
  active_groups <- as.integer(which(sapply(active_list, function(a) any(a != 0))))
  stopped_before <- setdiff(seq_len(K), active_groups)
  tab <- .decision_table(design, fit$pr_eff, fit$pr_safe, active_list, active_groups, n)
  if (!final) {
    res <- .interim_decide(fit$pr_eff, fit$pr_safe, active_list, active_groups, as.integer(n), design$N,
                           design$lambda, stopped_before)
    tab$decision <- ifelse(is.na(tab$pr_eff) | is.na(tab$pr_safe), "continue (no estimate)",
                           ifelse(tab$pr_eff < tab$cutoff | tab$pr_safe < tab$cutoff, "drop", "continue"))
    out <- list(table = tab, active = unname(res$active_list), stopped = sort(res$stopped_groups))
  } else {
    sel <- .final_select(fit$pr_eff, fit$pr_safe, active_list, stopped_before, as.integer(n), design$N,
                         design$lambda, K)
    tab$decision <- ifelse(!is.na(tab$pr_eff) & !is.na(tab$pr_safe) & tab$pr_eff >= tab$cutoff &
                             tab$pr_safe >= tab$cutoff, "recommend", "not recommended")
    out <- list(table = tab, active = unname(active_list), stopped = sort(stopped_before),
                recommended = unname(sel))
  }
  out$final <- final
  structure(out, class = "boat_decision")
}

#' Fit the randomization model
#'
#' For each open experimental treatment and the control in each active subgroup, fits a probit
#' regression of the patients' clinical responses on their mean PDO responses to that treatment
#' (Albert-Chib Gibbs sampler). [boat_randomize()] uses these fits to randomize the patients
#' enrolled until the next analysis. Only used by the adaptive design.
#'
#' @inheritParams boat_fit
#' @param active List of the treatments still open in each subgroup, usually the `active`
#'   element of [boat_decide()]. `NULL` means all treatments of the design.
#' @param seed Seed set before sampling, or `NULL` to use the current random number stream.
#' @return An object of class `boat_rar`.
#' @examples
#' \donttest{
#' design <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)), N = 150, lambda = 0.70,
#'                       rar_mcmc = list(iter = 400, burn = 200))
#' truth <- boat_scenario(pR = rbind(c(0.50, 0.30, NA, 0.70),
#'                                   c(0.50, 0.70, 0.50, 0.30),
#'                                   c(0.30,   NA, 0.50, 0.30)),
#'                        piT = c(0.25, 0.35, 0.20, 0.15))
#' dat <- boat_generate_data(design, truth, n = 30, seed = 1)
#' rar <- boat_rar_model(design, dat$pat, dat$pdo, seed = 2)
#' # a new patient in subgroup 1 whose PDOs responded best to A3
#' boat_randomize(design, k = 1, pdo_mean = c("0" = -0.2, "1" = 0.1, "3" = 0.8), rar = rar)
#' }
#' @export
boat_rar_model <- function(design, pat, pdo, active = NULL, seed = NULL) {
  .check_design(design)
  d <- .internal_data(design, pat, pdo)
  if (is.null(active)) active <- design$arms
  active_list <- setNames(lapply(active, function(a) a), as.character(seq_len(design$K)))
  active_groups <- as.integer(which(sapply(active_list, function(a) any(a != 0))))
  if (!is.null(seed)) set.seed(seed)
  structure(list(draws = .fit_rar(design, d$pat, d$pdo, active_groups, active_list), active = unname(active_list)),
            class = "boat_rar")
}

#' Randomize a new patient
#'
#' Gives the randomization probabilities of a newly enrolled patient in subgroup `k` and draws
#' the patient's treatment. Under adaptive randomization, the probability of treatment
#' \eqn{A_j} is proportional to the posterior probability that \eqn{A_j} gives the highest
#' response probability for this patient, predicted from the patient's own mean PDO response to
#' each open treatment. The probabilities are then truncated to `[trunc, 1 - trunc]` and
#' renormalized (so they can end slightly below `trunc`). Before the first interim analysis
#' (`rar = NULL`) and under equal randomization all open treatments have the same probability.
#'
#' @param design A `boat_design` object.
#' @param k Subgroup of the new patient.
#' @param pdo_mean The patient's mean PDO response to each treatment open in subgroup `k`: a
#'   vector named by treatment (for example `c("0" = 0.1, "1" = 0.4, "3" = 0.2)`). Not needed
#'   before the first interim analysis or under equal randomization.
#' @param rar A `boat_rar` object from [boat_rar_model()], or `NULL` before the first interim
#'   analysis.
#' @param active List of the treatments open in each subgroup. Defaults to `rar$active` or to
#'   all treatments of the design.
#' @return A list with the randomization probabilities `prob` (named by treatment) and the
#'   assigned treatment `treatment`.
#' @examples
#' design <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)))
#' # before the first interim analysis all open treatments are equally likely
#' boat_randomize(design, k = 2)
#' @export
boat_randomize <- function(design, k, pdo_mean = NULL, rar = NULL, active = NULL) {
  .check_design(design)
  if (is.null(active)) active <- if (!is.null(rar)) rar$active else design$arms
  choices <- active[[k]]
  if (design$randomization == "equal") {
    prob <- setNames(rep(1 / length(choices), length(choices)), as.character(choices))
    return(list(prob = prob, treatment = .assign_equal(choices)))
  }
  if (!is.null(rar) && is.null(pdo_mean)) stop("pdo_mean is needed for adaptive randomization")
  barYO <- setNames(rep(NA_real_, length(choices)), as.character(choices))
  if (!is.null(pdo_mean)) barYO[intersect(names(pdo_mean), names(barYO))] <- pdo_mean[intersect(names(pdo_mean), names(barYO))]
  prob <- .alloc_probs(design, if (is.null(rar)) NULL else rar$draws, k, barYO, choices)
  list(prob = prob, treatment = .assign_adaptive(prob, choices))
}

#' @export
print.boat_fit <- function(x, digits = 3, ...) {
  cat(sprintf("BOAT model fit to %d patients\n", x$n))
  s <- x$summary
  s$treatment <- paste0("A", s$j); s$subgroup <- paste0("B", s$k)
  print(format(s[, c("subgroup", "treatment", "n", "pR", "pr_eff", "pT", "pr_safe")], digits = digits),
        row.names = FALSE)
  cat("pR, pT: posterior mean response and toxicity probabilities\n")
  cat("pr_eff: Pr(pi_R,j,k > pi_R,0,k + delta_R | data); pr_safe: Pr(pi_T,j < pi_T,0 + delta_T | data)\n")
  invisible(x)
}

#' @export
print.boat_decision <- function(x, digits = 3, ...) {
  cat(if (x$final) "Final analysis\n" else "Interim analysis\n")
  t <- x$table
  if (nrow(t)) {
    t$subgroup <- paste0("B", t$k); t$treatment <- paste0("A", t$j)
    print(format(t[, c("subgroup", "treatment", "n", "pr_eff", "pr_safe", "cutoff", "decision")], digits = digits),
          row.names = FALSE)
  }
  if (x$final) {
    for (k in seq_along(x$recommended))
      cat(sprintf("  B%d: %s\n", k, if (length(x$recommended[[k]])) paste0("recommend ", paste0("A", x$recommended[[k]], collapse = ", "))
                  else "no experimental treatment recommended"))
  } else {
    for (k in seq_along(x$active))
      cat(sprintf("  B%d: %s\n", k, if (k %in% x$stopped) "accrual stopped"
                  else paste("open:", paste0("A", x$active[[k]], collapse = " "))))
  }
  invisible(x)
}

#' @export
print.boat_rar <- function(x, ...) {
  n <- sum(sapply(x$draws, length))
  cat(sprintf("BOAT randomization model: probit fits for %d (subgroup, treatment) pairs\n", n))
  invisible(x)
}

# ---- internal steps shared with the trial simulation ----

# user data -> the internal format (patients numbered 1..N in row order; YE = response)
.internal_data <- function(design, pat, pdo) {
  need <- c("id", "k", "j", "YT", "YR")
  if (!all(need %in% names(pat))) stop("pat must have columns ", paste(need, collapse = ", "))
  if (!all(c("id", "k", "j", "YO") %in% names(pdo))) stop("pdo must have columns id, k, j, YO")
  # subgroups and treatments may come as factors or characters; convert them by their labels
  lab <- function(x) suppressWarnings(as.integer(as.character(x)))
  pat$k <- lab(pat$k); pat$j <- lab(pat$j); pdo$k <- lab(pdo$k); pdo$j <- lab(pdo$j)
  if (anyDuplicated(pat$id)) stop("patient ids must be unique")
  if (any(!pat$k %in% seq_len(design$K))) stop("pat$k must lie in 1:K")
  ok <- mapply(function(k, j) j %in% design$arms[[k]], pat$k, pat$j)
  if (!all(ok)) stop("some patients were assigned a treatment not evaluated in their subgroup")
  if (any(!pat$YT %in% 0:1) || any(!pat$YR %in% 0:1)) stop("YT and YR must be 0 or 1")
  idx <- match(pdo$id, pat$id)
  if (anyNA(idx)) stop("every PDO id must be a patient id")
  if (any(is.na(pdo$k) | pdo$k != pat$k[idx])) stop("PDO subgroups differ from the patients' subgroups")
  list(pat = data.frame(id = seq_len(nrow(pat)), k = as.integer(pat$k), j = pat$j, YT = pat$YT, YE = pat$YR),
       pdo = data.frame(id = idx, k = as.integer(pdo$k), j = pdo$j, YO = pdo$YO))
}

.fit_efficacy <- function(design, df_pat, df_pdo, seed) {
  pr <- design$prior
  args <- list(df_pat = df_pat, df_pdo = df_pdo, K = design$K, J = design$J, Sk_list = design$arms,
               n_iter = design$mcmc$iter, n_burn = design$mcmc$burn, n_thin = design$mcmc$thin,
               delta_E = design$delta_R, prior_muE_mean = pr$mu_mean, prior_muE_sd = pr$mu_sd,
               prior_sigma_b_alpha = pr$var_a, prior_sigma_b_beta = pr$var_b,
               prior_sigma_eps_alpha = pr$var_a, prior_sigma_eps_beta = pr$var_b,
               prior_sigma_mu_alpha = pr$var_a, prior_sigma_mu_beta = pr$var_b,
               prior_rho_mean = pr$rho_mean, prior_rho_sd = pr$rho_sd, verbose = FALSE)
  if (is.null(seed)) {
    # the sampler always calls set.seed(); continue the current stream instead
    seed <- sample.int(.Machine$integer.max, 1)
  }
  args$seed <- seed
  do.call(.boat_mcmc_efficacy, args)
}

.fit_toxicity <- function(design, df_pat) {
  .boat_tox_posterior(df_pat = df_pat, J = design$J, prior_piT_a = design$prior$tox_a,
                      prior_piT_b = design$prior$tox_b, delta_T = design$delta_T)
}

# the probit fits order the patients with merge(), whose order follows the locale's collation;
# the fits run in the C locale so that the results do not depend on the session's locale
.fit_rar <- function(design, df_pat, df_pdo, active_groups, active_list) {
  old <- Sys.getlocale("LC_COLLATE")
  if (old != "C") {
    invisible(Sys.setlocale("LC_COLLATE", "C"))
    on.exit(invisible(Sys.setlocale("LC_COLLATE", old)), add = TRUE)
  }
  .boat_ac_probit(df_pat = df_pat, df_pdo = df_pdo, K = design$K, active_groups = active_groups,
                  active_list = active_list, n_iter_pp = design$rar_mcmc$iter, n_burn_pp = design$rar_mcmc$burn,
                  n_thin_pp = design$rar_mcmc$thin, prior_mean = design$prior$rar_mean,
                  prior_sd = design$prior$rar_sd, standardize_x = design$rar_standardize, verbose = FALSE)
}

# interim decisions for the active subgroups (drop if ineffective or unsafe)
.interim_decide <- function(prob_eff, prob_safe, active_list, active_groups, n_by_k, max_n_k, lambda,
                            stopped_groups = integer(0)) {
  for (k in as.integer(active_groups)) {
    k_chr <- as.character(k)
    drop_list <- c()
    nk_now <- n_by_k[k]
    PL <- .cutoff(nk_now, max_n_k[k], lambda)
    for (j in active_list[[k_chr]]) {
      if (j == 0) next
      post_eff <- prob_eff[k, j + 1]
      post_safe <- prob_safe[j + 1]
      if (is.na(post_eff) || is.na(post_safe)) next
      if (post_eff < PL || post_safe < PL) drop_list <- c(drop_list, j)
    }
    if (length(drop_list) > 0) active_list[[k_chr]] <- setdiff(active_list[[k_chr]], drop_list)
    if (!(0 %in% active_list[[k_chr]])) active_list[[k_chr]] <- union(active_list[[k_chr]], 0)
    remaining_exp <- setdiff(active_list[[k_chr]], 0)
    if (length(remaining_exp) == 0) {
      active_groups <- as.integer(setdiff(active_groups, k))
      stopped_groups <- union(stopped_groups, k)
    }
  }
  list(active_list = active_list, active_groups = active_groups, stopped_groups = stopped_groups)
}

# final recommendations: open treatments that are effective and safe
.final_select <- function(prob_eff, prob_safe, active_list, stopped_groups, n_by_k, max_n_k, lambda, K) {
  final_selection <- setNames(rep(list(integer(0)), K), as.character(seq_len(K)))
  for (k in setdiff(seq_len(K), stopped_groups)) {
    k_chr <- as.character(k)
    open_exp <- setdiff(active_list[[k_chr]], 0)
    if (length(open_exp) == 0) next
    PL <- .cutoff(n_by_k[k], max_n_k[k], lambda)
    for (j in open_exp) {
      pe <- prob_eff[k, j + 1]
      ps <- prob_safe[j + 1]
      if (is.na(pe) || is.na(ps)) next
      if (pe >= PL && ps >= PL) final_selection[[k_chr]] <- union(final_selection[[k_chr]], j)
    }
  }
  for (k in seq_len(K)) {
    k_chr <- as.character(k)
    final_selection[[k_chr]] <- if (length(final_selection[[k_chr]]) > 0) sort(unique(final_selection[[k_chr]])) else integer(0)
  }
  final_selection
}

.decision_table <- function(design, prob_eff, prob_safe, active_list, active_groups, n) {
  rows <- lapply(active_groups, function(k) {
    j <- setdiff(active_list[[as.character(k)]], 0)
    if (!length(j)) return(NULL)
    data.frame(k = k, j = j, n = n[k], pr_eff = prob_eff[k, j + 1], pr_safe = prob_safe[j + 1],
               cutoff = .cutoff(n[k], design$N[k], design$lambda))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame(k = integer(0), j = integer(0), n = integer(0), pr_eff = numeric(0),
                                      pr_safe = numeric(0), cutoff = numeric(0))
  out
}

# randomization probabilities of a new patient in subgroup k (adaptive design)
.alloc_probs <- function(design, alpha_draws_current, k_i, barYO_map, choices) {
  m_active <- length(choices)
  if (m_active == 1) return(setNames(1, as.character(choices)))
  if (is.null(alpha_draws_current)) {
    probs <- rep(1 / m_active, m_active)
    names(probs) <- as.character(choices)
    return(probs)
  }
  k_chr <- as.character(k_i)
  arms_all <- choices
  B_eff_list <- integer(0)
  for (j in arms_all) {
    adkj <- alpha_draws_current[[k_chr]][[as.character(j)]]
    if (!is.null(adkj)) B_eff_list <- c(B_eff_list, min(length(adkj$a0), length(adkj$a1)))
  }
  B_eff <- if (length(B_eff_list) == 0) 0 else min(B_eff_list)
  if (!is.finite(B_eff) || B_eff <= 0) {
    probs <- rep(1 / length(arms_all), length(arms_all))
    names(probs) <- as.character(arms_all)
    return(probs)
  }
  L <- length(arms_all)
  Pmat <- matrix(NA_real_, nrow = B_eff, ncol = L)
  colnames(Pmat) <- as.character(arms_all)
  for (cc in seq_along(arms_all)) {
    j <- arms_all[cc]
    j_chr <- as.character(j)
    xj_raw <- barYO_map[j_chr]
    adkj <- alpha_draws_current[[k_chr]][[j_chr]]
    if (is.null(adkj) || is.na(xj_raw) ||
        (design$rar_standardize && (!is.finite(adkj$x_sd) || adkj$x_sd < 1e-8))) {
      Pmat[, cc] <- 0.5
    } else {
      xj_std <- if (design$rar_standardize) (xj_raw - adkj$x_mean) / adkj$x_sd else xj_raw
      Pmat[, cc] <- pnorm(adkj$a0[1:B_eff] + adkj$a1[1:B_eff] * xj_std)
    }
  }
  rj <- .rmax_prob(Pmat)
  names(rj) <- as.character(arms_all)
  p_raw <- (pmax(rj, 1e-12))^design$gamma
  p_raw <- p_raw / sum(p_raw)
  probs <- .truncate_probs(p_raw, design$trunc)
  names(probs) <- as.character(arms_all)
  probs
}

# posterior probability that each column is the largest (ties broken at random)
.rmax_prob <- function(Pmat) {
  B <- nrow(Pmat); L <- ncol(Pmat)
  win <- integer(L)
  for (b in 1:B) {
    row <- Pmat[b, ]
    mx <- max(row)
    idx <- which(row == mx)
    jwin <- if (length(idx) == 1) idx else sample(idx, 1)
    win[jwin] <- win[jwin] + 1L
  }
  win / B
}

.truncate_probs <- function(p, eps) {
  K <- length(p)
  if (eps < 0 || eps >= 1 / K) stop("trunc_eps must satisfy 0 <= trunc_eps < 1/K.")
  p <- p / sum(p)
  p <- pmin(pmax(p, eps), 1 - eps)
  p / sum(p)
}

.assign_adaptive <- function(probs, choices) {
  if (length(choices) == 1) return(choices[1])
  sample(as.integer(names(probs)), 1, prob = probs)
}

.assign_equal <- function(choices) {
  if (length(choices) == 1) return(choices[1])
  sample(choices, 1)
}
