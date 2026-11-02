#' Specify a BOAT design
#'
#' Collects the design parameters of a BOAT trial.
#'
#' @param K Number of subgroups.
#' @param J Number of experimental treatments. Treatment 0 is the control.
#' @param arms List of length `K`. Element `k` gives the treatments evaluated in subgroup `k`; it
#'   must contain the control 0. `NULL` evaluates all treatments `0:J` in every subgroup.
#' @param N Maximum sample size of each subgroup, a single number or a vector of length `K`.
#' @param interim Fractions of `N` at which interim analyses are done.
#' @param lambda Tuning parameter of the decision cutoffs (see [boat_boundary()]).
#' @param delta_T Toxicity margin \eqn{\delta_T}.
#' @param delta_R Minimum clinically meaningful response improvement \eqn{\delta_R}.
#' @param randomization `"adaptive"` for PDO-guided response-adaptive randomization or
#'   `"equal"` for equal randomization among the open treatments.
#' @param gamma Power applied to the randomization probabilities.
#' @param trunc Truncation of the randomization probabilities: they are truncated to
#'   `[trunc, 1 - trunc]` and renormalized. Must be below 1 over the number of open treatments.
#' @param mcmc Iterations, burn-in and thinning of the RJMCMC sampler of the response model.
#' @param rar_mcmc Iterations, burn-in and thinning of the Gibbs sampler of the probit models
#'   used for randomization.
#' @param prior List of hyperparameters: `mu_mean` and `mu_sd` (normal prior of the probit
#'   response means), `var_a` and `var_b` (inverse-gamma priors of the variances), `rho_mean`
#'   and `rho_sd` (truncated normal prior of \eqn{\rho > 0}), `tox_a` and `tox_b` (beta priors
#'   of the toxicity probabilities), `rar_mean` and `rar_sd` (normal priors of the intercept and
#'   slope of the randomization probit models).
#' @param rar_standardize Standardize the PDO responses in the randomization probit models.
#'
#' @return An object of class `boat_design`.
#' @seealso [boat_boundary()]
#' @examples
#' # three subgroups and three experimental treatments; A2 is not evaluated in B1 and A1 not in B3
#' design <- boat_design(K = 3, J = 3, arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)), N = 150,
#'                       interim = c(0.5, 0.75), lambda = 0.70)
#' design
#' @export
boat_design <- function(K = 3, J = 3, arms = NULL, N = 150, interim = c(0.5, 0.75),
                        lambda = 0.70, delta_T = 0.05, delta_R = 0.05,
                        randomization = c("adaptive", "equal"), gamma = 1, trunc = 0.10,
                        mcmc = list(iter = 10000, burn = 5000, thin = 1),
                        rar_mcmc = list(iter = 2000, burn = 1000, thin = 1),
                        prior = list(), rar_standardize = TRUE) {
  randomization <- match.arg(randomization)
  stopifnot(K >= 1, J >= 1, K == round(K), J == round(J))
  if (is.null(arms)) arms <- rep(list(0:J), K)
  if (!is.list(arms) || length(arms) != K) stop("arms must be a list with one element per subgroup")
  arms <- lapply(arms, function(a) sort(unique(as.integer(a))))
  for (k in seq_len(K)) {
    if (!(0L %in% arms[[k]])) stop("arms[[", k, "]] must contain the control 0")
    if (any(arms[[k]] < 0 | arms[[k]] > J)) stop("arms[[", k, "]] must be within 0:J")
    if (length(arms[[k]]) < 2) stop("arms[[", k, "]] must contain at least one experimental treatment")
  }
  if (length(N) == 1) N <- rep(N, K)
  if (length(N) != K || any(N < 1)) stop("N must be a positive number or a vector of length K")
  interim <- sort(unique(interim))
  if (any(interim <= 0 | interim >= 1)) stop("interim fractions must lie in (0, 1)")
  if (lambda <= 0 || lambda > 1) stop("lambda must lie in (0, 1]")
  if (trunc < 0 || trunc >= 1 / max(lengths(arms))) stop("trunc must be below 1 over the number of treatments in a subgroup")
  mcmc <- .fill(mcmc, list(iter = 10000, burn = 5000, thin = 1))
  rar_mcmc <- .fill(rar_mcmc, list(iter = 2000, burn = 1000, thin = 1))
  if (mcmc$burn >= mcmc$iter || rar_mcmc$burn >= rar_mcmc$iter) stop("burn-in must be smaller than the number of iterations")
  prior <- .fill(prior, list(mu_mean = 0, mu_sd = 10, var_a = 0.01, var_b = 0.01, rho_mean = 0, rho_sd = 10,
                             tox_a = 1, tox_b = 1, rar_mean = c(0, 0), rar_sd = c(100, 100)))
  structure(list(K = as.integer(K), J = as.integer(J), arms = arms, N = as.integer(N), interim = interim,
                 lambda = lambda, delta_T = delta_T, delta_R = delta_R, randomization = randomization,
                 gamma = gamma, trunc = trunc, mcmc = mcmc, rar_mcmc = rar_mcmc, prior = prior,
                 rar_standardize = rar_standardize),
            class = "boat_design")
}

#' Decision cutoffs of the interim and final analyses
#'
#' At an analysis with \eqn{n_k} patients in subgroup \eqn{B_k}, treatment \eqn{A_j} is dropped
#' from \eqn{B_k} if
#' \deqn{\Pr(\pi_{R,j,k} > \pi_{R,0,k} + \delta_R \mid \mathcal{D}) < c_k \quad \text{or} \quad
#' \Pr(\pi_{T,j} < \pi_{T,0} + \delta_T \mid \mathcal{D}) < c_k,
#' \qquad c_k = \lambda \sqrt{n_k / N_k}.}
#' At the final analysis (\eqn{n_k = N_k}, so \eqn{c_k = \lambda}) a treatment still open is
#' recommended if both probabilities are at least \eqn{c_k}.
#'
#' @param design A `boat_design` object.
#' @return A data frame with one row per analysis and subgroup: the number of patients and the
#'   cutoff.
#' @examples
#' design <- boat_design(arms = list(c(0, 1, 3), 0:3, c(0, 2, 3)), N = 150, lambda = 0.70)
#' boat_boundary(design)
#' @export
boat_boundary <- function(design) {
  .check_design(design)
  tg <- .interim_targets(design)
  out <- do.call(rbind, lapply(seq_len(ncol(tg) + 1), function(a) {
    n <- if (a <= ncol(tg)) tg[, a] else design$N
    data.frame(analysis = if (a <= ncol(tg)) paste("Interim", a) else "Final", subgroup = seq_len(design$K),
               n = as.integer(n), cutoff = .cutoff(n, design$N, design$lambda))
  }))
  rownames(out) <- NULL
  out
}

#' @export
print.boat_design <- function(x, ...) {
  cat(if (x$randomization == "adaptive") "BOAT design" else "BOAT design with equal randomization",
      sprintf("with %d subgroups and %d experimental treatments\n", x$K, x$J))
  for (k in seq_len(x$K))
    cat(sprintf("  B%d: treatments %s, at most %d patients\n", k, paste0("A", x$arms[[k]], collapse = " "), x$N[k]))
  cat(sprintf("  interim analyses after %s of the patients\n", paste0(100 * x$interim, "%", collapse = " and ")))
  cat(sprintf("  lambda = %.2f, delta_T = %.2f, delta_R = %.2f\n", x$lambda, x$delta_T, x$delta_R))
  if (x$randomization == "adaptive")
    cat(sprintf("  adaptive randomization: gamma = %g, truncation at %g\n", x$gamma, x$trunc))
  cat(sprintf("  MCMC: %d iterations (burn-in %d); randomization model: %d iterations (burn-in %d)\n",
              x$mcmc$iter, x$mcmc$burn, x$rar_mcmc$iter, x$rar_mcmc$burn))
  invisible(x)
}

# ---- internal helpers ----
.fill <- function(x, default) {
  bad <- setdiff(names(x), names(default))
  if (length(bad)) stop("unknown element(s): ", paste(bad, collapse = ", "))
  default[names(x)] <- x
  default
}

.check_design <- function(design) {
  if (!inherits(design, "boat_design")) stop("design must be created by boat_design()")
  invisible(design)
}

.cutoff <- function(n, N, lambda) lambda * sqrt(n / N)

# interim sample sizes (K x number of interims)
.interim_targets <- function(design) {
  K <- design$K; max_n_k <- design$N
  tg <- sapply(design$interim, function(fr) pmin(as.integer(ceiling(fr * max_n_k)), max_n_k))
  tg <- matrix(tg, nrow = K)
  colnames(tg) <- paste0("I", seq_len(ncol(tg)))
  rownames(tg) <- paste0("k", seq_len(K))
  tg
}
