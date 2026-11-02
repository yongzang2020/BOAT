#' Specify a simulation scenario
#'
#' @param pR Matrix of true response probabilities with one row per subgroup and one column per
#'   treatment (first column = control). `NA` marks a treatment not evaluated in a subgroup.
#' @param piT Vector of true toxicity probabilities of the treatments (first = control).
#' @return An object of class `boat_scenario`.
#' @examples
#' boat_scenario(pR = rbind(c(0.5, 0.3, NA, 0.7), c(0.5, 0.7, 0.5, 0.3), c(0.3, NA, 0.5, 0.3)),
#'               piT = c(0.25, 0.35, 0.20, 0.15))
#' @export
boat_scenario <- function(pR, piT) {
  pR <- as.matrix(pR)
  if (length(piT) != ncol(pR)) stop("piT must have one value per column of pR")
  if (any(is.na(pR[, 1]))) stop("the control (first column of pR) must be evaluated in every subgroup")
  if (any(pR < 0 | pR > 1, na.rm = TRUE) || any(piT < 0 | piT > 1)) stop("probabilities must lie in [0, 1]")
  dimnames(pR) <- list(paste0("k", seq_len(nrow(pR))), paste0("T", seq_len(ncol(pR)) - 1))
  structure(list(pR = pR, piT = as.numeric(piT)), class = "boat_scenario")
}

#' Data-generating parameters of a simulation
#'
#' @param sigma_E Standard deviation of the patient random effect shared by the patient's
#'   clinical response and PDO responses.
#' @param rho Association between the patient random effect and the PDO responses.
#' @param sigma_O Residual standard deviation of the PDO responses.
#' @param sigma_OR Standard deviation linking the PDO response mean to the clinical response
#'   mean on the probit scale.
#' @param M Number of PDO replicates per patient and treatment.
#' @return A list of class `boat_dgm`.
#' @examples
#' boat_dgm(rho = 0.3)
#' @export
boat_dgm <- function(sigma_E = 1, rho = 0.6, sigma_O = 0.1, sigma_OR = 0.1, M = 20) {
  stopifnot(sigma_E > 0, sigma_O > 0, sigma_OR > 0, M >= 1)
  structure(list(sigma_E = sigma_E, rho = rho, sigma_O = sigma_O, sigma_OR = sigma_OR, M = as.integer(M)),
            class = "boat_dgm")
}

#' @export
print.boat_scenario <- function(x, ...) {
  tab <- cbind(x$pR)
  tab <- rbind(piT = x$piT, tab)
  cat("True toxicity (piT) and response probabilities (rows k1, k2, ...):\n")
  print(tab, na.print = "-")
  invisible(x)
}

# subgroup k's admissible treatments under a scenario
.admissible <- function(scenario, k, delta_R, delta_T) {
  pR <- scenario$pR; piT <- scenario$piT
  J <- ncol(pR) - 1
  which(sapply(seq_len(J), function(j) !is.na(pR[k, j + 1]) && pR[k, j + 1] > pR[k, 1] + delta_R &&
                 piT[j + 1] < piT[1] + delta_T))
}

.check_scenario <- function(scenario, design) {
  if (!inherits(scenario, "boat_scenario")) stop("scenario must be created by boat_scenario()")
  if (nrow(scenario$pR) != design$K || ncol(scenario$pR) != design$J + 1)
    stop("scenario does not have the design's numbers of subgroups and treatments")
  for (k in seq_len(design$K))
    if (!identical(as.integer(which(!is.na(scenario$pR[k, ])) - 1L), design$arms[[k]]))
      stop("the treatments evaluated in subgroup ", k, " (non-NA entries of pR) differ from design$arms")
  invisible(scenario)
}
