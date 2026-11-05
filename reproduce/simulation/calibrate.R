# ============================================================
# Calibration of lambda (Section 3; values in Sections 4.1 and 5). For every lambda on a grid, the
# design is simulated under Scenarios 1-3, which contain the null subgroups (no admissible
# treatment): B2 in Scenario 1, B1-B3 in Scenario 2 and B2-B3 in Scenario 3. The chosen lambda is
# the grid value whose largest FWER_k over these subgroups is closest to ALPHA = 0.10. Ties go to
# the smallest value with FWER_k at most ALPHA. The lambda of the illustrative trial (Section 5) is
# calibrated in the same way under its null setting (Supplementary Figure S.8), with nominal FWER
# ALPHA_ILLUSTRATION = 0.05.
# This is Step D3 of Supplementary Figure S.7 for the illustrative trial.
#
# BOAT and BOAT-ER are calibrated with boat_calibrate() of the BOAT package and the comparators
# with comparator_calibrate() (functions.R), which applies the same rule.
#
# Usage (from this folder), with <name> one of the calibrations listed below:
#   Rscript calibrate.R <name> [cores]   e.g. Rscript calibrate.R BOAT 16
# The result is printed and saved in results/calibration/<name>.rds and <name>.csv (the largest
# FWER_k for every lambda). Each lambda uses 1000 replications with seeds 1001-2000.
# ============================================================
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.f) == 1) dirname(normalizePath(.f)) else getwd(), "functions.R"))

# design, data-generating values, true probabilities, delta_R, lambda grid and nominal FWER of each
# calibration
calibration <- function(design, grid, dgm = DGM_MAIN, truths = "main", delta_R = 0.05, alpha = ALPHA) {
  list(design = design, grid = grid, dgm = dgm, truths = truths, delta_R = delta_R, alpha = alpha)
}
CALIBRATIONS <- list(
  # main simulation (Supplementary Table S.1)
  "CD-SEP"               = calibration("CD-SEP",  seq(0.66, 0.86, by = 0.02)),
  "CD-BHM"               = calibration("CD-BHM",  seq(0.68, 0.88, by = 0.02)),
  "CD-AB"                = calibration("CD-AB",   seq(0.70, 0.90, by = 0.02)),
  "BOAT-ER"              = calibration("BOAT-ER", seq(0.60, 0.76, by = 0.02)),
  "BOAT"                 = calibration("BOAT",    seq(0.62, 0.78, by = 0.02)),
  # sensitivity analyses with a recalibrated lambda (Supplementary Tables S.4-S.6)
  "BOAT_sigmaOR0.05"     = calibration("BOAT", seq(0.62, 0.72, by = 0.01), dgm = boat_dgm(sigma_OR = 0.05)),
  "BOAT_sigmaOR0.20"     = calibration("BOAT", seq(0.71, 0.81, by = 0.01), dgm = boat_dgm(sigma_OR = 0.20)),
  "BOAT_deltaR0.20"      = calibration("BOAT", seq(0.66, 0.76, by = 0.01), truths = "deltaR020", delta_R = 0.20),
  "BOAT_rho0.90"         = calibration("BOAT", seq(0.66, 0.76, by = 0.01), dgm = boat_dgm(rho = 0.9)),
  "BOAT_rho0.00"         = calibration("BOAT", seq(0.66, 0.76, by = 0.01), dgm = boat_dgm(rho = 0)),
  "BOAT_rho-0.60"        = calibration("BOAT", seq(0.76, 0.88, by = 0.01), dgm = boat_dgm(rho = -0.6)),
  # illustrative trial (Section 5)
  "illustration_BOAT-ER" = calibration("BOAT-ER", seq(0.60, 0.76, by = 0.02), truths = "illustration_null",
                                       alpha = ALPHA_ILLUSTRATION),
  "illustration_BOAT"    = calibration("BOAT",    seq(0.64, 0.80, by = 0.02), truths = "illustration_null",
                                       alpha = ALPHA_ILLUSTRATION)
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !(args[1] %in% names(CALIBRATIONS)))
  stop("give one of the calibrations: ", paste(names(CALIBRATIONS), collapse = ", "))
name <- args[1]
cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
cal <- CALIBRATIONS[[name]]

scenarios <- lapply(if (cal$truths == "illustration_null") 1 else 1:3, function(s) scenario_truth(s, cal$truths))
st <- setting(name, cal$design, cal$grid[1], dgm = cal$dgm, truths = cal$truths, delta_R = cal$delta_R)
design <- make_design(st, scenarios[[1]])
res <- if (cal$design %in% c("BOAT", "BOAT-ER")) {
  boat_calibrate(design, scenarios, lambda = cal$grid, nsim = NSIM, alpha = cal$alpha, dgm = cal$dgm, cores = cores)
} else {
  comparator_calibrate(cal$design, design, scenarios, lambda = cal$grid, nsim = NSIM, alpha = cal$alpha, cores = cores)
}

out_dir <- file.path(CODE_DIR, "results", "calibration")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(res, file.path(out_dir, paste0(name, ".rds")))
write.csv(data.frame(lambda = res$max_fwer$lambda, max_fwer = res$max_fwer$fwer),
          file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
cat(sprintf("Calibration %s (%s)\n", name, cal$design))
print(res)
