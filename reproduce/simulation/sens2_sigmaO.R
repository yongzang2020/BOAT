# ============================================================
# Sensitivity analysis 2: residual standard deviation sigma_O of the PDO responses
# (Supplementary Tables S.6-S.9). BOAT with sigma_O = 0.05 and 0.20, lambda = 0.70 for both.
# The rows with sigma_O = 0.1 are the BOAT rows of the main simulation (main.R).
# ============================================================
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.f) == 1) dirname(normalizePath(.f)) else getwd(), "functions.R"))

ANALYSIS <- list(
  name = "sens2_sigmaO", tables = "Supplementary Tables S.6-S.9", scenarios = 1:8,
  settings = list(setting("BOAT_sigmaO0.05", "BOAT", 0.70, boat_dgm(sigma_O = 0.05)),
                  setting("BOAT_sigmaO0.20", "BOAT", 0.70, boat_dgm(sigma_O = 0.20)))
)
run_analysis(ANALYSIS)
