# ============================================================
# Sensitivity analysis 3: PDO-clinical linkage standard deviation sigma_OR
# (Supplementary Table S.4). BOAT with sigma_OR = 0.05 (lambda = 0.67) and 0.20
# (lambda = 0.76). The rows with sigma_OR = 0.1 are the BOAT rows of the main simulation (main.R).
# ============================================================
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.f) == 1) dirname(normalizePath(.f)) else getwd(), "functions.R"))

ANALYSIS <- list(
  name = "sens3_sigmaOR", tables = "Supplementary Table S.4", scenarios = 1:8,
  settings = list(setting("BOAT_sigmaOR0.05", "BOAT", 0.67, boat_dgm(sigma_OR = 0.05)),
                  setting("BOAT_sigmaOR0.20", "BOAT", 0.76, boat_dgm(sigma_OR = 0.20)))
)
run_analysis(ANALYSIS)
