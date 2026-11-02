# ============================================================
# Sensitivity analysis 1: number of PDO replicates M (Supplementary Tables S.2-S.5).
# BOAT with M = 10 and M = 30, lambda = 0.70 for both. The rows with M = 20 are the BOAT rows of
# the main simulation (main.R).
# ============================================================
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.f) == 1) dirname(normalizePath(.f)) else getwd(), "functions.R"))

ANALYSIS <- list(
  name = "sens1_M", tables = "Supplementary Tables S.2-S.5", scenarios = 1:8,
  settings = list(setting("BOAT_M10", "BOAT", 0.70, boat_dgm(M = 10)),
                  setting("BOAT_M30", "BOAT", 0.70, boat_dgm(M = 30)))
)
run_analysis(ANALYSIS)
