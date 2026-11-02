# ============================================================
# Sensitivity analysis 5: association rho between the patient random effect and the PDO responses
# (Supplementary Table S.18). BOAT with rho = 0.9, 0 and -0.6, with lambda recalibrated for each
# value: 0.71, 0.71 and 0.82. The rows with rho = 0.6 are the BOAT rows of the main simulation
# (main.R).
# ============================================================
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.f) == 1) dirname(normalizePath(.f)) else getwd(), "functions.R"))

ANALYSIS <- list(
  name = "sens5_rho", tables = "Supplementary Table S.18", scenarios = 1:8,
  settings = list(setting("BOAT_rho0.90", "BOAT", 0.71, boat_dgm(rho = 0.9)),
                  setting("BOAT_rho0.00", "BOAT", 0.71, boat_dgm(rho = 0)),
                  setting("BOAT_rho-0.60", "BOAT", 0.82, boat_dgm(rho = -0.6)))
)
run_analysis(ANALYSIS)
