# ============================================================
# Sensitivity analysis 4: minimum clinically meaningful response improvement delta_R
# (Supplementary Tables S.14-S.17). BOAT with delta_R = 0.20 and lambda = 0.71 under the true
# probabilities SCENARIOS_DELTAR020 (settings.R). The rows with delta_R = 0.05 are the BOAT rows
# of the main simulation (main.R).
# ============================================================
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.f) == 1) dirname(normalizePath(.f)) else getwd(), "functions.R"))

ANALYSIS <- list(
  name = "sens4_deltaR", tables = "Supplementary Tables S.14-S.17", scenarios = 1:8,
  settings = list(setting("BOAT_deltaR0.20", "BOAT", 0.71, truths = "deltaR020", delta_R = 0.20))
)
run_analysis(ANALYSIS)
