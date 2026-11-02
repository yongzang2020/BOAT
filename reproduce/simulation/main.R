# ============================================================
# Main simulation (Supplementary Table S.1): Scenarios 1-8 under CD-SEP, CD-BHM, CD-AB, BOAT-ER and BOAT.
# lambda: CD-SEP 0.76, CD-BHM 0.78, CD-AB 0.80, BOAT-ER 0.68, BOAT 0.70.
# Usage: Rscript main.R [cores], or Rscript main.R <setting> <scenario> [cores] for one setting and
# scenario, then Rscript main.R summary (see functions.R).
# ============================================================
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.f) == 1) dirname(normalizePath(.f)) else getwd(), "functions.R"))

ANALYSIS <- list(
  name = "main", tables = "Supplementary Table S.1", scenarios = 1:8,
  settings = lapply(names(LAMBDA_MAIN), function(d) setting(d, d, LAMBDA_MAIN[[d]]))
)
run_analysis(ANALYSIS)
