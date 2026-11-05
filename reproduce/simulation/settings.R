# ============================================================
# Settings of the simulation study (Section 4): the true probabilities of the eight scenarios,
# the data-generating values and the lambda of every design (calibrated with calibrate.R).
# ============================================================

N_MAX <- 150                 # maximum sample size of each subgroup
INTERIM <- c(0.5, 0.75)      # interim analyses after 50% and 75% of the patients
DELTA_T <- 0.05              # toxicity margin
ALPHA <- 0.10                # nominal FWER used to calibrate lambda in the simulation study
ALPHA_ILLUSTRATION <- 0.05   # nominal FWER used to calibrate lambda of the illustrative trial

# True response probabilities (rows: subgroups B1-B3; columns: A0-A3, A0 = control; NA = not
# evaluated) and true toxicity probabilities (A0-A3). In the simulation study A2 is not evaluated
# in B1 and A1 is not evaluated in B3.
make_scenario <- function(pR, piT) boat_scenario(pR = matrix(pR, nrow = 3, byrow = TRUE), piT = piT)

SCENARIOS <- list(
  make_scenario(c(0.50, 0.30,   NA, 0.70,
                  0.50, 0.70, 0.50, 0.30,
                  0.30,   NA, 0.50, 0.30), c(0.25, 0.35, 0.20, 0.15)),
  make_scenario(c(0.50, 0.70,   NA, 0.30,
                  0.50, 0.70, 0.50, 0.30,
                  0.50,   NA, 0.50, 0.30), c(0.25, 0.35, 0.20, 0.20)),
  make_scenario(c(0.50, 0.70,   NA, 0.70,
                  0.50, 0.70, 0.50, 0.30,
                  0.50,   NA, 0.50, 0.30), c(0.25, 0.35, 0.20, 0.20)),
  make_scenario(c(0.50, 0.20,   NA, 0.70,
                  0.30, 0.60, 0.30, 0.30,
                  0.30,   NA, 0.50, 0.30), c(0.25, 0.15, 0.20, 0.22)),
  make_scenario(c(0.50, 0.50,   NA, 0.70,
                  0.30, 0.50, 0.30, 0.30,
                  0.30,   NA, 0.50, 0.30), c(0.25, 0.15, 0.20, 0.22)),
  make_scenario(c(0.50, 0.45,   NA, 0.60,
                  0.30, 0.45, 0.30, 0.60,
                  0.30,   NA, 0.30, 0.60), c(0.25, 0.15, 0.20, 0.22)),
  make_scenario(c(0.50, 0.45,   NA, 0.60,
                  0.30, 0.45, 0.50, 0.60,
                  0.50,   NA, 0.50, 0.60), c(0.25, 0.15, 0.20, 0.22)),
  make_scenario(c(0.30, 0.60,   NA, 0.40,
                  0.30, 0.60, 0.50, 0.40,
                  0.30,   NA, 0.50, 0.40), c(0.25, 0.15, 0.20, 0.22))
)

# delta_R = 0.20 (Supplementary Tables S.14-S.17): the response probabilities of the experimental
# treatments are raised so that the same treatments stay admissible; toxicity is unchanged
SCENARIOS_DELTAR020 <- list(
  make_scenario(c(0.50, 0.45,   NA, 0.85,
                  0.50, 0.85, 0.65, 0.45,
                  0.30,   NA, 0.65, 0.45), c(0.25, 0.35, 0.20, 0.15)),
  make_scenario(c(0.50, 0.85,   NA, 0.45,
                  0.50, 0.85, 0.65, 0.45,
                  0.50,   NA, 0.65, 0.45), c(0.25, 0.35, 0.20, 0.20)),
  make_scenario(c(0.50, 0.85,   NA, 0.85,
                  0.50, 0.85, 0.65, 0.45,
                  0.50,   NA, 0.65, 0.45), c(0.25, 0.35, 0.20, 0.20)),
  make_scenario(c(0.50, 0.35,   NA, 0.85,
                  0.30, 0.75, 0.45, 0.45,
                  0.30,   NA, 0.65, 0.45), c(0.25, 0.15, 0.20, 0.22)),
  make_scenario(c(0.50, 0.65,   NA, 0.85,
                  0.30, 0.65, 0.45, 0.45,
                  0.30,   NA, 0.65, 0.45), c(0.25, 0.15, 0.20, 0.22)),
  make_scenario(c(0.50, 0.60,   NA, 0.75,
                  0.30, 0.60, 0.45, 0.75,
                  0.30,   NA, 0.45, 0.75), c(0.25, 0.15, 0.20, 0.22)),
  make_scenario(c(0.50, 0.60,   NA, 0.75,
                  0.30, 0.60, 0.65, 0.75,
                  0.50,   NA, 0.65, 0.75), c(0.25, 0.15, 0.20, 0.22)),
  make_scenario(c(0.30, 0.75,   NA, 0.55,
                  0.30, 0.75, 0.65, 0.55,
                  0.30,   NA, 0.65, 0.55), c(0.25, 0.15, 0.20, 0.22))
)

# null setting of the illustrative trial (Supplementary Figure S.8), used to calibrate its lambda:
# all four treatments in every subgroup and no admissible treatment
ILLUSTRATION_NULL <- make_scenario(c(0.50, 0.70, 0.50, 0.30,
                                     0.50, 0.70, 0.50, 0.30,
                                     0.50, 0.70, 0.50, 0.30), c(0.25, 0.35, 0.20, 0.20))

SCENARIO_SETS <- list(main = SCENARIOS, deltaR020 = SCENARIOS_DELTAR020, illustration_null = list(ILLUSTRATION_NULL))

# data-generating values of the main simulation: sigma_E = 1, rho = 0.6, sigma_O = 0.1,
# sigma_OR = 0.1 and M = 20 (the defaults of boat_dgm())
DGM_MAIN <- boat_dgm()

# lambda of the designs in the main simulation (Supplementary Table S.1)
LAMBDA_MAIN <- c("CD-SEP" = 0.76, "CD-BHM" = 0.78, "CD-AB" = 0.80, "BOAT-ER" = 0.68, "BOAT" = 0.70)
