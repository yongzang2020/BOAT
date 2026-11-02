# ============================================================
# Illustrative trial (Section 5)
#
# A hypothetical umbrella trial motivated by BATTLE-1 with subgroups B1 (EGFR), B2 (KRAS/BRAF)
# and B3 (VEGF/VEGFR-2), the control A0 (docetaxel) and the experimental treatments A1 (erlotinib),
# A2 (sorafenib) and A3 (vandetanib), all four evaluated in every subgroup.
#
# The data are simulated. No patient data of BATTLE-1 are used, and BATTLE-1 did not collect PDO
# data. Part 1 simulates the trial under BOAT and keeps its data. Part 2 analyzes these data step
# by step as during a real trial, following Steps C1-C5 of Supplementary Figure S.9. Part 3
# simulates the same trial under BOAT-ER for comparison.
# Parts 2 and 3 produce the decisions in Figure 4 and the numbers of patients in
# Supplementary Figure S.11.
#
# Usage: Rscript illustration.R   (about 6 minutes)
# ============================================================
library(BOAT)

cat("\n==================== Part 1: simulated data ====================\n")

# true response probabilities (rows B1-B3, columns A0-A3) and toxicity probabilities (Figure 3)
truth <- boat_scenario(pR = rbind(c(0.50, 0.70, 0.20, 0.50),
                                  c(0.30, 0.30, 0.60, 0.30),
                                  c(0.30, 0.30, 0.20, 0.50)),
                       piT = c(0.25, 0.22, 0.15, 0.20))
print(truth)

# Data-generating model. Patient i has a random effect b_i ~ N(0, sigma_E^2). The patient's PDOs
# are tested with every treatment open in the patient's subgroup, M replicates per treatment, with
# responses YO ~ N(mu_O,jk + rho * b_i, sigma_O^2), where mu_O,jk ~ N(mu_R,jk, sigma_OR^2). The
# patient is then randomized to treatment A_j. The clinical response is YR = I(mu_R,jk + b_i + e_i > 0)
# with e_i ~ N(0, 1), and the toxicity is YT ~ Bernoulli(piT_j). mu_R,jk is set so that the
# response probability equals the true one above.
dgm <- boat_dgm(sigma_E = 1, rho = 0.6, sigma_O = 0.1, sigma_OR = 0.1, M = 20)

# 150 patients per subgroup with interim analyses after 75 and 113 of them. lambda = 0.72 was
# calibrated with nominal FWER 0.05 under the null setting of Supplementary Figure S.10
# (Step D3 of Supplementary Figure S.9, simulation/calibrate.R).
design <- boat_design(N = 150, interim = c(0.5, 0.75), lambda = 0.72)
print(design)
bounds <- boat_boundary(design)
print(bounds)

# One trial under BOAT. keep_data = TRUE keeps its data. pat has one row per patient in enrollment
# order: id, subgroup k, assigned treatment j (0 = control), toxicity YT and response YR. pdo has
# one row per PDO replicate: patient id, subgroup k, treatment tested j and PDO response YO.
trial <- boat_simulate_trial(design, truth, dgm = dgm, seed = 1002, keep_data = TRUE)
pat <- trial$pat
pdo <- trial$pdo
print(head(pat))
print(head(pdo))

cat("\n=============== Part 2: step-by-step analysis of the data ===============\n")

# The MCMC of analysis l uses seed 1002 + 1000 l.

# data available when every subgroup has enrolled n patients
data_at <- function(n) {
  ids <- pat$id[ave(pat$id, pat$k, FUN = seq_along) <= n]
  list(pat = pat[pat$id %in% ids, ], pdo = pdo[pdo$id %in% ids, ], next_id = max(ids) + 1)
}
# randomization probabilities of a newly enrolled patient from the patient's mean PDO responses
randomize_next <- function(id, rar = NULL) {
  p <- pat[pat$id == id, ]
  x <- pdo[pdo$id == id, ]
  pdo_mean <- tapply(x$YO, x$j, mean)
  cat(sprintf("\nPatient %d (subgroup B%d), mean PDO response to each open treatment:\n", id, p$k))
  print(round(pdo_mean, 3))
  r <- boat_randomize(design, k = p$k, pdo_mean = pdo_mean, rar = rar)
  cat("Randomization probabilities:\n")
  print(round(r$prob, 3))
  cat(sprintf("Treatment received in the trial: A%d\n", p$j))
  invisible(r$prob)
}

# Step C1. Before the first interim analysis all treatments are equally likely.
randomize_next(1)

# Steps C2 and C3 for interim analysis 1 with 75 patients per subgroup. C2 fits the model, applies
# the decision rule and fits the randomization model used until the next analysis. C3 randomizes
# the next patient.
n1 <- bounds$n[bounds$analysis == "Interim 1"]
d1 <- data_at(n1[1])
fit1 <- boat_fit(design, d1$pat, d1$pdo, seed = 2002)
print(fit1)
dec1 <- boat_decide(design, fit1, n = n1)
print(dec1)
rar1 <- boat_rar_model(design, d1$pat, d1$pdo, active = dec1$active)
randomize_next(d1$next_id, rar1)

# Steps C2 and C3 for interim analysis 2 with 113 patients per subgroup.
n2 <- bounds$n[bounds$analysis == "Interim 2"]
d2 <- data_at(n2[1])
fit2 <- boat_fit(design, d2$pat, d2$pdo, seed = 3002)
dec2 <- boat_decide(design, fit2, n = n2, active = dec1$active)
print(dec2)
rar2 <- boat_rar_model(design, d2$pat, d2$pdo, active = dec2$active)
randomize_next(d2$next_id, rar2)

# Step C5. Final analysis with 150 patients per subgroup.
fit3 <- boat_fit(design, pat, pdo, seed = 4002)
print(fit3)
final <- boat_decide(design, fit3, n = design$N, active = dec2$active, final = TRUE)
print(final)

cat("\n=============== Part 3: the same trial under BOAT-ER ===============\n")

# BOAT-ER randomizes equally among the open treatments. lambda = 0.68 was calibrated in the same
# way as for BOAT.
design_er <- boat_design(N = 150, interim = c(0.5, 0.75), lambda = 0.68, randomization = "equal")
print(boat_boundary(design_er))
trial_er <- boat_simulate_trial(design_er, truth, dgm = dgm, seed = 1002)
for (l in seq_along(trial_er$interims)) {
  cat(sprintf("\nInterim analysis %d\n", l))
  t <- trial_er$interims[[l]]
  print(data.frame(subgroup = paste0("B", t$k), treatment = paste0("A", t$j), n = t$n,
                   pr_eff = round(t$pr_eff, 3), pr_safe = round(t$pr_safe, 3),
                   cutoff = round(t$cutoff, 3), decision = t$decision), row.names = FALSE)
}
cat("\nFinal analysis\n")
t <- trial_er$final
print(data.frame(subgroup = paste0("B", t$k), treatment = paste0("A", t$j), n = t$n,
                 pr_eff = round(t$pr_eff, 3), pr_safe = round(t$pr_safe, 3), cutoff = round(t$cutoff, 3),
                 decision = ifelse(t$pr_eff >= t$cutoff & t$pr_safe >= t$cutoff, "recommend", "not recommended")),
      row.names = FALSE)
print(trial_er)

# numbers of patients by subgroup and treatment under the two designs (Supplementary Figure S.11)
n_boat <- trial$n_by_trt_by_k
n_er <- trial_er$n_by_trt_by_k
dimnames(n_boat) <- dimnames(n_er) <- list(paste0("B", 1:3), paste0("A", 0:3))
cat("\nNumber of patients, BOAT-ER\n"); print(n_er)
cat("\nNumber of patients, BOAT\n"); print(n_boat)
