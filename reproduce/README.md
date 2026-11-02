# Supplementary code for "A Bayesian Precision Umbrella Trial Design Incorporating Patient-Derived Organoid Response Data"

This archive contains the R package **BOAT**, which implements the proposed design, the code of
the illustrative trial (Section 5) and the code of the simulation study (Section 4). BOAT and its
version with equal randomization (BOAT-ER) are run with the functions of the package. The
comparators CD-SEP, CD-BHM and CD-AB have their own code.

| Result | Code |
|---|---|
| Illustrative trial (Section 5, Figure 4 and Supplementary Figure S.11) | `illustration/illustration.R` |
| Calibration of lambda (procedure in Section 3, values in Sections 4.1 and 5) | `simulation/calibrate.R` |
| Main simulation (Supplementary Table S.1) | `simulation/main.R` |
| Sensitivity to the number of PDO replicates M (Supplementary Tables S.2-S.5) | `simulation/sens1_M.R` |
| Sensitivity to sigma_O (Supplementary Tables S.6-S.9) | `simulation/sens2_sigmaO.R` |
| Sensitivity to sigma_OR (Supplementary Tables S.10-S.13) | `simulation/sens3_sigmaOR.R` |
| Sensitivity to delta_R (Supplementary Tables S.14-S.17) | `simulation/sens4_deltaR.R` |
| Sensitivity to rho (Supplementary Table S.18) | `simulation/sens5_rho.R` |

## Contents

```
README.md
BOAT_0.1.0.tar.gz          R package BOAT (source package)
illustration/
  illustration.R           illustrative trial: simulated data and their step-by-step analysis
simulation/
  settings.R               scenarios, data-generating values and lambda of every design
  functions.R              runs the settings and scenarios of an analysis and summarizes them
  main.R                   main simulation
  sens1_M.R ... sens5_rho.R  sensitivity analyses 1-5
  calibrate.R              calibration of lambda
  comparators/
    trial_noPDO.R          trial flow shared by CD-SEP, CD-BHM and CD-AB
    CD_SEP.R               efficacy analysis of CD-SEP (separate beta posteriors by subgroup)
    CD_BHM.R               efficacy analysis of CD-BHM (Bayesian hierarchical model)
    CD_AB.R                efficacy analysis of CD-AB (subgroup clustering without PDO data)
```

## Requirements and installation

R (version 4.1 or later) with the package `msm`. We used R 4.5.3. In R:

```r
install.packages("msm")
install.packages("BOAT_0.1.0.tar.gz", repos = NULL, type = "source")
```

The simulation scripts run replications in parallel by forking (Linux and macOS). On Windows they
run on one core.

## 1. The BOAT package

| Step | Function | Purpose |
|---|---|---|
| Design | `boat_design()` | specify a design: subgroups, treatments, sample sizes, interim analyses, lambda, margins, randomization |
| | `boat_boundary()` | cutoffs of the interim and final analyses |
| Conduct | `boat_fit()` | fit the joint model of clinical and PDO responses (with adaptive subgroup clustering) and the toxicity model to the trial data |
| | `boat_decide()` | interim decisions (drop treatments, stop subgroups) and final recommendations |
| | `boat_rar_model()` | fit the model used for PDO-guided randomization |
| | `boat_randomize()` | randomization probabilities and treatment assignment of a new patient |
| Simulation | `boat_scenario()`, `boat_dgm()` | true probabilities and data-generating values |
| | `boat_generate_data()`, `boat_simulate_trial()` | trial data and one simulated trial |
| | `boat_get_oc()` | operating characteristics |
| | `boat_calibrate()` | calibrate lambda |

Every function has a help page with examples (for example `?boat_fit`). The vignette
(`vignette("BOAT")`) goes through the design, the conduct and the simulation of a trial.

## 2. Illustrative trial (Section 5)

Run from the folder `illustration/`:

```
Rscript illustration.R      # about 6 minutes
```

The illustrative trial is a hypothetical umbrella trial motivated by BATTLE-1, with subgroups
EGFR, KRAS/BRAF and VEGF/VEGFR-2, the control docetaxel and the experimental treatments erlotinib,
sorafenib and vandetanib. **Its data are simulated.** No patient data of BATTLE-1 are used, and
BATTLE-1 did not collect PDO data. `illustration.R` has three parts.

1. **Simulated data.** The script states the true toxicity and response probabilities (Figure 3)
   and the data-generating model with sigma_E = 1, rho = 0.6, sigma_O = 0.1, sigma_OR = 0.1 and
   M = 20 PDO replicates per treatment, and simulates one trial under BOAT (lambda = 0.72,
   seed 1002). The data are two data frames in the form that `boat_fit()` takes. `pat` has one
   row per patient in enrollment order with the subgroup `k`, the assigned treatment `j`
   (0 = control), the toxicity `YT` and the response `YR`. `pdo` has one row per PDO replicate
   with the patient `id`, `k`, the treatment tested `j` and the PDO response `YO`.
2. **Step-by-step analysis.** The data are analyzed as during the trial. At each interim analysis
   the script selects the patients enrolled so far, fits the model (`boat_fit()`), applies the
   decision rule (`boat_decide()`), fits the randomization model (`boat_rar_model()`) and computes
   the randomization probabilities of the next patient from the patient's PDO responses
   (`boat_randomize()`). The final analysis gives the recommended treatments.
3. **BOAT-ER.** The same trial is simulated under BOAT-ER (lambda = 0.68) for comparison.

Part 2 follows Steps C1-C5 of the flowchart in Supplementary Figure S.9. Step D3, the
calibration of lambda, is done by the calibrations `illustration_BOAT` and
`illustration_BOAT-ER` of `simulation/calibrate.R`.

Parts 2 and 3 produce the decisions in Figure 4 and the numbers of patients in
Supplementary Figure S.11.

## 3. Simulation study (Section 4)

| Script | Designs and settings | lambda |
|---|---|---|
| `main.R` | CD-SEP, CD-BHM, CD-AB, BOAT-ER, BOAT | 0.76, 0.78, 0.80, 0.68, 0.70 |
| `sens1_M.R` | BOAT with M = 10 and 30 | 0.70, 0.70 |
| `sens2_sigmaO.R` | BOAT with sigma_O = 0.05 and 0.20 | 0.70, 0.70 |
| `sens3_sigmaOR.R` | BOAT with sigma_OR = 0.05 and 0.20 | 0.67, 0.76 |
| `sens4_deltaR.R` | BOAT with delta_R = 0.20 | 0.71 |
| `sens5_rho.R` | BOAT with rho = 0.9, 0 and -0.6 | 0.71, 0.71, 0.82 |

Every analysis uses the eight scenarios of `settings.R` with 1000 replications, and replication r
uses seed 1000 + r. In each sensitivity analysis the rows of the main setting (M = 20,
sigma_O = 0.1, sigma_OR = 0.1, delta_R = 0.05, rho = 0.6) are the BOAT rows of the main
simulation. BOAT and BOAT-ER are simulated with `boat_get_oc()` of the package. The three
comparators share one trial flow (equal randomization, the interim and final decision rules of
Section 3, and the toxicity model) and differ only in the efficacy analysis. They are simulated
with `comparator_get_oc()` in `functions.R`, which returns its summary in the same format as
`boat_get_oc()`.

Run from the folder `simulation/`, giving the number of cores:

```
Rscript main.R 16            # every setting and scenario in turn
Rscript main.R BOAT 3 16     # one setting and scenario, for example as one job on a cluster
Rscript main.R summary
```

The result of each setting and scenario is saved in `results/<analysis>/` and is skipped when the
script is run again. `summary` writes `results/<analysis>/summary_<analysis>.csv` in the layout of
the supplementary tables: for each setting, scenario and subgroup the rows % Final Rec,
# Patients and % Patients for the treatments A0-A3, together with (P_k1, P_k2) for a subgroup
with admissible treatments or FWER_k for a subgroup without.

### Calibration of lambda

For every lambda on a grid, `calibrate.R` simulates the design under Scenarios 1-3, which contain
the null subgroups (B2 in Scenario 1, B1-B3 in Scenario 2 and B2-B3 in Scenario 3). The chosen
lambda is the grid value whose largest FWER_k over these subgroups is closest to 0.10, ties going
to the smallest value with FWER_k at most 0.10. BOAT and BOAT-ER are calibrated with
`boat_calibrate()` of the package and the comparators with `comparator_calibrate()` in
`functions.R`, which applies the same rule. The lambda of the illustrative trial is calibrated in
the same way under the null setting of Supplementary Figure S.10, in which no treatment is
admissible in any subgroup, with nominal FWER 0.05. The sensitivity analyses on M and sigma_O use
the lambda of the main simulation.

```
Rscript calibrate.R BOAT 16
```

The calibrations are `CD-SEP`, `CD-BHM`, `CD-AB`, `BOAT-ER` and `BOAT` (grid step 0.02);
`BOAT_sigmaOR0.05`, `BOAT_sigmaOR0.20`, `BOAT_deltaR0.20`, `BOAT_rho0.90`,
`BOAT_rho0.00` and `BOAT_rho-0.60` (grid step 0.01); and `illustration_BOAT-ER` and
`illustration_BOAT` (grid step 0.02). The result is printed and saved in `results/calibration/`.

### Computing time

On one core of a desktop computer, one simulated trial takes about 3 minutes for BOAT, 2 minutes
for BOAT-ER, 20 seconds for CD-AB, 10 seconds for CD-BHM and less than a second for CD-SEP. The main
simulation (8 scenarios, 5 designs, 1000 replications) thus takes about 800 core-hours, each BOAT
setting of a sensitivity analysis about 450 core-hours and each lambda of a BOAT calibration about
150 core-hours. A computing cluster is therefore recommended. We ran the simulations on a Linux
computing cluster.
