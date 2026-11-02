# ============================================================
# Shared code of the simulation study. BOAT and BOAT-ER use boat_get_oc() and boat_calibrate() of
# the BOAT package. The comparators use the code in comparators/, and comparator_get_oc() and
# comparator_calibrate() below return their output in the same format.
#
# Usage of an analysis script (main.R, sens1_M.R, ...):
#   Rscript main.R [cores]                        every setting and scenario
#   Rscript main.R <setting> <scenario> [cores]   one setting and scenario, e.g. Rscript main.R BOAT 3 16
#   Rscript main.R summary                        writes results/<analysis>/summary_<analysis>.csv
# Each setting and scenario uses 1000 replications, replication r with seed 1000 + r. Its result is
# saved in results/<analysis>/ and skipped when the script is run again.
# ============================================================
suppressMessages({
  library(BOAT)
  library(parallel)
})

get_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
CODE_DIR <- get_script_dir()
source(file.path(CODE_DIR, "settings.R"))
NSIM <- 1000

# ---- settings ----
# A setting is one row block of a table: a design with its lambda, data-generating values,
# true probabilities and delta_R. key names the setting in the result files.
setting <- function(key, design, lambda, dgm = DGM_MAIN, truths = "main", delta_R = 0.05) {
  list(key = key, design = design, lambda = lambda, dgm = dgm, truths = truths, delta_R = delta_R)
}
scenario_truth <- function(s, truths) SCENARIO_SETS[[truths]][[s]]

# the design of a setting; the treatments evaluated in each subgroup are the non-NA entries of
# the scenario's response probabilities
make_design <- function(st, sc, lambda = st$lambda) {
  boat_design(arms = lapply(seq_len(nrow(sc$pR)), function(k) which(!is.na(sc$pR[k, ])) - 1),
              N = N_MAX, interim = INTERIM, lambda = lambda, delta_T = DELTA_T, delta_R = st$delta_R,
              randomization = if (st$design == "BOAT") "adaptive" else "equal")
}

# ---- comparators ----
comparators <- new.env()
for (f in c("trial_noPDO.R", "CD_SEP.R", "CD_BHM.R", "CD_AB.R"))
  sys.source(file.path(CODE_DIR, "comparators", f), envir = comparators)

# one simulated trial of a comparator, with the sample sizes, interim analyses, margins and lambda
# of a boat_design object
comparator_trial <- function(name, design, sc, seed) {
  args <- list(pE_mat = sc$pR, piT_vec = sc$piT, max_n_k = design$N, interim_fracs = design$interim,
               delta_E = design$delta_R, delta_T = design$delta_T, lambda = design$lambda,
               prior_piT_a = 1, prior_piT_b = 1, seed = seed)
  mcmc <- list(mcmc_iter = 10000, mcmc_burn = 5000, mcmc_thin = 1)
  switch(name,
    "CD-SEP" = do.call(comparators$simulate_trial_CD_SEP, c(args, list(prior_pE_a = 1, prior_pE_b = 1))),
    "CD-BHM" = do.call(comparators$simulate_trial_CD_BHM, c(args, mcmc, list(prior_tau2_a = 0.01, prior_tau2_b = 0.01))),
    "CD-AB"  = do.call(comparators$simulate_trial_CD_AB, c(args, mcmc)),
    stop("unknown design ", name))
}

# admissible treatments of subgroup k
admissible <- function(sc, k, delta_R, delta_T) {
  J <- ncol(sc$pR) - 1
  which(sapply(seq_len(J), function(j) !is.na(sc$pR[k, j + 1]) && sc$pR[k, j + 1] > sc$pR[k, 1] + delta_R &&
                 sc$piT[j + 1] < sc$piT[1] + delta_T))
}

# operating characteristics of simulated trials, in the format of boat_get_oc()$table
oc_table <- function(trials, sc, design) {
  R <- length(trials); K <- nrow(sc$pR); J <- ncol(sc$pR) - 1
  rec <- matrix(0, K, J + 1); n_sum <- matrix(0, K, J + 1); hit1 <- hit2 <- numeric(K)
  for (tr in trials) {
    for (k in 1:K) {
      ch <- tr$final_selection[[as.character(k)]]; ch <- sort(unique(ch[ch %in% 1:J]))
      if (length(ch) == 0) rec[k, 1] <- rec[k, 1] + 1 else rec[k, ch + 1] <- rec[k, ch + 1] + 1
      A <- admissible(sc, k, design$delta_R, design$delta_T)
      if (length(ch) > 0 && length(A) > 0) {
        hit1[k] <- hit1[k] + all(ch %in% A)
        hit2[k] <- hit2[k] + any(ch %in% A)
      }
    }
    n_sum <- n_sum + tr$n_by_trt_by_k
  }
  mean_n <- n_sum / R; pct_n <- 100 * mean_n / rowSums(mean_n)
  out <- list()
  for (k in 1:K) {
    for (a in 0:J) {
      if (is.na(sc$pR[k, a + 1])) next
      out[[length(out) + 1]] <- data.frame(subgroup = k, treatment = a,
                                           measure = c("final_rec", "n_patients", "pct_patients"),
                                           value = round(c(100 * rec[k, a + 1] / R, mean_n[k, a + 1], pct_n[k, a + 1]), 1))
    }
    if (length(admissible(sc, k, design$delta_R, design$delta_T)) > 0) {
      out[[length(out) + 1]] <- data.frame(subgroup = k, treatment = NA, measure = c("power1", "power2"),
                                           value = round(100 * c(hit1[k], hit2[k]) / R, 1))
    } else {
      out[[length(out) + 1]] <- data.frame(subgroup = k, treatment = NA, measure = "fwer",
                                           value = round(100 * (1 - rec[k, 1] / R), 1))
    }
  }
  do.call(rbind, out)
}

# counterpart of boat_get_oc() for a comparator
comparator_get_oc <- function(name, design, sc, nsim = NSIM, seeds = 1000 + seq_len(nsim), cores = 1) {
  run <- function(s) comparator_trial(name, design, sc, s)
  trials <- if (cores > 1 && .Platform$OS.type != "windows")
    mclapply(seeds, run, mc.cores = cores, mc.preschedule = FALSE) else lapply(seeds, run)
  list(table = oc_table(trials, sc, design), trials = trials)
}

# counterpart of boat_calibrate() for a comparator: the lambda whose largest FWER over the null
# subgroups is closest to alpha, ties going to the smallest lambda with FWER at most alpha
comparator_calibrate <- function(name, design, scenarios, lambda, nsim = NSIM, seeds = 1000 + seq_len(nsim),
                                 alpha = ALPHA, cores = 1) {
  rows <- list()
  for (lam in lambda) {
    d <- design; d$lambda <- lam
    for (s in seq_along(scenarios)) {
      tab <- comparator_get_oc(name, d, scenarios[[s]], nsim, seeds, cores)$table
      f <- tab[tab$measure == "fwer", ]
      if (nrow(f)) rows[[length(rows) + 1]] <- data.frame(lambda = lam, scenario = s, subgroup = f$subgroup,
                                                          fwer = f$value / 100)
    }
  }
  tab <- do.call(rbind, rows)
  mx <- aggregate(fwer ~ lambda, data = tab, FUN = max)
  d <- abs(mx$fwer - alpha)
  best <- mx[abs(d - min(d)) < 1e-12, ]
  chosen <- if (any(best$fwer <= alpha)) min(best$lambda[best$fwer <= alpha]) else min(best$lambda)
  structure(list(table = tab, max_fwer = mx, lambda = chosen, alpha = alpha), class = "boat_calibration")
}

# ---- running an analysis ----
result_dir <- function(analysis) file.path(CODE_DIR, "results", analysis$name)
job_file <- function(analysis, st, s) file.path(result_dir(analysis), sprintf("%s_S%d.rds", st$key, s))

# operating characteristics of one setting under one scenario
run_job <- function(analysis, st, s, cores) {
  f <- job_file(analysis, st, s)
  if (file.exists(f)) return(invisible(NULL))
  sc <- scenario_truth(s, st$truths)
  design <- make_design(st, sc)
  t0 <- proc.time()[["elapsed"]]
  tab <- if (st$design %in% c("BOAT", "BOAT-ER")) boat_get_oc(design, sc, nsim = NSIM, dgm = st$dgm, cores = cores)$table
         else comparator_get_oc(st$design, design, sc, cores = cores)$table
  dir.create(result_dir(analysis), showWarnings = FALSE, recursive = TRUE)
  saveRDS(tab, f)
  cat(sprintf("%s, Scenario %d: done in %.1f minutes\n", st$key, s, (proc.time()[["elapsed"]] - t0) / 60))
}

# the results in the layout of the supplementary tables
summarize_analysis <- function(analysis) {
  rows <- list()
  for (st in analysis$settings) for (s in analysis$scenarios) {
    f <- job_file(analysis, st, s)
    if (!file.exists(f)) next
    tab <- readRDS(f)
    for (k in sort(unique(tab$subgroup))) {
      t <- tab[tab$subgroup == k, ]
      get <- function(ms) sapply(0:3, function(a) { v <- t$value[t$treatment %in% a & t$measure == ms]; if (length(v)) v else NA })
      blk <- data.frame(setting = st$key, design = st$design, lambda = st$lambda, scenario = s, subgroup = k,
                        measure = c("% Final Rec", "# Patients", "% Patients"),
                        rbind(get("final_rec"), get("n_patients"), get("pct_patients")),
                        P_k1 = NA, P_k2 = NA, FWER_k = NA)
      names(blk)[7:10] <- paste0("A", 0:3)
      p <- t[is.na(t$treatment), ]
      if (any(p$measure == "fwer")) {
        blk$FWER_k[1] <- p$value[p$measure == "fwer"]
      } else {
        blk$P_k1[1] <- p$value[p$measure == "power1"]; blk$P_k2[1] <- p$value[p$measure == "power2"]
      }
      rows[[length(rows) + 1]] <- blk
    }
  }
  if (length(rows) == 0) { cat("No results yet in", result_dir(analysis), "\n"); return(invisible(NULL)) }
  out <- do.call(rbind, rows)
  f <- file.path(result_dir(analysis), sprintf("summary_%s.csv", analysis$name))
  write.csv(out, f, row.names = FALSE, na = "")
  cat(sprintf("%s (%s): summary of %d settings and scenarios written to %s\n", analysis$name, analysis$tables,
              length(unique(paste(out$setting, out$scenario))), f))
  invisible(out)
}

run_analysis <- function(analysis, args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) >= 1 && args[1] == "summary") return(invisible(summarize_analysis(analysis)))
  keys <- sapply(analysis$settings, `[[`, "key")
  if (length(args) >= 2 && args[1] %in% keys) {
    jobs <- list(list(st = analysis$settings[[match(args[1], keys)]], s = as.integer(args[2])))
    cores <- if (length(args) >= 3) as.integer(args[3]) else 1L
  } else {
    jobs <- list()
    for (st in analysis$settings) for (s in analysis$scenarios) jobs[[length(jobs) + 1]] <- list(st = st, s = s)
    cores <- if (length(args) >= 1) as.integer(args[1]) else 1L
  }
  for (j in jobs) run_job(analysis, j$st, j$s, cores)
  summarize_analysis(analysis)
}
