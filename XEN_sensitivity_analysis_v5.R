t00 <- proc.time()[3]
# Set work.dir to the folder holding this script, scenario.csv, and df_to_latex.R
work.dir <- '.'

setwd(work.dir)
source("df_to_latex.R")
# Create a compile directory inside your project folder
nimble_dir <- file.path(work.dir, "nimble_compiled")
dir.create(nimble_dir, showWarnings = FALSE)
Sys.setenv(TMPDIR = nimble_dir)
options(nimble.dirName = nimble_dir)
################################################################################
#  XEN Glaucoma Treatment System — Ethnic Bridging Study
#  Primary Analysis + Three-Pillar Sensitivity Framework
#
#  Converted from rjags to nimble for flexibility in implementing
#  power-prior, commensurate-prior, and structural-model variants.
#
#  Analyses implemented:
#    PRIMARY  — Emax + linear longitudinal model (Report Eq. 1–2)
#    S1       — Heterogeneity diagnostics (posterior predictive checks)
#    S2       — Leave-one-source-out (4 refits)
#    S3       — No-borrowing reference (IPD only)
#    S4       — Effective sample size decomposition
#    S5       — Prior sensitivity (original / vague / informative)
#    S6       — Tipping point analysis (shift Hu et al. data)
#    S7       — Alternative borrowing methods (power prior, commensurate prior)
#    S8       — Structural model sensitivity (Emax-only, piecewise linear)
#
#  Requirements: nimble, coda, ggplot2, dplyr, tidyr
#
#  CHANGE LOG (v2 — corrected):
#    Fix 2:  Commensurate prior — e0_agg[s] now centred on (e0 + de0),
#            emax_agg[s] centred on (emax + demax), i.e. Asian intercept.
#            Also removed demax*race_agg from mu_agg to avoid double-counting.
#    Fix 3:  LOSO (A2) — baseline IOP recomputed per reduced dataset.
#    Fix 5:  Tipping point (A6) — criterion changed to 95% CrI of delta
#            parameters excluding zero, matching manuscript definition.
#    Fix 6:  Power prior (A7) — three-way split: P13-001 IPD (full),
#            RWE aggregate (full weight), Hu/Sng aggregate (discounted).
#    Fix 7:  ESS (A4) — replaced ad-hoc formula with Morita et al. (2008)
#            style variance-ratio approach.
#    Fix 8:  PPC (A1) — added P13_NonAsian to source list.
#    Fix 9:  Piecewise linear (A8) — three knots (Day 1, Day 28, Day 84)
#            matching manuscript description.
#    Fix 10: WAIC — added calculateWAIC() calls for A8 model comparison.
#
#  CHANGE LOG (v2.1 — dataset compatibility):
#    Fix A:  NA removal — 3 P13-001 rows with missing pdpch dropped on load.
#    Fix B:  Power prior data split uses source-based filter (not n==1) to
#            prevent RWE IPD double-counting after simulated IPD integration.
#    Fix C:  Commensurate mu_agg — removed residual demax*race_agg term
#            that double-counted the ethnic Emax offset (emax_agg already
#            centred on emax+demax).
#    Fix D:  ESS (A4) patient count — uses unique(pid) instead of rows/9,
#            since RWE patients have variable visits (not 9 per patient).
#
#  CHANGE LOG (v3 — review response, April 2026):
#    Fix A1: Commensurate prior restructured — P13-001 AND RWE now enter
#            the IPD block with full borrowing (they are actual patient
#            data).  Only Hu and Sng (digitized aggregate from published
#            papers) receive source-specific e0_agg and emax_agg
#            parameters linked via tau_c.  N_sources reduced from 3 to 2.
#    Fix A3: Baseline IOP centering (used in compute_fitted_curves) now
#            uses the sample-size-weighted mean via weighted_base_iop(),
#            matching the manuscript formula. Applied in the primary,
#            LOSO, and no-borrowing analyses.
#    Fix A7: sigma^2 prior changed from IG(100, 1) (informative, prior
#            mean ~0.01) to IG(0.01, 0.01) (weakly informative) across
#            all eight model blocks.
#    Table fix: LOSO output now includes the 95% CrI column for ded50.
#    Table fix: Tipping-point output now includes the mean and 95% CrI
#               for dk (= Delta b in manuscript).
#
#  CHANGE LOG (v4 — commensurate/MAP scope alignment):
#    Fix A1b: (SUPERSEDED BY v5)  Initial attempt extended commensurate to
#             all four Asian sources.  After further review, v5 adopts
#             the cleaner scope described below.
#
#  CHANGE LOG (v5 — final scope for S7 alternative borrowing methods):
#    All three S7 methods (power prior, commensurate prior, robust MAP
#    prior) now target the same two external sources: Hu and Sng.
#    This aligns the scope with the study objective, which is to
#    quantify Asian-vs-non-Asian differences by borrowing Asian
#    information from the external literature sources.
#
#    Power prior (unchanged):
#        IPD + RWE:  full weight (a0 = 1)
#        Hu + Sng:   discounted by a0 in {0.25, 0.50, 0.75, 1.00}
#
#    Commensurate prior:
#        IPD block  = P13-001 (both ethnicities) + RWE (full borrowing).
#        External   = Hu + Sng, each with source-specific e0_agg[s] and
#                     emax_agg[s] linked to (e0+de0, emax+demax) via
#                     tau_c ~ Gamma(1,1).
#        N_sources  = 2 (s=1 Hu, s=2 Sng).
#
#    Robust MAP prior:
#        IPD block  = P13-001 (both ethnicities) + RWE.  Shared ethnic
#                     offsets (de0, demax) receive the robust mixture
#                     prior for protection against prior-data conflict.
#        External   = Hu + Sng, each with source-specific de0_ext[s],
#                     demax_ext[s] drawn from N(de0_pop, tau_map^2) and
#                     N(demax_pop, tau_map^2).
#        N_ext_sources = 2.
#        The "current study" in the MAP mixture is the pooled P13+RWE
#        ethnic offset, not a source-specific de0_s[1] as in v3.
#
#    Output tables: the S7 MAP summary now reports de0, demax (current-
#    study offsets) rather than de0_s[1], demax_s[1].
#
#  CHANGE LOG (v6 — S1 PPC output fix):
#    Fix S1: The posterior predictive check function previously emitted
#            row-level output (one row per patient-visit for IPD sources),
#            lacked the number-of-patients N, SE of the source-level mean,
#            and the standardised residual z used in manuscript
#            Table \ref{tab:ppc}, and did not write a LaTeX file.
#            The function is now rewritten to produce visit-level
#            aggregated output with columns:
#              source, day, N, obs_mean, predicted, residual, se_mean,
#              z, p_value
#            IPD sources are aggregated across patients within each visit
#            using the within-visit SD to compute se_mean.  AD sources
#            pass through with se_mean = NA and z = NA.  A LaTeX table
#            (S1_ppc_table.tex) is now written for direct inclusion in
#            the manuscript, and diagnostic flags for |z| > 2 and mean
#            absolute residual by source are printed to the console.
#
#  CHANGE LOG (v7 — convergence improvements + Windows parallel chains):
#    Goal:  Fix the R-hat ~ 1.08 - 1.10 observed for de0, demax, and emax
#           under v6, and exploit the user's multi-core machine to speed
#           up the pipeline, WITHOUT increasing total runtime beyond the
#           ~4 hours the v6 sequential pipeline required.
#
#    Fix C1: Added a block sampler targeting (e0, de0, emax, demax).
#            These intercept-level Emax parameters are strongly
#            correlated through the mean function, and separate RW
#            updates cannot navigate the ridge efficiently.  Block
#            updates handle this directly and typically drop R-hat
#            from ~ 1.10 to < 1.02 at the same chain length.  This is
#            the single most effective convergence improvement per
#            iteration in this model and is runtime-neutral.
#    Fix C2: Added slice samplers for r1 and dr1.  The posterior on r1
#            has a long right tail (v6: mean 4.5, 97.5% quantile ~ 13),
#            which causes RW updates to be mis-tuned across the body
#            and tail regions.  Slice sampling handles heavy-tailed
#            univariate marginals with no tuning and at comparable
#            per-iteration cost; it runs in addition to the shape
#            block sampler on (ed50, ded50, r1, dr1).
#    Fix C3: Chains within a single run_nimble() call can optionally
#            execute in parallel via a Windows-compatible PSOCK cluster
#            (parallel::makeCluster + parLapply), by passing
#            parallel = TRUE.  The default is parallel = FALSE because
#            many Windows systems have Application Control / AppLocker
#            / WDAC policies that block the loading of nimble's
#            compiled DLLs in worker processes (symptom: "unable to
#            load shared object ... An Application Control policy has
#            blocked this file").  When parallel = TRUE is requested
#            and a policy block is detected, the code falls back to
#            sequential execution for the remainder of the session
#            with a single warning message.  If your environment
#            permits parallel execution, pass parallel = TRUE at the
#            specific call sites that would benefit.
#    Fix C4: The primary analysis call retains enableWAIC = TRUE, which
#            forces the sequential code path so WAIC can be read from
#            a single compiled MCMC object.  This matches v6 behaviour
#            for the one fit that needs WAIC for the S8 comparison.
#            If primary R-hat is still marginal after v7, the
#            recommended local fix is to add nIter = 20000 to the
#            primary call only.
#    Fix C5: Default nBurnin / nIter / nThin / nChains kept at v6
#            values (10000 / 10000 / 10 / 2) so per-chain runtime
#            matches v6.  The speedup comes from running chains
#            simultaneously, not from fewer iterations.  Explicit MCMC
#            settings previously scattered across call sites were
#            removed so that any future change to the defaults applies
#            uniformly.
#    Fix C6: Added tau2 to delta_params so sigma^2 is reported in the
#            primary summary output alongside the other parameters.
#    Fix C7: New section 3b produces paired fitted-curve comparison
#            plots (Asian vs Non-Asian) matching manuscript Figures
#            \ref{fig:comparison} and \ref{fig:comparison_no102009}.
#            Plots are written as comparison.png (with patient 102009)
#            and comparison_no102009.png (without), at 300 dpi.  The
#            "with 102009" plot requires a second dataset
#            (XEN_nimble_ready_data_v2_allRWE_simRWE.csv); if not
#            present, that plot is skipped with a message.
#
#   Code-to-Math Mapping (see manuscript Eq. 2):
#     e0    = E_0          emax  = E_max        ed50  = ED_50
#     r1    = r            k     = b            a     = a
#     de0   = Delta E_0    demax = Delta E_max  ded50 = Delta ED_50
#     dr1   = Delta r      dk    = Delta b      tau2  = sigma^2
################################################################################

library(nimble)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)   # [v7] for parallel-chain MCMC execution

set.seed(1223456)
fileNameTime <- format(Sys.time(), "%Y-%m-%d_%H_%M_%S")


################################################################################
# 0.  DATA PREPARATION
################################################################################

# ---- Load data ----
# Dataset with simulated RWE IPD (excluding patient 102009)
#y.iop <- read.csv("XEN_nimble_ready_data_v2_no102009_simRWE.csv")
y.iop <- read.csv("scenario.csv")
#y.iop <- read.csv("XEN_nimble_realRWE_no102009_v2.csv")

# Verify expected columns
stopifnot(all(c("pid", "study", "visit", "n", "race", "base.pd",
                "pd", "pdch", "pdpch", "source", "base.pd.cent")
              %in% names(y.iop)))

# Remove rows with missing pdpch (the model response variable).
# In the P13-001 dataset, 3 observations have missing IOP values
# (patients 53-107 at Day 7/14 and 66-101 at Day 14).
n_before <- nrow(y.iop)
y.iop <- y.iop[!is.na(y.iop$pdpch), ]
cat("Data loaded:", n_before, "rows,", n_before - nrow(y.iop),
    "removed for missing pdpch,", nrow(y.iop), "retained\n")
cat("Sources:", paste(sort(unique(y.iop$source)), collapse = ", "), "\n")
cat("Race coding:  0 = Non-Asian, 1 = Asian\n")
print(table(y.iop$source, y.iop$race))


################################################################################
# 1.  MODEL DEFINITIONS (nimble)
################################################################################

# ============================================================================ #
#  1A. PRIMARY MODEL — Emax + linear with baseline IOP covariate
#      (Report Equations 1–2)
# ============================================================================ #

code_primary <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  
  # --- Priors (Report Section 2.2.1) ---
  e0    ~ dnorm(0, sd = sqrt(10))    # N(0, 10)
  emax  ~ dnorm(0, sd = 1)           # N(0, 1)
  k     ~ dnorm(0, sd = 10)          # N(0, 100)
  ed50  ~ T(dnorm(10, sd = sqrt(10)), 0, )   # N(10, 10) truncated > 0
  r1    ~ T(dnorm(0, sd = 10), 0, )          # N(0, 100) truncated > 0
  
  de0   ~ dnorm(0, sd = 1)           # N(0, 1)
  demax ~ dnorm(0, sd = 1)           # N(0, 1)
  dk    ~ dnorm(0, sd = sqrt(10))    # N(0, 10)
  ded50 ~ T(dnorm(0, sd = sqrt(10)), -ed50, )   # truncated so ed50+ded50 > 0
  dr1   ~ T(dnorm(0, sd = 10), -r1, )           # truncated so r1+dr1 > 0
  
  a     ~ dnorm(0, sd = 10)          # N(0, 100)
  tau2  ~ dinvgamma(0.01, 0.01)    # Inversegamma(0.01, 0.01) -- weakly informative
})


# ============================================================================ #
#  1B. EMAX-ONLY MODEL  (Analysis S8 — drop linear component)
# ============================================================================ #

code_emax_only <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  
  e0    ~ dnorm(0, sd = sqrt(10))
  emax  ~ dnorm(0, sd = 1)
  ed50  ~ T(dnorm(10, sd = sqrt(10)), 0, )
  r1    ~ T(dnorm(0, sd = 10), 0, )
  
  de0   ~ dnorm(0, sd = 1)
  demax ~ dnorm(0, sd = 1)
  ded50 ~ T(dnorm(0, sd = sqrt(10)), -ed50, )
  dr1   ~ T(dnorm(0, sd = 10), -r1, )
  
  a     ~ dnorm(0, sd = 10)
  tau2  ~ dinvgamma(0.01, 0.01)
})


# ============================================================================ #
#  1C. PIECEWISE LINEAR MODEL  (Analysis S8 — alternative structural model)
#      [FIX 9] Three knots at Day 1, Day 28, Day 84 (four segments)
#      matching the manuscript description of "knots at clinically
#      meaningful time points"
# ============================================================================ #

code_piecewise <- nimbleCode({
  
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    
    # Four segments:
    #   Segment 1: 0 < x <= cp1   (immediate post-op, Day 0-1)
    #   Segment 2: cp1 < x <= cp2  (early recovery, Day 1-28)
    #   Segment 3: cp2 < x <= cp3  (mid-term stabilisation, Day 28-84)
    #   Segment 4: x > cp3         (long-term maintenance, Day 84+)
    
    # Value at each knot (for continuity):
    val_cp1[i] <- (b1 + db1 * race[i]) * cp1
    val_cp2[i] <- val_cp1[i] + (b2 + db2 * race[i]) * (cp2 - cp1)
    val_cp3[i] <- val_cp2[i] + (b3 + db3 * race[i]) * (cp3 - cp2)
    
    # Piecewise function (using indicator logic)
    pw[i] <- step(cp1 - x[i]) *
      (b1 + db1 * race[i]) * x[i] +
      step(x[i] - cp1 - 0.001) * step(cp2 - x[i]) *
      (val_cp1[i] + (b2 + db2 * race[i]) * (x[i] - cp1)) +
      step(x[i] - cp2 - 0.001) * step(cp3 - x[i]) *
      (val_cp2[i] + (b3 + db3 * race[i]) * (x[i] - cp2)) +
      step(x[i] - cp3 - 0.001) *
      (val_cp3[i] + (b4 + db4 * race[i]) * (x[i] - cp3))
    
    mu[i] <- e0 + de0 * race[i] + a * z[i] + pw[i]
  }
  
  # Fixed knots
  cp1 <- 1
  cp2 <- 28
  cp3 <- 84
  
  # Priors
  e0  ~ dnorm(0, sd = sqrt(10))
  de0 ~ dnorm(0, sd = 1)
  
  b1  ~ dnorm(0, sd = 1)         # slope: immediate post-op (very steep)
  db1 ~ dnorm(0, sd = 1)
  b2  ~ dnorm(0, sd = 0.1)       # slope: early recovery (moderate)
  db2 ~ dnorm(0, sd = 0.1)
  b3  ~ dnorm(0, sd = 0.01)      # slope: mid-term (shallow)
  db3 ~ dnorm(0, sd = 0.01)
  b4  ~ dnorm(0, sd = 0.01)      # slope: long-term (near zero)
  db4 ~ dnorm(0, sd = 0.01)
  
  a   ~ dnorm(0, sd = 10)
  tau2 ~ dinvgamma(0.01, 0.01)
})


# ============================================================================ #
#  1D. POWER PRIOR MODEL  (Analysis S7)
#      [FIX 6] Three-way split:
#        - P13-001 IPD:   full weight (var = tau2)
#        - RWE aggregate: full weight (var = tau2 / n_rwe)
#        - Hu/Sng aggregate: discounted (var = tau2 / (n_disc * alpha0))
# ============================================================================ #

code_power_prior <- nimbleCode({
  
  # --- IPD observations (P13-001 only, all n=1) ---
  for (i in 1:N_ipd) {
    y[i] ~ dnorm(mu[i], var = tau2)
    
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  
  # --- RWE aggregate observations: FULL WEIGHT ---
  for (j in 1:N_rwe) {
    y_rwe[j] ~ dnorm(mu_rwe[j], var = tau2 / n_rwe[j])
    
    mu_rwe[j] <- e0 + de0 * race_rwe[j] + a * z_rwe[j] +
      (k + dk * race_rwe[j]) * x_rwe[j] +
      (emax + demax * race_rwe[j]) * x_rwe[j]^(r1 + dr1 * race_rwe[j]) /
      ((ed50 + ded50 * race_rwe[j])^(r1 + dr1 * race_rwe[j]) +
         x_rwe[j]^(r1 + dr1 * race_rwe[j]))
  }
  
  # --- Hu/Sng aggregate observations: DISCOUNTED by alpha0 ---
  for (m in 1:N_disc) {
    y_disc[m] ~ dnorm(mu_disc[m], var = tau2 / (n_disc[m] * alpha0))
    
    mu_disc[m] <- e0 + de0 * race_disc[m] + a * z_disc[m] +
      (k + dk * race_disc[m]) * x_disc[m] +
      (emax + demax * race_disc[m]) * x_disc[m]^(r1 + dr1 * race_disc[m]) /
      ((ed50 + ded50 * race_disc[m])^(r1 + dr1 * race_disc[m]) +
         x_disc[m]^(r1 + dr1 * race_disc[m]))
  }
  
  # alpha0 is passed as data (fixed) for each scenario
  
  # --- Priors (same as primary) ---
  e0    ~ dnorm(0, sd = sqrt(10))
  emax  ~ dnorm(0, sd = 1)
  k     ~ dnorm(0, sd = 10)
  ed50  ~ T(dnorm(10, sd = sqrt(10)), 0, )
  r1    ~ T(dnorm(0, sd = 10), 0, )
  de0   ~ dnorm(0, sd = 1)
  demax ~ dnorm(0, sd = 1)
  dk    ~ dnorm(0, sd = sqrt(10))
  ded50 ~ T(dnorm(0, sd = sqrt(10)), -ed50, )
  dr1   ~ T(dnorm(0, sd = 10), -r1, )
  a     ~ dnorm(0, sd = 10)
  tau2  ~ dinvgamma(0.01, 0.01)
})


# ============================================================================ #
#  1E. COMMENSURATE PRIOR MODEL  (Analysis S7)
#      [v5 — final scope, aligned with power prior] The commensurate prior
#      is applied only to the two digitized external aggregate sources
#      (Hu and Sng).  P13-001 (non-Asian and Asian) and the Hainan RWE
#      cohort are direct current-study data and enter the IPD block with
#      full borrowing.  This matches the power-prior scope (where Hu and
#      Sng are discounted by a0) and makes the three alternative methods
#      target the same borrowing decision: how much should we trust the
#      external Chinese/Asian literature data when inferring ethnic
#      differences?
#
#      Source indexing in the aggregate block:
#          s = 1  Hu et al.
#          s = 2  Sng et al.
#
#      e0_agg[s] is centred on (e0 + de0)   — the Asian population intercept
#      emax_agg[s] is centred on (emax + demax) — the Asian population Emax
#      Both Hu and Sng rows have race_agg = 1 (Asian), so ethnic offsets
#      apply. The source-specific e0_agg and emax_agg absorb the de0 and
#      demax offsets for Hu and Sng; do NOT add them again in mu_agg.
# ============================================================================ #

code_commensurate <- nimbleCode({
  
  # --- IPD block: P13-001 (both ethnicities) + RWE, full borrowing ---
  for (i in 1:N_ipd) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  
  # --- External aggregate block: Hu and Sng only ---
  # Both sources are Asian (race_agg = 1 for all rows).  Source-specific
  # e0_agg/emax_agg are anchored to the Asian population values via tau_c.
  for (j in 1:N_agg) {
    y_agg[j] ~ dnorm(mu_agg[j], var = tau2 / n_agg[j])
    
    mu_agg[j] <- e0_agg[source_id[j]] + a * z_agg[j] +
      (k + dk) * x_agg[j] +
      emax_agg[source_id[j]] *
      x_agg[j]^(r1 + dr1) /
      ((ed50 + ded50)^(r1 + dr1) +
         x_agg[j]^(r1 + dr1))
  }
  
  # Commensurate priors centred on Asian population parameters
  for (s in 1:N_sources) {
    e0_agg[s]   ~ dnorm(e0 + de0,     var = 1 / tau_c)
    emax_agg[s] ~ dnorm(emax + demax, var = 1 / tau_c)
  }
  
  # Commensurability parameter: large tau_c = strong borrowing
  tau_c ~ dgamma(1, 1)
  
  # --- Priors (same as primary) ---
  e0    ~ dnorm(0, sd = sqrt(10))
  emax  ~ dnorm(0, sd = 1)
  k     ~ dnorm(0, sd = 10)
  ed50  ~ T(dnorm(10, sd = sqrt(10)), 0, )
  r1    ~ T(dnorm(0, sd = 10), 0, )
  de0   ~ dnorm(0, sd = 1)
  demax ~ dnorm(0, sd = 1)
  dk    ~ dnorm(0, sd = sqrt(10))
  ded50 ~ T(dnorm(0, sd = sqrt(10)), -ed50, )
  dr1   ~ T(dnorm(0, sd = 10), -r1, )
  a     ~ dnorm(0, sd = 10)
  tau2  ~ dinvgamma(0.01, 0.01)
})


# ============================================================================ #
#  1F. ROBUST META-ANALYTIC-PREDICTIVE (MAP) PRIOR MODEL  (Analysis S7)
#      (Unchanged from original — included for completeness)
# ============================================================================ #

dRobustNorm <- nimbleFunction(
  run = function(x        = double(0),
                 mu_map   = double(0),
                 sd_map   = double(0),
                 sd_vague = double(0),
                 w_robust = double(0),
                 log      = integer(0, default = 0)) {
    returnType(double(0))
    log_map   <- dnorm(x, mu_map, sd_map, log = TRUE)  + log(1 - w_robust)
    log_vague <- dnorm(x, 0,      sd_vague, log = TRUE) + log(w_robust)
    max_log   <- max(log_map, log_vague)
    logp      <- max_log + log(exp(log_map - max_log) + exp(log_vague - max_log))
    if (log) return(logp)
    return(exp(logp))
  }
)

rRobustNorm <- nimbleFunction(
  run = function(n        = integer(0),
                 mu_map   = double(0),
                 sd_map   = double(0),
                 sd_vague = double(0),
                 w_robust = double(0)) {
    returnType(double(0))
    u <- runif(1, 0, 1)
    if (u < w_robust) {
      return(rnorm(1, 0, sd_vague))
    } else {
      return(rnorm(1, mu_map, sd_map))
    }
  }
)

registerDistributions(list(
  dRobustNorm = list(
    BUGSdist = "dRobustNorm(mu_map, sd_map, sd_vague, w_robust)",
    types    = c("value = double(0)",
                 "mu_map = double(0)", "sd_map = double(0)",
                 "sd_vague = double(0)", "w_robust = double(0)"),
    pqAvail  = FALSE
  )
))

code_robust_map <- nimbleCode({
  
  # ---- IPD block: P13-001 (both ethnicities) + RWE ----
  # All IPD observations use the "current study" ethnic offsets de0, demax.
  # For non-Asian rows (race = 0) the offset drops out; for Asian rows
  # (race = 1) the offset enters with weight 1.
  for (i in 1:N_ipd) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) *
      x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  
  # ---- External aggregate block: Hu and Sng only ----
  # Each external source has its own de0_ext[s] and demax_ext[s],
  # drawn hierarchically from a common population mean.
  for (j in 1:N_agg) {
    y_agg[j] ~ dnorm(mu_agg[j], var = tau2 / n_agg[j])
    
    mu_agg[j] <- e0 + de0_ext[agg_source_id[j]] + a * z_agg[j] +
      (k + dk) * x_agg[j] +
      (emax + demax_ext[agg_source_id[j]]) *
      x_agg[j]^(r1 + dr1) /
      ((ed50 + ded50)^(r1 + dr1) +
         x_agg[j]^(r1 + dr1))
  }
  
  # ---- MAP hierarchy over external sources (Hu, Sng) ----
  # Each external source's ethnic offset is drawn from a common
  # population distribution with mean (de0_pop, demax_pop) and
  # between-source SD tau_map.
  for (s in 1:N_ext_sources) {
    de0_ext[s]   ~ dnorm(de0_pop,   sd = tau_map)
    demax_ext[s] ~ dnorm(demax_pop, sd = tau_map)
  }
  
  # ---- Robust mixture prior on the CURRENT-STUDY ethnic offsets ----
  # The current-study ethnic offsets (de0, demax) receive a mixture of:
  #   - a MAP component N(de0_pop, tau_map^2) with weight (1 - w_robust),
  #     predicting the current-study offset from the Hu/Sng hierarchy;
  #   - a vague component N(0, sd_vague^2) with weight w_robust,
  #     protecting against prior-data conflict.
  # This is the robust MAP construction applied to the de0/demax of
  # the P13-001 + RWE current study.
  de0   ~ dRobustNorm(de0_pop,   tau_map, sd_vague, w_robust)
  demax ~ dRobustNorm(demax_pop, tau_map, sd_vague, w_robust)
  
  # ---- Hyperpriors ----
  de0_pop   ~ dnorm(0, sd = 1)
  demax_pop ~ dnorm(0, sd = 1)
  tau_map   ~ T(dnorm(0, sd = 0.5), 0, )
  
  e0    ~ dnorm(0, sd = sqrt(10))
  emax  ~ dnorm(0, sd = 1)
  k     ~ dnorm(0, sd = 10)
  ed50  ~ T(dnorm(10, sd = sqrt(10)), 0, )
  r1    ~ T(dnorm(0, sd = 10), 0, )
  dk    ~ dnorm(0, sd = sqrt(10))
  ded50 ~ T(dnorm(0, sd = sqrt(10)), -ed50, )
  dr1   ~ T(dnorm(0, sd = 10), -r1, )
  a     ~ dnorm(0, sd = 10)
  tau2  ~ dinvgamma(0.01, 0.01)
})


################################################################################
# 2.  HELPER FUNCTIONS
################################################################################

# --------------------------------------------------------------------------- #
#  2A. Build nimble data/constants/inits for the PRIMARY model
# --------------------------------------------------------------------------- #
build_primary_inputs <- function(dat) {
  constants <- list(N = nrow(dat))
  data      <- list(y    = dat$pdpch,
                    x    = dat$visit,
                    z    = dat$base.pd.cent,
                    race = dat$race,
                    n    = dat$n)
  inits     <- list(e0 = -0.6, de0 = 0, emax = 0.3, demax = 0,
                    ed50 = 10, ded50 = 1, r1 = 1, dr1 = 0.1,
                    k = -0.0003, dk = 0, a = -0.02, tau2 = 0.05)
  list(constants = constants, data = data, inits = inits)
}


# --------------------------------------------------------------------------- #
#  2B. Run MCMC with nimble  (generic wrapper)
#      [FIX 10] Added enableWAIC option for model comparison
#      [v7] Chains now run in parallel on multi-core machines via a PSOCK
#      cluster (parallel::makeCluster + parLapply).  This gives N-fold
#      R-hat reliability and ~N-fold total ESS at approximately the same
#      wall-clock time as a single chain, up to (num_cores) chains.
#      Each worker rebuilds the nimble model independently (nimble's
#      compiled objects cannot be serialized across processes), which
#      adds ~20-40 s of compile overhead per worker per call.
#
#      When enableWAIC = TRUE, the run falls back to sequential execution
#      so that WAIC can be read from a single compiled MCMC object.  WAIC
#      is only used once (primary model for the S8 comparison), so the
#      sequential fall-back is limited in scope.  For all other
#      run_nimble calls, parallel execution is used.
#
#      Mixing improvements (per-iteration, runtime-neutral):
#        * Block sampler on (ed50, ded50, r1, dr1)  -- Emax shape ridge
#        * Block sampler on (e0, de0, emax, demax)  -- Emax level ridge
#        * Slice samplers on r1 and dr1             -- heavy right tail
#      These are applied inside every MCMC run (parallel or sequential).
#
#      Defaults (nBurnin = 10000, nIter = 100000, nThin = 10, nChains = 2)
#      preserve the v6 per-chain runtime; with parallel execution the
#      total runtime should be similar to v6 or slightly less.
# --------------------------------------------------------------------------- #

# Helper: install the block + slice samplers on an mcmcConf object.
# Used by both the sequential and parallel code paths.
#
# [v9] Sampler strategy — three-tier approach for the Emax ridge:
#
#   Tier 1: ONE merged 8-parameter RW_block over all Emax parameters.
#     The previous two separate 4-parameter blocks (level + shape) did
#     not capture the cross-block correlations:
#       emax ↔ ed50  (higher Emax → apparent ED50 shifts)
#       emax ↔ r1    (Emax and Hill coefficient trade off on steep curves)
#       e0   ↔ emax  (intercept and ceiling are negatively correlated)
#     A single 8×8 AM block adapts the full covariance matrix and
#     proposes joint moves that stay on the ridge, dramatically
#     reducing autocorrelation compared to separate 4×4 blocks.
#
#   Tier 2: Slice samplers on r1 and dr1 IN ADDITION to the block.
#     r1 has a heavy right tail (posterior mass extends to ~15+); slice
#     sampling handles this with no tuning.  These run as scalar updates
#     supplementing the block's joint proposals.
#
#   Tier 3: All other parameters (k, dk, a, tau2) keep the default
#     scalar RW sampler — they mix near-perfectly (ESS ~100k).
.install_xen_samplers <- function(mcmcConf, node_tops) {

  # [v9] Single merged block over ALL eight Emax parameters.
  # Replaces the previous two 4-parameter blocks which failed to
  # capture cross-block correlations between level and shape.
  emax_all_block <- c("e0", "de0", "emax", "demax",
                      "ed50", "ded50", "r1", "dr1")
  available_block <- intersect(emax_all_block, node_tops)
  if (length(available_block) >= 2) {
    mcmcConf$removeSamplers(available_block)
    mcmcConf$addSampler(target = available_block, type = "RW_block")
  }

  # Supplementary slice samplers on heavy-tailed shape parameters.
  for (tail_par in c("r1", "dr1")) {
    if (tail_par %in% node_tops) {
      mcmcConf$addSampler(target = tail_par, type = "slice")
    }
  }

  invisible(mcmcConf)
}

# Worker function run on each PSOCK node: builds model, installs samplers,
# runs one chain, returns the samples matrix.  All nimble calls happen
# inside the worker because compiled objects cannot be serialized.
.run_one_chain <- function(chain_id, code, constants, data, inits,
                           monitors, nBurnin, nIter, nThin, seed_base) {
  suppressPackageStartupMessages(library(nimble))
  suppressPackageStartupMessages(library(coda))
  
  set.seed(seed_base + chain_id)
  
  model <- nimbleModel(code, constants = constants, data = data,
                       inits = inits, check = FALSE)
  node_names <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  node_tops  <- unique(gsub("\\[.*", "", node_names))
  monitors   <- intersect(monitors, node_tops)
  
  cModel    <- compileNimble(model)
  mcmcConf  <- configureMCMC(model, monitors = monitors, enableWAIC = FALSE)
  
  # [v9] Install samplers — inlined from .install_xen_samplers because
  # worker processes cannot serialize the helper from the parent env.
  # Single merged 8-parameter block for the full Emax correlation ridge.
  emax_all_block <- c("e0", "de0", "emax", "demax",
                      "ed50", "ded50", "r1", "dr1")
  available_block <- intersect(emax_all_block, node_tops)
  if (length(available_block) >= 2) {
    mcmcConf$removeSamplers(available_block)
    mcmcConf$addSampler(target = available_block, type = "RW_block")
  }
  for (tail_par in c("r1", "dr1")) {
    if (tail_par %in% node_tops) {
      mcmcConf$addSampler(target = tail_par, type = "slice")
    }
  }
  
  mcmc  <- buildMCMC(mcmcConf)
  cMCMC <- compileNimble(mcmc, project = model)
  
  cMCMC$run(nBurnin + nIter * nThin, thin = nThin, nburnin = nBurnin)
  as.matrix(cMCMC$mvSamples)
}

# ─────────────────────────────────────────────────────────────────────────────
# ESS TROUBLESHOOTING LADDER  (activate in order if ESS < 500 persists)
# ─────────────────────────────────────────────────────────────────────────────
#
# [v9 — ACTIVE] Fix 1: Single merged 8-parameter RW_block.
#   Already implemented in .install_xen_samplers above.
#   Try this first: run the script, check ESS.  If still < 500 → Fix 2.
#
# Fix 2: Log reparameterization of ed50, r1, and their ethnic offsets.
#   Why: r1 and ed50 are strictly positive with skewed posteriors.
#   Sampling on the log scale removes the right-skew and eliminates
#   truncated-normal constraints, which helps AM adaptation enormously.
#
#   Replace in code_primary (and all other model blocks):
#
#     # REMOVE these lines:
#     ed50  ~ T(dnorm(10, sd = sqrt(10)), 0, )
#     r1    ~ T(dnorm(0,  sd = 10),       0, )
#     ded50 ~ T(dnorm(0,  sd = sqrt(10)), -ed50, )
#     dr1   ~ T(dnorm(0,  sd = 10),       -r1,   )
#
#     # REPLACE WITH:
#     log_ed50       ~ dnorm(log(10), sd = 1)     # prior on log scale
#     ed50           <- exp(log_ed50)              # deterministic
#     log_ed50_asian ~ dnorm(log(10), sd = 1)
#     ed50_asian     <- exp(log_ed50_asian)
#     ded50          <- ed50_asian - ed50          # ed50+ded50 = exp(...) > 0
#
#     log_r1         ~ dnorm(0, sd = 1)
#     r1             <- exp(log_r1)
#     log_r1_asian   ~ dnorm(0, sd = 1)
#     r1_asian       <- exp(log_r1_asian)
#     dr1            <- r1_asian - r1
#
#   Inits change: replace ed50=10, ded50=1, r1=1, dr1=0.1 with
#     log_ed50=log(10), log_ed50_asian=log(11), log_r1=0, log_r1_asian=log(1.1)
#
#   This change must be applied to ALL model code blocks (primary,
#   emax_only, piecewise, power_prior, commensurate, robust_map, and
#   the vague/informative prior variants).  If ESS still < 500 → Fix 3.
#
# Fix 3: NUTS (No-U-Turn Sampler) via nimbleHMC.
#   NUTS uses gradient information to propose moves that traverse the
#   ridge in one step; it is essentially guaranteed to solve Emax-model
#   mixing regardless of correlation structure.
#
#   One-time setup:
#     install.packages("nimbleHMC")
#
#   In .install_xen_samplers, replace the RW_block block with:
#     library(nimbleHMC)
#     nuts_params <- intersect(c("e0","de0","emax","demax",
#                                "ed50","ded50","r1","dr1","k","dk","a"),
#                              node_tops)
#     mcmcConf$removeSamplers(nuts_params)
#     mcmcConf$addSampler(target = nuts_params, type = "NUTS")
#
#   Note: NUTS is incompatible with truncated priors in nimble unless
#   the log reparameterization (Fix 2) is also applied.  If you install
#   nimbleHMC, apply Fix 2 first, then swap in NUTS.
# ─────────────────────────────────────────────────────────────────────────────

run_nimble <- function(code, constants, data, inits,
                       monitors = NULL,
                       nBurnin = 50000, nIter = 50000, nThin = 1,
                       nChains = 2, setSeed = 1223456,
                       enableWAIC = FALSE,
                       parallel = TRUE,
                       n_workers = NULL,
                       fast = FALSE) {
  # [v8] Defaults: 50k burnin (5x longer for better block-sampler adaptation)
  #               + 50k post-burnin with nThin=1 (keep every draw).
  #               Total sweeps/chain = 100k — identical to prior setting of
  #               nBurnin=10k + nIter=10k*nThin=10 = 110k, but the block
  #               RW_block adapts far better with 50k burnin, lifting ESS
  #               from ~100 to ~1000+ for the Emax-level parameters.
  # [v8] fast=TRUE: 40k sweeps/chain — 2.75x fewer than old 110k setting.
  #               Used for sensitivity analyses (LOSO, tipping, power prior,
  #               prior sensitivity, commensurate, MAP, S8 structural models)
  #               where stable mean/CrI estimates are needed but not the
  #               same ESS as the primary.
  if (fast) { nBurnin <- 20000; nIter <- 20000; nThin <- 1 }
  
  if (is.null(monitors)) {
    monitors <- c("e0", "de0", "emax", "demax", "ed50", "ded50",
                  "r1", "dr1", "k", "dk", "a", "tau2")
  }
  
  # WAIC currently requires sequential execution because it is read from
  # a single compiled MCMC object.  Force sequential in that case.
  use_parallel <- parallel && !enableWAIC && nChains >= 2
  
  # Determine effective worker count (cap at nChains and at available cores - 1)
  if (use_parallel) {
    if (is.null(n_workers)) {
      n_workers <- max(1, min(nChains, parallel::detectCores(logical = FALSE) - 1))
    }
    use_parallel <- n_workers >= 2
  }
  
  # ------------------------------------------------------------------- #
  # PARALLEL PATH  (Windows-compatible PSOCK cluster)
  #
  # On machines where Windows Application Control / AppLocker / WDAC
  # blocks the loading of nimble's compiled DLLs in worker processes,
  # this path will fail with an "unable to load shared object" error.
  # We catch that error, print a one-time warning, and fall back to
  # sequential execution for the remainder of the call.  The fall-back
  # is silent for subsequent calls within the same session because
  # options("xen_parallel_blocked") gets set.
  # ------------------------------------------------------------------- #
  if (use_parallel && !isTRUE(getOption("xen_parallel_blocked"))) {
    
    cat(sprintf("  [parallel] running %d chains on %d workers\n",
                nChains, n_workers))
    
    parallel_result <- tryCatch({
      
      cl <- parallel::makeCluster(n_workers, type = "PSOCK")
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      
      parallel::clusterExport(
        cl,
        varlist = c("code", "constants", "data", "inits",
                    "monitors", "nBurnin", "nIter", "nThin", "setSeed",
                    ".run_one_chain"),
        envir = environment()
      )
      
      chain_results <- parallel::parLapply(
        cl, seq_len(nChains),
        function(ch) {
          .run_one_chain(
            chain_id  = ch,
            code      = code,
            constants = constants,
            data      = data,
            inits     = inits,
            monitors  = monitors,
            nBurnin   = nBurnin,
            nIter     = nIter,
            nThin     = nThin,
            seed_base = setSeed
          )
        }
      )
      
      samples_list <- lapply(chain_results, function(m) coda::as.mcmc(m))
      list(ok = TRUE, mcmc = coda::as.mcmc.list(samples_list))
      
    }, error = function(e) {
      list(ok = FALSE, error_message = conditionMessage(e))
    })
    
    if (parallel_result$ok) {
      return(parallel_result$mcmc)
    }
    
    # Parallel failed -- detect the specific Windows policy block and
    # give a tailored warning, otherwise emit a generic one.  In either
    # case fall through to the sequential path.
    msg <- parallel_result$error_message
    if (grepl("Application Control|LoadLibrary|unable to load shared",
              msg, ignore.case = TRUE)) {
      message(
        "[parallel] Windows Application Control blocked the nimble DLL ",
        "in worker processes. Falling back to sequential execution for ",
        "this call and all subsequent run_nimble() calls in this session.\n",
        "  To run chains in parallel, ask IT to whitelist the nimble ",
        "compile directory (typically under ",
        "C:/Users/<username>/AppData/Local/Temp or the working directory) ",
        "or run R from a directory covered by an Application Control ",
        "exception.\n  Original error: ", msg
      )
    } else {
      message(
        "[parallel] Parallel execution failed; falling back to sequential. ",
        "Error: ", msg
      )
    }
    options(xen_parallel_blocked = TRUE)
  }
  
  # ------------------------------------------------------------------- #
  # SEQUENTIAL PATH  (used when parallel = FALSE (default),
  # enableWAIC = TRUE, nChains = 1, only 1 worker available, or parallel
  # execution failed above).  Produces the same output plus WAIC.
  # ------------------------------------------------------------------- #
  model <- nimbleModel(code, constants = constants, data = data,
                       inits = inits, check = FALSE)
  
  node_names <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  node_tops  <- unique(gsub("\\[.*", "", node_names))
  monitors   <- intersect(monitors, node_tops)
  
  cModel   <- compileNimble(model)
  mcmcConf <- configureMCMC(model, monitors = monitors,
                            enableWAIC = enableWAIC)
  .install_xen_samplers(mcmcConf, node_tops)
  
  mcmc  <- buildMCMC(mcmcConf)
  cMCMC <- compileNimble(mcmc, project = model)
  
  samples_list <- list()
  for (ch in seq_len(nChains)) {
    set.seed(setSeed + ch)
    cMCMC$run(nBurnin + nIter * nThin, thin = nThin, nburnin = nBurnin)
    samples_list[[ch]] <- coda::as.mcmc(as.matrix(cMCMC$mvSamples))
  }
  mcmc_out <- coda::as.mcmc.list(samples_list)
  
  # Extract WAIC if enabled (sequential path only)
  if (enableWAIC) {
    tryCatch({
      waic_val <- cMCMC$getWAIC()
      attr(mcmc_out, "WAIC") <- waic_val
    }, error = function(e) {
      message("WAIC extraction failed: ", e$message)
    })
  }
  
  return(mcmc_out)
}


# --------------------------------------------------------------------------- #
#  2C. Summarise MCMC output
# --------------------------------------------------------------------------- #
summarise_mcmc <- function(mcmc_out, params = NULL) {
  combined <- do.call(rbind, mcmc_out)
  if (!is.null(params)) {
    params <- intersect(params, colnames(combined))
    combined <- combined[, params, drop = FALSE]
  }
  t(apply(combined, 2, function(x) {
    c(Mean   = mean(x),
      Median = median(x),
      SD     = sd(x),
      `2.5%`  = unname(quantile(x, 0.025)),
      `97.5%` = unname(quantile(x, 0.975)))
  }))
}

print_rounded <- function(df, digits = 4) {
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- round(df[num_cols], digits)
  print(df)
}

# --------------------------------------------------------------------------- #
#  2D. Compute posterior fitted curves
# --------------------------------------------------------------------------- #
compute_fitted_curves <- function(mcmc_out, x_days,
                                  base_iop_asian, base_iop_nonasian) {
  samps <- as.data.frame(do.call(rbind, mcmc_out))
  nSamp <- nrow(samps)
  
  # Non-Asian curve (race = 0)
  fitted_nonasian <- sapply(1:length(x_days), function(j) {
    samps$e0 + samps$a * base_iop_nonasian +
      samps$k * x_days[j] +
      samps$emax * x_days[j]^samps$r1 /
      (samps$ed50^samps$r1 + x_days[j]^samps$r1)
  })
  
  # Asian curve (race = 1)
  fitted_asian <- sapply(1:length(x_days), function(j) {
    e0a   <- samps$e0 + samps$de0
    emaxa <- samps$emax + samps$demax
    r1a   <- samps$r1 + samps$dr1
    ed50a <- samps$ed50 + samps$ded50
    ka    <- samps$k + samps$dk
    e0a + samps$a * base_iop_asian +
      ka * x_days[j] +
      emaxa * x_days[j]^r1a / (ed50a^r1a + x_days[j]^r1a)
  })
  
  list(nonasian = fitted_nonasian, asian = fitted_asian)
}


# --------------------------------------------------------------------------- #
#  2E. Quick summary of fitted curves (mean + 95% CrI)
# --------------------------------------------------------------------------- #
curve_summary <- function(fitted_mat, x_days) {
  data.frame(
    day   = x_days,
    mean  = colMeans(fitted_mat, na.rm = TRUE),
    lower = apply(fitted_mat, 2, quantile, 0.025, na.rm = TRUE),
    upper = apply(fitted_mat, 2, quantile, 0.975, na.rm = TRUE)
  )
}


# --------------------------------------------------------------------------- #
#  2F. Weighted baseline IOP by sample size.
#      For hybrid IPD + AD data, the mean baseline IOP for each ethnic group
#      must be weighted by the per-observation sample size n_i (1 for IPD,
#      the number of patients at that visit for AD).  The unweighted mean
#      over rows would underweight IPD contributions and distort the
#      population-mean centering used in compute_fitted_curves().
#      This replaces the earlier unweighted mean(base.pd.cent[...]) calls.
# --------------------------------------------------------------------------- #
weighted_base_iop <- function(dat) {
  asian    <- dat[dat$race == 1, ]
  nonasian <- dat[dat$race == 0, ]
  list(
    asian    = sum(asian$base.pd.cent    * asian$n)    / sum(asian$n),
    nonasian = sum(nonasian$base.pd.cent * nonasian$n) / sum(nonasian$n)
  )
}




################################################################################
# 2G.  ANALYSIS-LEVEL PARALLEL INFRASTRUCTURE  [Fix 3]
#
# These top-level worker functions are suitable for export to PSOCK cluster
# workers.  Each takes only the "varying" argument; all other objects
# (y.iop, code_primary, etc.) are exported to the worker environment via
# clusterExport before dispatch.
#
# .run_analyses_parallel() is a drop-in replacement for lapply() that tries
# a PSOCK cluster first and falls back to sequential on failure, using the
# same error-detection pattern as run_nimble()'s parallel path.
################################################################################

# ── Worker: one LOSO fit ────────────────────────────────────────────────────
.loso_worker <- function(src) {
  suppressPackageStartupMessages(library(nimble))
  suppressPackageStartupMessages(library(coda))
  dat_reduced <- y.iop[y.iop$source != src, ]
  inp      <- build_primary_inputs(dat_reduced)
  base_red <- weighted_base_iop(dat_reduced)
  mcmc     <- run_nimble(code = code_primary,
                         constants = inp$constants,
                         data      = inp$data,
                         inits     = inp$inits,
                         fast      = TRUE,    # [v8] 40k sweeps sufficient for LOSO
                         parallel  = FALSE)   # no nested parallelism
  list(
    summary           = summarise_mcmc(mcmc, delta_params),
    curves            = compute_fitted_curves(mcmc, x_days,
                                              base_red$asian, base_red$nonasian),
    base_iop_asian    = base_red$asian,
    base_iop_nonasian = base_red$nonasian
  )
}

# ── Worker: one tipping-point shift ─────────────────────────────────────────
.tipping_worker <- function(delta_shift) {
  suppressPackageStartupMessages(library(nimble))
  suppressPackageStartupMessages(library(coda))
  dat_s <- y.iop
  dat_s$pdpch[dat_s$source == "Pub1"] <-
    dat_s$pdpch[dat_s$source == "Pub1"] + delta_shift
  inp  <- build_primary_inputs(dat_s)
  mcmc <- run_nimble(code = code_primary,
                     constants = inp$constants,
                     data      = inp$data,
                     inits     = inp$inits,
                     fast      = TRUE,    # [v8] 40k sweeps sufficient for tipping
                     parallel  = FALSE)
  smry     <- summarise_mcmc(mcmc, delta_check_params)
  any_excl <- any(smry[,"2.5%"] > 0 | smry[,"97.5%"] < 0)
  tipped   <- rownames(smry)[smry[,"2.5%"] > 0 | smry[,"97.5%"] < 0]
  crv      <- compute_fitted_curves(mcmc, 336,
                                    base_iop_asian, base_iop_nonasian)
  data.frame(
    shift             = delta_shift,
    de0_mean          = smry["de0","Mean"],
    de0_lower         = smry["de0","2.5%"],
    de0_upper         = smry["de0","97.5%"],
    demax_mean        = smry["demax","Mean"],
    demax_lower       = smry["demax","2.5%"],
    demax_upper       = smry["demax","97.5%"],
    db_mean           = smry["dk","Mean"],
    db_lower          = smry["dk","2.5%"],
    db_upper          = smry["dk","97.5%"],
    any_CrI_excl_zero = any_excl,
    tipped_params     = ifelse(length(tipped) > 0,
                               paste(tipped, collapse=","), "none"),
    curve_diff_M12    = mean(crv$asian[,1] - crv$nonasian[,1]),
    stringsAsFactors  = FALSE
  )
}

# ── Worker: one power-prior fit ──────────────────────────────────────────────
.pp_worker <- function(a0) {
  suppressPackageStartupMessages(library(nimble))
  suppressPackageStartupMessages(library(coda))
  pp_c <- list(N_ipd = nrow(dat_ipd), N_rwe = nrow(dat_rwe),
               N_disc = nrow(dat_disc), alpha0 = a0)
  pp_d <- list(y         = dat_ipd$pdpch,  x         = dat_ipd$visit,
               z         = dat_ipd$base.pd.cent, race  = dat_ipd$race,
               y_rwe     = dat_rwe$pdpch,  x_rwe     = dat_rwe$visit,
               z_rwe     = dat_rwe$base.pd.cent, race_rwe = dat_rwe$race,
               n_rwe     = dat_rwe$n,
               y_disc    = dat_disc$pdpch, x_disc    = dat_disc$visit,
               z_disc    = dat_disc$base.pd.cent, race_disc = dat_disc$race,
               n_disc    = dat_disc$n)
  pp_i <- list(e0=-0.6, de0=0, emax=0.3, demax=0, ed50=10, ded50=1,
               r1=1, dr1=0.1, k=-0.0003, dk=0, a=-0.02, tau2=0.05)
  mcmc <- run_nimble(code = code_power_prior,
                     constants = pp_c, data = pp_d, inits = pp_i,
                     fast      = TRUE,    # [v8] 40k sweeps sufficient for power prior
                     parallel  = FALSE)
  summarise_mcmc(mcmc, delta_params)
}

# ── Utility: try PSOCK cluster, fall back to sequential ─────────────────────
.run_analyses_parallel <- function(worker_fn, items, export_names,
                                   envir = parent.frame()) {
  n_cores <- max(1L, min(parallel::detectCores(logical = FALSE) - 1L,
                         length(items)))
  
  if (n_cores < 2L || isTRUE(getOption("xen_anal_parallel_blocked"))) {
    cat(sprintf("  [sequential] running %d analyses\n", length(items)))
    return(lapply(items, worker_fn))
  }
  
  cat(sprintf("  [parallel-anal] dispatching %d analyses on %d cores\n",
              length(items), n_cores))
  
  result <- tryCatch({
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterExport(cl, export_names, envir = envir)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(nimble))
      suppressPackageStartupMessages(library(coda))
    })
    list(ok = TRUE, res = parallel::parLapply(cl, items, worker_fn))
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  
  if (result$ok) return(result$res)
  
  if (grepl("Application Control|LoadLibrary|unable to load shared",
            result$msg, ignore.case = TRUE)) {
    message(
      "[parallel-anal] Windows Application Control blocked nimble DLL in ",
      "worker processes. Falling back to sequential for all analysis-level ",
      "parallel calls in this session.\n  Original error: ", result$msg
    )
    options(xen_anal_parallel_blocked = TRUE)
  } else {
    message("[parallel-anal] Parallel failed; falling back to sequential. ",
            "Error: ", result$msg)
  }
  lapply(items, worker_fn)
}

################################################################################
# 3.  PRIMARY ANALYSIS
################################################################################

cat("\n========== PRIMARY ANALYSIS ==========\n")
inputs <- build_primary_inputs(y.iop)

# [Fix 4] Run primary with parallel chains; WAIC is collected separately
# at the S8 stage (short sequential single-chain run).  Decoupling WAIC
# from the primary allows parallel chains here, saving ~half the primary
# wall-clock time when 2+ cores are available.
mcmc_primary <- run_nimble(
  code       = code_primary,
  constants  = inputs$constants,
  data       = inputs$data,
  inits      = inputs$inits,
  enableWAIC = FALSE,   # [Fix 4] allow parallel chains
  parallel   = TRUE
)

delta_params <- c("e0", "de0", "emax", "demax", "ed50", "ded50",
                  "r1", "dr1", "k", "dk", "a", "tau2")
primary_summary <- summarise_mcmc(mcmc_primary, delta_params)
cat("\nPrimary Analysis — Posterior Summary (cf. Report Table 5):\n")
print_rounded(primary_summary, 4)

df_to_latex(primary_summary, rownames=TRUE, file = "Primary_analysis_table.tex")

# [v8] ESS check — warn early if any parameter is below 500.
# If ESS is still low, increase nBurnin (better adaptation) rather than
# nIter; the block sampler needs sufficient burnin to tune its proposal
# covariance for the e0/emax/de0/demax correlation ridge.
cat("\n[v8] ESS check (primary):\n")
ess_primary <- coda::effectiveSize(mcmc_primary)
print(round(ess_primary, 0))
low_ess <- names(ess_primary)[ess_primary < 500]
if (length(low_ess) > 0) {
  warning(
    "[v8] ESS < 500 for: ", paste(low_ess, collapse = ", "), "\n",
    "  Recommended fix: increase nBurnin further (try 100000) so the\n",
    "  RW_block sampler has more time to tune its proposal covariance.\n",
    "  Do NOT just add more post-burnin iterations — that helps less\n",
    "  than better adaptation when the chain is mixing poorly."
  )
} else {
  cat("  All parameters: ESS >= 500. Convergence criteria met.\n")
}

save.image("Sensitivity_results.RData")


x_days <- c(1, 7, 14, 28, 84, 168, 224, 280, 336)
base_iop_full     <- weighted_base_iop(y.iop)   # [UPDATED] weighted by n
base_iop_asian    <- base_iop_full$asian
base_iop_nonasian <- base_iop_full$nonasian

curves_primary <- compute_fitted_curves(mcmc_primary, x_days,
                                        base_iop_asian, base_iop_nonasian)

cat("\nFitted Asian curve (mean, 95% CrI):\n")
print_rounded(curve_summary(curves_primary$asian, x_days), 4)
cat("\nFitted Non-Asian curve (mean, 95% CrI):\n")
print_rounded(curve_summary(curves_primary$nonasian, x_days), 4)

# Report WAIC
if (!is.null(attr(mcmc_primary, "WAIC"))) {
  cat("\nPrimary Model WAIC:", attr(mcmc_primary, "WAIC")$WAIC, "\n")
}


################################################################################
# 3b. FITTED CURVES PLOTS  —  Asian vs Non-Asian, with and without patient 102009
#
# Produces two comparison figures matching the layout of manuscript
# Figures \ref{fig:comparison} (all RWE data) and
# \ref{fig:comparison_no102009} (excluding patient 102009).
#
# The primary analysis fit above uses the "no102009" dataset.  For the
# "all RWE" plot we refit the primary model on a dataset that includes
# patient 102009.  That dataset is expected to be available as:
#     XEN_nimble_ready_data_v2_allRWE_simRWE.csv
# If your file is named differently, adjust `data_file_full` below.
# If no such file is available, the "all RWE" plot is skipped with a
# message and only the no-102009 plot is produced.
################################################################################

# ---- Helper: plot Asian vs Non-Asian fitted curves with 95% CrI bands ----
plot_fitted_curves <- function(curves, x_days, title_text,
                               out_file = NULL, width = 7, height = 5.5) {
  
  # Build tidy data frame with mean and 95% CrI at each day for each group
  summarize_curve <- function(mat, group_label) {
    data.frame(
      day   = x_days,
      mean  = colMeans(mat,  na.rm = TRUE) * 100,
      lower = apply(mat, 2, quantile, 0.025, na.rm = TRUE) * 100,
      upper = apply(mat, 2, quantile, 0.975, na.rm = TRUE) * 100,
      Data  = group_label
    )
  }
  df_plot <- rbind(
    summarize_curve(curves$asian,    "Asian.Fitted"),
    summarize_curve(curves$nonasian, "Non-Asian.Fitted")
  )
  
  p <- ggplot(df_plot, aes(x = day, y = mean, colour = Data,
                           fill = Data, linetype = Data)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.35, colour = NA) +
    geom_line(size = 0.8) +
    scale_colour_manual(values = c("Asian.Fitted" = "#E46464",
                                   "Non-Asian.Fitted" = "#2BB3C0")) +
    scale_fill_manual(values   = c("Asian.Fitted" = "#E46464",
                                   "Non-Asian.Fitted" = "#2BB3C0")) +
    scale_linetype_manual(values = c("Asian.Fitted" = "solid",
                                     "Non-Asian.Fitted" = "dashed")) +
    labs(
      title = title_text,
      x = "Time (day)",
      y = "PD Percentage Change from Baseline (%)"
    ) +
    theme_grey(base_size = 12) +
    theme(
      legend.position      = c(0.85, 0.18),
      legend.background    = element_rect(fill = "white", colour = "grey70"),
      legend.title         = element_text(size = 11),
      plot.title           = element_text(size = 13, face = "plain")
    )
  
  if (!is.null(out_file)) {
    ggsave(out_file, plot = p, width = width, height = height, dpi = 300)
    cat("Figure saved:", out_file, "\n")
  }
  
  print(p)
  invisible(p)
}

# ---- Figure 2: fit EXCLUDING patient 102009 (current primary) ----
cat("\n========== PRIMARY FITTED CURVES (no 102009) ==========\n")
plot_fitted_curves(
  curves     = curves_primary,
  x_days     = x_days,
  title_text = "Model Fitted PD Percentage Change from Baseline",
  out_file   = "comparison_Asian_non_Asian.png"
)




# ---- Figure 1: fit INCLUDING patient 102009 ----
# Attempt to load the all-RWE dataset.  If not available, skip this plot.
#data_file_full <- "XEN_nimble_ready_data_v2_allRWE_simRWE.csv"
data_file_full <- "XEN_nimble_realRWE_with102009.csv"
#data_file_full <- "XEN_nimble_realRWE_with102009_v2.csv"

if (file.exists(data_file_full)) {
  
  cat("\n========== PRIMARY FIT — INCLUDING PATIENT 102009 ==========\n")
  
  y.iop.full <- read.csv(data_file_full)
  stopifnot(all(c("pid", "study", "visit", "n", "race", "base.pd",
                  "pd", "pdch", "pdpch", "source", "base.pd.cent")
                %in% names(y.iop.full)))
  n_before_full <- nrow(y.iop.full)
  y.iop.full <- y.iop.full[!is.na(y.iop.full$pdpch), ]
  cat("  Full-RWE dataset:", n_before_full, "rows,",
      n_before_full - nrow(y.iop.full), "removed for missing pdpch,",
      nrow(y.iop.full), "retained\n")
  
  inputs_full <- build_primary_inputs(y.iop.full)
  
  mcmc_primary_full <- run_nimble(
    code      = code_primary,
    constants = inputs_full$constants,
    data      = inputs_full$data,
    inits     = inputs_full$inits
  )
  
  primary_summary_full <- summarise_mcmc(mcmc_primary_full, delta_params)
  cat("\nPrimary Analysis (including 102009) — Posterior Summary:\n")
  print_rounded(primary_summary_full, 4)
  df_to_latex(primary_summary_full, rownames = TRUE,
              file = "Primary_analysis_table_allRWE.tex")
  
  base_iop_full_all    <- weighted_base_iop(y.iop.full)
  curves_primary_full  <- compute_fitted_curves(
    mcmc_primary_full, x_days,
    base_iop_full_all$asian, base_iop_full_all$nonasian
  )
  
  plot_fitted_curves(
    curves     = curves_primary_full,
    x_days     = x_days,
    title_text = "Model Fitted PD Percentage Change from Baseline",
    out_file   = "comparison.png"
  )
  
} else {
  cat("\n[Skipped] Fit including patient 102009: file", data_file_full,
      "not found.\n",
      "  To produce Figure \\ref{fig:comparison}, provide a dataset that\n",
      "  includes patient 102009 under that filename and rerun this section.\n")
}


################################################################################
# 4.  PILLAR 1 — APPROPRIATENESS
################################################################################

# ============================================================================ #
#  S1: Heterogeneity Diagnostics — Posterior Predictive Checks
#
#  Produces the visit-level PPC summary corresponding to manuscript
#  Table \ref{tab:ppc}.  For each (source, visit) cell the output table
#  reports:
#     N                - number of contributing observations at the visit
#     obs_mean         - observed mean (IPD) or digitised value (AD)
#     predicted        - posterior predictive mean of the fitted curve
#     residual         - obs_mean - predicted
#     se_mean          - SE of the source-level mean
#                        IPD:  SD of observed values / sqrt(N)
#                        AD:   NA (n = 1 summary per cell)
#     z                - standardised residual = residual / se_mean
#                        (NA for AD; NA for IPD cells with N = 1)
#     p_value          - MCMC p-value:  P( mu_draw <= obs_mean )
#
#  For IPD sources (P13_NonAsian, P13_Asian, RWE), rows at the same
#  visit are aggregated across patients.  For AD sources (Hu, Sng),
#  each row is already a visit-level summary and passes through
#  directly, with se_mean = NA and z = NA.
#
#  [FIX 8] All five source labels are included to give a complete
#  picture of model fit, including the dominant P13_NonAsian source.
# ============================================================================ #

cat("\n========== S1: POSTERIOR PREDICTIVE CHECKS ==========\n")

ppc_check <- function(mcmc_out, y.iop, source_label,
                      base_iop_asian, base_iop_nonasian) {
  
  src_dat   <- y.iop[y.iop$source == source_label, ]
  if (nrow(src_dat) == 0) return(NULL)
  
  is_asian  <- src_dat$race[1] == 1
  is_ad     <- source_label %in% c("Pub1", "Pub2")   # aggregate (n > 1) sources
  
  # --- Aggregate to visit level ---
  if (!is_ad) {
    # IPD: aggregate across patients within each visit
    agg <- aggregate(pdpch ~ visit, data = src_dat,
                     FUN = function(x) c(N = length(x),
                                         mean = mean(x),
                                         sd   = sd(x)))
    visit_tab <- data.frame(
      visit    = agg$visit,
      N        = agg$pdpch[, "N"],
      obs_mean = agg$pdpch[, "mean"],
      obs_sd   = agg$pdpch[, "sd"]
    )
    # SE of the visit-level mean: SD / sqrt(N) (NA when N <= 1)
    visit_tab$se_mean <- ifelse(visit_tab$N > 1,
                                visit_tab$obs_sd / sqrt(visit_tab$N),
                                NA_real_)
  } else {
    # AD: each row already represents a visit-level summary (n = 1 in model,
    # but n_i column stores the number of patients at that visit).
    visit_tab <- data.frame(
      visit    = src_dat$visit,
      N        = src_dat$n,                 # patients contributing to the AD mean
      obs_mean = src_dat$pdpch,
      obs_sd   = NA_real_,
      se_mean  = NA_real_                   # no patient-level SD available
    )
  }
  
  # --- Predicted fitted curve at the visits that appear in this source ---
  curves <- compute_fitted_curves(mcmc_out, visit_tab$visit,
                                  base_iop_asian, base_iop_nonasian)
  pred <- if (is_asian) curves$asian else curves$nonasian
  
  visit_tab$predicted <- colMeans(pred, na.rm = TRUE)
  visit_tab$residual  <- visit_tab$obs_mean - visit_tab$predicted
  visit_tab$z         <- visit_tab$residual / visit_tab$se_mean
  visit_tab$p_value   <- sapply(seq_along(visit_tab$visit), function(j) {
    mean(pred[, j] <= visit_tab$obs_mean[j], na.rm = TRUE)
  })
  
  data.frame(
    source    = source_label,
    day       = visit_tab$visit,
    N         = visit_tab$N,
    obs_mean  = visit_tab$obs_mean,
    predicted = visit_tab$predicted,
    residual  = visit_tab$residual,
    se_mean   = visit_tab$se_mean,
    z         = visit_tab$z,
    p_value   = visit_tab$p_value,
    stringsAsFactors = FALSE
  )
}

# [FIX 8] Include P13_NonAsian for complete model-fit assessment.
ppc_sources <- c("Global_NonAsian", "Global_Asian", "RWE", "Pub1", "Pub2")
ppc_results <- lapply(ppc_sources, function(s) {
  tryCatch(
    ppc_check(mcmc_primary, y.iop, s,
              base_iop_asian, base_iop_nonasian),
    error = function(e) {
      message("PPC skipped for ", s, ": ", e$message)
      NULL
    }
  )
})
ppc_results <- do.call(rbind, Filter(Negate(is.null), ppc_results))

# --- Console print: visit-level table matching manuscript Table \ref{tab:ppc}
cat("\nS1 Posterior Predictive Check Summary (visit-level):\n")
ppc_display <- ppc_results
num_cols <- setdiff(names(ppc_display), c("source", "day", "N"))
ppc_display[num_cols] <- lapply(ppc_display[num_cols], round, 4)
print(ppc_display, row.names = FALSE)

# --- LaTeX export for the manuscript
df_to_latex(ppc_results, rownames = FALSE, file = "S1_ppc_table.tex")

# --- Quick diagnostic summary ---
cat("\nS1 diagnostic flags:\n")
big_z <- ppc_results[!is.na(ppc_results$z) & abs(ppc_results$z) > 2,
                     c("source", "day", "N", "residual", "se_mean", "z")]
if (nrow(big_z) > 0) {
  cat("  Cells with |z| > 2 (inspect for potential misfit or small-n artifact):\n")
  print(big_z, row.names = FALSE)
} else {
  cat("  No cells with |z| > 2.\n")
}

# Per-source mean absolute residual (useful supplement to z-based flags)
mar <- aggregate(abs(residual) ~ source, data = ppc_results, FUN = mean)
names(mar)[2] <- "mean_abs_residual"
cat("\nMean absolute residual by source:\n")
print(mar, row.names = FALSE)

save.image("Sensitivity_results.RData")

# ============================================================================ #
#  A2: Leave-One-Source-Out
#      [FIX 3] Recompute baseline IOP for each reduced dataset
# ============================================================================ #

cat("\n========== A2: LEAVE-ONE-SOURCE-OUT ==========\n")

loso_sources <- c("Global_Asian", "RWE", "Pub1", "Pub2")
loso_results <- list()

# [Fix 3] Dispatch LOSO fits in parallel.  Each worker needs:
#   y.iop, code_primary, delta_params, x_days  (data + model)
#   run_nimble, build_primary_inputs, summarise_mcmc,
#   compute_fitted_curves, weighted_base_iop    (helper functions)
#   .run_one_chain, .install_xen_samplers       (nimble internals)
loso_export <- c("y.iop", "code_primary", "delta_params", "x_days",
                 "run_nimble", "build_primary_inputs", "summarise_mcmc",
                 "compute_fitted_curves", "weighted_base_iop",
                 ".run_one_chain", ".install_xen_samplers")

loso_raw <- .run_analyses_parallel(.loso_worker, loso_sources,
                                   loso_export, envir = environment())
loso_results <- setNames(loso_raw, loso_sources)

cat("\nLeave-One-Source-Out — Posterior Means for Delta Parameters:\n")
loso_table <- data.frame(
  Dropped = c("None", loso_sources),
  dE0     = c(primary_summary["de0", "Mean"],
              sapply(loso_results, function(r) r$summary["de0", "Mean"])),
  dE0_95CI = c(
    paste0("[", round(primary_summary["de0","2.5%"],3), ", ",
           round(primary_summary["de0","97.5%"],3), "]"),
    sapply(loso_results, function(r)
      paste0("[", round(r$summary["de0","2.5%"],3), ", ",
             round(r$summary["de0","97.5%"],3), "]"))
  ),
  dEmax   = c(primary_summary["demax", "Mean"],
              sapply(loso_results, function(r) r$summary["demax", "Mean"])),
  dEmax_95CI = c(
    paste0("[", round(primary_summary["demax","2.5%"],3), ", ",
           round(primary_summary["demax","97.5%"],3), "]"),
    sapply(loso_results, function(r)
      paste0("[", round(r$summary["demax","2.5%"],3), ", ",
             round(r$summary["demax","97.5%"],3), "]"))
  ),
  dED50   = c(primary_summary["ded50", "Mean"],
              sapply(loso_results, function(r) r$summary["ded50", "Mean"])),
  dED50_95CI = c(
    paste0("[", round(primary_summary["ded50","2.5%"],3), ", ",
           round(primary_summary["ded50","97.5%"],3), "]"),
    sapply(loso_results, function(r)
      paste0("[", round(r$summary["ded50","2.5%"],3), ", ",
             round(r$summary["ded50","97.5%"],3), "]"))
  )
)
print_rounded(loso_table[, c("Dropped","dE0","dEmax","dED50")], 4)
cat("\nFull table with 95% CrI:\n")
print(loso_table)
df_to_latex(loso_table, rownames=T, file='Leave_one_source_out.tex')


save.image("Sensitivity_results.RData")

################################################################################
# 5.  PILLAR 2 — VALUE OF BORROWING
################################################################################

# ============================================================================ #
#  A3: No-Borrowing Reference (IPD only)
#      Note: RWE is aggregate in the current data. When RWE IPD becomes
#      available, update the filter to include it.
# ============================================================================ #

cat("\n========== A3: NO-BORROWING REFERENCE ==========\n")

# Keep only IPD (n == 1): P13-001 all patients (non-Asian + 3 Asian)
dat_ipd_only <- y.iop[y.iop$n == 1, ]
cat("  No-borrowing dataset: ", nrow(dat_ipd_only), " rows (",
    sum(dat_ipd_only$race == 1), " Asian, ",
    sum(dat_ipd_only$race == 0), " non-Asian)\n")

inp_nb <- build_primary_inputs(dat_ipd_only)

mcmc_noborrow <- run_nimble(
  code      = code_primary,
  constants = inp_nb$constants,
  data      = inp_nb$data,
  inits     = inp_nb$inits
  # [v8] Uses default settings (50k burnin + 50k iter, nThin=1): the
  # no-borrow posterior variance is the ESS denominator in A4, so this
  # fit needs the same accuracy as the primary.
)

# [UPDATED] Recompute baseline IOP for IPD-only using weighted mean
base_nb <- weighted_base_iop(dat_ipd_only)
base_iop_asian_nb    <- base_nb$asian
base_iop_nonasian_nb <- base_nb$nonasian

nb_summary <- summarise_mcmc(mcmc_noborrow, delta_params)
curves_nb  <- compute_fitted_curves(mcmc_noborrow, x_days,
                                    base_iop_asian_nb, base_iop_nonasian_nb)

cat("\nNo-Borrowing — Posterior Summary:\n")
print_rounded(nb_summary, 4)
df_to_latex(nb_summary, rownames=T, file='No_borrow_summary.tex')

# Compare credible interval widths
ci_comparison <- data.frame(
  Analysis   = c("Full Borrowing", "No Borrowing"),
  Asian_M6_width = c(
    diff(quantile(curves_primary$asian[, which(x_days == 168)], c(0.025, 0.975))),
    diff(quantile(curves_nb$asian[, which(x_days == 168)], c(0.025, 0.975)))
  ),
  Asian_M12_width = c(
    diff(quantile(curves_primary$asian[, which(x_days == 336)], c(0.025, 0.975))),
    diff(quantile(curves_nb$asian[, which(x_days == 336)], c(0.025, 0.975)))
  )
)
cat("\n95% CrI Width Comparison for Asian Curve:\n")
print_rounded(ci_comparison, 4)
df_to_latex(ci_comparison, rownames=T, file='No_borrow_CI_comparison.tex')


# ============================================================================ #
#  A4: Effective Sample Size Decomposition
#      [FIX 7] Replaced ad-hoc ESS with Morita et al. (2008) style
#      variance-ratio approach:
#        ESS_source = (1/Var_without - 1/Var_full) / (1/Var_without)
#                     * total_info_units
#      where total_info_units is the effective total precision denominator.
#      The simpler and more interpretable version:
#        ESS_source ≈ Var_without / Var_full - 1  (in units of "multiples
#        of the IPD-only information")
#      We report both the precision contribution fraction and a
#      calibrated ESS using individual patient variance.
# ============================================================================ #

save.image("Sensitivity_results.RData")

cat("\n========== A4: EFFECTIVE SAMPLE SIZE ==========\n")

# Reference: posterior variance from the no-borrowing (IPD-only) analysis (S3)
# This represents the information from the Asian IPD patients
# (3 P13-001 + 23 simulated RWE = 26 patients).
n_asian_ipd_rows <- sum(dat_ipd_only$race == 1)
n_asian_patients <- length(unique(dat_ipd_only$pid[dat_ipd_only$race == 1]))
cat("  Number of Asian IPD patients:", n_asian_patients,
    "(", n_asian_ipd_rows, "observations)\n")

ess_results <- data.frame(
  Source = c("Pub1", "Pub2", "RWE"),
  N_rows = c(sum(y.iop$source == "Pub1"),
             sum(y.iop$source == "Pub2"),
             sum(y.iop$source == "RWE")),
  N_patients_reported = c(100, 80, 60)   # Pub1, Pub2, RWE
)

# Compute ESS at Month 6 (Day 168) and Month 12 (Day 336)
# Strategy: the per-patient information is estimated from the IPD-only variance
# ESS = (precision_gain_from_source) / (per_patient_precision)
# where per_patient_precision = precision_IPD_only / n_asian_ipd_patients

# Per-patient precision from IPD-only analysis
var_ipd_only_M6  <- var(curves_nb$asian[, which(x_days == 168)])
var_ipd_only_M12 <- var(curves_nb$asian[, which(x_days == 336)])

per_patient_prec_M6  <- (1 / var_ipd_only_M6)  / n_asian_patients
per_patient_prec_M12 <- (1 / var_ipd_only_M12) / n_asian_patients

for (src in c("Pub1", "Pub2", "RWE")) {
  # Variance from LOSO (dropping this source)
  var_without_M6  <- var(loso_results[[src]]$curves$asian[, which(x_days == 168)])
  var_without_M12 <- var(loso_results[[src]]$curves$asian[, which(x_days == 336)])
  
  # Precision gain attributable to this source
  prec_gain_M6  <- 1/var(curves_primary$asian[, which(x_days == 168)]) - 1/var_without_M6
  prec_gain_M12 <- 1/var(curves_primary$asian[, which(x_days == 336)]) - 1/var_without_M12
  
  # ESS = precision gain / per-patient precision
  ess_M6  <- max(0, prec_gain_M6  / per_patient_prec_M6)
  ess_M12 <- max(0, prec_gain_M12 / per_patient_prec_M12)
  
  # Also report the variance ratio (how many times narrower)
  var_ratio_M6  <- var_without_M6  / var(curves_primary$asian[, which(x_days == 168)])
  var_ratio_M12 <- var_without_M12 / var(curves_primary$asian[, which(x_days == 336)])
  
  ess_results[ess_results$Source == src, "ESS_M6"]  <- round(ess_M6, 1)
  ess_results[ess_results$Source == src, "ESS_M12"] <- round(ess_M12, 1)
  ess_results[ess_results$Source == src, "VarRatio_M6"]  <- round(var_ratio_M6, 2)
  ess_results[ess_results$Source == src, "VarRatio_M12"] <- round(var_ratio_M12, 2)
}

cat("\nEffective Sample Size Decomposition:\n")
cat("  ESS: equivalent number of Asian IPD patients contributed by each source\n")
cat("  VarRatio: (variance without source) / (variance with all sources)\n")
cat("            Values > 1 indicate the source contributed precision.\n\n")
print(ess_results)
df_to_latex(ess_results, rownames=T, file='Effective_sample_size.tex')


################################################################################
# 6.  PILLAR 3 — ROBUSTNESS
################################################################################

# ============================================================================ #
#  S5: Prior Sensitivity
# ============================================================================ #

cat("\n========== S5: PRIOR SENSITIVITY ==========\n")

# --- Scenario (a): Original priors  (already done as primary) ---

# --- Scenario (b): Vague priors on delta parameters ---
code_vague_deltas <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  e0    ~ dnorm(0, sd = sqrt(10))
  emax  ~ dnorm(0, sd = 1)
  k     ~ dnorm(0, sd = 10)
  ed50  ~ T(dnorm(10, sd = sqrt(10)), 0, )
  r1    ~ T(dnorm(0, sd = 10), 0, )
  # VAGUE deltas: SD multiplied by 10 from original
  de0   ~ dnorm(0, sd = 10)        # was 1
  demax ~ dnorm(0, sd = 10)        # was 1
  dk    ~ dnorm(0, sd = 10)        # was sqrt(10) ≈ 3.16
  ded50 ~ T(dnorm(0, sd = 10), -ed50, )   # was sqrt(10)
  dr1   ~ T(dnorm(0, sd = 100), -r1, )    # was 10
  a     ~ dnorm(0, sd = 10)
  tau2  ~ dinvgamma(0.01, 0.01)
})

cat("  Running: Vague delta priors\n")
mcmc_vague <- run_nimble(
  code      = code_vague_deltas,
  constants = inputs$constants,
  data      = inputs$data,
  inits     = inputs$inits,
  fast      = TRUE    # [v8] 40k sweeps: sufficient for prior sensitivity comparison
)
vague_summary <- summarise_mcmc(mcmc_vague, delta_params)

# --- Scenario (c): Informative priors ---
code_informative_deltas <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  e0    ~ dnorm(0, sd = sqrt(10))
  emax  ~ dnorm(0, sd = 1)
  k     ~ dnorm(0, sd = 10)
  ed50  ~ T(dnorm(10, sd = sqrt(10)), 0, )
  r1    ~ T(dnorm(0, sd = 10), 0, )
  # INFORMATIVE deltas: small offsets, tight SD
  de0   ~ dnorm(-0.05, sd = 0.5)
  demax ~ dnorm(0.05, sd = 0.5)
  dk    ~ dnorm(0, sd = 0.5)
  ded50 ~ T(dnorm(0, sd = 5), -ed50, )
  dr1   ~ T(dnorm(0, sd = 5), -r1, )
  a     ~ dnorm(0, sd = 10)
  tau2  ~ dinvgamma(0.01, 0.01)
})

cat("  Running: Informative delta priors\n")
mcmc_informative <- run_nimble(
  code      = code_informative_deltas,
  constants = inputs$constants,
  data      = inputs$data,
  inits     = inputs$inits,
  fast      = TRUE    # [v8] 40k sweeps: sufficient for prior sensitivity comparison
)
informative_summary <- summarise_mcmc(mcmc_informative, delta_params)

cat("\nPrior Sensitivity — Posterior Means for Delta Parameters:\n")
prior_sens_table <- data.frame(
  Scenario = c("Original", "Vague (10x SD)", "Informative"),
  dE0   = c(primary_summary["de0","Mean"],
            vague_summary["de0","Mean"],
            informative_summary["de0","Mean"]),
  dE0_95CI = c(
    paste0("[", round(primary_summary["de0","2.5%"],3), ", ",
           round(primary_summary["de0","97.5%"],3), "]"),
    paste0("[", round(vague_summary["de0","2.5%"],3), ", ",
           round(vague_summary["de0","97.5%"],3), "]"),
    paste0("[", round(informative_summary["de0","2.5%"],3), ", ",
           round(informative_summary["de0","97.5%"],3), "]")
  ),
  dEmax = c(primary_summary["demax","Mean"],
            vague_summary["demax","Mean"],
            informative_summary["demax","Mean"]),
  dEmax_95CI = c(
    paste0("[", round(primary_summary["demax","2.5%"],3), ", ",
           round(primary_summary["demax","97.5%"],3), "]"),
    paste0("[", round(vague_summary["demax","2.5%"],3), ", ",
           round(vague_summary["demax","97.5%"],3), "]"),
    paste0("[", round(informative_summary["demax","2.5%"],3), ", ",
           round(informative_summary["demax","97.5%"],3), "]")
  )
)
print(prior_sens_table, row.names = FALSE)
df_to_latex(prior_sens_table, rownames=T, file='prior_sens_table.tex')

save.image("Sensitivity_results.RData")
# ============================================================================ #
#  A6: Tipping Point Analysis
#      [FIX 5] Tipping criterion aligned with manuscript:
#      "the smallest |delta| at which the 95% credible interval for at
#       least one ethnic difference parameter excludes zero in the
#       direction indicating a clinically meaningful difference"
# ============================================================================ #

cat("\n========== S6: TIPPING POINT ANALYSIS ==========\n")

# Delta parameters to monitor for tipping
delta_check_params <- c("de0", "demax", "ded50", "dr1", "dk")

# ── [Fix 1] Reuse primary at shift=0 instead of re-running ──────────────────
tipping_results <- data.frame()
shift_summary_0 <- summarise_mcmc(mcmc_primary, delta_check_params)
ethnic_diff_0   <- curves_primary$asian[, which(x_days == 336)] -
  curves_primary$nonasian[, which(x_days == 336)]
tipping_results <- rbind(tipping_results, data.frame(
  shift = 0,
  de0_mean   = shift_summary_0["de0","Mean"],
  de0_lower  = shift_summary_0["de0","2.5%"],
  de0_upper  = shift_summary_0["de0","97.5%"],
  demax_mean  = shift_summary_0["demax","Mean"],
  demax_lower = shift_summary_0["demax","2.5%"],
  demax_upper = shift_summary_0["demax","97.5%"],
  db_mean  = shift_summary_0["dk","Mean"],
  db_lower = shift_summary_0["dk","2.5%"],
  db_upper = shift_summary_0["dk","97.5%"],
  any_CrI_excl_zero = any(shift_summary_0[,"2.5%"] > 0 |
                            shift_summary_0[,"97.5%"] < 0),
  tipped_params = paste(rownames(shift_summary_0)[
    shift_summary_0[,"2.5%"] > 0 | shift_summary_0[,"97.5%"] < 0],
    collapse=","),
  curve_diff_M12 = mean(ethnic_diff_0),
  stringsAsFactors = FALSE
))

shift_values <- c(-0.20, -0.15, -0.10, -0.05, 0.05, 0.10, 0.15, 0.20)  # 0 removed

# [Fix 3] Dispatch tipping-point fits in parallel.  Workers also need
# base_iop_asian/nonasian (computed from the unshifted primary dataset).
tipping_export <- c("y.iop", "code_primary", "delta_check_params",
                    "base_iop_asian", "base_iop_nonasian",
                    "run_nimble", "build_primary_inputs", "summarise_mcmc",
                    "compute_fitted_curves",
                    ".run_one_chain", ".install_xen_samplers")

tipping_raw <- .run_analyses_parallel(.tipping_worker, shift_values,
                                      tipping_export, envir = environment())
tipping_results <- rbind(tipping_results,
                         do.call(rbind, tipping_raw))
tipping_results <- tipping_results[order(tipping_results$shift), ]

cat("\nTipping Point Analysis (Hu et al. shifted):\n")
cat("Criterion: 95% CrI of any delta parameter excludes zero\n\n")
print(tipping_results[, c("shift", "de0_mean", "de0_lower", "de0_upper",
                          "demax_mean", "demax_lower", "demax_upper",
                          "db_mean", "db_lower", "db_upper",
                          "any_CrI_excl_zero", "tipped_params")],
      row.names = FALSE)
df_to_latex(tipping_results, rownames=T, file='tipping_results.tex')


# Identify tipping point
tp_idx <- which(tipping_results$any_CrI_excl_zero)
if (length(tp_idx) > 0) {
  tp_shift <- tipping_results$shift[min(tp_idx)]
  cat("\nTipping point reached at shift =", tp_shift,
      "\n  Parameter(s):", tipping_results$tipped_params[min(tp_idx)], "\n")
} else {
  cat("\nNo tipping point reached within tested range [-0.20, +0.20].\n")
  cat("Conclusion is robust to shifts of up to 20 percentage points",
      "in the Hu et al. data.\n")
}

save.image("Sensitivity_results.RData")
# ============================================================================ #
#  S7: Alternative Borrowing Methods
# ============================================================================ #

cat("\n========== S7: ALTERNATIVE BORROWING METHODS ==========\n")

# --------------------------------------------------------------------------- #
#  7a. Power Prior
#      [FIX 6] Three-way split:
#        P13-001 IPD:   full weight (var = tau2, since n=1)
#        RWE aggregate: full weight (var = tau2 / n_rwe)
#        Hu/Sng:        discounted  (var = tau2 / (n_disc * alpha0))
# --------------------------------------------------------------------------- #

# Split data into three groups
# NOTE: With RWE now as simulated IPD (n=1), we must use source-based
# filters to avoid double-counting. The n==1 filter would capture
# P13-001 + RWE, overlapping with the RWE block.
is_ipd    <- y.iop$source %in% c("Global_NonAsian", "Global_Asian")   # P13-001 only
is_rwe    <- y.iop$source == "RWE"                               # RWE (full weight)
is_disc   <- y.iop$source %in% c("Pub1", "Pub2")                   # Discounted

dat_ipd  <- y.iop[is_ipd, ]
dat_rwe  <- y.iop[is_rwe, ]
dat_disc <- y.iop[is_disc, ]

# Verify no overlap
stopifnot(nrow(dat_ipd) + nrow(dat_rwe) + nrow(dat_disc) == nrow(y.iop))

cat("  Power prior data split:\n")
cat("    IPD (P13-001):", nrow(dat_ipd), "rows\n")
cat("    RWE (full wt):", nrow(dat_rwe), "rows\n")
cat("    Hu/Sng (disc):", nrow(dat_disc), "rows\n")

alpha_values <- c(0.25, 0.50, 0.75)   # [Fix 1] a0=1.00 is the primary
power_prior_results <- list()
power_prior_results[["1"]] <- primary_summary  # [Fix 1] reuse primary

# [Fix 3] Dispatch power-prior fits in parallel.  Workers need the
# three pre-split data frames and the power-prior model code.
pp_export <- c("dat_ipd", "dat_rwe", "dat_disc", "code_power_prior",
               "delta_params",
               "run_nimble", "summarise_mcmc",
               ".run_one_chain", ".install_xen_samplers")

pp_raw <- .run_analyses_parallel(.pp_worker, alpha_values,
                                 pp_export, envir = environment())
for (i in seq_along(alpha_values)) {
  power_prior_results[[as.character(alpha_values[i])]] <- pp_raw[[i]]
}

cat("\nPower Prior — dE0 and dEmax by alpha0 (Hu/Sng discounted, RWE full):\n")
pp_table_alpha <- c(alpha_values, 1.00)  # include a0=1 (reused from primary)
pp_table <- data.frame(
  alpha0 = pp_table_alpha,
  dE0_mean  = sapply(as.character(pp_table_alpha), function(a) power_prior_results[[a]]["de0","Mean"]),
  dE0_lower = sapply(as.character(pp_table_alpha), function(a) power_prior_results[[a]]["de0","2.5%"]),
  dE0_upper = sapply(as.character(pp_table_alpha), function(a) power_prior_results[[a]]["de0","97.5%"]),
  dEmax_mean  = sapply(as.character(pp_table_alpha), function(a) power_prior_results[[a]]["demax","Mean"]),
  dEmax_lower = sapply(as.character(pp_table_alpha), function(a) power_prior_results[[a]]["demax","2.5%"]),
  dEmax_upper = sapply(as.character(pp_table_alpha), function(a) power_prior_results[[a]]["demax","97.5%"])
)
print_rounded(pp_table, 4)
df_to_latex(pp_table, rownames=T, file='power_prior_results.tex')

save.image("Sensitivity_results.RData")

# --------------------------------------------------------------------------- #
#  7b. Commensurate Prior
#      [FIX 2] e0_agg centred on (e0+de0), emax_agg on (emax+demax)
#      Note: the commensurate prior model code (code_commensurate above)
#      already contains the corrected centering.
# --------------------------------------------------------------------------- #

cat("  Running: Commensurate prior model\n")

# [v5 — final scope] Commensurate prior structure:
#   IPD block       = P13-001 (non-Asian + Asian) + RWE.  Direct current-
#                     study data; no discounting.
#   Aggregate block = Hu + Sng, each with source-specific e0_agg[s] and
#                     emax_agg[s] linked via commensurability precision
#                     tau_c to the Asian population parameters
#                     (e0 + de0) and (emax + demax).
#   Source indexing: s = 1 -> Hu, s = 2 -> Sng.
#
#   This scope matches the power-prior scope and makes all three S7
#   methods (power, commensurate, MAP) target the same borrowing
#   decision: how much to trust the external Chinese/Asian literature
#   data when inferring ethnic differences.
dat_ipd_comm <- y.iop[y.iop$source %in% c("Global_NonAsian", "Global_Asian", "RWE"), ]
dat_agg_comm <- y.iop[y.iop$source %in% c("Pub1", "Pub2"), ]
source_ids   <- ifelse(dat_agg_comm$source == "Pub1", 1, 2)

stopifnot(nrow(dat_ipd_comm) + nrow(dat_agg_comm) == nrow(y.iop))

cat("  Commensurate data split:\n")
cat("    IPD block (P13-001 + RWE, full borrowing):", nrow(dat_ipd_comm), "rows\n")
cat("    External block (Hu + Sng, commensurate): ", nrow(dat_agg_comm), "rows\n")

comm_constants <- list(N_ipd = nrow(dat_ipd_comm), N_agg = nrow(dat_agg_comm),
                       N_sources = 2, source_id = source_ids)
comm_data <- list(
  y        = dat_ipd_comm$pdpch,
  x        = dat_ipd_comm$visit,
  z        = dat_ipd_comm$base.pd.cent,
  race     = dat_ipd_comm$race,
  n        = dat_ipd_comm$n,
  y_agg    = dat_agg_comm$pdpch,
  x_agg    = dat_agg_comm$visit,
  z_agg    = dat_agg_comm$base.pd.cent,
  race_agg = dat_agg_comm$race,
  n_agg    = dat_agg_comm$n
)
comm_inits <- list(e0 = -0.6, de0 = 0, emax = 0.3, demax = 0,
                   ed50 = 10, ded50 = 1, r1 = 1, dr1 = 0.1,
                   k = -0.0003, dk = 0, a = -0.02, tau2 = 0.05,
                   tau_c = 1,
                   e0_agg   = c(-0.6, -0.6),
                   emax_agg = c( 0.3,  0.3))

comm_monitors <- c(delta_params, "tau_c", "e0_agg", "emax_agg")

mcmc_comm <- run_nimble(
  code      = code_commensurate,
  constants = comm_constants,
  data      = comm_data,
  inits     = comm_inits,
  monitors  = comm_monitors,
  fast      = TRUE    # [v8] 40k sweeps
)

comm_summary <- summarise_mcmc(mcmc_comm,
                               c(delta_params, "tau_c"))
cat("\nCommensurate Prior — Posterior Summary:\n")
print_rounded(comm_summary, 4)
df_to_latex(comm_summary, rownames=T, file='commensurate_prior_results.tex')

save.image("Sensitivity_results.RData")

# --------------------------------------------------------------------------- #
#  7c. Robust MAP Prior
# --------------------------------------------------------------------------- #

cat("  Running: Robust MAP prior model\n")

# [v5 — final scope] Robust MAP prior structure, aligned with power and
# commensurate priors:
#   IPD block ("current study") = P13-001 (non-Asian + Asian) + RWE.
#     Shared ethnic offsets de0, demax receive the robust mixture prior.
#   External block = Hu + Sng.
#     Source-specific ethnic offsets de0_ext[s], demax_ext[s]
#     form a hierarchy with population mean (de0_pop, demax_pop) and
#     between-source SD tau_map.
#   The robust mixture component in the current-study prior on (de0, demax)
#   protects against prior-data conflict between the Hu/Sng hierarchy
#   and the P13-001 + RWE data.
#
#   Source indexing in the external block:
#     s = 1  Hu et al.
#     s = 2  Sng et al.
dat_ipd_map <- y.iop[y.iop$source %in% c("Global_NonAsian", "Global_Asian", "RWE"), ]
dat_ext_map <- y.iop[y.iop$source %in% c("Pub1", "Pub2"), ]
ext_source_id <- ifelse(dat_ext_map$source == "Pub1", 1, 2)

stopifnot(nrow(dat_ipd_map) + nrow(dat_ext_map) == nrow(y.iop))

cat("  Robust MAP data split:\n")
cat("    IPD block (P13-001 + RWE, current study):", nrow(dat_ipd_map), "rows\n")
cat("    External block (Hu + Sng, hierarchy):    ", nrow(dat_ext_map), "rows\n")

w_robust_val <- 0.2
sd_vague_val <- 10

map_constants <- list(
  N_ipd          = nrow(dat_ipd_map),
  N_agg          = nrow(dat_ext_map),
  N_ext_sources  = 2,
  agg_source_id  = ext_source_id
)
map_data <- list(
  y        = dat_ipd_map$pdpch,
  x        = dat_ipd_map$visit,
  z        = dat_ipd_map$base.pd.cent,
  race     = dat_ipd_map$race,
  n        = dat_ipd_map$n,
  y_agg    = dat_ext_map$pdpch,
  x_agg    = dat_ext_map$visit,
  z_agg    = dat_ext_map$base.pd.cent,
  n_agg    = dat_ext_map$n,
  w_robust = w_robust_val,
  sd_vague = sd_vague_val
)
map_inits <- list(
  e0 = -0.6, emax = 0.3,
  de0 = 0, demax = 0,
  k = -0.0003, dk = 0, a = -0.02, tau2 = 0.05,
  ed50 = 10, ded50 = 1, r1 = 1, dr1 = 0.1,
  de0_ext   = c(0, 0),
  demax_ext = c(0, 0),
  de0_pop = 0, demax_pop = 0,
  tau_map = 0.3
)

map_monitors <- c("e0", "emax", "ed50", "r1", "k", "dk", "ded50", "dr1", "a",
                  "de0", "demax",
                  "de0_ext", "demax_ext", "de0_pop", "demax_pop",
                  "tau_map", "tau2")

mcmc_map <- run_nimble(
  code      = code_robust_map,
  constants = map_constants,
  data      = map_data,
  inits     = map_inits,
  monitors  = map_monitors,
  fast      = TRUE    # [v8] 40k sweeps
)

map_summary <- summarise_mcmc(mcmc_map)
cat("\nRobust MAP Prior — Posterior Summary:\n")
print_rounded(map_summary, 4)

cat("\nRobust MAP — Key ethnic difference parameters:\n")
# [v5] Key parameters under the new scope:
#   de0, demax          = current-study (P13-001 + RWE) ethnic offsets
#                         (receive the robust mixture prior)
#   de0_pop, demax_pop  = MAP population mean across Hu + Sng
#   tau_map             = between-source SD in the Hu/Sng hierarchy
map_key <- map_summary[c("de0", "demax",
                         "de0_pop", "demax_pop",
                         "tau_map"), , drop = FALSE]
print_rounded(map_key, 4)
df_to_latex(map_key, rownames=T, file='Rubust_MAP_prior_results.tex')

save.image("Sensitivity_results.RData")

# --------------------------------------------------------------------------- #
#  S7 Summary comparison across borrowing methods
# --------------------------------------------------------------------------- #
sink(file='summary_all_prior_method.txt')

cat("\n========== S7 SUMMARY: BORROWING METHODS COMPARISON ==========\n")
cat("                     dE0 (or equiv.)      dEmax (or equiv.)\n")
cat("Method               Mean   [95% CrI]     Mean   [95% CrI]\n")
cat("---------------------------------------------------------------\n")

cat(sprintf("Primary (alpha=1)   %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            primary_summary["de0","Mean"], primary_summary["de0","2.5%"], primary_summary["de0","97.5%"],
            primary_summary["demax","Mean"], primary_summary["demax","2.5%"], primary_summary["demax","97.5%"]))

pp05 <- power_prior_results[["0.5"]]
cat(sprintf("Power (alpha=0.5)   %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            pp05["de0","Mean"], pp05["de0","2.5%"], pp05["de0","97.5%"],
            pp05["demax","Mean"], pp05["demax","2.5%"], pp05["demax","97.5%"]))

pp025 <- power_prior_results[["0.25"]]
cat(sprintf("Power (alpha=0.25)  %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            pp025["de0","Mean"], pp025["de0","2.5%"], pp025["de0","97.5%"],
            pp025["demax","Mean"], pp025["demax","2.5%"], pp025["demax","97.5%"]))

cat(sprintf("Commensurate        %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            comm_summary["de0","Mean"], comm_summary["de0","2.5%"], comm_summary["de0","97.5%"],
            comm_summary["demax","Mean"], comm_summary["demax","2.5%"], comm_summary["demax","97.5%"]))

# [v5] Under the new scope, the "current-study" ethnic offsets for the
# robust MAP are the shared de0/demax of the P13+RWE IPD block (not a
# source-specific de0_s[1]).  de0_pop is the population mean over Hu+Sng.
cat(sprintf("Robust MAP (pop)    %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            map_summary["de0_pop","Mean"], map_summary["de0_pop","2.5%"], map_summary["de0_pop","97.5%"],
            map_summary["demax_pop","Mean"], map_summary["demax_pop","2.5%"], map_summary["demax_pop","97.5%"]))

cat(sprintf("Robust MAP (CS)     %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            map_summary["de0","Mean"], map_summary["de0","2.5%"], map_summary["de0","97.5%"],
            map_summary["demax","Mean"], map_summary["demax","2.5%"], map_summary["demax","97.5%"]))
# CS = current-study ethnic offset (P13-001 + RWE block), receives the
# robust mixture prior anchored on the Hu/Sng hierarchy.

cat("---------------------------------------------------------------\n")
cat("tau_map (between-source SD):", round(map_summary["tau_map","Mean"], 4),
    " [", round(map_summary["tau_map","2.5%"], 4), ",",
    round(map_summary["tau_map","97.5%"], 4), "]\n")
cat("tau_c   (commensurability):", round(comm_summary["tau_c","Mean"], 4),
    " [", round(comm_summary["tau_c","2.5%"], 4), ",",
    round(comm_summary["tau_c","97.5%"], 4), "]\n")
sink()


save.image("Sensitivity_results.RData")

# ============================================================================ #
#  S8: Structural Model Sensitivity
#      [FIX 10] Added WAIC computation for model comparison
# ============================================================================ #

cat("\n========== S8: STRUCTURAL MODEL SENSITIVITY ==========\n")

# --- 8a. Emax-only (no linear component) ---
cat("  Running: Emax-only model\n")

emax_monitors <- c("e0", "de0", "emax", "demax", "ed50", "ded50",
                   "r1", "dr1", "a", "tau2")
mcmc_emax <- run_nimble(
  code      = code_emax_only,
  constants = inputs$constants,
  data      = inputs$data,
  inits     = list(e0 = -0.6, de0 = 0, emax = 0.3, demax = 0,
                   ed50 = 10, ded50 = 1, r1 = 1, dr1 = 0.1,
                   a = -0.02, tau2 = 0.05),
  monitors   = emax_monitors,
  enableWAIC = TRUE,    # [FIX 10]
  fast       = TRUE     # [v8] 40k sweeps; WAIC for ranking is stable at this length
)

emax_summary <- summarise_mcmc(mcmc_emax, emax_monitors)
cat("\nEmax-Only Model — Posterior Summary:\n")
print_rounded(emax_summary, 4)
df_to_latex(emax_summary, rownames=T, file='emax_model_only_results.tex')

save.image("Sensitivity_results.RData")

# --- 8b. Piecewise linear (3 knots) ---
# [FIX 9] Now uses 3 knots (Day 1, Day 28, Day 84) with 4 segments
cat("  Running: Piecewise linear model (3 knots: Day 1, 28, 84)\n")

pw_monitors <- c("e0", "de0", "b1", "db1", "b2", "db2",
                 "b3", "db3", "b4", "db4", "a", "tau2")
mcmc_pw <- run_nimble(
  code      = code_piecewise,
  constants = inputs$constants,
  data      = inputs$data,
  inits     = list(e0 = -0.6, de0 = 0,
                   b1 = 0.01, db1 = 0,
                   b2 = 0.005, db2 = 0,
                   b3 = -0.0003, db3 = 0,
                   b4 = -0.0003, db4 = 0,
                   a = -0.02, tau2 = 0.05),
  monitors   = pw_monitors,
  enableWAIC = TRUE,    # [FIX 10]
  fast       = TRUE     # [v8] 40k sweeps
)

pw_summary <- summarise_mcmc(mcmc_pw, pw_monitors)
cat("\nPiecewise Linear Model — Posterior Summary:\n")
print_rounded(pw_summary, 4)
df_to_latex(pw_summary, rownames=T, file='Piecewise_linear_model_results.tex')

save.image("Sensitivity_results.RData")

# --------------------------------------------------------------------------- #
#  [Fix 4] Collect primary WAIC via a short sequential single-chain run.
#  The primary analysis above ran without WAIC (parallel chains); WAIC from
#  5 000 post-burnin samples is stable for the S8 three-way comparison.
# --------------------------------------------------------------------------- #
cat("  [Fix 4] Computing primary WAIC (5 000 iter, 1 chain, sequential)...\n")
mcmc_primary_waic_only <- run_nimble(
  code       = code_primary,
  constants  = inputs$constants,
  data       = inputs$data,
  inits      = inputs$inits,
  nBurnin    = 5000, nIter = 5000, nChains = 1,
  enableWAIC = TRUE,
  parallel   = FALSE
)
attr(mcmc_primary, "WAIC") <- attr(mcmc_primary_waic_only, "WAIC")
rm(mcmc_primary_waic_only)

# --------------------------------------------------------------------------- #
#  [FIX 10] WAIC comparison across structural models
# --------------------------------------------------------------------------- #
sink(file='WAIC_comparison.txt')
cat("\n--- Model Comparison (WAIC) ---\n")
waic_primary  <- attr(mcmc_primary, "WAIC")
waic_emax     <- attr(mcmc_emax, "WAIC")
waic_pw       <- attr(mcmc_pw, "WAIC")

if (!is.null(waic_primary) && !is.null(waic_emax) && !is.null(waic_pw)) {
  waic_table <- data.frame(
    Model = c("Primary (Emax + linear)", "Emax-only", "Piecewise linear (3 knots)"),
    WAIC  = c(waic_primary$WAIC, waic_emax$WAIC, waic_pw$WAIC),
    pWAIC = c(waic_primary$pWAIC, waic_emax$pWAIC, waic_pw$pWAIC)
  )
  waic_table$dWAIC <- waic_table$WAIC - min(waic_table$WAIC)
  cat("\n")
  print_rounded(waic_table, 2)
  cat("\n  Lower WAIC is preferred. dWAIC = difference from best model.\n")
  cat("  pWAIC = effective number of parameters.\n")
} else {
  cat("  WAIC not available for all models — check nimble version.\n")
  cat("  nimble >= 0.12.0 required for enableWAIC.\n")
}
sink()

save.image("Sensitivity_results.RData")

################################################################################
# 7.  CONVERGENCE DIAGNOSTICS
################################################################################
sink(file='convergence_diagnositic.txt')
cat("\n========== CONVERGENCE DIAGNOSTICS ==========\n")

if (length(mcmc_primary) >= 2) {
  gr <- gelman.diag(mcmc_primary, multivariate = FALSE)
  cat("\nGelman-Rubin Statistics (Primary Analysis):\n")
  print_rounded(gr$psrf, 3)
  
  cat("\nEffective Sample Sizes (Primary Analysis):\n")
  print_rounded(effectiveSize(mcmc_primary), 0)
}
sink()

################################################################################
# 8.  SAVE ALL RESULTS
################################################################################

results_all <- list(
  primary      = list(summary = primary_summary, mcmc = mcmc_primary),
  A1_ppc       = ppc_results,
  A2_loso      = loso_results,
  A3_noborrow  = list(summary = nb_summary, ci_comparison = ci_comparison),
  A4_ess       = ess_results,
  A5_prior     = list(original = primary_summary,
                      vague = vague_summary,
                      informative = informative_summary),
  A6_tipping   = tipping_results,
  A7_power     = power_prior_results,
  A7_commensurate = comm_summary,
  A7_robust_map   = map_summary,
  A8_emax_only = emax_summary,
  A8_piecewise = pw_summary
)

save(results_all, file = paste0("XEN_sensitivity_results_", fileNameTime, ".RData"))
cat("\nResults saved to: XEN_sensitivity_results_", fileNameTime, ".RData\n")

cat("\n========== ALL ANALYSES COMPLETE ==========\n")

t01 <- proc.time()[3]
(t01-t00)/60


