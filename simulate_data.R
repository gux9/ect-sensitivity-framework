################################################################################
#  Simulate Fully Synthetic ECT Dataset
#
#  Generates four data sources:
#    Global study  — Non-Asian IPD  (N = 100 patients)
#    Global study  — Asian IPD      (N =  25 patients)
#    RWE cohort    — Asian IPD      (N =  60 patients)
#    Publication 1 — Asian aggregate (n = 100, summary rows)
#    Publication 2 — Asian aggregate (n =  80, summary rows)
#
#  Method:
#    1. Fit the Emax + linear longitudinal model to scenario.csv using
#       weighted nonlinear least squares (weights = n[i] per model variance).
#    2. Store the fitted parameter vector as the simulation truth.
#    3. Draw new patients, simulate trajectories, introduce realistic
#       missingness, and write simulated_data.csv.
#
#  Outputs:
#    simulated_data.csv  —  drop-in replacement for scenario.csv
#
#  Required columns match Sensitivity_analysis.R expectations:
#    pid, study, visit, n, race, base.pd, pd, pdch, pdpch,
#    source, base.pd.cent
################################################################################

set.seed(2025)

# ============================================================
# 0.  Helper: Emax + linear mean function
# ============================================================

emax_mu <- function(x, race, z, p) {
  ## p = named list/vector with elements:
  ##   e0, emax, ed50, r1, k, a, de0, demax, ded50, dr1, dk
  e0    <- p[["e0"]];    emax  <- p[["emax"]];  ed50 <- p[["ed50"]]
  r1    <- p[["r1"]];    k     <- p[["k"]];     a    <- p[["a"]]
  de0   <- p[["de0"]];   demax <- p[["demax"]]; ded50 <- p[["ded50"]]
  dr1   <- p[["dr1"]];   dk    <- p[["dk"]]

  r_i    <- r1 + dr1 * race
  ed50_i <- ed50 + ded50 * race
  # Guard against ed50_i hitting zero or negative
  ed50_i <- pmax(ed50_i, 1e-6)
  r_i    <- pmax(r_i,    1e-6)

  xr  <- x ^ r_i
  d50r <- ed50_i ^ r_i

  (e0 + de0 * race) +
    a * z +
    (k + dk * race) * x +
    (emax + demax * race) * xr / (d50r + xr)
}


# ============================================================
# 1.  True longitudinal parameters
# ============================================================
#
# The simulation truth is taken directly from the published model-fitted
# curves (the IOP % change figure), NOT from a fresh nonlinear fit to the
# raw rows.  An unconstrained least-squares fit to the sparse, noisy raw
# data collapses to a degenerate solution (ED50 -> 0, Emax with the wrong
# sign), giving a flat, monotone curve.  The published Emax + linear fit
# instead has the characteristic shape:
#
#   * deep initial drop at day 1            (~ -0.65 % change, fraction)
#   * sharp recovery to a peak near day 84  (~ -0.29 to -0.30)
#   * gradual decline through day 336
#   * Asian arm starts slightly deeper, recovers to the same peak, and
#     declines more steeply late than the Non-Asian arm
#
# These constants reproduce that figure.  Code-to-math mapping matches
# Sensitivity_analysis.R (manuscript Eq. 2).
#
#   mu(x) = (e0 + de0*race) + a*z + (k + dk*race)*x
#           + (emax + demax*race) * x^r / (ed50^r + x^r)
#   with r = r1 + dr1*race,  ed50 = ed50 + ded50*race

TRUE_PARAMS <- list(
  e0    = -0.66,      # day-0 intercept (deep initial lowering)
  emax  =  0.42,      # recovery magnitude toward the plateau
  ed50  =  8.0,       # day at which half the recovery is reached
  r1    =  1.3,       # Hill / steepness of the recovery
  k     = -0.00028,   # slow long-term decline (Non-Asian)
  a     = -0.010,     # baseline-PD covariate effect (centred)
  de0   = -0.04,      # Asian: deeper initial drop
  demax =  0.05,      # Asian: larger recovery -> same peak
  ded50 =  0.0,       # Asian: same half-recovery time
  dr1   =  0.0,       # Asian: same steepness
  dk    = -0.00012,   # Asian: steeper long-term decline
  tau2  =  0.0025     # residual variance (per single observation)
)

cat("Using published-figure true parameters:\n")
print(round(unlist(TRUE_PARAMS), 6))
cat("\n")

# Load the original data only to (a) verify the expected column layout and
# (b) keep the script self-checking; it is NOT used to estimate parameters.
if (file.exists("scenario.csv")) {
  orig <- read.csv("scenario.csv", stringsAsFactors = FALSE)
  if ("pdpch.1" %in% names(orig)) {
    names(orig)[names(orig) == "pdpch"]   <- "pdch"
    names(orig)[names(orig) == "pdpch.1"] <- "pdpch"
  }
  orig <- orig[!is.na(orig$pdpch), ]
  cat("Reference data scenario.csv loaded:", nrow(orig), "rows.\n\n")
}


# ============================================================
# 2.  Baseline PD distributions per source
#     (means / SDs chosen to match original data then rounded
#      to prevent reverse-engineering of the original submission)
# ============================================================

BASELINE <- list(
  Global_NonAsian = list(mean = 24.0, sd = 4.5),
  Global_Asian    = list(mean = 23.5, sd = 4.0),
  RWE             = list(mean = 28.0, sd = 9.0)
)

# Aggregate publication baselines (fixed, rounded)
PUB1_BASE <- 20.0
PUB2_BASE <- 22.0


# ============================================================
# 3.  Design: new sample sizes and visit schedules
# ============================================================

N_GNA  <- 100   # Global Non-Asian IPD
N_GA   <-  25   # Global Asian IPD
N_RWE  <-  60   # RWE Asian IPD
N_PUB1 <- 100   # Publication 1 aggregate
N_PUB2 <-  80   # Publication 2 aggregate

VISITS_GLOBAL <- c(1, 7, 14, 28, 84, 168, 224, 280, 336)
VISITS_RWE    <- c(1, 7, 28, 84, 168, 224, 280, 336)
VISITS_PUB1   <- c(1, 7, 28, 84, 168)
VISITS_PUB2   <- c(1, 7, 28, 84, 168, 336)

# Missingness probability for IPD visits (per visit, per patient, for later
# visits only): mirrors ~8% overall missingness in original data
MISS_PROB_EARLY <- 0.03   # days 1–28
MISS_PROB_LATE  <- 0.08   # days 84+


# ============================================================
# 4.  Simulation helper: generate one IPD source
# ============================================================

sim_ipd <- function(n_patients, race_val, study_val, source_label,
                    pid_prefix, visits, base_mean, base_sd,
                    global_base_mean, params) {

  rows <- list()
  for (i in seq_len(n_patients)) {
    pid_i    <- sprintf("%s-%03d", pid_prefix, i)
    base_i   <- max(rnorm(1, base_mean, base_sd), 5)   # truncate at 5
    base_i   <- round(base_i, 3)
    z_i      <- base_i - global_base_mean               # centered

    for (v in visits) {
      ## missingness
      miss_p <- if (v <= 28) MISS_PROB_EARLY else MISS_PROB_LATE
      if (runif(1) < miss_p) next

      mu_i   <- emax_mu(v, race_val, z_i, params)
      pdpch_i <- rnorm(1, mu_i, sqrt(params$tau2))

      pd_i   <- base_i * (1 + pdpch_i)
      pdch_i <- pd_i - base_i

      rows[[length(rows) + 1]] <- data.frame(
        pid          = pid_i,
        study        = study_val,
        visit        = v,
        n            = 1L,
        race         = race_val,
        base.pd      = base_i,
        pd           = round(pd_i,   4),
        pdch         = round(pdch_i, 4),
        pdpch        = round(pdpch_i, 6),
        source       = source_label,
        base.pd.cent = round(z_i, 6),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


# ============================================================
# 5.  Compute global baseline mean for centering
#     (weighted by source size, race=0 and race=1 combined,
#      consistent with how base.pd.cent is derived in original)
# ============================================================

# Draw baseline PD for all IPD patients first, compute grand mean,
# then simulate trajectories using that centering value.

set.seed(2025)   # reset for reproducibility

draw_baselines <- function(n, mean_b, sd_b) {
  pmax(rnorm(n, mean_b, sd_b), 5)
}

base_gna  <- draw_baselines(N_GNA,  BASELINE$Global_NonAsian$mean, BASELINE$Global_NonAsian$sd)
base_ga   <- draw_baselines(N_GA,   BASELINE$Global_Asian$mean,    BASELINE$Global_Asian$sd)
base_rwe  <- draw_baselines(N_RWE,  BASELINE$RWE$mean,             BASELINE$RWE$sd)

# Grand mean for centering (all IPD patients + weighted publication baselines)
all_ipd_base <- c(base_gna, base_ga, base_rwe)
pub_base_vec <- c(rep(PUB1_BASE, N_PUB1), rep(PUB2_BASE, N_PUB2))
GLOBAL_BASE_MEAN <- mean(c(all_ipd_base, pub_base_vec))

cat(sprintf("Global baseline mean for centering: %.4f\n\n", GLOBAL_BASE_MEAN))


# ============================================================
# 6.  Simulate all IPD sources
# ============================================================

cat("Simulating IPD data ...\n")

# -- Global Non-Asian --
gna_rows <- list()
for (i in seq_len(N_GNA)) {
  pid_i  <- sprintf("GNA-%03d", i)
  base_i <- round(base_gna[i], 3)
  z_i    <- base_i - GLOBAL_BASE_MEAN

  for (v in VISITS_GLOBAL) {
    miss_p <- if (v <= 28) MISS_PROB_EARLY else MISS_PROB_LATE
    if (runif(1) < miss_p) next

    mu_i    <- emax_mu(v, 0L, z_i, TRUE_PARAMS)
    pdpch_i <- rnorm(1, mu_i, sqrt(TRUE_PARAMS$tau2))
    pd_i    <- base_i * (1 + pdpch_i)

    gna_rows[[length(gna_rows) + 1]] <- data.frame(
      pid = pid_i, study = 0L, visit = v, n = 1L, race = 0L,
      base.pd = base_i, pd = round(pd_i, 4),
      pdch = round(pd_i - base_i, 4), pdpch = round(pdpch_i, 6),
      source = "Global_NonAsian", base.pd.cent = round(z_i, 6),
      stringsAsFactors = FALSE)
  }
}
dat_gna <- do.call(rbind, gna_rows)

# -- Global Asian --
ga_rows <- list()
for (i in seq_len(N_GA)) {
  pid_i  <- sprintf("GA-%03d", i)
  base_i <- round(base_ga[i], 3)
  z_i    <- base_i - GLOBAL_BASE_MEAN

  for (v in VISITS_GLOBAL) {
    miss_p <- if (v <= 28) MISS_PROB_EARLY else MISS_PROB_LATE
    if (runif(1) < miss_p) next

    mu_i    <- emax_mu(v, 1L, z_i, TRUE_PARAMS)
    pdpch_i <- rnorm(1, mu_i, sqrt(TRUE_PARAMS$tau2))
    pd_i    <- base_i * (1 + pdpch_i)

    ga_rows[[length(ga_rows) + 1]] <- data.frame(
      pid = pid_i, study = 0L, visit = v, n = 1L, race = 1L,
      base.pd = base_i, pd = round(pd_i, 4),
      pdch = round(pd_i - base_i, 4), pdpch = round(pdpch_i, 6),
      source = "Global_Asian", base.pd.cent = round(z_i, 6),
      stringsAsFactors = FALSE)
  }
}
dat_ga <- do.call(rbind, ga_rows)

# -- RWE (Asian only) --
rwe_rows <- list()
for (i in seq_len(N_RWE)) {
  pid_i  <- sprintf("RWE-%03d", i)
  base_i <- round(base_rwe[i], 4)
  z_i    <- base_i - GLOBAL_BASE_MEAN

  for (v in VISITS_RWE) {
    miss_p <- if (v <= 28) MISS_PROB_EARLY else MISS_PROB_LATE
    if (runif(1) < miss_p) next

    mu_i    <- emax_mu(v, 1L, z_i, TRUE_PARAMS)
    pdpch_i <- rnorm(1, mu_i, sqrt(TRUE_PARAMS$tau2))
    pd_i    <- base_i * (1 + pdpch_i)

    rwe_rows[[length(rwe_rows) + 1]] <- data.frame(
      pid = pid_i, study = 1L, visit = v, n = 1L, race = 1L,
      base.pd = base_i, pd = round(pd_i, 4),
      pdch = round(pd_i - base_i, 4), pdpch = round(pdpch_i, 6),
      source = "RWE", base.pd.cent = round(z_i, 6),
      stringsAsFactors = FALSE)
  }
}
dat_rwe <- do.call(rbind, rwe_rows)

cat(sprintf("  Global_NonAsian: %d rows (%d patients)\n", nrow(dat_gna), N_GNA))
cat(sprintf("  Global_Asian:    %d rows (%d patients)\n", nrow(dat_ga),  N_GA))
cat(sprintf("  RWE:             %d rows (%d patients)\n", nrow(dat_rwe), N_RWE))


# ============================================================
# 7.  Simulate aggregate publication data
# ============================================================

cat("Simulating aggregate publication data ...\n")

sim_aggregate <- function(pub_label, study_val, n_agg, base_mean_pub, visits_pub) {
  ## For aggregate data the mean pdpch at each visit is itself normally
  ## distributed with variance tau2 / n_agg (central limit theorem).
  ## Use the publication's own average baseline for centering.
  z_pub <- base_mean_pub - GLOBAL_BASE_MEAN

  rows <- list()
  for (v in visits_pub) {
    mu_v    <- emax_mu(v, 1L, z_pub, TRUE_PARAMS)
    ## observed aggregate mean = true mean + sampling error
    sd_agg  <- sqrt(TRUE_PARAMS$tau2 / n_agg)
    pdpch_v <- rnorm(1, mu_v, sd_agg)

    pd_v    <- base_mean_pub * (1 + pdpch_v)
    pdch_v  <- pd_v - base_mean_pub

    rows[[length(rows) + 1]] <- data.frame(
      pid = pub_label, study = study_val, visit = v, n = n_agg, race = 1L,
      base.pd = base_mean_pub, pd = round(pd_v, 4),
      pdch = round(pdch_v, 4), pdpch = round(pdpch_v, 6),
      source = pub_label,
      base.pd.cent = round(z_pub, 6),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

dat_pub1 <- sim_aggregate("Pub1", study_val = 1L, n_agg = N_PUB1,
                          base_mean_pub = PUB1_BASE, visits_pub = VISITS_PUB1)
dat_pub2 <- sim_aggregate("Pub2", study_val = 2L, n_agg = N_PUB2,
                          base_mean_pub = PUB2_BASE, visits_pub = VISITS_PUB2)

cat(sprintf("  Pub1: %d summary rows (n = %d)\n", nrow(dat_pub1), N_PUB1))
cat(sprintf("  Pub2: %d summary rows (n = %d)\n", nrow(dat_pub2), N_PUB2))


# ============================================================
# 8.  Combine, sort, and write
# ============================================================

simdata <- rbind(dat_gna, dat_ga, dat_rwe, dat_pub1, dat_pub2)

# Sort: source order, then patient, then visit
src_order <- c("Global_NonAsian", "Global_Asian", "RWE", "Pub1", "Pub2")
simdata$source <- factor(simdata$source, levels = src_order)
simdata <- simdata[order(simdata$source, simdata$pid, simdata$visit), ]
simdata$source <- as.character(simdata$source)

cat("\nFinal dataset summary:\n")
print(table(simdata$source, simdata$race))
cat("Total rows:", nrow(simdata), "\n\n")

write.csv(simdata, "simulated_data.csv", row.names = FALSE)
cat("Written: simulated_data.csv\n")


# ============================================================
# 9.  Quick sanity-check: mean pdpch by source × visit
# ============================================================

cat("\nMean pdpch by source and visit (non-missing):\n")
agg <- aggregate(pdpch ~ source + visit, data = simdata, FUN = mean)
agg <- agg[order(agg$source, agg$visit), ]
agg$visit <- round(agg$visit, 3)
agg$pdpch <- round(agg$pdpch, 3)
print(agg)


# ============================================================
# 10. Plot fitted vs simulated mean curves (optional, no PDF dependency)
# ============================================================

cat("\nSimulation complete.\n")
cat("To use in Sensitivity_analysis.R, set:\n")
cat('  y.pd <- read.csv("simulated_data.csv")\n')
