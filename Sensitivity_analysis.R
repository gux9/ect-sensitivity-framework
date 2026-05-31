################################################################################
#  Ethnic Bridging Study
#  Bayesian Hybrid Evidence Synthesis: Primary Analysis + Sensitivity Framework
#
#  Analyses produced by this script:
#    PRIMARY  Emax + linear longitudinal model
#    S1       Posterior predictive checks (heterogeneity diagnostics)
#    S2       Leave-one-source-out (LOSO)
#    S3       No-borrowing reference (IPD only)
#    S4       Effective sample size decomposition
#    S5       Prior sensitivity (vague / informative delta priors)
#    S6       Tipping point on Pub1 . data
#    S7       Alternative borrowing methods
#               7a. Power prior      (Pub1/Pub2 discounted by alpha0)
#               7b. Commensurate prior
#               7c. Robust MAP prior
#    S8       Structural model sensitivity (Emax-only, piecewise linear)
#               with WAIC comparison
#    Convergence diagnostics (Gelman-Rubin, ESS)
#
#  Input data file:  sensitivity.csv
#  Required columns:
#    pid, study, visit, n, race, base.pd, pd, pdch, pdpch,
#    source, base.pd.cent
#
#  Race coding:  0 = Non-Asian, 1 = Asian
#  Sources:      Global_NonAsian, Global_Asian, RWE, Pub1, Pub2
#
#  Required R packages: nimble, coda
#
#  Code-to-math mapping (manuscript Eq. 2):
#    e0    = E_0          emax  = E_max        ed50  = ED_50
#    r1    = r            k     = b            a     = a
#    de0   = Delta E_0    demax = Delta E_max  ded50 = Delta ED_50
#    dr1   = Delta r      dk    = Delta b      tau2  = sigma^2
################################################################################

library(nimble)
library(coda)

set.seed(1223456)


################################################################################
# 0. LOAD DATA
################################################################################

y.pd <- read.csv("sensitivity.csv")

stopifnot(all(c("pid", "study", "visit", "n", "race", "base.pd",
                "pd", "pdch", "pdpch", "source", "base.pd.cent")
              %in% names(y.pd)))

# Drop rows with missing response
y.pd <- y.pd[!is.na(y.pd$pdpch), ]

cat("Data loaded:", nrow(y.pd), "rows\n")
cat("Sources:", paste(sort(unique(y.pd$source)), collapse = ", "), "\n")
cat("Race coding: 0 = Non-Asian, 1 = Asian\n")
print(table(y.pd$source, y.pd$race))


################################################################################
# 1. MODEL DEFINITIONS (nimble)
################################################################################

# ---------------------------------------------------------------------------- #
# 1A. PRIMARY MODEL — Emax + linear with baseline pd covariate
# ---------------------------------------------------------------------------- #
code_primary <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  # Priors
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


# ---------------------------------------------------------------------------- #
# 1B. EMAX-ONLY MODEL (S8) — drops linear component
# ---------------------------------------------------------------------------- #
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


# ---------------------------------------------------------------------------- #
# 1C. PIECEWISE LINEAR MODEL (S8) — three knots at Day 1, 28, 84
# ---------------------------------------------------------------------------- #
code_piecewise <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    # Values at each knot (continuity)
    val_cp1[i] <- (b1 + db1 * race[i]) * cp1
    val_cp2[i] <- val_cp1[i] + (b2 + db2 * race[i]) * (cp2 - cp1)
    val_cp3[i] <- val_cp2[i] + (b3 + db3 * race[i]) * (cp3 - cp2)
    # Piecewise function via step indicators
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
  b1  ~ dnorm(0, sd = 1);    db1 ~ dnorm(0, sd = 1)
  b2  ~ dnorm(0, sd = 0.1);  db2 ~ dnorm(0, sd = 0.1)
  b3  ~ dnorm(0, sd = 0.01); db3 ~ dnorm(0, sd = 0.01)
  b4  ~ dnorm(0, sd = 0.01); db4 ~ dnorm(0, sd = 0.01)
  a    ~ dnorm(0, sd = 10)
  tau2 ~ dinvgamma(0.01, 0.01)
})


# ---------------------------------------------------------------------------- #
# 1D. POWER PRIOR MODEL (S7) — Pub1/Pub2 discounted by alpha0
#     IPD (Global-001): full weight (var = tau2)
#     RWE aggregate: full weight (var = tau2 / n_rwe)
#     Pub1/Pub2:        discounted  (var = tau2 / (n_disc * alpha0))
# ---------------------------------------------------------------------------- #
code_power_prior <- nimbleCode({
  # IPD block
  for (i in 1:N_ipd) {
    y[i] ~ dnorm(mu[i], var = tau2)
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  # RWE block (full weight)
  for (j in 1:N_rwe) {
    y_rwe[j] ~ dnorm(mu_rwe[j], var = tau2 / n_rwe[j])
    mu_rwe[j] <- e0 + de0 * race_rwe[j] + a * z_rwe[j] +
      (k + dk * race_rwe[j]) * x_rwe[j] +
      (emax + demax * race_rwe[j]) * x_rwe[j]^(r1 + dr1 * race_rwe[j]) /
      ((ed50 + ded50 * race_rwe[j])^(r1 + dr1 * race_rwe[j]) +
         x_rwe[j]^(r1 + dr1 * race_rwe[j]))
  }
  # Pub1/Pub2 block (discounted by alpha0)
  for (m in 1:N_disc) {
    y_disc[m] ~ dnorm(mu_disc[m], var = tau2 / (n_disc[m] * alpha0))
    mu_disc[m] <- e0 + de0 * race_disc[m] + a * z_disc[m] +
      (k + dk * race_disc[m]) * x_disc[m] +
      (emax + demax * race_disc[m]) * x_disc[m]^(r1 + dr1 * race_disc[m]) /
      ((ed50 + ded50 * race_disc[m])^(r1 + dr1 * race_disc[m]) +
         x_disc[m]^(r1 + dr1 * race_disc[m]))
  }
  # Priors (same as primary)
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


# ---------------------------------------------------------------------------- #
# 1E. COMMENSURATE PRIOR MODEL (S7)
#     IPD block (Global-001 + RWE): full borrowing.
#     External block (Pub1, Pub2): source-specific e0_agg[s] and emax_agg[s]
#       linked to Asian population (e0 + de0, emax + demax) via tau_c.
# ---------------------------------------------------------------------------- #
code_commensurate <- nimbleCode({
  # IPD block
  for (i in 1:N_ipd) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) * x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  # External aggregate block (Pub1, Pub2) - both Asian
  for (j in 1:N_agg) {
    y_agg[j] ~ dnorm(mu_agg[j], var = tau2 / n_agg[j])
    mu_agg[j] <- e0_agg[source_id[j]] + a * z_agg[j] +
      (k + dk) * x_agg[j] +
      emax_agg[source_id[j]] *
      x_agg[j]^(r1 + dr1) /
      ((ed50 + ded50)^(r1 + dr1) +
         x_agg[j]^(r1 + dr1))
  }
  # Commensurate priors anchored on Asian population parameters
  for (s in 1:N_sources) {
    e0_agg[s]   ~ dnorm(e0   + de0,   var = 1 / tau_c)
    emax_agg[s] ~ dnorm(emax + demax, var = 1 / tau_c)
  }
  tau_c ~ dgamma(1, 1)
  # Priors (same as primary)
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


# ---------------------------------------------------------------------------- #
# 1F. ROBUST MAP PRIOR MODEL (S7)
#     Mixture prior on current-study ethnic offsets (de0, demax):
#       (1 - w_robust) * N(de0_pop, tau_map^2)  +  w_robust * N(0, sd_vague^2)
#     Pub1 and Pub2 form the MAP hierarchy.
# ---------------------------------------------------------------------------- #

# Custom robust-mixture distribution
dRobustNorm <- nimbleFunction(
  run = function(x        = double(0),
                 mu_map   = double(0),
                 sd_map   = double(0),
                 sd_vague = double(0),
                 w_robust = double(0),
                 log      = integer(0, default = 0)) {
    returnType(double(0))
    log_map   <- dnorm(x, mu_map, sd_map,   log = TRUE) + log(1 - w_robust)
    log_vague <- dnorm(x, 0,      sd_vague, log = TRUE) + log(w_robust)
    max_log <- max(log_map, log_vague)
    logp    <- max_log + log(exp(log_map - max_log) + exp(log_vague - max_log))
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
    if (u < w_robust) return(rnorm(1, 0, sd_vague))
    return(rnorm(1, mu_map, sd_map))
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
  # IPD block (current study)
  for (i in 1:N_ipd) {
    y[i] ~ dnorm(mu[i], var = tau2 / n[i])
    mu[i] <- e0 + de0 * race[i] + a * z[i] +
      (k + dk * race[i]) * x[i] +
      (emax + demax * race[i]) *
      x[i]^(r1 + dr1 * race[i]) /
      ((ed50 + ded50 * race[i])^(r1 + dr1 * race[i]) +
         x[i]^(r1 + dr1 * race[i]))
  }
  # External block (Pub1, Pub2) - source-specific ethnic offsets
  for (j in 1:N_agg) {
    y_agg[j] ~ dnorm(mu_agg[j], var = tau2 / n_agg[j])
    mu_agg[j] <- e0 + de0_ext[agg_source_id[j]] + a * z_agg[j] +
      (k + dk) * x_agg[j] +
      (emax + demax_ext[agg_source_id[j]]) *
      x_agg[j]^(r1 + dr1) /
      ((ed50 + ded50)^(r1 + dr1) +
         x_agg[j]^(r1 + dr1))
  }
  # MAP hierarchy across external sources
  for (s in 1:N_ext_sources) {
    de0_ext[s]   ~ dnorm(de0_pop,   sd = tau_map)
    demax_ext[s] ~ dnorm(demax_pop, sd = tau_map)
  }
  # Robust mixture prior on current-study ethnic offsets
  de0   ~ dRobustNorm(de0_pop,   tau_map, sd_vague, w_robust)
  demax ~ dRobustNorm(demax_pop, tau_map, sd_vague, w_robust)
  # Hyperpriors
  de0_pop   ~ dnorm(0, sd = 1)
  demax_pop ~ dnorm(0, sd = 1)
  tau_map   ~ T(dnorm(0, sd = 0.5), 0, )
  # Priors (same as primary)
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
# 2. HELPER FUNCTIONS
################################################################################

# Build nimble data / constants / inits for the primary model
build_primary_inputs <- function(dat) {
  list(
    constants = list(N = nrow(dat)),
    data      = list(y    = dat$pdpch,
                     x    = dat$visit,
                     z    = dat$base.pd.cent,
                     race = dat$race,
                     n    = dat$n),
    inits     = list(e0 = -0.6, de0 = 0, emax = 0.3, demax = 0,
                     ed50 = 10, ded50 = 1, r1 = 1, dr1 = 0.1,
                     k = -0.0003, dk = 0, a = -0.02, tau2 = 0.05)
  )
}


# Generic nimble MCMC wrapper.
# Block samplers on (ed50, ded50, r1, dr1) and (e0, de0, emax, demax) and
# slice samplers on r1 / dr1 are installed because these intercept- and
# shape-level Emax parameters are strongly correlated; block / slice
# updates navigate the posterior ridge far more efficiently than RW.
run_nimble <- function(code, constants, data, inits,
                       monitors   = NULL,
                       nBurnin    = 10000,
                       nIter      = 10000,
                       nThin      = 1,
                       nChains    = 2,
                       setSeed    = 1223456,
                       enableWAIC = FALSE) {
  if (is.null(monitors)) {
    monitors <- c("e0", "de0", "emax", "demax", "ed50", "ded50",
                  "r1", "dr1", "k", "dk", "a", "tau2")
  }
  model <- nimbleModel(code, constants = constants, data = data,
                       inits = inits, check = FALSE)
  node_names <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  node_tops  <- unique(gsub("\\[.*", "", node_names))
  monitors   <- intersect(monitors, node_tops)

  cModel   <- compileNimble(model)
  mcmcConf <- configureMCMC(model, monitors = monitors,
                            enableWAIC = enableWAIC)

  # Block samplers on correlated Emax parameters
  shape_block <- c("ed50", "ded50", "r1", "dr1")
  if (all(shape_block %in% node_tops)) {
    mcmcConf$removeSamplers(shape_block)
    mcmcConf$addSampler(target = shape_block, type = "RW_block")
  }
  level_block <- c("e0", "de0", "emax", "demax")
  if (all(level_block %in% node_tops)) {
    mcmcConf$removeSamplers(level_block)
    mcmcConf$addSampler(target = level_block, type = "RW_block")
  }
  # Slice samplers handle the heavy right tail on r1, dr1
  for (tail_par in c("r1", "dr1")) {
    if (tail_par %in% node_tops) {
      mcmcConf$addSampler(target = tail_par, type = "slice")
    }
  }

  mcmc  <- buildMCMC(mcmcConf)
  cMCMC <- compileNimble(mcmc, project = model)

  samples_list <- list()
  for (ch in seq_len(nChains)) {
    set.seed(setSeed + ch)
    cMCMC$run(nBurnin + nIter * nThin, thin = nThin, nburnin = nBurnin)
    samples_list[[ch]] <- coda::as.mcmc(as.matrix(cMCMC$mvSamples))
  }
  mcmc_out <- coda::as.mcmc.list(samples_list)

  if (enableWAIC) {
    tryCatch(
      attr(mcmc_out, "WAIC") <- cMCMC$getWAIC(),
      error = function(e) message("WAIC extraction failed: ", e$message)
    )
  }
  mcmc_out
}


# Summarise MCMC output: Mean, Median, SD, 95% CrI
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


# Round numeric columns and print
print_rounded <- function(df, digits = 4) {
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- round(df[num_cols], digits)
  print(df)
}


# Posterior fitted curves for Asian (race = 1) and Non-Asian (race = 0)
compute_fitted_curves <- function(mcmc_out, x_days,
                                  base_pd_asian, base_pd_nonasian) {
  samps <- as.data.frame(do.call(rbind, mcmc_out))

  fitted_nonasian <- sapply(seq_along(x_days), function(j) {
    samps$e0 + samps$a * base_pd_nonasian +
      samps$k * x_days[j] +
      samps$emax * x_days[j]^samps$r1 /
      (samps$ed50^samps$r1 + x_days[j]^samps$r1)
  })
  fitted_asian <- sapply(seq_along(x_days), function(j) {
    e0a   <- samps$e0   + samps$de0
    emaxa <- samps$emax + samps$demax
    r1a   <- samps$r1   + samps$dr1
    ed50a <- samps$ed50 + samps$ded50
    ka    <- samps$k    + samps$dk
    e0a + samps$a * base_pd_asian +
      ka * x_days[j] +
      emaxa * x_days[j]^r1a / (ed50a^r1a + x_days[j]^r1a)
  })
  list(nonasian = fitted_nonasian, asian = fitted_asian)
}


# Quick summary of fitted curves (mean and 95% CrI by day)
curve_summary <- function(fitted_mat, x_days) {
  data.frame(
    day   = x_days,
    mean  = colMeans(fitted_mat, na.rm = TRUE),
    lower = apply(fitted_mat, 2, quantile, 0.025, na.rm = TRUE),
    upper = apply(fitted_mat, 2, quantile, 0.975, na.rm = TRUE)
  )
}


# Sample-size-weighted baseline pd by ethnic group
# (needed for hybrid IPD + AD data so AD rows are not down-weighted)
weighted_base_pd <- function(dat) {
  asian    <- dat[dat$race == 1, ]
  nonasian <- dat[dat$race == 0, ]
  list(
    asian    = sum(asian$base.pd.cent    * asian$n)    / sum(asian$n),
    nonasian = sum(nonasian$base.pd.cent * nonasian$n) / sum(nonasian$n)
  )
}


################################################################################
# 3. PRIMARY ANALYSIS
################################################################################

cat("\n========== PRIMARY ANALYSIS ==========\n")

inputs <- build_primary_inputs(y.pd)

# WAIC enabled here so the S8 model-comparison block below has it available
mcmc_primary <- run_nimble(
  code       = code_primary,
  constants  = inputs$constants,
  data       = inputs$data,
  inits      = inputs$inits,
  enableWAIC = TRUE
)

delta_params <- c("e0", "de0", "emax", "demax", "ed50", "ded50",
                  "r1", "dr1", "k", "dk", "a", "tau2")
primary_summary <- summarise_mcmc(mcmc_primary, delta_params)

cat("\nPrimary Analysis - Posterior Summary:\n")
print_rounded(primary_summary, 4)

# Fitted curves at canonical visit days (Day 1 through Month 12)
x_days <- c(1, 7, 14, 28, 84, 168, 224, 280, 336)
base_pd_full     <- weighted_base_pd(y.pd)
base_pd_asian    <- base_pd_full$asian
base_pd_nonasian <- base_pd_full$nonasian

curves_primary <- compute_fitted_curves(mcmc_primary, x_days,
                                        base_pd_asian, base_pd_nonasian)

cat("\nFitted Asian curve (mean, 95% CrI):\n")
print_rounded(curve_summary(curves_primary$asian, x_days), 4)
cat("\nFitted Non-Asian curve (mean, 95% CrI):\n")
print_rounded(curve_summary(curves_primary$nonasian, x_days), 4)


################################################################################
# 4. PILLAR 1 - APPROPRIATENESS
################################################################################

# ---------------------------------------------------------------------------- #
# S1. Posterior predictive checks (visit-level)
# ---------------------------------------------------------------------------- #

cat("\n========== S1: POSTERIOR PREDICTIVE CHECKS ==========\n")

ppc_check <- function(mcmc_out, y.pd, source_label,
                      base_pd_asian, base_pd_nonasian) {
  src_dat <- y.pd[y.pd$source == source_label, ]
  if (nrow(src_dat) == 0) return(NULL)

  is_asian <- src_dat$race[1] == 1
  is_ad    <- source_label %in% c("Pub1", "Pub2")   # aggregate-data sources

  if (!is_ad) {
    # IPD: aggregate across patients within each visit
    agg <- aggregate(pdpch ~ visit, data = src_dat,
                     FUN = function(x) c(N = length(x),
                                         mean = mean(x),
                                         sd   = sd(x)))
    visit_tab <- data.frame(visit    = agg$visit,
                            N        = agg$pdpch[, "N"],
                            obs_mean = agg$pdpch[, "mean"],
                            obs_sd   = agg$pdpch[, "sd"])
    visit_tab$se_mean <- ifelse(visit_tab$N > 1,
                                visit_tab$obs_sd / sqrt(visit_tab$N),
                                NA_real_)
  } else {
    # AD: each row is already a visit-level summary
    visit_tab <- data.frame(visit    = src_dat$visit,
                            N        = src_dat$n,
                            obs_mean = src_dat$pdpch,
                            obs_sd   = NA_real_,
                            se_mean  = NA_real_)
  }

  curves <- compute_fitted_curves(mcmc_out, visit_tab$visit,
                                  base_pd_asian, base_pd_nonasian)
  pred <- if (is_asian) curves$asian else curves$nonasian

  visit_tab$predicted <- colMeans(pred, na.rm = TRUE)
  visit_tab$residual  <- visit_tab$obs_mean - visit_tab$predicted
  visit_tab$z         <- visit_tab$residual / visit_tab$se_mean
  visit_tab$p_value   <- sapply(seq_along(visit_tab$visit), function(j) {
    mean(pred[, j] <= visit_tab$obs_mean[j], na.rm = TRUE)
  })

  data.frame(source    = source_label,
             day       = visit_tab$visit,
             N         = visit_tab$N,
             obs_mean  = visit_tab$obs_mean,
             predicted = visit_tab$predicted,
             residual  = visit_tab$residual,
             se_mean   = visit_tab$se_mean,
             z         = visit_tab$z,
             p_value   = visit_tab$p_value,
             stringsAsFactors = FALSE)
}

ppc_sources <- c("Global_NonAsian", "Global_Asian", "RWE", "Pub1", "Pub2")
ppc_results <- do.call(rbind, Filter(Negate(is.null), lapply(ppc_sources, function(s) {
  tryCatch(
    ppc_check(mcmc_primary, y.pd, s, base_pd_asian, base_pd_nonasian),
    error = function(e) {
      message("PPC skipped for ", s, ": ", e$message); NULL
    }
  )
})))

cat("\nS1 Posterior Predictive Check Summary (visit-level):\n")
ppc_display <- ppc_results
num_cols    <- setdiff(names(ppc_display), c("source", "day", "N"))
ppc_display[num_cols] <- lapply(ppc_display[num_cols], round, 4)
print(ppc_display, row.names = FALSE)

# Diagnostic flags
cat("\nS1 diagnostic flags:\n")
big_z <- ppc_results[!is.na(ppc_results$z) & abs(ppc_results$z) > 2,
                     c("source", "day", "N", "residual", "se_mean", "z")]
if (nrow(big_z) > 0) {
  cat("  Cells with |z| > 2:\n")
  print(big_z, row.names = FALSE)
} else {
  cat("  No cells with |z| > 2.\n")
}
mar <- aggregate(abs(residual) ~ source, data = ppc_results, FUN = mean)
names(mar)[2] <- "mean_abs_residual"
cat("\nMean absolute residual by source:\n")
print(mar, row.names = FALSE)


# ---------------------------------------------------------------------------- #
# S2. Leave-one-source-out (LOSO)
# ---------------------------------------------------------------------------- #

cat("\n========== S2: LEAVE-ONE-SOURCE-OUT ==========\n")

loso_sources <- c("Global_Asian", "RWE", "Pub1", "Pub2")
loso_results <- list()

for (src in loso_sources) {
  cat("  Dropping:", src, "\n")
  dat_reduced <- y.pd[y.pd$source != src, ]
  inp         <- build_primary_inputs(dat_reduced)
  base_red    <- weighted_base_pd(dat_reduced)
  mcmc_red    <- run_nimble(code      = code_primary,
                            constants = inp$constants,
                            data      = inp$data,
                            inits     = inp$inits)
  loso_results[[src]] <- list(
    summary           = summarise_mcmc(mcmc_red, delta_params),
    curves            = compute_fitted_curves(mcmc_red, x_days,
                                              base_red$asian, base_red$nonasian),
    base_pd_asian    = base_red$asian,
    base_pd_nonasian = base_red$nonasian
  )
}

loso_table <- data.frame(
  Dropped = c("None", loso_sources),
  dE0     = c(primary_summary["de0", "Mean"],
              sapply(loso_results, function(r) r$summary["de0", "Mean"])),
  dE0_95CI = c(
    paste0("[", round(primary_summary["de0", "2.5%"], 3), ", ",
                round(primary_summary["de0", "97.5%"], 3), "]"),
    sapply(loso_results, function(r)
      paste0("[", round(r$summary["de0", "2.5%"], 3), ", ",
                  round(r$summary["de0", "97.5%"], 3), "]"))
  ),
  dEmax   = c(primary_summary["demax", "Mean"],
              sapply(loso_results, function(r) r$summary["demax", "Mean"])),
  dEmax_95CI = c(
    paste0("[", round(primary_summary["demax", "2.5%"], 3), ", ",
                round(primary_summary["demax", "97.5%"], 3), "]"),
    sapply(loso_results, function(r)
      paste0("[", round(r$summary["demax", "2.5%"], 3), ", ",
                  round(r$summary["demax", "97.5%"], 3), "]"))
  ),
  dED50   = c(primary_summary["ded50", "Mean"],
              sapply(loso_results, function(r) r$summary["ded50", "Mean"])),
  dED50_95CI = c(
    paste0("[", round(primary_summary["ded50", "2.5%"], 3), ", ",
                round(primary_summary["ded50", "97.5%"], 3), "]"),
    sapply(loso_results, function(r)
      paste0("[", round(r$summary["ded50", "2.5%"], 3), ", ",
                  round(r$summary["ded50", "97.5%"], 3), "]"))
  )
)

cat("\nLOSO - Posterior means for delta parameters:\n")
print_rounded(loso_table[, c("Dropped", "dE0", "dEmax", "dED50")], 4)
cat("\nLOSO - Full table with 95% CrI:\n")
print(loso_table, row.names = FALSE)


################################################################################
# 5. PILLAR 2 - VALUE OF BORROWING
################################################################################

# ---------------------------------------------------------------------------- #
# S3. No-borrowing reference (IPD only)
# ---------------------------------------------------------------------------- #

cat("\n========== S3: NO-BORROWING REFERENCE ==========\n")

dat_ipd_only <- y.pd[y.pd$n == 1, ]
cat("  No-borrowing dataset:", nrow(dat_ipd_only), "rows (",
    sum(dat_ipd_only$race == 1), "Asian,",
    sum(dat_ipd_only$race == 0), "non-Asian)\n")

inp_nb <- build_primary_inputs(dat_ipd_only)
mcmc_noborrow <- run_nimble(
  code      = code_primary,
  constants = inp_nb$constants,
  data      = inp_nb$data,
  inits     = inp_nb$inits
)

base_nb    <- weighted_base_pd(dat_ipd_only)
nb_summary <- summarise_mcmc(mcmc_noborrow, delta_params)
curves_nb  <- compute_fitted_curves(mcmc_noborrow, x_days,
                                    base_nb$asian, base_nb$nonasian)

cat("\nNo-Borrowing - Posterior Summary:\n")
print_rounded(nb_summary, 4)

# CrI width comparison vs. full borrowing
ci_comparison <- data.frame(
  Analysis        = c("Full Borrowing", "No Borrowing"),
  Asian_M6_width  = c(
    diff(quantile(curves_primary$asian[, which(x_days == 168)], c(0.025, 0.975))),
    diff(quantile(curves_nb$asian[,      which(x_days == 168)], c(0.025, 0.975)))
  ),
  Asian_M12_width = c(
    diff(quantile(curves_primary$asian[, which(x_days == 336)], c(0.025, 0.975))),
    diff(quantile(curves_nb$asian[,      which(x_days == 336)], c(0.025, 0.975)))
  )
)
cat("\n95% CrI Width Comparison (Asian curve):\n")
print_rounded(ci_comparison, 4)


# ---------------------------------------------------------------------------- #
# S4. Effective sample size decomposition
#     ESS = (precision gain from source) / (per-patient precision)
#     VarRatio = Var(without source) / Var(with all sources)
# ---------------------------------------------------------------------------- #

cat("\n========== S4: EFFECTIVE SAMPLE SIZE ==========\n")

n_asian_patients <- length(unique(dat_ipd_only$pid[dat_ipd_only$race == 1]))
cat("  Asian IPD patients:", n_asian_patients, "\n")

ess_results <- data.frame(
  Source              = c("Pub1", "Pub2", "RWE"),
  N_rows              = c(sum(y.pd$source == "Pub1"),
                          sum(y.pd$source == "Pub2"),
                          sum(y.pd$source == "RWE")),
  N_patients_reported = c(63, 31, 23)
)

var_ipd_only_M6  <- var(curves_nb$asian[, which(x_days == 168)])
var_ipd_only_M12 <- var(curves_nb$asian[, which(x_days == 336)])
per_patient_prec_M6  <- (1 / var_ipd_only_M6)  / n_asian_patients
per_patient_prec_M12 <- (1 / var_ipd_only_M12) / n_asian_patients

for (src in c("Pub1", "Pub2", "RWE")) {
  var_without_M6  <- var(loso_results[[src]]$curves$asian[, which(x_days == 168)])
  var_without_M12 <- var(loso_results[[src]]$curves$asian[, which(x_days == 336)])
  prec_gain_M6   <- 1/var(curves_primary$asian[, which(x_days == 168)]) - 1/var_without_M6
  prec_gain_M12  <- 1/var(curves_primary$asian[, which(x_days == 336)]) - 1/var_without_M12
  ess_M6         <- max(0, prec_gain_M6  / per_patient_prec_M6)
  ess_M12        <- max(0, prec_gain_M12 / per_patient_prec_M12)
  var_ratio_M6   <- var_without_M6  / var(curves_primary$asian[, which(x_days == 168)])
  var_ratio_M12  <- var_without_M12 / var(curves_primary$asian[, which(x_days == 336)])

  ess_results[ess_results$Source == src, "ESS_M6"]       <- round(ess_M6, 1)
  ess_results[ess_results$Source == src, "ESS_M12"]      <- round(ess_M12, 1)
  ess_results[ess_results$Source == src, "VarRatio_M6"]  <- round(var_ratio_M6, 2)
  ess_results[ess_results$Source == src, "VarRatio_M12"] <- round(var_ratio_M12, 2)
}

cat("\nEffective Sample Size Decomposition:\n")
cat("  ESS: equivalent Asian IPD patients contributed by each source\n")
cat("  VarRatio: (variance without source) / (variance with all sources)\n\n")
print(ess_results, row.names = FALSE)


################################################################################
# 6. PILLAR 3 - ROBUSTNESS
################################################################################

# ---------------------------------------------------------------------------- #
# S5. Prior sensitivity (vague vs. informative delta priors)
# ---------------------------------------------------------------------------- #

cat("\n========== S5: PRIOR SENSITIVITY ==========\n")

# Scenario (b): VAGUE delta priors (SD multiplied by 10 vs. primary)
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
  de0   ~ dnorm(0, sd = 10)
  demax ~ dnorm(0, sd = 10)
  dk    ~ dnorm(0, sd = 10)
  ded50 ~ T(dnorm(0, sd = 10), -ed50, )
  dr1   ~ T(dnorm(0, sd = 100), -r1, )
  a     ~ dnorm(0, sd = 10)
  tau2  ~ dinvgamma(0.01, 0.01)
})

cat("  Running: Vague delta priors\n")
mcmc_vague <- run_nimble(code_vague_deltas,
                         inputs$constants, inputs$data, inputs$inits)
vague_summary <- summarise_mcmc(mcmc_vague, delta_params)

# Scenario (c): INFORMATIVE delta priors (small offsets, tight SD)
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
  de0   ~ dnorm(-0.05, sd = 0.5)
  demax ~ dnorm( 0.05, sd = 0.5)
  dk    ~ dnorm(0, sd = 0.5)
  ded50 ~ T(dnorm(0, sd = 5), -ed50, )
  dr1   ~ T(dnorm(0, sd = 5), -r1, )
  a     ~ dnorm(0, sd = 10)
  tau2  ~ dinvgamma(0.01, 0.01)
})

cat("  Running: Informative delta priors\n")
mcmc_informative <- run_nimble(code_informative_deltas,
                               inputs$constants, inputs$data, inputs$inits)
informative_summary <- summarise_mcmc(mcmc_informative, delta_params)

prior_sens_table <- data.frame(
  Scenario = c("Original", "Vague (10x SD)", "Informative"),
  dE0   = c(primary_summary["de0", "Mean"],
            vague_summary["de0", "Mean"],
            informative_summary["de0", "Mean"]),
  dE0_95CI = c(
    paste0("[", round(primary_summary["de0", "2.5%"], 3), ", ",
                round(primary_summary["de0", "97.5%"], 3), "]"),
    paste0("[", round(vague_summary["de0", "2.5%"], 3), ", ",
                round(vague_summary["de0", "97.5%"], 3), "]"),
    paste0("[", round(informative_summary["de0", "2.5%"], 3), ", ",
                round(informative_summary["de0", "97.5%"], 3), "]")
  ),
  dEmax = c(primary_summary["demax", "Mean"],
            vague_summary["demax", "Mean"],
            informative_summary["demax", "Mean"]),
  dEmax_95CI = c(
    paste0("[", round(primary_summary["demax", "2.5%"], 3), ", ",
                round(primary_summary["demax", "97.5%"], 3), "]"),
    paste0("[", round(vague_summary["demax", "2.5%"], 3), ", ",
                round(vague_summary["demax", "97.5%"], 3), "]"),
    paste0("[", round(informative_summary["demax", "2.5%"], 3), ", ",
                round(informative_summary["demax", "97.5%"], 3), "]")
  )
)
cat("\nPrior Sensitivity - Posterior means for delta parameters:\n")
print(prior_sens_table, row.names = FALSE)


# ---------------------------------------------------------------------------- #
# S6. Tipping point analysis on Pub1 . data
#     Shift Pub1 . pdpch by delta in {-0.20, ..., +0.20} and
#     find the smallest |delta| at which any delta parameter's
#     95% CrI excludes zero.
# ---------------------------------------------------------------------------- #

cat("\n========== S6: TIPPING POINT ANALYSIS ==========\n")

delta_check_params <- c("de0", "demax", "ded50", "dr1", "dk")

# Shift = 0: reuse the primary fit
tipping_results <- data.frame()
shift_summary_0 <- summarise_mcmc(mcmc_primary, delta_check_params)
ethnic_diff_0   <- curves_primary$asian[,    which(x_days == 336)] -
                   curves_primary$nonasian[, which(x_days == 336)]
tipping_results <- rbind(tipping_results, data.frame(
  shift             = 0,
  de0_mean          = shift_summary_0["de0", "Mean"],
  de0_lower         = shift_summary_0["de0", "2.5%"],
  de0_upper         = shift_summary_0["de0", "97.5%"],
  demax_mean        = shift_summary_0["demax", "Mean"],
  demax_lower       = shift_summary_0["demax", "2.5%"],
  demax_upper       = shift_summary_0["demax", "97.5%"],
  db_mean           = shift_summary_0["dk", "Mean"],
  db_lower          = shift_summary_0["dk", "2.5%"],
  db_upper          = shift_summary_0["dk", "97.5%"],
  any_CrI_excl_zero = any(shift_summary_0[, "2.5%"] > 0 |
                          shift_summary_0[, "97.5%"] < 0),
  tipped_params     = paste(rownames(shift_summary_0)[
                              shift_summary_0[, "2.5%"] > 0 |
                              shift_summary_0[, "97.5%"] < 0],
                            collapse = ","),
  curve_diff_M12    = mean(ethnic_diff_0),
  stringsAsFactors  = FALSE
))

shift_values <- c(-0.20, -0.15, -0.10, -0.05, 0.05, 0.10, 0.15, 0.20)
for (delta_shift in shift_values) {
  cat("  Shift on Pub1 .:", delta_shift, "\n")
  dat_s <- y.pd
  dat_s$pdpch[dat_s$source == "Pub1"] <-
    dat_s$pdpch[dat_s$source == "Pub1"] + delta_shift
  inp_s    <- build_primary_inputs(dat_s)
  mcmc_s   <- run_nimble(code_primary,
                         inp_s$constants, inp_s$data, inp_s$inits)
  smry_s   <- summarise_mcmc(mcmc_s, delta_check_params)
  any_excl <- any(smry_s[, "2.5%"] > 0 | smry_s[, "97.5%"] < 0)
  tipped   <- rownames(smry_s)[smry_s[, "2.5%"] > 0 | smry_s[, "97.5%"] < 0]
  crv      <- compute_fitted_curves(mcmc_s, 336,
                                    base_pd_asian, base_pd_nonasian)

  tipping_results <- rbind(tipping_results, data.frame(
    shift             = delta_shift,
    de0_mean          = smry_s["de0", "Mean"],
    de0_lower         = smry_s["de0", "2.5%"],
    de0_upper         = smry_s["de0", "97.5%"],
    demax_mean        = smry_s["demax", "Mean"],
    demax_lower       = smry_s["demax", "2.5%"],
    demax_upper       = smry_s["demax", "97.5%"],
    db_mean           = smry_s["dk", "Mean"],
    db_lower          = smry_s["dk", "2.5%"],
    db_upper          = smry_s["dk", "97.5%"],
    any_CrI_excl_zero = any_excl,
    tipped_params     = ifelse(length(tipped) > 0,
                               paste(tipped, collapse = ","), "none"),
    curve_diff_M12    = mean(crv$asian[, 1] - crv$nonasian[, 1]),
    stringsAsFactors  = FALSE
  ))
}
tipping_results <- tipping_results[order(tipping_results$shift), ]

cat("\nTipping Point Analysis (Pub1 . shifted):\n")
cat("Criterion: 95% CrI of any delta parameter excludes zero\n\n")
print(tipping_results[, c("shift",
                          "de0_mean", "de0_lower", "de0_upper",
                          "demax_mean", "demax_lower", "demax_upper",
                          "db_mean", "db_lower", "db_upper",
                          "any_CrI_excl_zero", "tipped_params")],
      row.names = FALSE)

tp_idx <- which(tipping_results$any_CrI_excl_zero)
if (length(tp_idx) > 0) {
  tp_shift <- tipping_results$shift[min(tp_idx)]
  cat("\nTipping point reached at shift =", tp_shift, "\n",
      " Parameter(s):", tipping_results$tipped_params[min(tp_idx)], "\n")
} else {
  cat("\nNo tipping point reached within tested range [-0.20, +0.20].\n")
}


# ---------------------------------------------------------------------------- #
# S7. Alternative borrowing methods
#     7a. Power prior         (Pub1/Pub2 discounted by alpha0)
#     7b. Commensurate prior  (Pub1, Pub2 linked via tau_c)
#     7c. Robust MAP prior    (mixture protects against prior-data conflict)
# ---------------------------------------------------------------------------- #

cat("\n========== S7: ALTERNATIVE BORROWING METHODS ==========\n")

# ---- 7a. Power prior ------------------------------------------------------- #
# Split: IPD = Global-001 (full); RWE = full weight; Pub1/Pub2 = discounted
is_ipd   <- y.pd$source %in% c("Global_NonAsian", "Global_Asian")
is_rwe   <- y.pd$source == "RWE"
is_disc  <- y.pd$source %in% c("Pub1", "Pub2")
dat_ipd  <- y.pd[is_ipd,  ]
dat_rwe  <- y.pd[is_rwe,  ]
dat_disc <- y.pd[is_disc, ]
stopifnot(nrow(dat_ipd) + nrow(dat_rwe) + nrow(dat_disc) == nrow(y.pd))

cat("  Power-prior data split:\n")
cat("    IPD (Global-001):", nrow(dat_ipd), "rows\n")
cat("    RWE (full wt):", nrow(dat_rwe), "rows\n")
cat("    Pub1/Pub2 (disc):", nrow(dat_disc), "rows\n")

alpha_values <- c(0.25, 0.50, 0.75)
power_prior_results <- list()
power_prior_results[["1"]] <- primary_summary   # alpha0 = 1 reuses primary

for (a0 in alpha_values) {
  cat("  Power prior alpha0 =", a0, "\n")
  pp_c <- list(N_ipd = nrow(dat_ipd), N_rwe = nrow(dat_rwe),
               N_disc = nrow(dat_disc), alpha0 = a0)
  pp_d <- list(y         = dat_ipd$pdpch,    x         = dat_ipd$visit,
               z         = dat_ipd$base.pd.cent, race  = dat_ipd$race,
               y_rwe     = dat_rwe$pdpch,    x_rwe     = dat_rwe$visit,
               z_rwe     = dat_rwe$base.pd.cent, race_rwe = dat_rwe$race,
               n_rwe     = dat_rwe$n,
               y_disc    = dat_disc$pdpch,   x_disc    = dat_disc$visit,
               z_disc    = dat_disc$base.pd.cent, race_disc = dat_disc$race,
               n_disc    = dat_disc$n)
  pp_i <- list(e0 = -0.6, de0 = 0, emax = 0.3, demax = 0,
               ed50 = 10, ded50 = 1, r1 = 1, dr1 = 0.1,
               k = -0.0003, dk = 0, a = -0.02, tau2 = 0.05)
  mcmc_pp <- run_nimble(code_power_prior, pp_c, pp_d, pp_i)
  power_prior_results[[as.character(a0)]] <- summarise_mcmc(mcmc_pp, delta_params)
}

cat("\nPower Prior - dE0 and dEmax by alpha0:\n")
pp_table_alpha <- c(alpha_values, 1.00)
pp_table <- data.frame(
  alpha0     = pp_table_alpha,
  dE0_mean   = sapply(as.character(pp_table_alpha),
                      function(a) power_prior_results[[a]]["de0", "Mean"]),
  dE0_lower  = sapply(as.character(pp_table_alpha),
                      function(a) power_prior_results[[a]]["de0", "2.5%"]),
  dE0_upper  = sapply(as.character(pp_table_alpha),
                      function(a) power_prior_results[[a]]["de0", "97.5%"]),
  dEmax_mean = sapply(as.character(pp_table_alpha),
                      function(a) power_prior_results[[a]]["demax", "Mean"]),
  dEmax_lower= sapply(as.character(pp_table_alpha),
                      function(a) power_prior_results[[a]]["demax", "2.5%"]),
  dEmax_upper= sapply(as.character(pp_table_alpha),
                      function(a) power_prior_results[[a]]["demax", "97.5%"])
)
print_rounded(pp_table, 4)


# ---- 7b. Commensurate prior ------------------------------------------------ #
cat("  Running: Commensurate prior\n")

dat_ipd_comm <- y.pd[y.pd$source %in% c("Global_NonAsian", "Global_Asian", "RWE"), ]
dat_agg_comm <- y.pd[y.pd$source %in% c("Pub1", "Pub2"), ]
source_ids   <- ifelse(dat_agg_comm$source == "Pub1", 1, 2)
stopifnot(nrow(dat_ipd_comm) + nrow(dat_agg_comm) == nrow(y.pd))

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
                   tau_c    = 1,
                   e0_agg   = c(-0.6, -0.6),
                   emax_agg = c( 0.3,  0.3))
comm_monitors <- c(delta_params, "tau_c", "e0_agg", "emax_agg")

mcmc_comm <- run_nimble(code      = code_commensurate,
                        constants = comm_constants,
                        data      = comm_data,
                        inits     = comm_inits,
                        monitors  = comm_monitors)

comm_summary <- summarise_mcmc(mcmc_comm, c(delta_params, "tau_c"))
cat("\nCommensurate Prior - Posterior Summary:\n")
print_rounded(comm_summary, 4)


# ---- 7c. Robust MAP prior -------------------------------------------------- #
cat("  Running: Robust MAP prior\n")

dat_ipd_map   <- y.pd[y.pd$source %in% c("Global_NonAsian", "Global_Asian", "RWE"), ]
dat_ext_map   <- y.pd[y.pd$source %in% c("Pub1", "Pub2"), ]
ext_source_id <- ifelse(dat_ext_map$source == "Pub1", 1, 2)

w_robust_val <- 0.2
sd_vague_val <- 10

map_constants <- list(N_ipd         = nrow(dat_ipd_map),
                      N_agg         = nrow(dat_ext_map),
                      N_ext_sources = 2,
                      agg_source_id = ext_source_id)
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
map_inits <- list(e0 = -0.6, emax = 0.3, de0 = 0, demax = 0,
                  k = -0.0003, dk = 0, a = -0.02, tau2 = 0.05,
                  ed50 = 10, ded50 = 1, r1 = 1, dr1 = 0.1,
                  de0_ext   = c(0, 0),
                  demax_ext = c(0, 0),
                  de0_pop = 0, demax_pop = 0, tau_map = 0.3)
map_monitors <- c("e0", "emax", "ed50", "r1", "k", "dk", "ded50", "dr1", "a",
                  "de0", "demax", "de0_ext", "demax_ext",
                  "de0_pop", "demax_pop", "tau_map", "tau2")

mcmc_map <- run_nimble(code      = code_robust_map,
                       constants = map_constants,
                       data      = map_data,
                       inits     = map_inits,
                       monitors  = map_monitors)

map_summary <- summarise_mcmc(mcmc_map)
cat("\nRobust MAP - Posterior Summary:\n")
print_rounded(map_summary, 4)

cat("\nRobust MAP - Key ethnic difference parameters:\n")
print_rounded(map_summary[c("de0", "demax",
                            "de0_pop", "demax_pop", "tau_map"), ,
                          drop = FALSE], 4)


# ---- S7 summary comparison ------------------------------------------------- #
cat("\n========== S7 SUMMARY: BORROWING METHODS COMPARISON ==========\n")
cat("                     dE0 (or equiv.)      dEmax (or equiv.)\n")
cat("Method               Mean   [95% CrI]     Mean   [95% CrI]\n")
cat("---------------------------------------------------------------\n")

cat(sprintf("Primary (alpha=1)   %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            primary_summary["de0",   "Mean"],
            primary_summary["de0",   "2.5%"],
            primary_summary["de0",   "97.5%"],
            primary_summary["demax", "Mean"],
            primary_summary["demax", "2.5%"],
            primary_summary["demax", "97.5%"]))

pp05 <- power_prior_results[["0.5"]]
cat(sprintf("Power (alpha=0.5)   %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            pp05["de0", "Mean"], pp05["de0", "2.5%"], pp05["de0", "97.5%"],
            pp05["demax", "Mean"], pp05["demax", "2.5%"], pp05["demax", "97.5%"]))

pp025 <- power_prior_results[["0.25"]]
cat(sprintf("Power (alpha=0.25)  %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            pp025["de0", "Mean"], pp025["de0", "2.5%"], pp025["de0", "97.5%"],
            pp025["demax", "Mean"], pp025["demax", "2.5%"], pp025["demax", "97.5%"]))

cat(sprintf("Commensurate        %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            comm_summary["de0",   "Mean"],
            comm_summary["de0",   "2.5%"],
            comm_summary["de0",   "97.5%"],
            comm_summary["demax", "Mean"],
            comm_summary["demax", "2.5%"],
            comm_summary["demax", "97.5%"]))

# Pop = MAP population mean across Pub1+Pub2;
# CS  = current-study (Global-001 + RWE) ethnic offset under the robust mixture
cat(sprintf("Robust MAP (pop)    %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            map_summary["de0_pop",   "Mean"],
            map_summary["de0_pop",   "2.5%"],
            map_summary["de0_pop",   "97.5%"],
            map_summary["demax_pop", "Mean"],
            map_summary["demax_pop", "2.5%"],
            map_summary["demax_pop", "97.5%"]))

cat(sprintf("Robust MAP (CS)     %6.3f [%6.3f,%6.3f]  %6.3f [%6.3f,%6.3f]\n",
            map_summary["de0",   "Mean"],
            map_summary["de0",   "2.5%"],
            map_summary["de0",   "97.5%"],
            map_summary["demax", "Mean"],
            map_summary["demax", "2.5%"],
            map_summary["demax", "97.5%"]))

cat("---------------------------------------------------------------\n")
cat("tau_map (between-source SD):", round(map_summary["tau_map", "Mean"], 4),
    "[", round(map_summary["tau_map", "2.5%"], 4), ",",
    round(map_summary["tau_map", "97.5%"], 4), "]\n")
cat("tau_c   (commensurability):", round(comm_summary["tau_c", "Mean"], 4),
    "[", round(comm_summary["tau_c", "2.5%"], 4), ",",
    round(comm_summary["tau_c", "97.5%"], 4), "]\n")


# ---------------------------------------------------------------------------- #
# S8. Structural model sensitivity + WAIC comparison
# ---------------------------------------------------------------------------- #

cat("\n========== S8: STRUCTURAL MODEL SENSITIVITY ==========\n")

# ---- S8a. Emax-only model ------------------------------------------------- #
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
  monitors  = emax_monitors,
  enableWAIC = TRUE
)
emax_summary <- summarise_mcmc(mcmc_emax, emax_monitors)
cat("\nEmax-Only - Posterior Summary:\n")
print_rounded(emax_summary, 4)


# ---- S8b. Piecewise linear model (3 knots: Day 1, 28, 84) ----------------- #
cat("  Running: Piecewise linear model (3 knots at Day 1, 28, 84)\n")
pw_monitors <- c("e0", "de0", "b1", "db1", "b2", "db2",
                 "b3", "db3", "b4", "db4", "a", "tau2")
mcmc_pw <- run_nimble(
  code      = code_piecewise,
  constants = inputs$constants,
  data      = inputs$data,
  inits     = list(e0 = -0.6, de0 = 0,
                   b1 = 0.01,   db1 = 0,
                   b2 = 0.005,  db2 = 0,
                   b3 = -0.0003, db3 = 0,
                   b4 = -0.0003, db4 = 0,
                   a = -0.02, tau2 = 0.05),
  monitors  = pw_monitors,
  enableWAIC = TRUE
)
pw_summary <- summarise_mcmc(mcmc_pw, pw_monitors)
cat("\nPiecewise Linear - Posterior Summary:\n")
print_rounded(pw_summary, 4)


# ---- WAIC comparison ------------------------------------------------------- #
cat("\n--- Model Comparison (WAIC) ---\n")
waic_primary <- attr(mcmc_primary, "WAIC")
waic_emax    <- attr(mcmc_emax,    "WAIC")
waic_pw      <- attr(mcmc_pw,      "WAIC")

if (!is.null(waic_primary) && !is.null(waic_emax) && !is.null(waic_pw)) {
  waic_table <- data.frame(
    Model = c("Primary (Emax + linear)", "Emax-only", "Piecewise linear (3 knots)"),
    WAIC  = c(waic_primary$WAIC,  waic_emax$WAIC,  waic_pw$WAIC),
    pWAIC = c(waic_primary$pWAIC, waic_emax$pWAIC, waic_pw$pWAIC)
  )
  waic_table$dWAIC <- waic_table$WAIC - min(waic_table$WAIC)
  print_rounded(waic_table, 2)
  cat("  Lower WAIC is preferred. dWAIC = difference from best model.\n")
  cat("  pWAIC = effective number of parameters.\n")
} else {
  cat("  WAIC not available for all models (requires nimble >= 0.12.0).\n")
}


################################################################################
# 7. CONVERGENCE DIAGNOSTICS
################################################################################

cat("\n========== CONVERGENCE DIAGNOSTICS ==========\n")

if (length(mcmc_primary) >= 2) {
  gr <- gelman.diag(mcmc_primary, multivariate = FALSE)
  cat("\nGelman-Rubin statistics (Primary Analysis):\n")
  print_rounded(gr$psrf, 3)

  cat("\nEffective sample sizes (Primary Analysis):\n")
  print_rounded(effectiveSize(mcmc_primary), 0)
}


################################################################################
# 8. SAVE ALL RESULTS
################################################################################

results_all <- list(
  primary         = list(summary = primary_summary, mcmc = mcmc_primary),
  S1_ppc          = ppc_results,
  S2_loso         = loso_results,
  S3_noborrow     = list(summary = nb_summary, ci_comparison = ci_comparison),
  S4_ess          = ess_results,
  S5_prior        = list(original    = primary_summary,
                         vague       = vague_summary,
                         informative = informative_summary),
  S6_tipping      = tipping_results,
  S7_power        = power_prior_results,
  S7_commensurate = comm_summary,
  S7_robust_map   = map_summary,
  S8_emax_only    = emax_summary,
  S8_piecewise    = pw_summary
)

save(results_all, file = "Sensitivity_results.RData")
cat("\nResults saved to: Sensitivity_results.RData\n")
cat("\n========== ALL ANALYSES COMPLETE ==========\n")
