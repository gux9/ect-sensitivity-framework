# Sensitivity Analysis Framework for Externally Controlled Trials

R code and simulated data accompanying the manuscript:

> **A Practical Framework for Sensitivity Analysis in Externally Controlled Trials:
> An Illustration with a Bayesian Hybrid Evidence Synthesis Case Study**  
> *Submitted to [Statistics in Biopharmaceutical Research], Special Issue on Externally Controlled Trials, 2026*

---

## Overview

This repository provides a fully reproducible implementation of a three-pillar
sensitivity analysis framework for evaluating information borrowing in externally
controlled trials (ECTs) and hybrid evidence synthesis. The framework is organized
around three questions:

- **Pillar 1 — Appropriateness**: Was the borrowing justified given the available data?
- **Pillar 2 — Value**: Did the external data contribute meaningful precision?
- **Pillar 3 — Robustness**: Are the conclusions stable under reasonable perturbations?

Eight modular sensitivity analyses (S1–S8) are implemented:

| Analysis | Pillar | Description |
|----------|--------|-------------|
| S1 | Appropriateness | Heterogeneity diagnostics (posterior predictive checks) |
| S2 | Appropriateness | Leave-one-source-out analysis |
| S6 | Appropriateness | Tipping point analysis |
| S3 | Value | No-borrowing reference (IPD only) |
| S4 | Value | Effective sample size decomposition |
| S5 | Robustness | Prior sensitivity (original / vague / informative) |
| S7 | Robustness | Alternative borrowing methods (power prior, commensurate prior, robust MAP prior) |
| S8 | Robustness | Structural model sensitivity (Emax-only, piecewise linear) |

The worked example uses **simulated data** that mimic a hybrid evidence synthesis
combining individual patient data (IPD) from a pivotal study and a real-world
evidence (RWE) cohort, together with aggregate data from two published cohorts.
The primary model is a Bayesian longitudinal Emax-plus-linear model with
ethnic-difference parameters (ΔE₀, ΔEmax, ΔED₅₀, Δr, Δb).

---

## Repository Structure

```
.
├── Sensitivity_analysis.R   # Main analysis script (primary + S1–S8)
├── scenario.csv       # Simulated dataset (primary analysis)
└── README.md
```

### File descriptions

**`Sensitivity_analysis.R`**  
Main script implementing all eight sensitivity analyses plus the primary analysis.
Each analysis block is self-contained and labeled (PRIMARY, S1–S8). The script
exposes NIMBLE model code as named objects at the top of each block, and all
analyses return a standardized named list containing posterior summaries, MCMC
samples, convergence diagnostics, and LaTeX-ready table output.

**`XEN_nimble_realRWE_no102009.csv`**  
Simulated dataset used in the primary analysis and most sensitivity analyses.
Each row is a patient–visit observation with columns:
`pid`, `study`, `visit`, `n`, `race`, `base.pd`, `pd`, `pdch`, `pdpch`,
`source`, `base.pd.cent`.


---

## Requirements

### R version
R ≥ 4.0.0

### R packages

| Package | Purpose |
|---------|---------|
| `nimble` | MCMC via NIMBLE; implements all Bayesian models |
| `coda` | Convergence diagnostics (R-hat, ESS) |
| `ggplot2` | Posterior and fitted-curve plots |
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |

Install all dependencies with:

```r
install.packages(c("nimble", "coda", "ggplot2", "dplyr", "tidyr"))
```

---

## Usage

1. Clone or download this repository.

2. Open `Sensitivity_analysis.R` and update the working directory
   at the top of the script:

```r
work.dir <- 'path/to/your/local/folder'
```

3. By default, the script loads `senario.csv`. To use the
   sensitivity dataset, update the `read.csv()` call near line 213:

```r
y.iop <- read.csv("senario.csv")
```

4. Source the script:

```r
source("Sensitivity_analysis.R")
```

Each analysis block runs sequentially to the working directory. Console output 
includes convergence summaries (R-hat, ESS) for all monitored parameters.


### MCMC settings

The default chain configuration uses 10⁵ post-warmup iterations per chain
(consistent with the manuscript). Runtime on a standard desktop is approximately
4 hours for the full pipeline. A quick-test block at the bottom of the script
reduces iterations for development purposes.


---

## Notes

- All datasets in this repository are **simulated** and contain no real patient data.
- Source labels in the data files (`P13_NonAsian`, `P13_Asian`, `RWE`, `Hu`, `Sng`)
  are internal identifiers used for script compatibility. In the manuscript, these
  are masked as: Global non-Asian, Global Asian, RWE, Publication 1, Publication 2.
- The NIMBLE compiled object directory defaults to `nimble_compiled/` inside the
  working directory. This folder is created automatically on first run.

---

## Citation

If you use this code or framework, please cite:

> [Authors]. A Practical Framework for Sensitivity Analysis in Externally Controlled
> Trials: An Illustration with a Bayesian Hybrid Evidence Synthesis Case Study.
> *[Journal Name]*, 2026. [DOI to be added upon publication]

---

## License

[MIT License / or as appropriate]
