"""
simulate_data.py
================
Simulate a fully synthetic ECT longitudinal dataset and write it to
simulated_data.csv.  This script does NOT produce any plots.

Run:
    python3 simulate_data.py

Output (kept local, NOT committed):
    simulated_data.csv
"""

import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

rng = np.random.default_rng(2025)

# ──────────────────────────────────────────────────────────────────────────────
# 0. Emax + linear mean function  (manuscript Eq. 2)
# ──────────────────────────────────────────────────────────────────────────────

def emax_mu(x, race, z, p):
    r_i    = np.maximum(p["r1"]  + p["dr1"]  * race, 1e-6)
    ed50_i = np.maximum(p["ed50"] + p["ded50"] * race, 1e-6)
    xr     = np.power(x, r_i)
    d50r   = np.power(ed50_i, r_i)
    return (
        (p["e0"] + p["de0"] * race)
        + p["a"] * z
        + (p["k"] + p["dk"] * race) * x
        + (p["emax"] + p["demax"] * race) * xr / (d50r + xr)
    )


# ──────────────────────────────────────────────────────────────────────────────
# 1. True longitudinal parameters (from the published fitted curves)
# ──────────────────────────────────────────────────────────────────────────────

TRUE_PARAMS = {
    "e0"   : -0.66,
    "emax" :  0.42,
    "ed50" :  8.0,
    "r1"   :  1.3,
    "k"    : -0.00028,
    "a"    : -0.010,
    "de0"  : -0.04,
    "demax":  0.05,
    "ded50":  0.0,
    "dr1"  :  0.0,
    "dk"   : -0.00012,
    "tau2" :  0.0025,
}

print("True parameters:")
for k, v in TRUE_PARAMS.items():
    print(f"  {k:8s} = {v:.6f}")
print()

if os.path.exists("scenario.csv"):
    raw = pd.read_csv("scenario.csv", header=0)
    if "pdpch.1" in raw.columns:
        raw = raw.rename(columns={"pdpch": "pdch", "pdpch.1": "pdpch"})
    orig = raw.dropna(subset=["pdpch"]).copy()
    print(f"Reference scenario.csv loaded: {len(orig)} rows.\n")


# ──────────────────────────────────────────────────────────────────────────────
# 2. Simulation design
# ──────────────────────────────────────────────────────────────────────────────

N_GNA  = 100
N_GA   =  25
N_RWE  =  60
N_PUB1 = 100
N_PUB2 =  80

VISITS_GLOBAL = [1, 7, 14, 28, 84, 168, 224, 280, 336]
VISITS_RWE    = [1, 7, 28, 84, 168, 224, 280, 336]
VISITS_PUB1   = [1, 7, 28, 84, 168]
VISITS_PUB2   = [1, 7, 28, 84, 168, 336]

MISS_EARLY = 0.03
MISS_LATE  = 0.08

PUB1_BASE = 20.0
PUB2_BASE = 22.0

BASE_DISTS = {
    "Global_NonAsian": (24.0, 4.5),
    "Global_Asian"   : (23.5, 4.0),
    "RWE"            : (28.0, 9.0),
}


# ──────────────────────────────────────────────────────────────────────────────
# 3. Draw baselines and compute grand centering mean
# ──────────────────────────────────────────────────────────────────────────────

def draw_baselines(n, mean_b, sd_b):
    return np.maximum(rng.normal(mean_b, sd_b, n), 5.0)

rng = np.random.default_rng(2025)
base_gna = draw_baselines(N_GNA, *BASE_DISTS["Global_NonAsian"])
base_ga  = draw_baselines(N_GA,  *BASE_DISTS["Global_Asian"])
base_rwe = draw_baselines(N_RWE, *BASE_DISTS["RWE"])

all_ipd   = np.concatenate([base_gna, base_ga, base_rwe])
pub_bases = np.concatenate([np.full(N_PUB1, PUB1_BASE),
                            np.full(N_PUB2, PUB2_BASE)])
GLOBAL_BASE_MEAN = np.mean(np.concatenate([all_ipd, pub_bases]))
print(f"Global baseline mean for centering: {GLOBAL_BASE_MEAN:.4f}\n")


# ──────────────────────────────────────────────────────────────────────────────
# 4. IPD simulation helper
# ──────────────────────────────────────────────────────────────────────────────

def sim_ipd(n_patients, baselines, race_val, study_val, source_label,
            pid_prefix, visits):
    rows = []
    sd = np.sqrt(TRUE_PARAMS["tau2"])
    for i in range(n_patients):
        pid_i  = f"{pid_prefix}-{i+1:03d}"
        base_i = round(float(baselines[i]), 3)
        z_i    = base_i - GLOBAL_BASE_MEAN
        for v in visits:
            miss_p = MISS_EARLY if v <= 28 else MISS_LATE
            if rng.random() < miss_p:
                continue
            mu_i    = emax_mu(v, race_val, z_i, TRUE_PARAMS)
            pdpch_i = float(rng.normal(mu_i, sd))
            pd_i    = base_i * (1 + pdpch_i)
            rows.append({
                "pid"         : pid_i,
                "study"       : study_val,
                "visit"       : v,
                "n"           : 1,
                "race"        : race_val,
                "base.pd"     : base_i,
                "pd"          : round(pd_i,          4),
                "pdch"        : round(pd_i - base_i, 4),
                "pdpch"       : round(pdpch_i,        6),
                "source"      : source_label,
                "base.pd.cent": round(z_i,            6),
            })
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────────────────
# 5. Simulate IPD sources
# ──────────────────────────────────────────────────────────────────────────────

print("Simulating IPD data …")
dat_gna = sim_ipd(N_GNA, base_gna, 0, 0, "Global_NonAsian", "GNA", VISITS_GLOBAL)
dat_ga  = sim_ipd(N_GA,  base_ga,  1, 0, "Global_Asian",    "GA",  VISITS_GLOBAL)
dat_rwe = sim_ipd(N_RWE, base_rwe, 1, 1, "RWE",             "RWE", VISITS_RWE)

print(f"  Global_NonAsian : {len(dat_gna)} rows  ({N_GNA} patients)")
print(f"  Global_Asian    : {len(dat_ga)} rows  ({N_GA} patients)")
print(f"  RWE             : {len(dat_rwe)} rows  ({N_RWE} patients)")


# ──────────────────────────────────────────────────────────────────────────────
# 6. Simulate aggregate publication data
# ──────────────────────────────────────────────────────────────────────────────

def sim_aggregate(pub_label, study_val, n_agg, base_mean_pub, visits_pub):
    z_pub  = base_mean_pub - GLOBAL_BASE_MEAN
    sd_agg = np.sqrt(TRUE_PARAMS["tau2"] / n_agg)
    rows   = []
    for v in visits_pub:
        mu_v    = emax_mu(v, 1, z_pub, TRUE_PARAMS)
        pdpch_v = float(rng.normal(mu_v, sd_agg))
        pd_v    = base_mean_pub * (1 + pdpch_v)
        rows.append({
            "pid"         : pub_label,
            "study"       : study_val,
            "visit"       : v,
            "n"           : n_agg,
            "race"        : 1,
            "base.pd"     : base_mean_pub,
            "pd"          : round(pd_v,                  4),
            "pdch"        : round(pd_v - base_mean_pub,  4),
            "pdpch"       : round(pdpch_v,                6),
            "source"      : pub_label,
            "base.pd.cent": round(z_pub,                  6),
        })
    return pd.DataFrame(rows)

print("Simulating aggregate publication data …")
dat_pub1 = sim_aggregate("Pub1", 1, N_PUB1, PUB1_BASE, VISITS_PUB1)
dat_pub2 = sim_aggregate("Pub2", 2, N_PUB2, PUB2_BASE, VISITS_PUB2)

print(f"  Pub1 : {len(dat_pub1)} summary rows  (n = {N_PUB1})")
print(f"  Pub2 : {len(dat_pub2)} summary rows  (n = {N_PUB2})")


# ──────────────────────────────────────────────────────────────────────────────
# 7. Combine and write CSV
# ──────────────────────────────────────────────────────────────────────────────

SOURCE_ORDER = ["Global_NonAsian", "Global_Asian", "RWE", "Pub1", "Pub2"]
simdata = pd.concat([dat_gna, dat_ga, dat_rwe, dat_pub1, dat_pub2],
                    ignore_index=True)
simdata["source"] = pd.Categorical(simdata["source"], categories=SOURCE_ORDER)
simdata = simdata.sort_values(["source", "pid", "visit"]).reset_index(drop=True)
simdata["source"] = simdata["source"].astype(str)

print(f"\nFinal dataset: {len(simdata)} rows")
print(simdata.groupby(["source", "race"]).size().to_string())

simdata.to_csv("simulated_data.csv", index=False)
print("\nWritten: simulated_data.csv")
