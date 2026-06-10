"""
simulate_and_plot.py
====================
1. Fit the Emax + linear longitudinal model to scenario.csv via weighted
   least-squares (mirrors simulate_data.R).
2. Simulate fully synthetic data at increased sample sizes and write
   simulated_data.csv.
3. Produce a publication-quality line plot (mean ± 1 SD band) for all
   five data sources and save as ect_source_trajectories.png.

Output files (kept local, NOT committed to the repo):
    simulated_data.csv
    ect_source_trajectories.png
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")

rng = np.random.default_rng(2025)

# ──────────────────────────────────────────────────────────────────────────────
# 0. Emax + linear mean function  (manuscript Eq. 2)
# ──────────────────────────────────────────────────────────────────────────────

def emax_mu(x, race, z, p):
    """
    x    : visit day (scalar or array)
    race : 0 = Non-Asian, 1 = Asian
    z    : centred baseline PD
    p    : dict of model parameters
    """
    r_i    = np.maximum(p["r1"]  + p["dr1"]  * race, 1e-6)
    ed50_i = np.maximum(p["ed50"] + p["ded50"] * race, 1e-6)

    xr   = np.power(x, r_i)
    d50r = np.power(ed50_i, r_i)

    return (
        (p["e0"] + p["de0"] * race)
        + p["a"] * z
        + (p["k"] + p["dk"] * race) * x
        + (p["emax"] + p["demax"] * race) * xr / (d50r + xr)
    )


# ──────────────────────────────────────────────────────────────────────────────
# 1. Load original data and fit model
# ──────────────────────────────────────────────────────────────────────────────

print("Loading scenario.csv …")
raw = pd.read_csv("scenario.csv", header=0)

# The CSV has duplicate header names "pdpch,pdpch"; pandas appends ".1"
# Column 8 (pdpch)   = pdch  (absolute change)
# Column 9 (pdpch.1) = pdpch (% change) — the model outcome
rename = {}
if "pdpch.1" in raw.columns:
    rename["pdpch"]   = "pdch"
    rename["pdpch.1"] = "pdpch"
raw = raw.rename(columns=rename)

orig = raw.dropna(subset=["pdpch"]).copy()
print(f"  {len(orig)} rows after dropping NAs")
print(f"  Columns: {list(orig.columns)}")

# Weighted objective: sum n[i] * (y[i] - mu[i])^2
PARAM_NAMES = ["e0","emax","log_ed50","log_r1","k","a",
               "de0","demax","ded50","dr1","dk"]

def obj(par_vec):
    p = {
        "e0"   : par_vec[0],
        "emax" : par_vec[1],
        "ed50" : np.exp(par_vec[2]),
        "r1"   : np.exp(par_vec[3]),
        "k"    : par_vec[4],
        "a"    : par_vec[5],
        "de0"  : par_vec[6],
        "demax": par_vec[7],
        "ded50": par_vec[8],
        "dr1"  : par_vec[9],
        "dk"   : par_vec[10],
    }
    mu_hat = emax_mu(orig["visit"].values,
                     orig["race"].values,
                     orig["base.pd.cent"].values, p)
    resid = orig["pdpch"].values - mu_hat
    return float(np.sum(orig["n"].values * resid**2))

par0 = np.array([0.00, -0.45, np.log(0.5), np.log(1.2),
                 0.0008, -0.012, -0.03, -0.08, 0.0, 0.0, -0.0001])

print("Fitting model …")
res = minimize(obj, par0, method="Nelder-Mead",
               options={"maxiter": 200_000, "xatol": 1e-9, "fatol": 1e-9})

tp = res.x
TRUE_PARAMS = {
    "e0"   : tp[0],
    "emax" : tp[1],
    "ed50" : np.exp(tp[2]),
    "r1"   : np.exp(tp[3]),
    "k"    : tp[4],
    "a"    : tp[5],
    "de0"  : tp[6],
    "demax": tp[7],
    "ded50": tp[8],
    "dr1"  : tp[9],
    "dk"   : tp[10],
}

# Estimate tau2 from weighted residuals
mu_fit = emax_mu(orig["visit"].values, orig["race"].values,
                 orig["base.pd.cent"].values, TRUE_PARAMS)
resid  = orig["pdpch"].values - mu_fit
TRUE_PARAMS["tau2"] = float(
    np.sum(orig["n"].values * resid**2) / (orig["n"].sum() - 11)
)

print("Fitted parameters:")
for k, v in TRUE_PARAMS.items():
    print(f"  {k:8s} = {v:.6f}")
print(f"Convergence: {res.success}  (message: {res.message})\n")


# ──────────────────────────────────────────────────────────────────────────────
# 2. Simulation design
# ──────────────────────────────────────────────────────────────────────────────

N_GNA  = 100   # Global Non-Asian  IPD
N_GA   =  25   # Global Asian      IPD
N_RWE  =  60   # RWE Asian         IPD
N_PUB1 = 100   # Publication 1     aggregate
N_PUB2 =  80   # Publication 2     aggregate

VISITS_GLOBAL = [1, 7, 14, 28, 84, 168, 224, 280, 336]
VISITS_RWE    = [1, 7, 28, 84, 168, 224, 280, 336]
VISITS_PUB1   = [1, 7, 28, 84, 168]
VISITS_PUB2   = [1, 7, 28, 84, 168, 336]

MISS_EARLY = 0.03   # days 1–28
MISS_LATE  = 0.08   # days 84+

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

rng = np.random.default_rng(2025)   # reset for reproducibility
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

            mu_i     = emax_mu(v, race_val, z_i, TRUE_PARAMS)
            pdpch_i  = float(rng.normal(mu_i, sd))
            pd_i     = base_i * (1 + pdpch_i)
            pdch_i   = pd_i - base_i

            rows.append({
                "pid"         : pid_i,
                "study"       : study_val,
                "visit"       : v,
                "n"           : 1,
                "race"        : race_val,
                "base.pd"     : base_i,
                "pd"          : round(pd_i,    4),
                "pdch"        : round(pdch_i,  4),
                "pdpch"       : round(pdpch_i, 6),
                "source"      : source_label,
                "base.pd.cent": round(z_i,     6),
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
print(f"  Global_Asian    : {len(dat_ga)}  rows  ({N_GA} patients)")
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
        pdch_v  = pd_v - base_mean_pub
        rows.append({
            "pid"         : pub_label,
            "study"       : study_val,
            "visit"       : v,
            "n"           : n_agg,
            "race"        : 1,
            "base.pd"     : base_mean_pub,
            "pd"          : round(pd_v,    4),
            "pdch"        : round(pdch_v,  4),
            "pdpch"       : round(pdpch_v, 6),
            "source"      : pub_label,
            "base.pd.cent": round(z_pub,   6),
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
print(simdata.groupby(["source","race"]).size().to_string())

simdata.to_csv("simulated_data.csv", index=False)
print("\nWritten: simulated_data.csv")


# ──────────────────────────────────────────────────────────────────────────────
# 8. Plotting
# ──────────────────────────────────────────────────────────────────────────────

print("\nBuilding plot …")

# ── Colour / style palette ──────────────────────────────────────────────────
PALETTE = {
    "Global_NonAsian": "#2166AC",   # blue
    "Global_Asian"   : "#D6604D",   # red-orange
    "RWE"            : "#4DAC26",   # green
    "Pub1"           : "#8B5CF6",   # purple
    "Pub2"           : "#F59E0B",   # amber
}
LABELS = {
    "Global_NonAsian": "Global – Non-Asian  (N = 100 IPD)",
    "Global_Asian"   : "Global – Asian  (N = 25 IPD)",
    "RWE"            : "RWE – Asian  (N = 60 IPD)",
    "Pub1"           : "Publication 1 – Asian  (n = 100 aggregate)",
    "Pub2"           : "Publication 2 – Asian  (n = 80 aggregate)",
}

# ── Compute summary stats per source × visit ────────────────────────────────
def source_summary(df, source):
    sub = df[df["source"] == source]
    grp = sub.groupby("visit")["pdpch"]

    # For IPD sources: mean ± 1 SD across patients at each visit
    # For aggregate sources: the single summary value; SD from model SE
    mean_v = grp.mean()
    std_v  = grp.std()

    # For aggregate rows n > 1: the single row IS the mean; SD = sqrt(tau2/n)
    if sub["n"].max() > 1:
        n_val = sub["n"].iloc[0]
        std_v = np.sqrt(TRUE_PARAMS["tau2"] / n_val)
        std_v = pd.Series(std_v, index=mean_v.index)

    return mean_v, std_v

fig, ax = plt.subplots(figsize=(10, 6))

for src in SOURCE_ORDER:
    mean_v, std_v = source_summary(simdata, src)
    visits = mean_v.index.values
    m      = mean_v.values
    s      = std_v.values if hasattr(std_v, "values") else np.full_like(m, std_v)
    col    = PALETTE[src]

    # Band (±1 SD)
    ax.fill_between(visits, m - s, m + s,
                    alpha=0.15, color=col, linewidth=0)

    # Mean line
    is_aggregate = simdata.loc[simdata["source"] == src, "n"].max() > 1
    ls  = "--" if is_aggregate else "-"
    mk  = "s" if is_aggregate else "o"
    mks = 5   if is_aggregate else 4

    ax.plot(visits, m, color=col, linewidth=2,
            linestyle=ls, marker=mk, markersize=mks,
            label=LABELS[src], zorder=3)

# ── Axes formatting ─────────────────────────────────────────────────────────
ax.axhline(0, color="black", linewidth=0.8, linestyle=":", alpha=0.6)

ax.set_xlabel("Visit (day)", fontsize=12)
ax.set_ylabel("Mean % change from baseline (± 1 SD)", fontsize=12)
ax.set_title("Simulated longitudinal response by data source\n"
             "ECT sensitivity analysis framework — synthetic data",
             fontsize=13, fontweight="bold")

ax.set_xticks(sorted(set(VISITS_GLOBAL + VISITS_RWE + VISITS_PUB1 + VISITS_PUB2)))
ax.set_xticklabels(sorted(set(VISITS_GLOBAL + VISITS_RWE + VISITS_PUB1 + VISITS_PUB2)),
                   rotation=45, ha="right")
ax.set_xlim(-5, 350)
ax.tick_params(labelsize=10)
ax.grid(axis="y", linestyle="--", alpha=0.35)

# ── Legend (two groups: IPD vs aggregate) ───────────────────────────────────
ipd_handles  = []
agg_handles  = []
for src in SOURCE_ORDER:
    is_agg = simdata.loc[simdata["source"] == src, "n"].max() > 1
    col = PALETTE[src]
    ls  = "--" if is_agg else "-"
    mk  = "s" if is_agg else "o"
    h = Line2D([0], [0], color=col, linewidth=2, linestyle=ls,
               marker=mk, markersize=6, label=LABELS[src])
    (agg_handles if is_agg else ipd_handles).append(h)

band_note = mpatches.Patch(facecolor="grey", alpha=0.25,
                           label="±1 SD band")
all_handles = ipd_handles + agg_handles + [band_note]

ax.legend(handles=all_handles,
          title="Data source",
          title_fontsize=10,
          fontsize=9,
          loc="lower right",
          framealpha=0.9,
          edgecolor="#cccccc")

plt.tight_layout()
plt.savefig("ect_source_trajectories.png", dpi=150, bbox_inches="tight")
print("Written: ect_source_trajectories.png")
plt.close()

print("\nDone.")
