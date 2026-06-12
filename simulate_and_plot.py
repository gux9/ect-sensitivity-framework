"""
simulate_and_plot.py
====================
Read simulated_data.csv (produced by simulate_data.py) and produce a
publication-quality line plot (mean ± 1 SD band per data source).

Run:
    python3 simulate_and_plot.py

Input  (must exist):  simulated_data.csv
Output (kept local):  ect_source_trajectories.png
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ──────────────────────────────────────────────────────────────────────────────
# 1. Load data
# ──────────────────────────────────────────────────────────────────────────────

CSV_PATH = "simulated_data.csv"

if not pd.io.common.file_exists(CSV_PATH):
    sys.exit(f"ERROR: {CSV_PATH} not found. Run simulate_data.py first.")

simdata = pd.read_csv(CSV_PATH)
print(f"Loaded {CSV_PATH}: {len(simdata)} rows")
print(simdata.groupby(["source", "race"]).size().to_string())
print()

SOURCE_ORDER = ["Global_NonAsian", "Global_Asian", "RWE", "Pub1", "Pub2"]

# ──────────────────────────────────────────────────────────────────────────────
# 2. Colour / style palette
# ──────────────────────────────────────────────────────────────────────────────

PALETTE = {
    "Global_NonAsian": "#2166AC",   # blue
    "Global_Asian"   : "#D6604D",   # red-orange
    "RWE"            : "#4DAC26",   # green
    "Pub1"           : "#8B5CF6",   # purple
    "Pub2"           : "#F59E0B",   # amber
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Build legend labels dynamically from the data
# ──────────────────────────────────────────────────────────────────────────────

PREFIX = {
    "Global_NonAsian": "Global – Non-Asian",
    "Global_Asian"   : "Global – Asian",
    "RWE"            : "RWE – Asian",
    "Pub1"           : "Publication 1 – Asian",
    "Pub2"           : "Publication 2 – Asian",
}

def make_label(df, source):
    sub    = df[df["source"] == source]
    is_agg = sub["n"].max() > 1
    if is_agg:
        n = int(sub["n"].iloc[0])
        return f"{PREFIX[source]}  (n = {n} aggregate)"
    else:
        n = sub["pid"].nunique()
        return f"{PREFIX[source]}  (N = {n} IPD)"

LABELS = {src: make_label(simdata, src) for src in SOURCE_ORDER}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Summary statistics per source × visit
# ──────────────────────────────────────────────────────────────────────────────

def source_summary(df, source):
    """Return (mean_series, std_series) indexed by visit day."""
    sub    = df[df["source"] == source]
    grp    = sub.groupby("visit")["pdpch"]
    mean_v = grp.mean()
    is_agg = sub["n"].max() > 1
    if is_agg:
        # Aggregate row: SD = sqrt(tau2 / n) derived from within-visit spread
        std_v = grp.std().fillna(sub["pdpch"].std())
    else:
        std_v = grp.std()
    return mean_v, std_v

# ──────────────────────────────────────────────────────────────────────────────
# 5. Draw the plot
# ──────────────────────────────────────────────────────────────────────────────

all_visits = sorted(simdata["visit"].unique())

fig, ax = plt.subplots(figsize=(10, 6))

for src in SOURCE_ORDER:
    mean_v, std_v = source_summary(simdata, src)
    visits = mean_v.index.values
    m      = mean_v.values
    s      = std_v.values
    col    = PALETTE[src]
    is_agg = simdata.loc[simdata["source"] == src, "n"].max() > 1
    ls     = "--" if is_agg else "-"
    mk     = "s"  if is_agg else "o"
    mks    = 5    if is_agg else 4

    ax.fill_between(visits, m - s, m + s,
                    alpha=0.15, color=col, linewidth=0)
    ax.plot(visits, m, color=col, linewidth=2,
            linestyle=ls, marker=mk, markersize=mks,
            label=LABELS[src], zorder=3)

# Axes
ax.axhline(0, color="black", linewidth=0.8, linestyle=":", alpha=0.6)
ax.set_xlabel("Visit (day)", fontsize=12)
ax.set_ylabel("Mean % change from baseline (± 1 SD)", fontsize=12)
ax.set_title("Simulated longitudinal response by data source\n"
             "ECT sensitivity analysis framework — synthetic data",
             fontsize=13, fontweight="bold")
ax.set_xticks(all_visits)
ax.set_xticklabels(all_visits, rotation=45, ha="right")
ax.set_xlim(-5, 350)
ax.tick_params(labelsize=10)
ax.grid(axis="y", linestyle="--", alpha=0.35)

# Legend — IPD sources first, then aggregate
ipd_handles, agg_handles = [], []
for src in SOURCE_ORDER:
    is_agg = simdata.loc[simdata["source"] == src, "n"].max() > 1
    col = PALETTE[src]
    h = Line2D([0], [0], color=col, linewidth=2,
               linestyle="--" if is_agg else "-",
               marker="s" if is_agg else "o",
               markersize=6, label=LABELS[src])
    (agg_handles if is_agg else ipd_handles).append(h)

band_note = mpatches.Patch(facecolor="grey", alpha=0.25, label="±1 SD band")

ax.legend(handles=ipd_handles + agg_handles + [band_note],
          title="Data source", title_fontsize=10,
          fontsize=9, loc="lower right",
          framealpha=0.9, edgecolor="#cccccc")

plt.tight_layout()
plt.savefig("ect_source_trajectories.png", dpi=150, bbox_inches="tight")
print("Written: ect_source_trajectories.png")
plt.close()
