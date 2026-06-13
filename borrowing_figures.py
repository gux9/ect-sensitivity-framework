"""
Bayesian Borrowing Methods — Illustration Figures
Generates four PNG files:
  borrowing_power_prior.png
  borrowing_commensurate.png
  borrowing_predictive.png
  borrowing_comparison.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from scipy.stats import norm

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
try:
    plt.style.use("seaborn-v0_8-whitegrid")
except OSError:
    try:
        plt.style.use("seaborn-whitegrid")
    except OSError:
        plt.style.use("ggplot")

COLORS = {
    "hist": "#4C72B0",
    "alpha025": "#DD8452",
    "alpha05": "#55A868",
    "alpha10": "#C44E52",
    "flat": "#8172B2",
    "tau_low": "#937860",
    "tau_med": "#DA8BC3",
    "tau_high": "#8C8C8C",
    "comp1": "#4C72B0",
    "comp2": "#DD8452",
    "comp3": "#55A868",
    "sum_": "#C44E52",
}

theta = np.linspace(-5, 5, 1000)


# ===========================================================================
# Figure 1 — Power Prior
# ===========================================================================
def make_power_prior():
    fig, ax = plt.subplots(figsize=(8, 5))

    # Historical likelihood: N(0,1)
    hist_lik = norm.pdf(theta, 0, 1)
    hist_lik /= hist_lik.max()  # normalise to 1 for display

    # Flat prior (constant)
    flat = np.full_like(theta, 0.05)

    ax.plot(theta, hist_lik, color=COLORS["hist"], lw=2.5,
            label=r"Historical likelihood $L(\theta\,|\,D_0)$", zorder=5)

    alphas = [0.25, 0.5, 1.0]
    alpha_colors = [COLORS["alpha025"], COLORS["alpha05"], COLORS["alpha10"]]
    alpha_labels = [r"$\alpha=0.25$", r"$\alpha=0.50$", r"$\alpha=1.00$"]

    for a, c, lbl in zip(alphas, alpha_colors, alpha_labels):
        powered = hist_lik ** a
        powered /= powered.max()
        ax.plot(theta, powered, color=c, lw=2, linestyle="--", label=lbl)

    ax.plot(theta, flat, color=COLORS["flat"], lw=1.5, linestyle=":",
            label=r"Flat prior $\pi_0(\theta)$")

    # Annotations
    ax.annotate(r"$\alpha=0$ $\rightarrow$ no borrowing",
                xy=(-4.5, 0.05), fontsize=9, color="gray",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
    ax.annotate(r"$\alpha=1$ $\rightarrow$ full borrowing",
                xy=(1.5, 0.85), fontsize=9, color=COLORS["alpha10"],
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COLORS["alpha10"], alpha=0.8))

    ax.set_xlabel(r"$\theta$ (treatment effect)", fontsize=11)
    ax.set_ylabel("Density (unnormalized)", fontsize=11)
    ax.set_title(r"Power Prior:  $\pi(\theta\,|\,D_0,\alpha) \propto L(\theta\,|\,D_0)^\alpha \cdot \pi_0(\theta)$",
                 fontsize=12, pad=12)
    ax.legend(fontsize=9, loc="upper right")
    ax.set_xlim(-5, 5)
    ax.set_ylim(-0.02, 1.15)

    fig.tight_layout()
    fname = "borrowing_power_prior.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Done: {fname}")


# ===========================================================================
# Figure 2 — Commensurate Prior
# ===========================================================================
def make_commensurate():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ---- Left panel: commensurate prior for different τ ----
    ax = axes[0]
    theta_hist = 0.0
    tau_vals = [0.5, 2.0, 10.0]
    tau_colors = [COLORS["tau_low"], COLORS["tau_med"], COLORS["tau_high"]]
    tau_labels = [
        r"$\tau=0.5$ (low commensurability)",
        r"$\tau=2$",
        r"$\tau=10$ (high commensurability)",
    ]

    for tau, c, lbl in zip(tau_vals, tau_colors, tau_labels):
        sd = 1.0 / np.sqrt(tau)
        dens = norm.pdf(theta, theta_hist, sd)
        ax.plot(theta, dens, color=c, lw=2, label=lbl)
        ax.fill_between(theta, dens, alpha=0.08, color=c)

    ax.axvline(theta_hist, color="black", lw=1, linestyle="--", alpha=0.5)
    ax.set_xlabel(r"$\theta_{\mathrm{current}}$", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title(r"Commensurate prior on $\theta_{\mathrm{current}}\,|\,\theta_{\mathrm{hist}},\tau$",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.set_xlim(-5, 5)

    # ---- Right panel: θ_hist vs θ_current under high/low τ ----
    ax2 = axes[1]

    # θ_hist distribution (gray)
    theta_hist_val = 0.0
    hist_dist = norm.pdf(theta, theta_hist_val, 0.6)
    ax2.plot(theta, hist_dist, color="gray", lw=2, linestyle="-",
             label=r"$\theta_{\mathrm{hist}}$ distribution")
    ax2.fill_between(theta, hist_dist, alpha=0.15, color="gray")

    # θ_current under high τ (close to θ_hist) — e.g. mean shifted slightly
    high_tau_dist = norm.pdf(theta, 0.15, 0.7)
    ax2.plot(theta, high_tau_dist, color=COLORS["tau_high"], lw=2.5,
             label=r"$\theta_{\mathrm{current}}$, high $\tau$ (strong borrowing)")
    ax2.fill_between(theta, high_tau_dist, alpha=0.12, color=COLORS["tau_high"])

    # θ_current under low τ (diffuse)
    low_tau_dist = norm.pdf(theta, 0.5, 2.0)
    ax2.plot(theta, low_tau_dist, color=COLORS["tau_low"], lw=2.5,
             label=r"$\theta_{\mathrm{current}}$, low $\tau$ (weak borrowing)")
    ax2.fill_between(theta, low_tau_dist, alpha=0.12, color=COLORS["tau_low"])

    ax2.annotate("Large $\\tau$ → strong borrowing\nSmall $\\tau$ → weak borrowing",
                 xy=(2.0, 0.45), fontsize=9,
                 bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", ec="orange", alpha=0.9))

    ax2.set_xlabel(r"$\theta$", fontsize=11)
    ax2.set_ylabel("Density", fontsize=11)
    ax2.set_title("Effect of τ on borrowing strength", fontsize=11)
    ax2.legend(fontsize=8, loc="upper left")
    ax2.set_xlim(-5, 5)

    fig.suptitle(
        r"Commensurate Prior:  $\theta_{\mathrm{current}}\,|\,\theta_{\mathrm{hist}},\tau \;\sim\; \mathcal{N}(\theta_{\mathrm{hist}},\,1/\tau)$",
        fontsize=13, y=1.02)
    fig.tight_layout()
    fname = "borrowing_commensurate.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Done: {fname}")


# ===========================================================================
# Figure 3 — MAP / Bayesian Predictive Prior
# ===========================================================================
def make_predictive():
    fig = plt.figure(figsize=(10, 7))

    # ---- Top portion: flow diagram ----
    ax_flow = fig.add_axes([0.0, 0.52, 1.0, 0.45])
    ax_flow.set_xlim(0, 10)
    ax_flow.set_ylim(0, 2)
    ax_flow.axis("off")

    box_style = dict(boxstyle="round,pad=0.5", fc="#EAF2FB", ec="#2E86C1", lw=2)
    arrow_kw = dict(arrowstyle="-|>", color="#2E86C1", lw=2,
                    mutation_scale=18)

    # Box 1
    ax_flow.text(1.5, 1.0, "Historical studies\n$D_1, \\ldots, D_k$",
                 ha="center", va="center", fontsize=11, bbox=box_style)
    # Arrow 1→2
    ax_flow.annotate("", xy=(3.6, 1.0), xytext=(2.7, 1.0),
                     arrowprops=arrow_kw)
    # Box 2
    ax_flow.text(4.85, 1.0, "Hierarchical\nmodel",
                 ha="center", va="center", fontsize=11, bbox=box_style)
    # Arrow 2→3
    ax_flow.annotate("", xy=(6.5, 1.0), xytext=(6.1, 1.0),
                     arrowprops=arrow_kw)
    # Box 3
    ax_flow.text(8.0, 1.0, "Predictive prior\n$\\pi(\\theta_{\\mathrm{new}})$",
                 ha="center", va="center", fontsize=11, bbox=box_style)

    # Marginalisation annotation
    ax_flow.text(5.0, 0.15,
                 r"$\pi(\theta_{\mathrm{new}}) = \int p(\theta_{\mathrm{new}}\,|\,\mu)\;p(\mu\,|\,D_{\mathrm{hist}})\;d\mu$",
                 ha="center", va="center", fontsize=12,
                 bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", ec="darkorange", lw=1.5))

    ax_flow.set_title("Meta-Analytic Predictive (MAP) Prior", fontsize=14, pad=8)

    # ---- Bottom portion: mixture components ----
    ax_mix = fig.add_axes([0.08, 0.05, 0.88, 0.42])

    means = [-0.8, 0.3, 1.2]
    sds = [0.5, 0.6, 0.4]
    weights = [0.3, 0.45, 0.25]
    comp_colors = [COLORS["comp1"], COLORS["comp2"], COLORS["comp3"]]
    comp_labels = [f"Component {i+1}: $w_{i+1}={w:.2f}$, $\\mu_{i+1}={m}$"
                   for i, (w, m) in enumerate(zip(weights, means))]

    mixture = np.zeros_like(theta)
    for w, m, s, c, lbl in zip(weights, means, sds, comp_colors, comp_labels):
        comp = w * norm.pdf(theta, m, s)
        ax_mix.plot(theta, comp, color=c, lw=2, linestyle="--", label=lbl)
        ax_mix.fill_between(theta, comp, alpha=0.12, color=c)
        mixture += comp

    ax_mix.plot(theta, mixture, color=COLORS["sum_"], lw=3,
                label="MAP prior (mixture sum)")
    ax_mix.fill_between(theta, mixture, alpha=0.10, color=COLORS["sum_"])

    ax_mix.set_xlabel(r"$\theta_{\mathrm{new}}$ (treatment effect)", fontsize=11)
    ax_mix.set_ylabel("Density", fontsize=11)
    ax_mix.set_title("Gaussian mixture approximation to MAP prior", fontsize=11)
    ax_mix.legend(fontsize=9, loc="upper right")
    ax_mix.set_xlim(-3, 3)

    fname = "borrowing_predictive.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Done: {fname}")


# ===========================================================================
# Figure 4 — Comparison (3 panels)
# ===========================================================================
def make_comparison():
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle("Comparison of Bayesian Borrowing Methods", fontsize=14, y=1.01)

    # ---- Panel 1: Power Prior ----
    ax = axes[0]
    hist_lik = norm.pdf(theta, 0, 1)
    hist_lik /= hist_lik.max()
    ax.plot(theta, hist_lik, color=COLORS["hist"], lw=2,
            label=r"$L(\theta|D_0)$", zorder=5)
    for a, c, lbl in zip([0.1, 0.5, 1.0],
                          [COLORS["alpha025"], COLORS["alpha05"], COLORS["alpha10"]],
                          [r"$\alpha=0.1$", r"$\alpha=0.5$", r"$\alpha=1.0$"]):
        p = hist_lik ** a
        p /= p.max()
        ax.plot(theta, p, color=c, lw=1.8, linestyle="--", label=lbl)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-0.05, 1.2)
    ax.set_xlabel(r"$\theta$", fontsize=10)
    ax.set_ylabel("Density (unnorm.)", fontsize=10)
    ax.set_title("Power Prior", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.text(0.5, -0.18,
            r"$\pi(\theta|D_0,\alpha)\propto L(\theta|D_0)^\alpha\cdot\pi_0(\theta)$",
            ha="center", transform=ax.transAxes, fontsize=8.5,
            bbox=dict(fc="lightyellow", ec="orange", boxstyle="round,pad=0.3"))

    # ---- Panel 2: Commensurate ----
    ax2 = axes[1]
    theta_hist_val = 0.0
    hist_d = norm.pdf(theta, theta_hist_val, 0.6)
    ax2.plot(theta, hist_d / hist_d.max(), color="gray", lw=2, linestyle="-",
             label=r"$\theta_{\mathrm{hist}}$")
    ax2.fill_between(theta, hist_d / hist_d.max(), alpha=0.10, color="gray")

    for tau, c, lbl in zip([0.5, 10.0],
                             [COLORS["tau_low"], COLORS["tau_high"]],
                             [r"Low $\tau$ (weak)", r"High $\tau$ (strong)"]):
        sd = 1.0 / np.sqrt(tau)
        d = norm.pdf(theta, theta_hist_val + 0.1, sd)
        d /= d.max()
        ax2.plot(theta, d, color=c, lw=2, label=lbl)
        ax2.fill_between(theta, d, alpha=0.10, color=c)

    ax2.set_xlim(-4, 4)
    ax2.set_ylim(-0.05, 1.2)
    ax2.set_xlabel(r"$\theta_{\mathrm{current}}$", fontsize=10)
    ax2.set_ylabel("Density (norm.)", fontsize=10)
    ax2.set_title("Commensurate Prior", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=7.5)
    ax2.text(0.5, -0.18,
             r"$\theta_{\mathrm{cur}}|\theta_{\mathrm{hist}},\tau\sim\mathcal{N}(\theta_{\mathrm{hist}},1/\tau)$",
             ha="center", transform=ax2.transAxes, fontsize=8.5,
             bbox=dict(fc="lightyellow", ec="orange", boxstyle="round,pad=0.3"))

    # ---- Panel 3: MAP Prior ----
    ax3 = axes[2]
    means = [-0.8, 0.3, 1.2]
    sds = [0.5, 0.6, 0.4]
    weights = [0.3, 0.45, 0.25]
    comp_colors = [COLORS["comp1"], COLORS["comp2"], COLORS["comp3"]]
    mixture = np.zeros_like(theta)
    for w, m, s, c in zip(weights, means, sds, comp_colors):
        comp = w * norm.pdf(theta, m, s)
        ax3.plot(theta, comp, color=c, lw=1.5, linestyle="--")
        ax3.fill_between(theta, comp, alpha=0.12, color=c)
        mixture += comp
    ax3.plot(theta, mixture, color=COLORS["sum_"], lw=2.5,
             label="MAP prior")
    ax3.fill_between(theta, mixture, alpha=0.12, color=COLORS["sum_"])

    # Illustrate resulting posterior (current data at 0.8, sd=0.4)
    current_lik = norm.pdf(theta, 0.8, 0.4)
    post_unnorm = mixture * current_lik
    post_unnorm /= np.trapezoid(post_unnorm, theta)
    ax3.plot(theta, post_unnorm, color="black", lw=2, linestyle="-.",
             label="Posterior")

    ax3.set_xlim(-3, 3)
    ax3.set_xlabel(r"$\theta_{\mathrm{new}}$", fontsize=10)
    ax3.set_ylabel("Density", fontsize=10)
    ax3.set_title("MAP Prior", fontsize=11, fontweight="bold")
    ax3.legend(fontsize=7.5)
    ax3.text(0.5, -0.18,
             r"$\pi(\theta_{\mathrm{new}})=\int p(\theta_{\mathrm{new}}|\mu)\,p(\mu|D_{\mathrm{hist}})\,d\mu$",
             ha="center", transform=ax3.transAxes, fontsize=8.5,
             bbox=dict(fc="lightyellow", ec="orange", boxstyle="round,pad=0.3"))

    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fname = "borrowing_comparison.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Done: {fname}")


# ===========================================================================
# Main
# ===========================================================================
if __name__ == "__main__":
    import os
    os.chdir("/home/user/ect-sensitivity-framework")
    make_power_prior()
    make_commensurate()
    make_predictive()
    make_comparison()
