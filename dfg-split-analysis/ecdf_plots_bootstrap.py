"""Empirical CDFs with bootstrapped uncertainty bands.

Same three panels as ``ecdf_plots.py`` (net / ligand / protein charge
change, split by DFG-in vs DFG-out) but with a bootstrap-based estimate
of the noise around each ECDF.

Method
------
For each conformation group we draw ``N_BOOTSTRAP`` samples with
replacement of the same size as the observed data, compute the ECDF of
each resample on a shared ``x``-grid, then take the 2.5 / 50 / 97.5
percentiles of ``F_hat(x)`` across resamples.  The median curve is
plotted as the line; the 95 % band is shown as a shaded region at 50 %
opacity.  (Equivalently, this is the pointwise Efron non-parametric
bootstrap CI for the empirical CDF.)
"""

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%%
# --- Global style (match ecdf_plots.py) ---
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["axes.labelsize"] = 18
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.linewidth"] = 1.2
sns.set_style("white")

N_BOOTSTRAP = 2000
CI = 95  # pointwise confidence interval width (percent)
RNG = np.random.default_rng(20260413)

palette = {
    "DFG-in":  "tab:blue",
    "DFG-out": "tab:orange",
}
hue_order = ["DFG-in", "DFG-out"]


#%%
# Load data (same handling as ecdf_plots.py)
df = pd.read_csv(
    "./broken-chain-net-charge-changes-dfg-label-table.csv",
    encoding="latin-1",
)
df.columns = [c.replace("\xc6", "_") for c in df.columns]
df["Conformation"] = df["Kinase conformation"].str.extract(r"(DFG-in|DFG-out)")
df["Net_crg"] = df["Inhib_crg"] + df["Prot_crg"]


#%%
def _ecdf_on_grid(values, grid):
    """Evaluate F_hat(x) for an array of samples on a dense x-grid."""
    v = np.sort(np.asarray(values, dtype=float))
    v = v[~np.isnan(v)]
    if v.size == 0:
        return np.zeros_like(grid, dtype=float)
    # searchsorted gives the number of samples <= x for each grid point.
    return np.searchsorted(v, grid, side="right") / v.size


def _bootstrap_ecdf(values, grid, n_boot=N_BOOTSTRAP, rng=RNG):
    """Return (median, lower, upper) ECDF curves on ``grid``.

    Each resample has the same size as the input and is drawn with
    replacement.  The bands are the pointwise 2.5/97.5 percentiles
    (for CI=95).
    """
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    n = v.size
    if n == 0:
        empty = np.zeros_like(grid, dtype=float)
        return empty, empty, empty

    alpha = (100 - CI) / 2.0
    boot = np.empty((n_boot, grid.size), dtype=float)
    for b in range(n_boot):
        sample = rng.choice(v, size=n, replace=True)
        boot[b] = _ecdf_on_grid(sample, grid)

    lo = np.percentile(boot, alpha, axis=0)
    hi = np.percentile(boot, 100 - alpha, axis=0)
    med = np.percentile(boot, 50, axis=0)
    return med, lo, hi


def _plot_group_with_band(ax, values, grid, color, label):
    med, lo, hi = _bootstrap_ecdf(values, grid)
    # Point estimate (the observed ECDF, not the bootstrap median — they
    # agree in expectation but the observed one is the natural line).
    obs = _ecdf_on_grid(values, grid)
    ax.fill_between(grid, lo, hi, color=color, alpha=0.5, linewidth=0)
    ax.plot(grid, obs, color=color, linewidth=2.5, label=label)


def _plot_panel(df, value_col, xlim, xlabel, outpath, title=None):
    grid = np.linspace(xlim[0], xlim[1], 600)

    fig, ax = plt.subplots(figsize=(8, 5))
    for conf in hue_order:
        sub = df.loc[df["Conformation"] == conf, value_col].values
        _plot_group_with_band(ax, sub, grid, palette[conf], conf)

    ax.set_xlim(*xlim)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel("Empirical CDF", fontweight="bold")
    ax.axvline(0, color="grey", linestyle=":", linewidth=1)
    sns.despine(ax=ax, top=True, right=True)
    ax.tick_params(
        bottom=True, left=True, top=False, right=False, length=5, width=1,
    )
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight("bold")

    leg = ax.legend(frameon=False, fontsize=16, loc="best")
    for txt in leg.get_texts():
        txt.set_fontweight("bold")

    if title:
        ax.set_title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.show()


#%%
# ---- Plot 1: overall net charge change
_plot_panel(
    df,
    value_col="Net_crg",
    xlim=(-1.2, 0.6),
    xlabel=r"$\Delta$ net charge (inhibitor + kinase)",
    outpath="./broken_chain_ecdf_bootstrap_net_change_per_kinaseconf.pdf",
)

#%%
# ---- Plot 2: ligand (inhibitor) charge change
_plot_panel(
    df,
    value_col="Inhib_crg",
    xlim=(-0.9, 0.5),
    xlabel=r"$\Delta$Charge$_{inhibitor}$",
    outpath="./broken_chain_ecdf_bootstrap_inhibitor_change_per_kinaseconf.pdf",
)

#%%
# ---- Plot 3: protein (kinase) charge change
_plot_panel(
    df,
    value_col="Prot_crg",
    xlim=(-0.7, 0.5),
    xlabel=r"$\Delta$Charge$_{kinase}$",
    outpath="./broken_chain_ecdf_bootstrap_prot_change_per_kinaseconf.pdf",
)

# %%
