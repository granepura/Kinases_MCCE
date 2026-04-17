#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns

#%%

# --- Global style ---
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["axes.labelsize"] = 18
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.linewidth"] = 1.2
sns.set_style("white")

#%%
# Load the broken-chain csv.
# The header in this file uses a stray non-ASCII byte (0xC6) where an
# underscore should be, so read as latin-1 and repair the column names.
df = pd.read_csv(
    "./broken-chain-net-charge-changes-dfg-label-table.csv",
    encoding="latin-1",
)
df.columns = [c.replace("\xc6", "_") for c in df.columns]

# Normalize conformation labels to just "DFG-in" / "DFG-out"
df["Conformation"] = df["Kinase conformation"].str.extract(r"(DFG-in|DFG-out)")

#%%
palette = {
    "DFG-in":  "tab:blue",
    "DFG-out": "tab:orange",
}
hue_order = ["DFG-in", "DFG-out"]


def _ecdf_xy(values, xlim):
    """Return (x, y) for a step-ECDF that spans the full xlim.

    The line enters the plot at y=0 on the left edge, steps up at each
    data point, and is extended horizontally at y=1 to the right edge so
    the curve doesn't just stop at the last observation.
    """
    v = np.sort(np.asarray(values, dtype=float))
    v = v[~np.isnan(v)]
    n = v.size
    if n == 0:
        return np.array([xlim[0], xlim[1]]), np.array([0.0, 0.0])
    y = np.arange(1, n + 1) / n
    x = np.concatenate(([xlim[0]], v, [xlim[1]]))
    y_full = np.concatenate(([0.0], y, [1.0]))
    return x, y_full


def _plot_ecdf(ax, data, value_col, xlim):
    for conf in hue_order:
        sub = data.loc[data["Conformation"] == conf, value_col]
        x, y = _ecdf_xy(sub.values, xlim)
        ax.step(
            x, y,
            where="post",
            color=palette[conf],
            linewidth=2.5,
            label=conf,
            solid_capstyle="round",
        )


def _style_ecdf(ax, xlabel, xlim):
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel("Empirical CDF", fontweight="bold")
    ax.axvline(0, color="grey", linestyle=":", linewidth=1)
    sns.despine(ax=ax, top=True, right=True)
    ax.tick_params(
        bottom=True,
        left=True,
        top=False,
        right=False,
        length=5,
        width=1,
    )
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight("bold")

    leg = ax.legend(
        frameon=False,
        fontsize=16,
        loc="best",
    )
    for txt in leg.get_texts():
        txt.set_fontweight("bold")


#%%
# ---- Plot 1: ECDF of overall net charge change per row
# (ligand change + protein change), split by DFG-in / DFG-out.
df["Net_crg"] = df["Inhib_crg"] + df["Prot_crg"]

xlim_net = (-1.2, 0.6)
fig, ax = plt.subplots(figsize=(8, 5))
_plot_ecdf(ax, df, "Net_crg", xlim_net)
_style_ecdf(ax, r"$\Delta$ net charge (inhibitor + kinase)", xlim_net)
plt.tight_layout()
plt.savefig("./broken_chain_ecdf_net_change_per_kinaseconf.pdf", dpi=300)
plt.show()

#%%
# ---- Plot 2: ECDF of ligand (inhibitor) net charge change,
# DFG-in and DFG-out on the same plot.
xlim_inh = (-0.9, 0.5)
fig, ax = plt.subplots(figsize=(8, 5))
_plot_ecdf(ax, df, "Inhib_crg", xlim_inh)
_style_ecdf(ax, r"$\Delta$ inhibitor charge", xlim_inh)
plt.tight_layout()
plt.savefig("./broken_chain_ecdf_inhibitor_change_per_kinaseconf.pdf", dpi=300)
plt.show()

#%%
# ---- Plot 3: ECDF of protein (kinase) net charge change,
# DFG-in and DFG-out on the same plot.
xlim_prot = (-0.7, 0.5)
fig, ax = plt.subplots(figsize=(8, 5))
_plot_ecdf(ax, df, "Prot_crg", xlim_prot)
_style_ecdf(ax, r"$\Delta$ kinase charge", xlim_prot)
plt.tight_layout()
plt.savefig("./broken_chain_ecdf_prot_change_per_kinaseconf.pdf", dpi=300)
plt.show()

# %%
