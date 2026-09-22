#!/usr/bin/env python
"""
Created on Aug 26 23:15:00 2025

@author: Gehan Ranepura
Updated: Oct 6, 2025 - Case-sensitive paths, verbose output, console-mirrored logging
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# === Configuration (across trials) ===
# Aggregates the per-trial Fig4A point CSVs written by
# TrialNN/3-plot_..._Fig4A.py.  Each point is one residue: x = apo charge,
# y = holo charge, averaged over the trials, with error bars.
#
# The holo-apo difference is computed WITHIN each trial and then averaged --
# apo and holo in one trial share a seed and the same conformers, so the paired
# difference cancels noise that combining two independent SEMs would not.
import argparse
import csv
import glob
import statistics as _st
from collections import defaultdict

_ap = argparse.ArgumentParser(
    description="Fig4A across trials: mean per-residue charge with error bars "
                "from the replicate trials.")
_ap.add_argument("--root", help="directory holding Trial*/ (default: this script's)")
_ap.add_argument("--glob", default="Trial*", help="which trials to include (default: Trial*)")
_ap.add_argument("--err", choices=["sem", "sd", "range"], default="sem",
                 help="error bars: sem = SD/sqrt(n) on the mean (default), "
                      "sd = spread of trials, range = min..max")
_ap.add_argument("--title", action="store_true",
                 help="draw the figure title on the PNGs (default: off)")
_args = _ap.parse_args()

ROOT = os.path.abspath(_args.root) if _args.root \
       else os.path.dirname(os.path.abspath(__file__))
LABEL = "Trials"

plot_dir = os.path.join(ROOT, "plots_Trials_Fig4A_holo_vs_apo_xts")
log_file = os.path.join(plot_dir, f"{LABEL}_Fig4A_holo_vs_apo.log")
lst_file = os.path.join(ROOT, "pdb_inhibitor.lst")
x_label = "run_apo"
y_label = "run_holo"

# === Outlier criteria (two tiers) ===
# Strong outlier:  |holo - apo| >= STRONG_THRESHOLD   (half-ionization shift, ~strong
#                  protonation-state response to inhibitor binding)
# Affected residue: AFFECTED_THRESHOLD <= |holo - apo| < STRONG_THRESHOLD
#                  (modest but real shift, above MCCE Boltzmann/numerical noise ~0.05–0.1 e)
STRONG_THRESHOLD = 0.5
AFFECTED_THRESHOLD = 0.2
# Small tolerance so that residues sitting exactly on a threshold (e.g. Δq=0.20)
# are not excluded by IEEE-float rounding (0.59 - 0.39 == 0.19999999999999996).
THRESHOLD_EPS = 1e-9
outlier_table_file = f"{LABEL}_Fig4A_outliers_holo_vs_apo.tsv"

os.makedirs(plot_dir, exist_ok=True)

# Axis labels and names from directory paths
# x_label / y_label are set in the configuration block above.

# === Terminal Colors ===
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
RED = '\033[0;31m'
CYAN = '\033[0;36m'
BLUE = '\033[0;34m'
RESET = '\033[0m'

# === Helper function to print and log simultaneously ===
log_handle = None

def print_log(message=""):
    """Print to console and write to log file"""
    print(message)
    if log_handle:
        # Strip color codes from log file
        import re
        clean_msg = re.sub(r'\033\[[0-9;]+m', '', message)
        log_handle.write(clean_msg + "\n")
        log_handle.flush()

# === Start Logging ===
os.makedirs(plot_dir, exist_ok=True)
log_handle = open(log_file, "w")
print_log("="*80)
print_log(f"Fig4A ACROSS TRIALS  |  mean +/- {_args.err.upper()}")
print_log("="*80)
print_log(f"  root    : {ROOT}")
print_log(f"  x       : {x_label} (apo)")
print_log(f"  y       : {y_label} (holo)")
print_log(f"  outputs : {plot_dir}")

# === Load Inhibitor Code Mapping (COLUMN-BASED, KEEP ORIGINAL CASE) ===
inhibitor_codes = {}
with open(lst_file) as f:
    header = f.readline()
    for line in f:
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if len(cols) >= 3:
            pdb_entry = cols[0]        # Keep ORIGINAL case
            inhibitor = cols[1]
            icode = cols[2]
            inhibitor_codes[pdb_entry] = (inhibitor, icode)

print_log("="*80)
print_log(f"LOADED {len(inhibitor_codes)} PDB ENTRIES from {lst_file}")
print_log("="*80)

# === Helpers: parsing and stats ===
def parse_sum_crg(filepath):
    """Parse xts_sum_crg.out into {residue: charge}, skipping headers and summary lines."""
    charges = {}
    with open(filepath) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            # skip headers and summary lines
            if line.lower().startswith("ph"):
                continue
            if line.startswith("-"):   # dashed separator line
                continue
            if "Net_Charge" in line or "Protons" in line or "Electrons" in line:
                continue
            # expect "RESNAME   value"
            parts = line.split()
            if len(parts) != 2:
                continue
            residue, val = parts
            try:
                charges[residue] = float(val)
            except ValueError:
                continue
    return charges

def rankdata_average(x):
    """Return ranks with average method for ties (like scipy.stats.rankdata(method='average'))."""
    x = np.asarray(x)
    sorter = np.argsort(x, kind="mergesort")  # stable
    inv = np.empty_like(sorter)
    inv[sorter] = np.arange(len(x))
    x_sorted = x[sorter]
    diffs = np.diff(x_sorted)
    idx = np.concatenate(([0], np.nonzero(diffs)[0] + 1, [len(x)]))
    ranks = np.empty(len(x), dtype=float)
    for i in range(len(idx) - 1):
        start, end = idx[i], idx[i + 1]
        avg = 0.5 * (start + 1 + end)
        ranks[start:end] = avg
    return ranks[inv]

def spearman_rho(x, y):
    """Compute Spearman rank correlation (ρ) without SciPy."""
    x = np.asarray(x); y = np.asarray(y)
    if x.size < 2 or y.size < 2:
        return np.nan
    rx = rankdata_average(x); ry = rankdata_average(y)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return np.corrcoef(rx, ry)[0, 1]

# === Custom Color Map ===
aa_color_map = {
    # Acids
    'ASP': '#FF00FF',     # magenta
    'GLU': '#FF0000',     # red
    # Bases
    'ARG': '#4169E1',     # royal blue
    'LYS': '#00FFFF',     # cyan
    'HIS': '#800080',     # purple
    # Others (avoid red/blue/green families)
    'CYS': '#000000',     # black
    'TYR': '#FFA500',     # orange
    'SER': '#DAA520',     # goldenrod
    'THR': '#A0522D',     # sienna
    'ASN': '#9400D3',     # dark violet
    'GLN': '#9370DB',     # medium purple
    'ALA': '#808080',     # gray
    'VAL': '#FF8C00',     # dark orange
    'LEU': '#D2691E',     # chocolate
    'ILE': '#CD853F',     # peru
    'MET': '#B8860B',     # dark goldenrod
    'PRO': '#696969',     # dimgray
    'TRP': '#4B0082',     # indigo
    'PHE': '#8B4513',     # saddle brown
    'GLY': '#708090',     # slate gray
}
def get_color(resname):
    return aa_color_map.get(resname, '#000000')  # default black

# === Residue key parsing (e.g., "ASP-A0761_" -> ("ASP","A","761","ASP761")) ===
_res_key_re = re.compile(r'^([A-Za-z]{3})[+\-0]?([A-Za-z0-9])?(\d+)')

def parse_res_key(key):
    """Return (resname, chain, resnum_str, short_label) from a sum_crg residue key."""
    m = _res_key_re.match(key)
    if not m:
        return key[:3], "", "", key.rstrip("_")
    resname = m.group(1).upper()
    chain = m.group(2) or ""
    resnum = str(int(m.group(3)))  # strip leading zeros
    return resname, chain, resnum, f"{resname}{resnum}"

# === Accumulators for the final combined plot ===
grouped_all = defaultdict(lambda: ([], []))  # aa -> (X_list, Y_list)
grouped_all_err = defaultdict(lambda: ([], []))  # aa -> (xerr_list, yerr_list)
all_x_all = []  # X across all PDBs (ALL points)
all_y_all = []  # Y across all PDBs (ALL points)

# Global outlier list of dicts: {pdb, inhibitor, key, resname, chain, resnum, label, apo, holo, delta}
outliers_all = []

# EVERY point, not just outliers -- this is the per-trial CSV that the
# cross-trial averaging reads.  apo and holo are paired within one trial, so a
# per-trial delta is the quantity to average across trials.
points_all = []

# === Statistics tracking ===
stats = {
    'total': len(inhibitor_codes),
    'files_found': 0,
    'files_missing': 0,
    'plots_created': 0,
    'no_residues': 0,
    'total_residues': 0
}

# === Main Plotting Loop (per PDB) ===
print_log()
print_log("="*80)
print_log("PROCESSING PDB DIRECTORIES")
print_log("="*80)
print_log()

# === Load every trial's per-point CSV and average ===
csv_glob = os.path.join(ROOT, _args.glob,
                        "plots_Fig4A_holo_vs_apo_xts", "*_holo_vs_apo_points.csv")
trial_files = sorted(glob.glob(csv_glob))
if not trial_files:
    sys.exit(f"ERROR: no per-trial CSVs matched {csv_glob}\n"
             f"       Run TrialNN/3-plot_sumcrg_comparison_xts_Fig4A.py first.")

raw = defaultdict(dict)        # (pdb, key) -> {trial: (apo, holo)}
meta = {}                      # (pdb, key) -> resname/chain/resnum/label/is_inhibitor
for f in trial_files:
    n = 0
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh):
            k = (r["pdb"], r["key"])
            raw[k][r["trial"]] = (float(r["apo"]), float(r["holo"]))
            meta[k] = r
            n += 1
    print_log(f"  {os.path.relpath(f, ROOT)}  ->  {n} rows")

TRIALS = sorted({tr for v in raw.values() for tr in v})
print_log(f"\n  trials found: {', '.join(TRIALS)}  (n = {len(TRIALS)})")
if len(TRIALS) < 2:
    print_log("  WARNING: fewer than 2 trials -- error bars will all be zero.")

def _stats(vals):
    m = _st.fmean(vals)
    sd = _st.stdev(vals) if len(vals) > 1 else 0.0
    return m, sd, (sd / len(vals) ** 0.5 if len(vals) > 1 else 0.0), min(vals), max(vals)

AGG = {}                       # (pdb, key) -> statistics
partial = []
for k, by_trial in raw.items():
    if len(by_trial) != len(TRIALS):
        partial.append(k)
        continue
    apos = [by_trial[tr][0] for tr in TRIALS]
    holos = [by_trial[tr][1] for tr in TRIALS]
    deltas = [h - a for a, h in zip(apos, holos)]      # paired within trial
    am, asd, asem, alo, ahi = _stats(apos)
    hm, hsd, hsem, hlo, hhi = _stats(holos)
    dm, dsd, dsem, dlo, dhi = _stats(deltas)
    AGG[k] = dict(n=len(TRIALS), apo=apos, holo=holos, delta=deltas,
                  am=am, asd=asd, asem=asem, alo=alo, ahi=ahi,
                  hm=hm, hsd=hsd, hsem=hsem, hlo=hlo, hhi=hhi,
                  dm=dm, dsd=dsd, dsem=dsem, dlo=dlo, dhi=dhi)
if partial:
    print_log(f"  {len(partial)} residue(s) absent from some trial -- excluded")

def _err(pdb, key, axis):
    a = AGG[(pdb, key)]
    if _args.err == "sd":
        return a[f"{axis}sd"]
    if _args.err == "range":
        return (a[f"{axis}hi"] - a[f"{axis}lo"]) / 2.0
    return a[f"{axis}sem"]

# residues per PDB, so the loop below can work off the averages
BY_PDB = defaultdict(list)
for (pdb, key) in AGG:
    BY_PDB[pdb].append(key)

for idx, (pdb, (inhibitor, icode)) in enumerate(sorted(inhibitor_codes.items()), 1):
    print_log(f"[{idx:2d}/{len(inhibitor_codes)}] {pdb:<25s} │ {inhibitor:<15s} │ Code: {icode}")

    if pdb not in BY_PDB:
        stats['files_missing'] += 1
        print_log(f"      └─ {YELLOW}⚠️  SKIPPED (not in the per-trial CSVs){RESET}")
        print_log()
        continue

    stats['files_found'] += 1

    # The averages stand in for the per-run charges the per-trial script read.
    charges1 = {k: AGG[(pdb, k)]["am"] for k in BY_PDB[pdb]}   # apo  mean
    charges2 = {k: AGG[(pdb, k)]["hm"] for k in BY_PDB[pdb]}   # holo mean

    shared_keys = sorted(charges1.keys())
    
    if not shared_keys:
        stats['no_residues'] += 1
        print_log(f"      ├─ {YELLOW}⚠️  No shared residues{RESET}")
        print_log(f"      └─ {YELLOW}⚠️  SKIPPED{RESET}")
        print_log()
        continue

    print_log(f"      ├─ ✅ Found {len(shared_keys)} shared residues")

    grouped_x = defaultdict(list)
    grouped_y = defaultdict(list)
    grouped_keys = defaultdict(list)      # same order as grouped_x/grouped_y
    all_x = []  # per-PDB regression X (ALL points)
    all_y = []  # per-PDB regression Y (ALL points)
    outliers_pdb = []  # per-PDB outliers for annotation

    inhibitor_count = 0
    aa_counts = defaultdict(int)

    for key in shared_keys:
        aa = key[:3]
        xval = charges1[key]
        yval = charges2[key]

        # group for plotting
        _g = 'INHIBITOR' if aa == icode else aa
        grouped_x[_g].append(xval)
        grouped_y[_g].append(yval)
        grouped_keys[_g].append(key)
        Xg, Yg = grouped_all[_g]
        Xg.append(xval); Yg.append(yval)
        _ex, _ey = grouped_all_err[_g]
        _ex.append(_err(pdb, key, "a")); _ey.append(_err(pdb, key, "h"))
        if aa == icode:
            inhibitor_count += 1
        else:
            aa_counts[aa] += 1

        # ALWAYS include in regression arrays (per-PDB and global)
        all_x.append(xval)
        all_y.append(yval)
        all_x_all.append(xval)
        all_y_all.append(yval)

        # Two-tier outlier classification (with FP tolerance on the boundary)
        delta = AGG[(pdb, key)]['dm']     # mean of the per-trial paired deltas
        adelta = abs(delta)

        _rn, _ch, _num, _lbl = parse_res_key(key)
        points_all.append({
            'trial': LABEL, 'pdb': pdb, 'inhibitor': inhibitor, 'key': key,
            'resname': _rn, 'chain': _ch, 'resnum': _num, 'label': _lbl,
            'is_inhibitor': int(aa == icode),
            'apo': xval, 'holo': yval, 'delta': delta,
        })
        tier = None
        if adelta >= STRONG_THRESHOLD - THRESHOLD_EPS:
            tier = 'strong'
        elif adelta >= AFFECTED_THRESHOLD - THRESHOLD_EPS:
            tier = 'affected'
        if tier is not None:
            resname, chain, resnum, label = parse_res_key(key)
            _a = AGG[(pdb, key)]
            rec = {
                'pdb': pdb,
                'inhibitor': inhibitor,
                'key': key,
                'resname': resname,
                'chain': chain,
                'resnum': resnum,
                'label': label,
                'apo': xval,
                'holo': yval,
                'delta': delta,
                'delta_sem': _a['dsem'],
                'delta_min': _a['dlo'],
                'delta_max': _a['dhi'],
                'tier': tier,
            }
            outliers_pdb.append(rec)
            outliers_all.append(rec)

    stats['total_residues'] += len(shared_keys)

    # Show breakdown
    if inhibitor_count > 0:
        print_log(f"      │    • Inhibitor ({icode}): {inhibitor_count} residue(s)")
    top_aas = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)[:3]
    if top_aas:
        aa_summary = ", ".join([f"{aa}:{cnt}" for aa, cnt in top_aas])
        print_log(f"      │    • Top residues: {aa_summary}")

    # === Plotting (per PDB) ===
    plt.figure(figsize=(5, 5))

    for aa, x_vals in grouped_x.items():
        y_vals = grouped_y[aa]
        if aa == 'INHIBITOR':
            color = 'green'
            label = f"Inhibitor ({icode})"
        else:
            color = get_color(aa)
            label = aa
        plt.errorbar(x_vals, y_vals,
                     xerr=[_err(pdb, k, "a") for k in grouped_keys[label]],
                     yerr=[_err(pdb, k, "h") for k in grouped_keys[label]],
                     fmt="none", ecolor="0.4", elinewidth=0.9, capsize=2,
                     capthick=0.9, zorder=1)
        plt.scatter(x_vals, y_vals, label=label, color=color, alpha=0.8, edgecolor='k', s=48)

    # Outlier overlay — two tiers:
    #   strong   (|Δq| ≥ STRONG_THRESHOLD)   → red circle + label
    #   affected (|Δq| ≥ AFFECTED_THRESHOLD) → orange circle + label
    strong_pdb = [r for r in outliers_pdb if r['tier'] == 'strong']
    affected_pdb = [r for r in outliers_pdb if r['tier'] == 'affected']
    if affected_pdb:
        plt.scatter([r['apo'] for r in affected_pdb],
                    [r['holo'] for r in affected_pdb],
                    s=100, facecolors='none', edgecolors='darkorange',
                    linewidths=1.4,
                    label=f"Affected ({AFFECTED_THRESHOLD}≤|Δq|<{STRONG_THRESHOLD})",
                    zorder=5)
        for r in affected_pdb:
            plt.annotate(r['label'], (r['apo'], r['holo']),
                         xytext=(6, 6), textcoords='offset points',
                         fontsize=6, color='darkorange', zorder=6)
    if strong_pdb:
        plt.scatter([r['apo'] for r in strong_pdb],
                    [r['holo'] for r in strong_pdb],
                    s=130, facecolors='none', edgecolors='red',
                    linewidths=1.8,
                    label=f"Outlier (|Δq|≥{STRONG_THRESHOLD})",
                    zorder=6)
        for r in strong_pdb:
            plt.annotate(r['label'], (r['apo'], r['holo']),
                         xytext=(6, 6), textcoords='offset points',
                         fontsize=7, color='red', zorder=7)
    if outliers_pdb:
        print_log(f"      │    • Outliers: strong={len(strong_pdb)} "
                  f"(|Δq|≥{STRONG_THRESHOLD}), "
                  f"affected={len(affected_pdb)} "
                  f"({AFFECTED_THRESHOLD}≤|Δq|<{STRONG_THRESHOLD})")

    # Identity line y = x (solid red)
    if all_x and all_y:
        xy_min = min(min(all_x), min(all_y))
        xy_max = max(max(all_x), max(all_y))
    else:
        xy_min, xy_max = -1.0, 1.0
    plt.plot([xy_min, xy_max], [xy_min, xy_max], color='red', linestyle='-', linewidth=1.5, label="y = x")

    # Linear regression (dashed black) on ALL points + Stats
    if all_x and all_y:
        x_arr = np.array(all_x)
        y_arr = np.array(all_y)

        m, b = np.polyfit(x_arr, y_arr, 1)
        y_fit = m * x_arr + b

        ss_res = np.sum((y_arr - y_fit) ** 2)
        ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 1.0

        r = np.nan if (np.std(x_arr) == 0 or np.std(y_arr) == 0) else np.corrcoef(x_arr, y_arr)[0, 1]
        rho = spearman_rho(x_arr, y_arr)

        diff = y_arr - x_arr
        rmse = float(np.sqrt(np.mean(diff**2)))
        mae = float(np.mean(np.abs(diff)))
        N = len(x_arr)

        x_line = np.linspace(xy_min, xy_max, 100)
        y_line = m * x_line + b
        #plt.plot(x_line, y_line, 'k--', linewidth=1.5, label="Linear fit")

        stats_text = (
            f"$y = {m:.2f}x + {b:.2f}$, $R^2 = {r2:.2f}$\n"
            f"$r = {r:.2f}$, $\\rho = {rho:.2f}$\n"
            f"RMSE = {rmse:.2f}, MAE = {mae:.2f}\n"
            f"N = {N}"
        )
        # Stats box pinned upper-left
        plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
                 fontsize=9, verticalalignment='top', horizontalalignment='left',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.70))

    # === Titles & Labels ===
    plt.xlabel(f"{x_label.replace('run_', '').title()} Charge", fontweight="bold")
    plt.ylabel(f"{y_label.replace('run_', '').title()} Charge", fontweight="bold")
    if _args.title:
        plt.title(f"MCCE Kinase Charge ({y_label.replace('run_', '').title()} vs {x_label.replace('run_', '').title()}) -- {LABEL}\n {inhibitor} (PDB: {pdb})", fontweight="bold")
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)

    # Legend pinned lower-right
    plt.legend(loc='lower right', fontsize='small', markerscale=1, framealpha=0.9)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    outpath = os.path.join(plot_dir, f"{LABEL}_Fig4A_{pdb}_{inhibitor}.png")
    plt.savefig(outpath, dpi=200)
    plt.close()

    stats['plots_created'] += 1

    print_log(f"      └─ {GREEN}✅ Plot saved:{RESET} {outpath}")
    print_log()

# === FINAL COMBINED PLOT (across all PDBs; ALL points) ===
print_log("="*80)
print_log("CREATING COMBINED PLOT")
print_log("="*80)
print_log()

if all_x_all and all_y_all:
    plt.figure(figsize=(6, 6))

    # scatter by residue groups (global)
    for aa, (Xg, Yg) in grouped_all.items():
        if aa == 'INHIBITOR':
            color = 'green'
            label = "Inhibitor (all)"
        else:
            color = get_color(aa)
            label = aa
        plt.errorbar(Xg, Yg,
                     xerr=grouped_all_err[label][0], yerr=grouped_all_err[label][1],
                     fmt="none", ecolor="0.55", elinewidth=0.7, capsize=1.5,
                     capthick=0.7, zorder=1)
        plt.scatter(Xg, Yg, label=label, color=color, alpha=0.6, edgecolor='k', s=30)

    xy_min = min(min(all_x_all), min(all_y_all))
    xy_max = max(max(all_x_all), max(all_y_all))
    plt.plot([xy_min, xy_max], [xy_min, xy_max], color='red', linestyle='-', linewidth=1.5, label="y = x")

    # Two-tier outlier overlay across all PDBs.
    # Affected (orange) is drawn first and unlabeled to avoid clutter; strong (red)
    # is drawn on top with PDB:residue labels.
    strong_all = [r for r in outliers_all if r['tier'] == 'strong']
    affected_all = [r for r in outliers_all if r['tier'] == 'affected']
    if affected_all:
        plt.scatter([r['apo'] for r in affected_all],
                    [r['holo'] for r in affected_all],
                    s=90, facecolors='none', edgecolors='darkorange',
                    linewidths=1.3,
                    label=f"Affected ({AFFECTED_THRESHOLD}≤|Δq|<{STRONG_THRESHOLD})",
                    zorder=5)
    if strong_all:
        plt.scatter([r['apo'] for r in strong_all],
                    [r['holo'] for r in strong_all],
                    s=120, facecolors='none', edgecolors='red',
                    linewidths=1.6,
                    label=f"Outlier (|Δq|≥{STRONG_THRESHOLD})",
                    zorder=6)
        for r in strong_all:
            plt.annotate(f"{r['pdb']}:{r['label']}", (r['apo'], r['holo']),
                         xytext=(5, 5), textcoords='offset points',
                         fontsize=6, color='red', zorder=7)

    # Global regression on ALL points
    x_arr = np.array(all_x_all)
    y_arr = np.array(all_y_all)
    m, b = np.polyfit(x_arr, y_arr, 1)
    y_fit = m * x_arr + b
    ss_res = np.sum((y_arr - y_fit) ** 2)
    ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 1.0
    r = np.nan if (np.std(x_arr) == 0 or np.std(y_arr) == 0) else np.corrcoef(x_arr, y_arr)[0, 1]
    rho = spearman_rho(x_arr, y_arr)
    diff = y_arr - x_arr
    rmse = float(np.sqrt(np.mean(diff**2)))
    mae = float(np.mean(np.abs(diff)))
    N = len(x_arr)

    x_line = np.linspace(xy_min, xy_max, 100)
    y_line = m * x_line + b
    #plt.plot(x_line, y_line, 'k--', linewidth=1.5, label="Linear fit")

    stats_text = (
        f"$y = {m:.2f}x + {b:.2f}$, $R^2 = {r2:.2f}$\n"
        f"$r = {r:.2f}$, $\\rho = {rho:.2f}$\n"
        f"RMSE = {rmse:.2f}, MAE = {mae:.2f}\n"
        f"N = {N}"
    )
    # Stats box upper-left
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
             fontsize=9, verticalalignment='top', horizontalalignment='left',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.70))

    # Labels & styling
    #plt.xlabel(f"Boltzmann Weighted {x_label.replace('run_', '').title()} Charge", fontweight="bold")
    #plt.ylabel(f"Boltzmann Weighted {y_label.replace('run_', '').title()} Charge", fontweight="bold")
    #plt.title(f"MCCE Kinase Protein Charge ({y_label.replace('run_', '').title()} vs {x_label.replace('run_', '').title()})", fontweight="bold")
    plt.xlabel(f"Apo-Protein Charge", fontweight="bold", fontsize=16)
    plt.ylabel(f"Holo-Protein Charge", fontweight="bold", fontsize=16)
    if _args.title:
        plt.title(f"MCCE Kinase Protein Charge (Holo vs Apo) -- {LABEL} Fig4A", fontweight="bold")
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.legend(loc='lower right', fontsize='small', markerscale=1, framealpha=0.9, ncol=2)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    all_out = os.path.join(plot_dir, f"{LABEL}_Fig4A_ALL_{y_label}_vs_{x_label}.png")
    plt.savefig(all_out, dpi=220)
    plt.close()

    print_log(f"{GREEN}✅ Combined plot saved:{RESET} {all_out}")
    print_log()
else:
    print_log(f"{YELLOW}⚠️  No data accumulated for combined plot — nothing to save.{RESET}")
    print_log()

# === OUTLIER TABLE (two tiers) ===
print_log("="*80)
print_log(f"OUTLIER TABLE  (strong: |Δq| ≥ {STRONG_THRESHOLD} e;  "
          f"affected: {AFFECTED_THRESHOLD} ≤ |Δq| < {STRONG_THRESHOLD} e)")
print_log("="*80)

n_strong = sum(1 for r in outliers_all if r['tier'] == 'strong')
n_affected = sum(1 for r in outliers_all if r['tier'] == 'affected')

outlier_path = os.path.join(plot_dir, outlier_table_file)
with open(outlier_path, "w") as of:
    header_cols = ["Tier", "PDB", "Inhibitor", "ResName", "Chain", "ResNum",
                   "Residue", "ApoCharge", "HoloCharge", "DeltaCharge",
                   "DeltaSEM", "DeltaMin", "DeltaMax", "Reproducible"]
    of.write("\t".join(header_cols) + "\n")
    # strong first, then affected; each sorted by |Δq| descending
    ordered = (
        sorted([r for r in outliers_all if r['tier'] == 'strong'],
               key=lambda d: abs(d['delta']), reverse=True)
        + sorted([r for r in outliers_all if r['tier'] == 'affected'],
                 key=lambda d: abs(d['delta']), reverse=True)
    )
    for r in ordered:
        of.write("\t".join([
            r['tier'], r['pdb'], r['inhibitor'], r['resname'], r['chain'], r['resnum'],
            r['key'].rstrip('_'),
            f"{r['apo']:.3f}", f"{r['holo']:.3f}", f"{r['delta']:+.3f}",
            f"{r['delta_sem']:.3f}", f"{r['delta_min']:+.3f}", f"{r['delta_max']:+.3f}",
            # an outlier only counts if every trial clears the tier threshold
            "yes" if min(abs(r['delta_min']), abs(r['delta_max'])) >=
                     (STRONG_THRESHOLD if r['tier'] == 'strong'
                      else AFFECTED_THRESHOLD) - THRESHOLD_EPS else "NO",
        ]) + "\n")

print_log(f"  Strong outliers (|Δq|≥{STRONG_THRESHOLD}):      {n_strong}")
print_log(f"  Affected ({AFFECTED_THRESHOLD}≤|Δq|<{STRONG_THRESHOLD}):         {n_affected}")
print_log(f"  Outlier table written to:    {outlier_path}")

# === Per-point CSV (input for the cross-trial averaging) ===
points_csv = os.path.join(plot_dir, f"{LABEL}_Fig4A_holo_vs_apo_points.csv")
with open(points_csv, "w", newline="") as fh:
    _cols = ["pdb", "inhibitor", "key", "resname", "chain", "resnum", "label",
             "is_inhibitor", "n_trials",
             "apo_mean", "apo_sd", "apo_sem", "apo_min", "apo_max",
             "holo_mean", "holo_sd", "holo_sem", "holo_min", "holo_max",
             "delta_mean", "delta_sd", "delta_sem", "delta_min", "delta_max",
             "per_trial_apo", "per_trial_holo", "per_trial_delta"]
    w = csv.DictWriter(fh, fieldnames=_cols)
    w.writeheader()
    for (pdb_u, key) in sorted(AGG):
        a = AGG[(pdb_u, key)]
        m = meta[(pdb_u, key)]
        w.writerow({
            "pdb": pdb_u, "inhibitor": m["inhibitor"], "key": key,
            "resname": m["resname"], "chain": m["chain"], "resnum": m["resnum"],
            "label": m["label"], "is_inhibitor": m["is_inhibitor"], "n_trials": a["n"],
            "apo_mean": f"{a['am']:.6f}", "apo_sd": f"{a['asd']:.6f}",
            "apo_sem": f"{a['asem']:.6f}", "apo_min": f"{a['alo']:.6f}",
            "apo_max": f"{a['ahi']:.6f}",
            "holo_mean": f"{a['hm']:.6f}", "holo_sd": f"{a['hsd']:.6f}",
            "holo_sem": f"{a['hsem']:.6f}", "holo_min": f"{a['hlo']:.6f}",
            "holo_max": f"{a['hhi']:.6f}",
            "delta_mean": f"{a['dm']:.6f}", "delta_sd": f"{a['dsd']:.6f}",
            "delta_sem": f"{a['dsem']:.6f}", "delta_min": f"{a['dlo']:.6f}",
            "delta_max": f"{a['dhi']:.6f}",
            "per_trial_apo": ";".join(f"{v:.2f}" for v in a["apo"]),
            "per_trial_holo": ";".join(f"{v:.2f}" for v in a["holo"]),
            "per_trial_delta": ";".join(f"{v:+.2f}" for v in a["delta"]),
        })
print_log(f"  Statistics CSV written to:   {points_csv}  ({len(AGG)} residues)")

# Mirror the full outlier table into the log, grouped by tier
if outliers_all:
    hdr = f"  {'Tier':<9s} {'PDB':<8s} {'Inhib':<12s} {'Residue':<14s} {'Apo':>8s} {'Holo':>8s} {'Δq':>8s}"
    print_log(hdr)
    print_log("  " + "-"*(len(hdr)-2))
    for tier_name in ('strong', 'affected'):
        for r in sorted((x for x in outliers_all if x['tier'] == tier_name),
                        key=lambda d: (d['pdb'], -abs(d['delta']))):
            print_log(f"  {tier_name:<9s} {r['pdb']:<8s} {r['inhibitor']:<12s} "
                      f"{r['label']:<14s} {r['apo']:>8.3f} {r['holo']:>8.3f} "
                      f"{r['delta']:>+8.3f}")
print_log()

# === FINAL SUMMARY ===
print_log("="*80)
print_log("SUMMARY")
print_log("="*80)
print_log(f"  Total PDB entries:           {stats['total']}")
print_log(f"  Files found (both dirs):     {stats['files_found']}")
print_log(f"  Files missing:               {stats['files_missing']}")
print_log(f"  No shared residues:          {stats['no_residues']}")
print_log(f"  Plots created:               {stats['plots_created']}")
print_log(f"  Total residues processed:    {stats['total_residues']}")
print_log(f"  Strong outliers:             {n_strong}  (|Δq| ≥ {STRONG_THRESHOLD})")
print_log(f"  Affected residues:           {n_affected}  ({AFFECTED_THRESHOLD} ≤ |Δq| < {STRONG_THRESHOLD})")
print_log("="*80)
print_log(f"{GREEN}✅ All done! Check '{plot_dir}' for plots and '{log_file}' for log.{RESET}")

# === Close log file ===
log_handle.close()
