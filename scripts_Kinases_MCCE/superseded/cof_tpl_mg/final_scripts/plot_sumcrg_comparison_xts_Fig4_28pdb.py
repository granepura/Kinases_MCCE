#!/usr/bin/env python
"""
Created on Aug 26 23:15:00 2025

@author: Gehan Ranepura
Updated: Oct 6, 2025 - Case-sensitive paths, verbose output, console-mirrored logging
Updated: Jun 6, 2026 - Filtered to 28 PDB entries (27 unique) from excel
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# === Configuration ===
dir1 = "/data/home/granepura/Mar11/5-Kinases/final/cof_tpl_mg/run_prot2"
dir2 = "/data/home/granepura/Mar11/5-Kinases/final/cof_tpl_mg/run_kin"
lst_file = "/data/home/granepura/Mar11/5-Kinases/final/cof_tpl_mg/pdb_inhibitor.lst"
plot_dir = "plots_kin_vs_prot2_xts_28pdb"
log_file = "plot_sumcrg_comparison_xts_28pdb.log"

# Only include the 27 unique PDBs present in the excel table
INCLUDE_PDBS = {
    "1XKK", "2HYY", "2ITO", "2ITY", "2WGJ", "2XP2", "2YFX",
    "3AOX", "3CS9", "3LXK", "3PYY", "3UE4", "3WZD", "3WZE", "3ZOS",
    "4AG8", "4AGC", "4AGD", "4AN2", "4ASD", "4G5J", "4G5P", "4LMN",
    "4MKC", "4ZAU", "5L2I", "5MO4",
}

# === Outlier criteria (two tiers) ===
STRONG_THRESHOLD = 0.5
AFFECTED_THRESHOLD = 0.2
THRESHOLD_EPS = 1e-9
outlier_table_file = "outliers_holo_vs_apo.tsv"

os.makedirs(plot_dir, exist_ok=True)

# Axis labels and names from directory paths
x_label = os.path.basename(dir1.rstrip("/"))  # e.g., run_apo
y_label = os.path.basename(dir2.rstrip("/"))  # e.g., run_holo

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
        import re
        clean_msg = re.sub(r'\033\[[0-9;]+m', '', message)
        log_handle.write(clean_msg + "\n")
        log_handle.flush()

# === Start Logging ===
log_handle = open(log_file, "w")

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
            pdb_entry = cols[0]
            inhibitor = cols[1]
            icode = cols[2]
            if pdb_entry in INCLUDE_PDBS:
                inhibitor_codes[pdb_entry] = (inhibitor, icode)

print_log("="*80)
print_log(f"LOADED {len(inhibitor_codes)} PDB ENTRIES from {lst_file} (filtered to {len(INCLUDE_PDBS)} excel PDBs)")
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
            if line.lower().startswith("ph"):
                continue
            if line.startswith("-"):
                continue
            if "Net_Charge" in line or "Protons" in line or "Electrons" in line:
                continue
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
    sorter = np.argsort(x, kind="mergesort")
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
    """Compute Spearman rank correlation without SciPy."""
    x = np.asarray(x); y = np.asarray(y)
    if x.size < 2 or y.size < 2:
        return np.nan
    rx = rankdata_average(x); ry = rankdata_average(y)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return np.corrcoef(rx, ry)[0, 1]

# === Custom Color Map ===
aa_color_map = {
    'ASP': '#FF00FF',
    'GLU': '#FF0000',
    'ARG': '#4169E1',
    'LYS': '#00FFFF',
    'HIS': '#800080',
    'CYS': '#000000',
    'TYR': '#FFA500',
    'SER': '#DAA520',
    'THR': '#A0522D',
    'ASN': '#9400D3',
    'GLN': '#9370DB',
    'ALA': '#808080',
    'VAL': '#FF8C00',
    'LEU': '#D2691E',
    'ILE': '#CD853F',
    'MET': '#B8860B',
    'PRO': '#696969',
    'TRP': '#4B0082',
    'PHE': '#8B4513',
    'GLY': '#708090',
}
def get_color(resname):
    return aa_color_map.get(resname, '#000000')

# === Residue key parsing ===
_res_key_re = re.compile(r'^([A-Za-z]{3})[+\-0]?([A-Za-z0-9])?(\d+)')

def parse_res_key(key):
    """Return (resname, chain, resnum_str, short_label) from a sum_crg residue key."""
    m = _res_key_re.match(key)
    if not m:
        return key[:3], "", "", key.rstrip("_")
    resname = m.group(1).upper()
    chain = m.group(2) or ""
    resnum = str(int(m.group(3)))
    return resname, chain, resnum, f"{resname}{resnum}"

# === Accumulators for the final combined plot ===
grouped_all = defaultdict(lambda: ([], []))
all_x_all = []
all_y_all = []

outliers_all = []

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

for idx, (pdb, (inhibitor, icode)) in enumerate(sorted(inhibitor_codes.items()), 1):
    file1 = os.path.join(dir1, pdb, "xts_sum_crg.out")
    file2 = os.path.join(dir2, pdb, "xts_sum_crg.out")

    print_log(f"[{idx:2d}/{len(inhibitor_codes)}] {pdb:<25s} │ {inhibitor:<15s} │ Code: {icode}")

    f1_exists = os.path.isfile(file1)
    f2_exists = os.path.isfile(file2)

    if not f1_exists:
        print_log(f"      ├─ ❌ Missing: {x_label}/xts_sum_crg.out")
    if not f2_exists:
        print_log(f"      ├─ ❌ Missing: {y_label}/xts_sum_crg.out")

    if not f1_exists or not f2_exists:
        stats['files_missing'] += 1
        print_log(f"      └─ {YELLOW}⚠️  SKIPPED (files missing){RESET}")
        print_log()
        continue

    stats['files_found'] += 1

    charges1 = parse_sum_crg(file1)
    charges2 = parse_sum_crg(file2)

    shared_keys = sorted(set(charges1.keys()) & set(charges2.keys()))

    if not shared_keys:
        stats['no_residues'] += 1
        print_log(f"      ├─ {YELLOW}⚠️  No shared residues{RESET}")
        print_log(f"      └─ {YELLOW}⚠️  SKIPPED{RESET}")
        print_log()
        continue

    print_log(f"      ├─ ✅ Found {len(shared_keys)} shared residues")

    grouped_x = defaultdict(list)
    grouped_y = defaultdict(list)
    all_x = []
    all_y = []
    outliers_pdb = []

    inhibitor_count = 0
    aa_counts = defaultdict(int)

    for key in shared_keys:
        aa = key[:3]
        xval = charges1[key]
        yval = charges2[key]

        if aa == icode:
            grouped_x['INHIBITOR'].append(xval)
            grouped_y['INHIBITOR'].append(yval)
            Xg, Yg = grouped_all['INHIBITOR']
            Xg.append(xval); Yg.append(yval)
            inhibitor_count += 1
        else:
            grouped_x[aa].append(xval)
            grouped_y[aa].append(yval)
            Xg, Yg = grouped_all[aa]
            Xg.append(xval); Yg.append(yval)
            aa_counts[aa] += 1

        all_x.append(xval)
        all_y.append(yval)
        all_x_all.append(xval)
        all_y_all.append(yval)

        delta = yval - xval
        adelta = abs(delta)
        tier = None
        if adelta >= STRONG_THRESHOLD - THRESHOLD_EPS:
            tier = 'strong'
        elif adelta >= AFFECTED_THRESHOLD - THRESHOLD_EPS:
            tier = 'affected'
        if tier is not None:
            resname, chain, resnum, label = parse_res_key(key)
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
                'tier': tier,
            }
            outliers_pdb.append(rec)
            outliers_all.append(rec)

    stats['total_residues'] += len(shared_keys)

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
        plt.scatter(x_vals, y_vals, label=label, color=color, alpha=0.8, edgecolor='k', s=48)

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

    if all_x and all_y:
        xy_min = min(min(all_x), min(all_y))
        xy_max = max(max(all_x), max(all_y))
    else:
        xy_min, xy_max = -1.0, 1.0
    plt.plot([xy_min, xy_max], [xy_min, xy_max], color='red', linestyle='-', linewidth=1.5, label="y = x")

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

        stats_text = (
            f"$y = {m:.2f}x + {b:.2f}$, $R^2 = {r2:.2f}$\n"
            f"$r = {r:.2f}$, $\\rho = {rho:.2f}$\n"
            f"RMSE = {rmse:.2f}, MAE = {mae:.2f}\n"
            f"N = {N}"
        )
        plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
                 fontsize=9, verticalalignment='top', horizontalalignment='left',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.70))

    # === Titles & Labels ===
    plt.xlabel(f"{x_label.replace('run_', '').title()} Charge", fontweight="bold")
    plt.ylabel(f"{y_label.replace('run_', '').title()} Charge", fontweight="bold")
    plt.title(f"MCCE Kinase Charge ({y_label.replace('run_', '').title()} vs {x_label.replace('run_', '').title()}) \n {inhibitor} (PDB: {pdb})", fontweight="bold")
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)

    plt.legend(loc='lower right', fontsize='small', markerscale=1, framealpha=0.9)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    outpath = os.path.join(plot_dir, f"{pdb}_{inhibitor}.png")
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

    for aa, (Xg, Yg) in grouped_all.items():
        if aa == 'INHIBITOR':
            color = 'green'
            label = "Inhibitor (all)"
        else:
            color = get_color(aa)
            label = aa
        plt.scatter(Xg, Yg, label=label, color=color, alpha=0.6, edgecolor='k', s=30)

    xy_min = min(min(all_x_all), min(all_y_all))
    xy_max = max(max(all_x_all), max(all_y_all))
    plt.plot([xy_min, xy_max], [xy_min, xy_max], color='red', linestyle='-', linewidth=1.5, label="y = x")

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

    stats_text = (
        f"$y = {m:.2f}x + {b:.2f}$, $R^2 = {r2:.2f}$\n"
        f"$r = {r:.2f}$, $\\rho = {rho:.2f}$\n"
        f"RMSE = {rmse:.2f}, MAE = {mae:.2f}\n"
        f"N = {N}"
    )
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
             fontsize=9, verticalalignment='top', horizontalalignment='left',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.70))

    plt.xlabel(f"Apo-Protein Charge", fontweight="bold", fontsize=16)
    plt.ylabel(f"Holo-Protein Charge", fontweight="bold", fontsize=16)
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.legend(loc='lower right', fontsize='small', markerscale=1, framealpha=0.9, ncol=2)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    all_out = os.path.join(plot_dir, f"ALL_{y_label}_vs_{x_label}.png")
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
                   "Residue", "ApoCharge", "HoloCharge", "DeltaCharge"]
    of.write("\t".join(header_cols) + "\n")
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
        ]) + "\n")

print_log(f"  Strong outliers (|Δq|≥{STRONG_THRESHOLD}):      {n_strong}")
print_log(f"  Affected ({AFFECTED_THRESHOLD}≤|Δq|<{STRONG_THRESHOLD}):         {n_affected}")
print_log(f"  Outlier table written to:    {outlier_path}")

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
