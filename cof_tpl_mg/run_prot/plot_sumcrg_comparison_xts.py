#!/usr/bin/env python
"""
Created on Aug 26 23:15:00 2025

@author: Gehan Ranepura
Updated: Oct 6, 2025 - Case-sensitive paths, verbose output, console-mirrored logging
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# === Configuration ===
dir1 = "/home/granepura/5-kinases/cof_tpl_mg/run_prot"
dir2 = "/home/granepura/5-kinases/cof_tpl_mg/run_kin"
lst_file = "/home/granepura/5-kinases/cof_tpl_mg/pdb_inhibitor.lst"
plot_dir = "plots_kin_vs_prot_xts"
log_file = "plot_sumcrg_comparison_xts.log"



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
        # Strip color codes from log file
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

# === Accumulators for the final combined plot ===
grouped_all = defaultdict(lambda: ([], []))  # aa -> (X_list, Y_list)
all_x_all = []  # X across all PDBs (ALL points)
all_y_all = []  # Y across all PDBs (ALL points)

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

    # Check file existence
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
    all_x = []  # per-PDB regression X (ALL points)
    all_y = []  # per-PDB regression Y (ALL points)

    inhibitor_count = 0
    aa_counts = defaultdict(int)

    for key in shared_keys:
        aa = key[:3]
        xval = charges1[key]
        yval = charges2[key]

        # group for plotting
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

        # ALWAYS include in regression arrays (per-PDB and global)
        all_x.append(xval)
        all_y.append(yval)
        all_x_all.append(xval)
        all_y_all.append(yval)

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
        plt.scatter(x_vals, y_vals, label=label, color=color, alpha=0.8, edgecolor='k', s=48)

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
        plt.plot(x_line, y_line, 'k--', linewidth=1.5, label="Linear fit")

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
    plt.xlabel(f"Boltzmann Weighted {x_label.replace('run_', '').title()} Charge", fontweight="bold")
    plt.ylabel(f"Boltzmann Weighted {y_label.replace('run_', '').title()} Charge", fontweight="bold")
    plt.title(f"MCCE Kinase Charge ({y_label.replace('run_', '').title()} vs {x_label.replace('run_', '').title()}) \n {inhibitor} (PDB: {pdb})", fontweight="bold")
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)

    # Legend pinned lower-right
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

    # scatter by residue groups (global)
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
    plt.plot(x_line, y_line, 'k--', linewidth=1.5, label="Linear fit")

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
    plt.xlabel(f"Boltzmann Weighted {x_label.replace('run_', '').title()} Charge", fontweight="bold")
    plt.ylabel(f"Boltzmann Weighted {y_label.replace('run_', '').title()} Charge", fontweight="bold")
    plt.title(f"MCCE Kinase Charge ({y_label.replace('run_', '').title()} vs {x_label.replace('run_', '').title()})", fontweight="bold")
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
print_log("="*80)
print_log(f"{GREEN}✅ All done! Check '{plot_dir}' for plots and '{log_file}' for log.{RESET}")

# === Close log file ===
log_handle.close()
