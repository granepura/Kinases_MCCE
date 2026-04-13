#!/usr/bin/env python
"""
Aggregate inhibitor scatter plot: xts_sum_crg.out dir1 vs dir2

- Reads pdb_inhibitor.lst with columns:
    PDB   Inhibitor   Inhibitor_Code
- For each PDB, extracts only residues with res[:3] == Inhibitor_Code
- Plots ALL inhibitors across PDBs onto one scatter
- Creates individual plots for each inhibitor
- Outlier detection with 15% threshold
- Symmetrical axes starting from 0 when appropriate

Author: Gehan Ranepura
Created: Aug 27, 2025
Updated: Oct 8, 2025 - Symmetrical axes, individual plots, outlier detection, correct N count
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Patch

# === Config ===
dir1 = "/data/home/granepura/Mar11/5-Kinases/final/cof_tpl_mg/run_cof2"
dir2 = "/data/home/granepura/Mar11/5-Kinases/final/cof_tpl_mg/run_kin"
lst_file = "/data/home/granepura/Mar11/5-Kinases/final/cof_tpl_mg/pdb_inhibitor.lst"
plot_dir = "plots_inhibitors_kin_vs_cof2_xts"
log_file = "plot_sumcrg_inhibitors_xts.log"

# Outlier threshold: points with |y-x| > threshold are outliers
outlier_threshold = 0.15  # 15% difference

os.makedirs(plot_dir, exist_ok=True)
fig_path = os.path.join(plot_dir, "inhibitors_all_runs.png")
csv_path = os.path.join(plot_dir, "inhibitors_all_runs.csv")

# Axis labels
x_label = os.path.basename(dir1.rstrip("/"))
y_label = os.path.basename(dir2.rstrip("/"))

# === Helper function to print and log simultaneously ===
log_handle = None

def print_log(message=""):
    """Print to console and write to log file"""
    print(message)
    if log_handle:
        log_handle.write(message + "\n")
        log_handle.flush()

# === Start Logging ===
log_handle = open(log_file, "w")

# === Load inhibitor mapping (COLUMN-BASED) ===
inhibitor_map = {}
with open(lst_file) as f:
    _ = f.readline()  # skip header
    for line in f:
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if len(cols) < 3:
            print_log(f"Warning: Expected 3 columns, got {len(cols)} in line: {line}")
            continue
        pdb_entry = cols[0]
        inhibitor_name = cols[1]
        inhibitor_code = cols[2]
        inhibitor_map[pdb_entry] = (inhibitor_name, inhibitor_code)

print_log("="*80)
print_log(f"LOADED {len(inhibitor_map)} PDB ENTRIES from {lst_file}")
print_log("="*80)

# === Helpers ===
def parse_sum_crg(filepath):
    charges = {}
    with open(filepath) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("-") or line.lower().startswith("ph"):
                continue
            if "Net_Charge" in line or "Protons" in line or "Electrons" in line:
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            res, val = parts
            try:
                charges[res] = float(val)
            except ValueError:
                pass
    return charges

def spearman_rho(x, y):
    def rankdata_average(a):
        a = np.asarray(a)
        sorter = np.argsort(a, kind="mergesort")
        inv = np.empty_like(sorter); inv[sorter] = np.arange(len(a))
        a_sorted = a[sorter]
        diffs = np.diff(a_sorted)
        idx = np.concatenate(([0], np.nonzero(diffs)[0]+1, [len(a)]))
        ranks = np.empty(len(a), float)
        for i in range(len(idx)-1):
            start, end = idx[i], idx[i+1]
            ranks[start:end] = 0.5*(start+1+end)
        return ranks[inv]
    if len(x) < 2:
        return np.nan
    rx, ry = rankdata_average(x), rankdata_average(y)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return np.corrcoef(rx, ry)[0, 1]

# === Collect data with verbose output ===
print_log()
print_log("="*80)
print_log("PROCESSING PDB DIRECTORIES")
print_log("="*80)
print_log()

points = []  # (x, y, pdb, res, inhibitor)
stats = {
    'total': len(inhibitor_map),
    'files_found': 0,
    'files_missing': 0,
    'residues_found': 0,
    'pdbs_with_data': 0,
    'pdbs_without_data': 0
}

for idx, (pdb, (inhibitor, icode)) in enumerate(sorted(inhibitor_map.items()), 1):
    print_log(f"[{idx:2d}/{len(inhibitor_map)}] {pdb:<25s} │ {inhibitor:<15s} │ Code: {icode}")

    f1 = os.path.join(dir1, pdb, "xts_sum_crg.out")
    f2 = os.path.join(dir2, pdb, "xts_sum_crg.out")

    # Check file existence
    f1_exists = os.path.isfile(f1)
    f2_exists = os.path.isfile(f2)

    if not f1_exists:
        print_log(f"      ├─ ❌ Missing: {x_label}/xts_sum_crg.out")
    if not f2_exists:
        print_log(f"      ├─ ❌ Missing: {y_label}/xts_sum_crg.out")

    if not (f1_exists and f2_exists):
        stats['files_missing'] += 1
        print_log(f"      └─ ⚠️  SKIPPED (files missing)")
        print_log()
        continue

    stats['files_found'] += 1

    # Parse files
    ch1, ch2 = parse_sum_crg(f1), parse_sum_crg(f2)
    shared = set(ch1.keys()) & set(ch2.keys())

    # Find matching residues
    matches = [(res, ch1[res], ch2[res]) for res in shared if res[:3] == icode]

    if matches:
        stats['pdbs_with_data'] += 1
        stats['residues_found'] += 1  # Count only 1 per PDB since we only plot 1 per PDB
        print_log(f"      ├─ ✅ Found {len(matches)} residue(s) matching '{icode}':")
        for res, val1, val2 in matches:
            print_log(f"      │    • {res}: {x_label}={val1:.4f}, {y_label}={val2:.4f}")
        
        # Only use the FIRST matching residue per PDB for plotting
        first_res, first_val1, first_val2 = matches[0]
        points.append((first_val1, first_val2, pdb, first_res, inhibitor))
        
        if len(matches) > 1:
            print_log(f"      │    ⚠️  Multiple matches found - using first residue: {first_res}")
        
        print_log(f"      └─ ✓ Data collected (1 point per PDB)")
        print_log()
    else:
        stats['pdbs_without_data'] += 1
        print_log(f"      ├─ ⚠️  No residues matching '{icode}' in shared residues")
        print_log(f"      │    Total shared residues: {len(shared)}")
        # Show a sample of what residues ARE present
        sample = sorted(list(shared))[:5]
        print_log(f"      │    Sample residues: {', '.join(sample)}")
        print_log(f"      └─ ⚠️  NO DATA COLLECTED")
        print_log()

# === Summary ===
print_log("="*80)
print_log("SUMMARY")
print_log("="*80)
print_log(f"  Total PDB entries:           {stats['total']}")
print_log(f"  Files found (both dirs):     {stats['files_found']}")
print_log(f"  Files missing:               {stats['files_missing']}")
print_log(f"  PDBs with matching data:     {stats['pdbs_with_data']}")
print_log(f"  PDBs without matching data:  {stats['pdbs_without_data']}")
print_log(f"  Total residue data points:   {stats['residues_found']}")
print_log(f"  Points actually plotted:     {len(points)}")
print_log()
print_log(f"  Outlier threshold:           ±{outlier_threshold} (15%)")
print_log("="*80)
print_log()

if not points:
    print_log("❌ No inhibitor residues found. Check pdb_inhibitor.lst and codes.")
    log_handle.close()
    raise SystemExit("❌ No inhibitor residues found. Check pdb_inhibitor.lst and codes.")

# === Assign colors by inhibitor (excluding red) ===
unique_inhibitors = sorted(set(p[4] for p in points))
# Use tab20 but exclude red-ish colors (indices 6, 7 are red shades)
palette = list(plt.cm.tab20.colors)
# Remove red colors from palette
non_red_palette = [palette[i] for i in range(len(palette)) if i not in [6, 7]]
if len(unique_inhibitors) > len(non_red_palette):
    times = (len(unique_inhibitors) // len(non_red_palette)) + 1
    non_red_palette = (non_red_palette * times)[:len(unique_inhibitors)]
inhibitor_colors = {inh: non_red_palette[i] for i, inh in enumerate(unique_inhibitors)}

# === Combined scatter plot (ALL INHIBITORS) ===
print_log()
print_log("="*80)
print_log("CREATING COMBINED PLOT (ALL INHIBITORS)")
print_log("="*80)

fig = plt.figure(figsize=(9.0, 9.0))
ax = plt.subplot(111)

# Calculate outliers (difference > threshold)
X = np.array([p[0] for p in points], dtype=float)
Y = np.array([p[1] for p in points], dtype=float)
outlier_mask = np.abs(Y - X) > outlier_threshold

# N is simply the number of points (1 per PDB)
N = len(points)
n_outliers = int(np.sum(outlier_mask))
n_fitted = N - n_outliers

for inh in unique_inhibitors:
    inh_idx = [i for i, p in enumerate(points) if p[4] == inh]
    xs = [points[i][0] for i in inh_idx]
    ys = [points[i][1] for i in inh_idx]
    
    # Add legend entry with standard black edge (no red for legend)
    ax.scatter([], [], label=f"{inh} (N={len(xs)})",
              color=inhibitor_colors[inh], s=80, edgecolor='black', 
              linewidth=0.8, alpha=0.85)
    
    # Plot actual points with outlier coloring
    for i, idx in enumerate(inh_idx):
        x_val, y_val = xs[i], ys[i]
        
        if outlier_mask[idx]:
            # Outliers: thin red edge
            edge_color = 'red'
            edge_width = 1.0
        else:
            # Normal points: standard black edge
            edge_color = 'black'
            edge_width = 0.8
        
        ax.scatter([x_val], [y_val],
                  color=inhibitor_colors[inh], s=80, edgecolor=edge_color, 
                  linewidth=edge_width, alpha=0.85, zorder=3)

# Compute padded independent limits first
x_min, x_max = float(X.min()), float(X.max())
y_min, y_max = float(Y.min()), float(Y.max())
x_pad = 0.05 * (x_max - x_min if x_max > x_min else 1.0)
y_pad = 0.05 * (y_max - y_min if y_max > y_min else 1.0)

# Apply padding to both ends - use negative padding if data doesn't go negative
x_lower = x_min - x_pad if x_min < 0 else -x_pad
y_lower = y_min - y_pad if y_min < 0 else -y_pad

ax.set_xlim(x_lower, x_max + x_pad)
ax.set_ylim(y_lower, y_max + y_pad)

# Enforce common symmetrical limits so y = x is fully visible across both axes
x0, x1 = ax.get_xlim()
y0, y1 = ax.get_ylim()
lims = [min(x0, y0), max(x1, y1)]
ax.set_xlim(lims)
ax.set_ylim(lims)

# Identity line (red)
ax.plot(lims, lims, 'r-', lw=1.2, label="y = x", zorder=2)

# Regression, stats, and best-fit line using actual plotted points
if N >= 2 and np.std(X) > 0 and np.std(Y) > 0:
    m, b = np.polyfit(X, Y, 1)

    # Draw the fit across the CURRENT x-limits
    x0, x1 = ax.get_xlim()
    x_line = np.linspace(x0, x1, 200)
    y_line = m * x_line + b
    #ax.plot(x_line, y_line, 'k--', lw=1.5, label="Best linear fit", zorder=2)

    ss_res = float(np.sum((Y - (m * X + b)) ** 2))
    ss_tot = float(np.sum((Y - Y.mean()) ** 2))
    R2 = 1 - ss_res / ss_tot if ss_tot != 0 else 1.0
    r = float(np.corrcoef(X, Y)[0, 1])
    rho = float(spearman_rho(X, Y))
    diff = Y - X
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    mae = float(np.mean(np.abs(diff)))
else:
    m = b = R2 = r = rho = rmse = mae = float('nan')

stats_text = (
    f"y = {m:.2f}x + {b:.2f}, R² = {R2:.2f}, N = {N}\n"
    f"r = {r:.2f}, ρ = {rho:.2f}\n"
    f"RMSE = {rmse:.2f}, MAE = {mae:.2f}\n"
    f"Fitted = {n_fitted}, Outliers = {n_outliers}"
)

# Labels & styling
ax.set_xlabel("Inhibitor Charge in Solution", fontweight="bold", fontsize=11)
ax.set_ylabel("Inhibitor Bound Charge", fontweight="bold", fontsize=11)
#ax.set_title(f"MCCE Inhibitor Charge (Protein Bound vs Solution)", fontweight="bold", fontsize=12)

# Grid styling (darker major gridlines)
ax.grid(which="major", color="gray", linestyle="-", linewidth=0.8, alpha=0.6, zorder=1)

# Corner annotations
ax.text(
    0.20, 0.90, "Bound State\nMore Charged",
    transform=ax.transAxes, fontsize=20, fontweight="bold",
    color="blue", va="top", ha="center"
)
ax.text(
    0.80, 0.10, "Bound State\nLess Charged",
    transform=ax.transAxes, fontsize=20, fontweight="bold",
    color="black", va="bottom", ha="center"
)

# Add stats to legend with bar chart icon
stats_handle = mlines.Line2D([], [], color='none', marker='|', markersize=15, 
                             markeredgewidth=2, markeredgecolor='gray', 
                             label=stats_text, linestyle='None')
handles, labels = ax.get_legend_handles_labels()
#handles.append(stats_handle)
#labels.append(stats_text)

# Legend centered below plot with stats
legend = ax.legend(
    handles, labels,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.08),
    ncol=4,
    fontsize=9,
    framealpha=0.9,
    handlelength=1.5,
    handleheight=2.5
)

plt.tight_layout()
plt.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.close()

print_log(f"✅ Combined plot saved: {fig_path}")
print_log(f"   N = {N} (1 point per PDB), Fitted: {n_fitted}, Outliers: {n_outliers}")

# === Save combined table ===
with open(csv_path, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["pdb", "residue", "inhibitor",
                f"x({x_label})", f"y({y_label})", "diff(y-x)"])
    for x, y, pdb_u, res, inh in points:
        w.writerow([pdb_u, res, inh, f"{x:.6f}", f"{y:.6f}", f"{y - x:.6f}"])

print_log(f"✅ Combined data saved: {csv_path}")

# === Individual inhibitor plots ===
print_log()
print_log("="*80)
print_log("CREATING INDIVIDUAL INHIBITOR PLOTS")
print_log("="*80)

for inh in unique_inhibitors:
    # Filter points for this inhibitor
    inh_points = [p for p in points if p[4] == inh]
    
    if not inh_points:
        continue
    
    print_log(f"  Plotting {inh} ({len(inh_points)} points)...")
    
    # Safe filename
    safe_inh_name = inh.replace("/", "_").replace(" ", "_")
    inh_fig_path = os.path.join(plot_dir, f"inhibitor_{safe_inh_name}.png")
    inh_csv_path = os.path.join(plot_dir, f"inhibitor_{safe_inh_name}.csv")
    
    # Extract data
    X_inh = np.array([p[0] for p in inh_points], dtype=float)
    Y_inh = np.array([p[1] for p in inh_points], dtype=float)
    pdb_names = [p[2] for p in inh_points]
    
    # Calculate outliers for this inhibitor
    outlier_mask_inh = np.abs(Y_inh - X_inh) > outlier_threshold
    
    # Create plot
    fig = plt.figure(figsize=(9.0, 9.0))
    ax = plt.subplot(111)
    
    # Group PDBs by base name (before first underscore)
    unique_pdbs = sorted(set(pdb_names))
    pdb_base_names = {}
    for pdb in unique_pdbs:
        base = pdb.split('_')[0]
        if base not in pdb_base_names:
            pdb_base_names[base] = []
        pdb_base_names[base].append(pdb)
    
    # Assign colors to base names - exclude red colors
    base_names_sorted = sorted(pdb_base_names.keys())
    base_colors_palette = list(plt.cm.tab10.colors)
    # Remove red color from tab10 (index 3 is red)
    non_red_base_palette = [base_colors_palette[i] for i in range(len(base_colors_palette)) if i != 3]
    if len(base_names_sorted) > len(non_red_base_palette):
        times = (len(base_names_sorted) // len(non_red_base_palette)) + 1
        non_red_base_palette = (non_red_base_palette * times)[:len(base_names_sorted)]
    
    base_to_color = {base: non_red_base_palette[i] for i, base in enumerate(base_names_sorted)}
    
    # Create color variations for each PDB
    def create_shade_variations(base_color, n_variants):
        """Create n highly discernible shades from a base color"""
        if n_variants == 1:
            return [base_color]
        
        # Convert RGB to HSV-like for better color control
        r, g, b = base_color[:3]
        shades = []
        
        for i in range(n_variants):
            # Create more dramatic variations using both brightness and saturation
            # Range from dark (0.4) to bright (1.0) with more spread
            brightness = 0.4 + (0.6 * i / max(1, n_variants - 1))
            
            # Also vary saturation for more distinction
            saturation = 0.6 + (0.4 * (n_variants - 1 - i) / max(1, n_variants - 1))
            
            # Apply transformations
            # Desaturate by mixing with gray, then adjust brightness
            gray = 0.5
            new_r = (r * saturation + gray * (1 - saturation)) * brightness
            new_g = (g * saturation + gray * (1 - saturation)) * brightness
            new_b = (b * saturation + gray * (1 - saturation)) * brightness
            
            # Ensure values are in [0, 1]
            new_r = max(0.0, min(1.0, new_r))
            new_g = max(0.0, min(1.0, new_g))
            new_b = max(0.0, min(1.0, new_b))
            
            shades.append((new_r, new_g, new_b))
        
        return shades
    
    # Assign specific colors to each PDB
    pdb_to_color = {}
    for base, pdbs in pdb_base_names.items():
        base_color = base_to_color[base]
        shades = create_shade_variations(base_color, len(pdbs))
        for i, pdb in enumerate(pdbs):
            pdb_to_color[pdb] = shades[i]
    
    # Plot each PDB with its assigned color and outlier edge coloring
    for pdb_idx, pdb in enumerate(unique_pdbs):
        pdb_point_indices = [i for i, p in enumerate(inh_points) if p[2] == pdb]
        xs = [inh_points[i][0] for i in pdb_point_indices]
        ys = [inh_points[i][1] for i in pdb_point_indices]
        
        # Add legend entry with standard black edge (no red for legend)
        ax.scatter([], [], label=pdb, 
                  color=pdb_to_color[pdb], s=90, edgecolor='black', 
                  linewidth=0.8, alpha=0.85)
        
        # Plot each point with appropriate edge color
        for pt_idx, (x_val, y_val) in enumerate(zip(xs, ys)):
            if outlier_mask_inh[pdb_point_indices[pt_idx]]:
                # Outliers: thin red edge
                edge_color = 'red'
                edge_width = 1.0
            else:
                # Normal points: standard black edge
                edge_color = 'black'
                edge_width = 0.8
            
            ax.scatter([x_val], [y_val],
                      color=pdb_to_color[pdb], s=90, edgecolor=edge_color, 
                      linewidth=edge_width, alpha=0.85, zorder=3)
            
            # Add PDB labels to points
            ax.annotate(pdb, (x_val, y_val), 
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7)
    
    # Compute padded limits - apply same padding to both ends
    x_min, x_max = float(X_inh.min()), float(X_inh.max())
    y_min, y_max = float(Y_inh.min()), float(Y_inh.max())
    x_pad = 0.08 * (x_max - x_min if x_max > x_min else 1.0)
    y_pad = 0.08 * (y_max - y_min if y_max > y_min else 1.0)
    
    # Apply padding to both ends - use negative padding if data doesn't go negative
    x_lower = x_min - x_pad if x_min < 0 else -x_pad
    y_lower = y_min - y_pad if y_min < 0 else -y_pad
    
    ax.set_xlim(x_lower, x_max + x_pad)
    ax.set_ylim(y_lower, y_max + y_pad)
    
    # Enforce symmetrical limits
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    lims = [min(x0, y0), max(x1, y1)]
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    # Identity line (red)
    ax.plot(lims, lims, 'r-', lw=1.2, label="y = x", zorder=2)
    
    # Statistics and regression
    # N is simply the number of points (1 per PDB)
    N_inh = len(inh_points)
    n_outliers_inh = int(np.sum(outlier_mask_inh))
    n_fitted_inh = N_inh - n_outliers_inh
    
    if N_inh >= 2 and np.std(X_inh) > 0 and np.std(Y_inh) > 0:
        m_inh, b_inh = np.polyfit(X_inh, Y_inh, 1)
        
        x0, x1 = ax.get_xlim()
        x_line = np.linspace(x0, x1, 200)
        y_line = m_inh * x_line + b_inh
        #ax.plot(x_line, y_line, 'k--', lw=1.5, label="Best linear fit", zorder=2)
        
        ss_res = float(np.sum((Y_inh - (m_inh * X_inh + b_inh)) ** 2))
        ss_tot = float(np.sum((Y_inh - Y_inh.mean()) ** 2))
        R2_inh = 1 - ss_res / ss_tot if ss_tot != 0 else 1.0
        r_inh = float(np.corrcoef(X_inh, Y_inh)[0, 1])
        rho_inh = float(spearman_rho(X_inh, Y_inh))
        diff_inh = Y_inh - X_inh
        rmse_inh = float(np.sqrt(np.mean(diff_inh ** 2)))
        mae_inh = float(np.mean(np.abs(diff_inh)))
    else:
        m_inh = b_inh = R2_inh = r_inh = rho_inh = rmse_inh = mae_inh = float('nan')
    
    stats_text_inh = (
        f"y = {m_inh:.2f}x + {b_inh:.2f}, R² = {R2_inh:.2f}, N = {N_inh}\n"
        f"r = {r_inh:.2f}, ρ = {rho_inh:.2f}\n"
        f"RMSE = {rmse_inh:.2f}, MAE = {mae_inh:.2f}\n"
        f"Fitted = {n_fitted_inh}, Outliers = {n_outliers_inh}"
    )
    
    # Labels
    ax.set_xlabel("Boltzmann Weighted Solution Charge", fontweight="bold", fontsize=14)
    ax.set_ylabel("Boltzmann Weighted Bound Charge", fontweight="bold", fontsize=14)
    #ax.set_title(f"MCCE {inh} Charge (Protein Bound vs Solution)", fontweight="bold", fontsize=12)
    
    # Grid
    ax.grid(which="major", color="gray", linestyle="-", linewidth=0.8, alpha=0.6, zorder=1)
    
    # Corner annotations
    ax.text(
        0.20, 0.85, "Bound State\nMore Charged",
        transform=ax.transAxes, fontsize=20, fontweight="bold",
        color="blue", va="top", ha="center"
    )
    ax.text(
        0.80, 0.10, "Bound State\nLess Charged",
        transform=ax.transAxes, fontsize=20, fontweight="bold",
        color="black", va="bottom", ha="center"
    )
    
    # Add stats to legend with bar chart icon
    stats_handle = mlines.Line2D([], [], color='none', marker='|', markersize=15, 
                                 markeredgewidth=2, markeredgecolor='gray', 
                                 label=stats_text_inh, linestyle='None')
    handles, labels = ax.get_legend_handles_labels()
    handles.append(stats_handle)
    labels.append(stats_text_inh)
    
    # Legend centered below plot with stats
    n_legend_items = len(handles)
    ncol = min(4, max(2, n_legend_items // 3))
    legend = ax.legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=ncol,
        fontsize=9,
        framealpha=0.9,
        handlelength=1.5,
        handleheight=2.5
    )
    
    plt.tight_layout()
    plt.savefig(inh_fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Save individual CSV
    with open(inh_csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pdb", "residue", "inhibitor",
                    f"x({x_label})", f"y({y_label})", "diff(y-x)"])
        for x, y, pdb_u, res, inh_name in inh_points:
            w.writerow([pdb_u, res, inh_name, f"{x:.6f}", f"{y:.6f}", f"{y - x:.6f}"])
    
    print_log(f"    ✅ Plot: {inh_fig_path}")
    print_log(f"       N = {N_inh} (1 point per PDB), Fitted: {n_fitted_inh}, Outliers: {n_outliers_inh}")
    print_log(f"    ✅ Data: {inh_csv_path}")

print_log()
print_log("="*80)
print_log("ALL PLOTS COMPLETED")
print_log("="*80)
print_log(f"✅ Log saved: {log_file}")

# === Close log file ===
log_handle.close()
