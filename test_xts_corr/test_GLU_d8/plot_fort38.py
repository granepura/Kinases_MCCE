#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# File mapping with user labels
file_mapping = {
    'no_xts': 'fort.38',
    'xts_jj': 'jj_fort.38',
    'xts_gr': 'xts_fort.38'
}

def read_fort38_file(filename):
    """
    Read fort.38 file and extract GLU-1A0292_005 data
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    pH_values = None
    glu_fraction = None
    
    for i, line in enumerate(lines):
        # Find the pH header line
        if line.strip().startswith('ph'):
            pH_values = [float(x) for x in line.split()[1:]]
        
        # Find the GLU-1 line (deprotonated form)
        if 'GLU-1A0292_005' in line:
            glu_fraction = [float(x) for x in line.split()[1:]]
            break
    
    if pH_values is None or glu_fraction is None:
        raise ValueError(f"Could not find pH values or GLU-1A0292_005 data in {filename}")
    
    return np.array(pH_values), np.array(glu_fraction)

# Henderson-Hasselbalch function for deprotonation
# fraction_deprotonated = 1 / (1 + 10^(pKa - pH))
def henderson_hasselbalch(pH, pKa):
    return 1 / (1 + 10**(pKa - pH))

# Read data from files
datasets = {}
for label, filename in file_mapping.items():
    try:
        pH, fraction = read_fort38_file(filename)
        datasets[label] = {'pH': pH, 'fraction': fraction}
        print(f"Successfully read {filename} as '{label}'")
    except Exception as e:
        print(f"Error reading {filename}: {e}")

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8))

colors = {'no_xts': 'blue', 'xts_jj': 'red', 'xts_gr': 'green'}
markers = {'no_xts': 'o', 'xts_jj': 's', 'xts_gr': '^'}

# For smooth curve plotting
pH_smooth = np.linspace(0, 14, 300)

print("\n" + "=" * 60)
print("GLU-1A0292_005 pKa Analysis")
print("=" * 60)

for label, data in datasets.items():
    pH = data['pH']
    fraction = data['fraction']
    
    # Filter out data points where fraction is 0 or 1 for better fitting
    # (these saturated points can skew the fit)
    mask = (fraction > 0.01) & (fraction < 0.99)
    pH_fit = pH[mask]
    fraction_fit = fraction[mask]
    
    # Fit the Henderson-Hasselbalch curve
    try:
        popt, pcov = curve_fit(henderson_hasselbalch, pH_fit, fraction_fit, p0=[5.0])
        pKa_fitted = popt[0]
        pKa_error = np.sqrt(np.diag(pcov))[0]
        
        # Calculate R-squared
        residuals = fraction_fit - henderson_hasselbalch(pH_fit, pKa_fitted)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((fraction_fit - np.mean(fraction_fit))**2)
        r_squared = 1 - (ss_res / ss_tot)
        
        print(f"\n{label} ({file_mapping[label]}):")
        print(f"  Fitted pKa: {pKa_fitted:.2f} ± {pKa_error:.2f}")
        print(f"  R²: {r_squared:.4f}")
        
        # Plot data points
        ax.scatter(pH, fraction, color=colors[label], marker=markers[label], s=100, 
                  label=f'{label} (data)', zorder=3, edgecolors='black', linewidth=1)
        
        # Plot fitted curve
        fraction_smooth = henderson_hasselbalch(pH_smooth, pKa_fitted)
        ax.plot(pH_smooth, fraction_smooth, color=colors[label], linewidth=2.5, 
               label=f'{label} (fit, pKa={pKa_fitted:.2f})', linestyle='--', alpha=0.8)
        
        # Mark the pKa position
        ax.axvline(x=pKa_fitted, color=colors[label], linestyle=':', linewidth=1.5, alpha=0.5)
        ax.plot(pKa_fitted, 0.5, 'D', color=colors[label], markersize=10, 
               markeredgecolor='black', markeredgewidth=1, zorder=4)
        
    except Exception as e:
        print(f"\n{label}: Fitting failed - {e}")
        ax.scatter(pH, fraction, color=colors[label], marker=markers[label], s=100,
                  label=f'{label} (data)', zorder=3)

# Format the plot
ax.set_xlabel('pH', fontsize=14, fontweight='bold')
ax.set_ylabel('Fraction Deprotonated (GLU⁻)', fontsize=14, fontweight='bold')
ax.set_title('GLU-1A0292_005: Henderson-Hasselbalch Titration Curves', 
            fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(-0.5, 14.5)
ax.set_ylim(-0.05, 1.05)
ax.legend(loc='upper left', fontsize=10, framealpha=0.9)

# Add horizontal line at 0.5 (pKa point)
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.text(14, 0.52, 'pKa (50% deprotonated)', fontsize=10, ha='right', va='bottom')

plt.tight_layout()
plt.savefig('GLU1_pKa_fit.png', dpi=300, bbox_inches='tight')
print("\n" + "=" * 60)
print("Plot saved as 'GLU1_pKa_fit.png'")
print("=" * 60)
plt.show()
