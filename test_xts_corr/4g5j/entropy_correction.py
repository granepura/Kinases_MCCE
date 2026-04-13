#!/usr/bin/env python3
"""
Entropy Correction for fort.38 File
Applies conformational entropy corrections to non-amino acid residues

Formula: E_TS = -RT Σ(P_i ln P_i)
where:
- P_i is the normalized probability of conformer i within its charge state group
- R = 1.987×10⁻³ kcal·K⁻¹·mol⁻¹
- T = 298 K
- RT = 0.592 kcal/mol
"""

import re
import math
from collections import defaultdict
from typing import Dict, List, Tuple

# Constants
R = 1.987e-3  # kcal·K⁻¹·mol⁻¹
T = 298.0     # K
RT = R * T    # = 0.592 kcal/mol

# Amino acids to exclude from entropy correction
AMINO_ACIDS = {
    'ACE', 'NME', 'CTR', 'NTR', 'CTG', 'NTG',
    'ALA', 'ARG', 'ASN', 'ASP', 'CYD', 'GLN', 'GLU',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL', 'CYS'
}


class Conformer:
    """Represents a single conformer"""
    def __init__(self, name: str, probabilities: List[float]):
        self.full_name = name
        self.probabilities = probabilities
        self.original_probs = probabilities.copy()
        
        # Parse the conformer name
        parsed = self._parse_name(name)
        if parsed:
            self.res_type = parsed['res_type']
            self.charge_state = parsed['charge_state']
            self.chain = parsed['chain']
            self.res_num = parsed['res_num']
            self.conf_num = parsed['conf_num']
            self.full_res_id = f"{self.res_type}{self.chain}{self.res_num}"
            self.charge_group = f"{self.res_type}{self.charge_state}"
        else:
            raise ValueError(f"Cannot parse conformer name: {name}")
    
    @staticmethod
    def _parse_name(name: str) -> Dict[str, str]:
        """Parse conformer name: e.g., '0WN+1A1101_001'"""
        parts = name.split('_')
        if len(parts) != 2:
            return None
        
        conf_num = parts[1]
        main_part = parts[0]
        
        # Pattern: RESTYPE(3 chars) + CHARGE_STATE + CHAIN(1 letter) + RESNUM(4 digits)
        match = re.match(r'^([A-Z0-9]{3})([\+\-]?\w+)([A-Z])(\d{4})$', main_part)
        if not match:
            return None
        
        return {
            'res_type': match.group(1),
            'charge_state': match.group(2),
            'chain': match.group(3),
            'res_num': match.group(4),
            'conf_num': conf_num
        }


def calculate_ets(normalized_probs: List[float]) -> float:
    """Calculate entropy correction E_TS = -RT Σ(P_i ln P_i)"""
    total = 0.0
    for p in normalized_probs:
        if p > 1e-10:
            total += p * math.log(p)
    return -RT * total


def apply_entropy_correction(conformers: List[Conformer], ph_index: int) -> Tuple[List[float], Dict]:
    """Apply entropy correction to a group of conformers at a given pH"""
    
    # Group conformers by charge state
    charge_groups = defaultdict(list)
    for conf in conformers:
        charge_groups[conf.charge_group].append(conf)
    
    # Calculate E_TS for each charge group
    ets_values = {}
    group_info = {}
    
    for charge_group, confs in charge_groups.items():
        probs = [c.probabilities[ph_index] for c in confs]
        total_prob = sum(probs)
        
        if total_prob < 1e-10:
            # Unpopulated group
            ets_values[charge_group] = 0.0
            group_info[charge_group] = {
                'total_prob': total_prob,
                'normalized_probs': probs,
                'ets': 0.0,
                'num_conformers': len(confs)
            }
            continue
        
        # Normalize within group
        normalized_probs = [p / total_prob for p in probs]
        ets = calculate_ets(normalized_probs)
        
        ets_values[charge_group] = ets
        group_info[charge_group] = {
            'total_prob': total_prob,
            'normalized_probs': normalized_probs,
            'ets': ets,
            'num_conformers': len(confs)
        }
    
    # Apply corrections: probs -> energies -> correct -> back to probs
    # Use numerically stable approach
    corrected_energies = []
    for conf in conformers:
        prob = conf.probabilities[ph_index]
        if prob > 1e-10:
            energy = -math.log(prob)
        else:
            # For very small probabilities, use a large but finite energy
            energy = 30.0  # Corresponds to prob ~ 1e-13
        
        # Apply entropy correction for this conformer's charge group
        ets = ets_values[conf.charge_group]
        corrected_energy = energy - ets
        corrected_energies.append(corrected_energy)
    
    # Shift energies to avoid overflow: subtract minimum
    min_energy = min(corrected_energies)
    shifted_energies = [e - min_energy for e in corrected_energies]
    
    # Convert back to probabilities (now numerically stable)
    exp_terms = []
    for e in shifted_energies:
        if e > 100:  # Still too large
            exp_terms.append(0.0)
        else:
            exp_terms.append(math.exp(-e))
    
    partition = sum(exp_terms)
    if partition < 1e-100:
        # All probabilities are essentially zero, return original
        return [c.probabilities[ph_index] for c in conformers], group_info
    
    new_probs = [exp_val / partition for exp_val in exp_terms]
    
    return new_probs, group_info


def process_fort38(input_file: str, output_file: str, log_file: str):
    """Process fort.38 file and apply entropy corrections"""
    
    # Read the file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    header = lines[0].strip().split()
    ph_values = [float(x) for x in header[1:]]
    num_ph = len(ph_values)
    
    print(f"Processing {len(lines)-1} conformers at {num_ph} pH values")
    
    # Parse all conformers
    all_conformers = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        
        parts = line.split()
        name = parts[0]
        probs = [float(x) for x in parts[1:]]
        
        try:
            conf = Conformer(name, probs)
            all_conformers.append(conf)
        except ValueError as e:
            print(f"Warning: {e}")
            continue
    
    # Group by residue
    residue_groups = defaultdict(list)
    for conf in all_conformers:
        residue_groups[conf.full_res_id].append(conf)
    
    # Process each residue
    log_entries = []
    num_corrected = 0
    
    for res_id, res_conformers in residue_groups.items():
        res_type = res_conformers[0].res_type
        is_amino_acid = res_type in AMINO_ACIDS
        
        if not is_amino_acid:
            num_corrected += 1
            print(f"Applying entropy correction to {res_id} (type: {res_type})")
            
            # Apply correction at each pH
            for ph_idx in range(num_ph):
                new_probs, group_info = apply_entropy_correction(res_conformers, ph_idx)
                
                # Update probabilities
                for i, conf in enumerate(res_conformers):
                    conf.probabilities[ph_idx] = new_probs[i]
                
                # Log significant corrections
                for charge_group, info in group_info.items():
                    if info['total_prob'] > 1e-10 and info['ets'] > 1e-6:
                        log_entries.append({
                            'residue': res_id,
                            'res_type': res_type,
                            'ph': ph_values[ph_idx],
                            'charge_group': charge_group,
                            'num_conformers': info['num_conformers'],
                            'ets': info['ets']
                        })
    
    print(f"\nCorrected {num_corrected} non-amino acid residues")
    
    # Write output file
    with open(output_file, 'w') as f:
        # Write header - match exact fort.38 format
        # ' ph' (3 chars) + 14 spaces + pH values right-aligned in 6-char fields
        f.write(' ph              ')
        f.write(''.join(f'{ph:6.1f}' for ph in ph_values))
        f.write('\n')
        
        # Write conformer data - match exact fort.38 format
        # Conformer name (14 chars) + 1 space + probabilities in 6-char fields
        for conf in all_conformers:
            f.write(conf.full_name)  # Name is exactly 14 characters
            f.write(' ')              # Single space separator
            f.write(''.join(f'{p:6.3f}' for p in conf.probabilities))
            f.write('\n')
    
    print(f"Wrote corrected file to: {output_file}")
    
    # Write log file
    with open(log_file, 'w') as f:
        f.write('ENTROPY CORRECTION LOG\n')
        f.write('=' * 80 + '\n')
        f.write('Applied to non-amino acid residues only\n')
        f.write('Formula: E_TS = -RT Σ(P_i ln P_i)\n')
        f.write(f'R = {R} kcal·K⁻¹·mol⁻¹\n')
        f.write(f'T = {T} K\n')
        f.write(f'RT = {RT:.3f} kcal/mol\n\n')
        
        f.write('Amino acids excluded from correction:\n')
        f.write(', '.join(sorted(AMINO_ACIDS)) + '\n\n')
        
        f.write('=' * 80 + '\n')
        f.write('NON-AMINO ACID RESIDUES PROCESSED:\n')
        f.write('=' * 80 + '\n\n')
        
        # List all non-amino acid residues
        for res_id, res_conformers in sorted(residue_groups.items()):
            res_type = res_conformers[0].res_type
            if res_type not in AMINO_ACIDS:
                f.write(f'Residue: {res_id} (Type: {res_type})\n')
                f.write(f'  Conformers:\n')
                for conf in res_conformers:
                    f.write(f'    - {conf.full_name} (charge group: {conf.charge_group})\n')
                f.write('\n')
        
        f.write('=' * 80 + '\n')
        f.write('ENTROPY CORRECTIONS APPLIED:\n')
        f.write('=' * 80 + '\n\n')
        
        if not log_entries:
            f.write('No significant entropy corrections applied.\n')
            f.write('(All charge state groups had single conformers or zero population)\n\n')
        else:
            f.write(f'Total corrections applied: {len(log_entries)}\n\n')
            
            # Group log entries by residue
            by_residue = defaultdict(list)
            for entry in log_entries:
                by_residue[entry['residue']].append(entry)
            
            for res_id, entries in sorted(by_residue.items()):
                f.write(f"Residue: {res_id} ({entries[0]['res_type']})\n")
                for entry in entries:
                    f.write(f"  pH {entry['ph']:.1f}, Charge Group: {entry['charge_group']}\n")
                    f.write(f"    Number of conformers: {entry['num_conformers']}\n")
                    f.write(f"    E_TS correction: {entry['ets']:.4f} kcal/mol\n\n")
        
        f.write('=' * 80 + '\n')
        f.write('SUMMARY:\n')
        f.write('=' * 80 + '\n\n')
        f.write(f"Total residues processed: {len(residue_groups)}\n")
        f.write(f"Non-amino acid residues: {num_corrected}\n")
        f.write(f"Amino acid residues (unchanged): {len(residue_groups) - num_corrected}\n")
        f.write(f"Significant entropy corrections: {len(log_entries)}\n")
    
    print(f"Wrote log file to: {log_file}")


if __name__ == '__main__':
    import sys
    import os
    
    # Default to fort.38 if no argument provided
    if len(sys.argv) == 1:
        input_file = 'fort.38'
    elif len(sys.argv) == 2:
        input_file = sys.argv[1]
    else:
        print("Usage: python entropy_correction.py [fort.38]")
        print("\nOptional argument: path to fort.38 file (default: fort.38)")
        print("\nThis script will:")
        print("  - Read the fort.38 file")
        print("  - Apply entropy corrections to non-amino acids only")
        print("  - Generate xts_fort.38 (corrected file)")
        print("  - Generate entropy_correction.log (detailed log)")
        sys.exit(1)
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found!")
        print("\nPlease ensure fort.38 exists in the current directory,")
        print("or specify the correct path as an argument.")
        sys.exit(1)
    
    output_file = 'xts_fort.38'
    log_file = 'entropy_correction.log'
    
    print("=" * 60)
    print("ENTROPY CORRECTION FOR fort.38")
    print("=" * 60)
    print(f"Input file: {input_file}")
    print()
    
    try:
        process_fort38(input_file, output_file, log_file)
        print("\n" + "=" * 60)
        print("COMPLETE!")
        print("=" * 60)
        print(f"\nGenerated files:")
        print(f"  - {output_file} (corrected probabilities)")
        print(f"  - {log_file} (detailed correction log)")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
