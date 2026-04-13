#!/usr/bin/env python3
"""
Created on Oct 15 2025
@author: Gehan Ranepura

Name: xts_corr.py
Applies conformational entropy corrections to non-amino acid residues for file fort.38

ENTROPY CORRECTION METHOD:
==========================
This script corrects for the bias that arises when different charge states have 
different numbers of conformers. Charge states with more conformers would otherwise 
appear more probable simply due to having more "counts" in the ensemble.

The correction penalizes conformers belonging to charge states with multiple conformers.

EXAMPLE with GLU conformers:
----------------------------
Original probabilities:
  GLU01A0292_001: 0.005  }
  GLU01A0292_002: 0.007  } charge state 01, N=2
  GLU02A0292_003: 0.016  }
  GLU02A0292_004: 0.012  } charge state 02, N=2
  GLU-1A0292_005: 0.960  } charge state -1, N=1

Step 1: Convert probabilities to relative free energies
    ΔG_i/RT = -ln(P_i/P_ref)
    
    Using the most probable conformer as reference (P_ref = 0.960):
    
    GLU01A0292_001: -ln(0.005/0.960) = 5.257
    GLU01A0292_002: -ln(0.007/0.960) = 4.921
    GLU02A0292_003: -ln(0.016/0.960) = 4.094
    GLU02A0292_004: -ln(0.012/0.960) = 4.382
    GLU-1A0292_005: -ln(0.960/0.960) = 0.000
    
Step 2: Apply entropy penalty based on conformational degeneracy
    ΔG_corrected/RT = ΔG_original/RT + ln(N)
    
    where N = number of conformers in that charge state group
    
    Penalties:
      01 (N=2): ln(2) = 0.693
      02 (N=2): ln(2) = 0.693
      -1 (N=1): ln(1) = 0.000
    
    Corrected energies:
      GLU01A0292_001: 5.257 + 0.693 = 5.950
      GLU01A0292_002: 4.921 + 0.693 = 5.614
      GLU02A0292_003: 4.094 + 0.693 = 4.787
      GLU02A0292_004: 4.382 + 0.693 = 5.075
      GLU-1A0292_005: 0.000 + 0.000 = 0.000
    
Step 3: Convert back to probabilities using Boltzmann distribution
    P_i = exp(-ΔG_corrected/RT) / Σ_j exp(-ΔG_corrected/RT)
    
    Boltzmann factors (exp(-ΔG)):
      GLU01A0292_001: exp(-5.950) = 0.002605
      GLU01A0292_002: exp(-5.614) = 0.003649
      GLU02A0292_003: exp(-4.787) = 0.008339
      GLU02A0292_004: exp(-5.075) = 0.006271
      GLU-1A0292_005: exp(-0.000) = 1.000000
    
    Partition function Z = 1.020864
    
    Final corrected probabilities:
      GLU01A0292_001: 0.002605/1.020864 = 0.0026 (0.26%)  [was 0.50%]
      GLU01A0292_002: 0.003649/1.020864 = 0.0036 (0.36%)  [was 0.70%]
      GLU02A0292_003: 0.008339/1.020864 = 0.0082 (0.82%)  [was 1.60%]
      GLU02A0292_004: 0.006271/1.020864 = 0.0061 (0.61%)  [was 1.20%]
      GLU-1A0292_005: 1.000000/1.020864 = 0.9796 (97.96%) [was 96.00%]

RESULT: Neutral conformers (01, 02) are penalized ~50% due to having 2 conformers 
each, while the ionized state (-1) with only 1 conformer receives no penalty and 
becomes even more dominant.

EFFECT:
- Charge states with more conformers (larger N) get penalized more
- Each conformer in a multi-conformer charge state has ln(N) added to its energy
- This lowers the probability of conformers in charge states with many conformers
- Charge states with single conformers (N=1) get no penalty (ln(1) = 0)

CONSTANTS:
- R = 1.987×10⁻³ kcal·K⁻¹·mol⁻¹
- T = 298 K
- RT = 0.592 kcal/mol
- Energy penalty = RT × ln(N) kcal/mol
"""

import re
import math
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple

# Constants
R = 1.987e-3  # kcal·K⁻¹·mol⁻¹
T = 298.0     # K
RT = R * T    # = 0.592 kcal/mol

# Amino acids to exclude from entropy correction (unless --doAAs flag is used)
AMINO_ACIDS_LIST = {
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


def parse_charge_state(charge_state: str) -> float:
    """Parse charge state string to numeric value

    Examples:
        '+1' -> 1.0
        '+2' -> 2.0
        '-'  -> -1.0
        '-2' -> -2.0
        '01' -> 0.0
        '02' -> 0.0
        '+a' -> 2.0 (doubly protonated)
        '+b' -> 3.0 (triply protonated)
    """
    if not charge_state:
        return 0.0

    # Handle numeric formats
    if charge_state[0] in '+-':
        # '+1', '-1', '+2', etc.
        if len(charge_state) > 1 and charge_state[1:].isdigit():
            return float(charge_state)
        # Just '+' means +1, just '-' means -1
        elif charge_state == '+':
            return 1.0
        elif charge_state == '-':
            return -1.0
        # Handle letter codes: '+a' = +2, '+b' = +3, etc.
        elif len(charge_state) > 1 and charge_state[1].isalpha():
            letter_charge = ord(charge_state[1].lower()) - ord('a') + 2
            return float(letter_charge) if charge_state[0] == '+' else -float(letter_charge)

    # Handle '01', '02' formats (neutral)
    if charge_state.startswith('0'):
        return 0.0

    # Default to 0
    return 0.0


def apply_entropy_correction(conformers: List[Conformer], ph_index: int) -> Tuple[List[float], Dict]:
    """Apply Boltzmann entropy correction to penalize charge states with more conformers
    
    This corrects for the bias where charge states with more conformers appear more
    probable simply because they have more "counts" in the ensemble.
    
    Method:
    1. Convert probabilities to relative free energies: ΔG/RT = -ln(P_i/P_ref)
    2. Add penalty: ΔG_corrected = ΔG_original + ln(N), where N = number of conformers
    3. Convert back to probabilities using Boltzmann distribution
    """

    # Group conformers by charge state
    charge_groups = defaultdict(list)
    for conf in conformers:
        charge_groups[conf.charge_group].append(conf)

    # Calculate number of conformers in each charge group
    group_info = {}
    for charge_group, confs in charge_groups.items():
        num_conformers = len(confs)
        entropy_penalty = math.log(num_conformers) if num_conformers > 0 else 0.0
        
        probs = [c.probabilities[ph_index] for c in confs]
        total_prob = sum(probs)
        
        group_info[charge_group] = {
            'total_prob': total_prob,
            'num_conformers': num_conformers,
            'entropy_penalty': entropy_penalty,
            'penalty_kcal': entropy_penalty * RT
        }

    # Step 1: Find reference (most probable conformer)
    probs = [c.probabilities[ph_index] for c in conformers]
    p_ref = max(probs)
    
    # Step 2: Calculate relative free energies using reference
    relative_energies = []
    for conf in conformers:
        prob = conf.probabilities[ph_index]
        
        # ΔG/RT = -ln(P_i/P_ref)
        if prob > 1e-10:
            relative_energy = -math.log(prob / p_ref)
        else:
            # For very small probabilities, use a large but finite energy
            relative_energy = 30.0
        
        relative_energies.append(relative_energy)
    
    # Step 3: Apply entropy penalty
    corrected_energies = []
    for i, conf in enumerate(conformers):
        penalty = group_info[conf.charge_group]['entropy_penalty']
        corrected_energy = relative_energies[i] + penalty  # ADD penalty
        corrected_energies.append(corrected_energy)

    # Step 4: Calculate Boltzmann factors
    boltzmann_factors = []
    for e in corrected_energies:
        if e > 100:  # Prevent overflow
            boltzmann_factors.append(0.0)
        else:
            boltzmann_factors.append(math.exp(-e))

    # Step 5: Calculate partition function and normalize
    partition = sum(boltzmann_factors)
    if partition < 1e-100:
        # All probabilities are essentially zero, return original
        return [c.probabilities[ph_index] for c in conformers], group_info

    new_probs = [bf / partition for bf in boltzmann_factors]

    return new_probs, group_info


def process_fort38(input_file: str, output_file: str, log_file: str, amino_acids: set):
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
        is_amino_acid = res_type in amino_acids

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
                    if info['total_prob'] > 1e-10 and info['num_conformers'] > 1:
                        log_entries.append({
                            'residue': res_id,
                            'res_type': res_type,
                            'ph': ph_values[ph_idx],
                            'charge_group': charge_group,
                            'num_conformers': info['num_conformers'],
                            'entropy_penalty': info['entropy_penalty'],
                            'penalty_kcal': info['penalty_kcal']
                        })

    if amino_acids:
        print(f"\nCorrected {num_corrected} non-amino acid residues")
    else:
        print(f"\nCorrected {num_corrected} residues (including amino acids)")

    # Write output file - match exact fort.38 format
    with open(output_file, 'w') as f:
        # Write header
        f.write(' ph            ')
        f.write(''.join(f'{ph:>5.1f} ' for ph in ph_values))
        f.write('\n')

        # Write conformer data
        for conf in all_conformers:
            f.write(conf.full_name)
            f.write(' ')
            f.write(''.join(f'{p:5.3f} ' for p in conf.probabilities))
            f.write('\n')

    print(f"Wrote corrected file to: {output_file}")

    # Write log file
    with open(log_file, 'w') as f:
        f.write('BOLTZMANN ENTROPY CORRECTION LOG\n')
        f.write('=' * 80 + '\n')
        if amino_acids:
            f.write('Applied to non-amino acid residues only\n\n')
        else:
            f.write('Applied to ALL residues (including amino acids)\n\n')
        
        f.write('METHOD:\n')
        f.write('-------\n')
        f.write('Corrects for bias when charge states have different numbers of conformers.\n')
        f.write('Charge states with more conformers are penalized to account for the\n')
        f.write('entropic cost of selecting a specific conformer.\n\n')
        
        f.write('PROCEDURE:\n')
        f.write('1. Convert probabilities to free energies: ΔG/RT = -ln(P)\n')
        f.write('2. Add entropy penalty: ΔG_corrected/RT = ΔG_original/RT + ln(N)\n')
        f.write('   where N = number of conformers in that charge state\n')
        f.write('3. Convert back to probabilities: P_i = exp(-ΔG_corrected/RT) / Z\n\n')
        
        f.write('CONSTANTS:\n')
        f.write(f'R = {R} kcal·K⁻¹·mol⁻¹\n')
        f.write(f'T = {T} K\n')
        f.write(f'RT = {RT:.3f} kcal/mol\n\n')
        
        f.write('PENALTIES:\n')
        f.write('N=1: ln(1) = 0.000 → No penalty (0.000 kcal/mol)\n')
        f.write('N=2: ln(2) = 0.693 → Penalty = 0.410 kcal/mol\n')
        f.write('N=3: ln(3) = 1.099 → Penalty = 0.651 kcal/mol\n')
        f.write('N=4: ln(4) = 1.386 → Penalty = 0.821 kcal/mol\n\n')

        f.write('Amino acids excluded from correction:\n')
        if amino_acids:
            f.write(', '.join(sorted(amino_acids)) + '\n\n')
        else:
            f.write('(None - all residues corrected)\n\n')

        f.write('=' * 80 + '\n')
        if amino_acids:
            f.write('NON-AMINO ACID RESIDUES PROCESSED:\n')
        else:
            f.write('ALL RESIDUES PROCESSED:\n')
        f.write('=' * 80 + '\n\n')

        # List all non-amino acid residues
        for res_id, res_conformers in sorted(residue_groups.items()):
            res_type = res_conformers[0].res_type
            if res_type not in amino_acids:
                f.write(f'Residue: {res_id} (Type: {res_type})\n')
                f.write(f'  Conformers:\n')
                for conf in res_conformers:
                    f.write(f'    - {conf.full_name} (charge group: {conf.charge_group})\n')
                f.write('\n')

        f.write('=' * 80 + '\n')
        f.write('ENTROPY PENALTIES APPLIED:\n')
        f.write('=' * 80 + '\n\n')

        if not log_entries:
            f.write('No penalties applied.\n')
            f.write('(All charge state groups had single conformers)\n\n')
        else:
            f.write(f'Total penalties applied: {len(log_entries)}\n\n')

            # Group log entries by residue
            by_residue = defaultdict(list)
            for entry in log_entries:
                by_residue[entry['residue']].append(entry)

            for res_id, entries in sorted(by_residue.items()):
                f.write(f"Residue: {res_id} ({entries[0]['res_type']})\n")
                for entry in entries:
                    f.write(f"  pH {entry['ph']:.1f}, Charge Group: {entry['charge_group']}\n")
                    f.write(f"    Number of conformers (N): {entry['num_conformers']}\n")
                    f.write(f"    Entropy penalty (ln N): {entry['entropy_penalty']:.4f}\n")
                    f.write(f"    Energy penalty: {entry['penalty_kcal']:.4f} kcal/mol\n\n")

        f.write('=' * 80 + '\n')
        f.write('SUMMARY:\n')
        f.write('=' * 80 + '\n\n')
        f.write(f"Total residues processed: {len(residue_groups)}\n")
        if amino_acids:
            f.write(f"Non-amino acid residues corrected: {num_corrected}\n")
            f.write(f"Amino acid residues (unchanged): {len(residue_groups) - num_corrected}\n")
        else:
            f.write(f"All residues corrected (including amino acids): {num_corrected}\n")
        f.write(f"Charge groups with penalties (N>1): {len(set(e['charge_group'] for e in log_entries))}\n")

    print(f"Wrote log file to: {log_file}")

    return all_conformers, ph_values, residue_groups


def generate_sum_crg(conformers: List[Conformer], ph_values: List[float],
                     residue_groups: Dict, output_file: str, original_sum_crg: str = 'sum_crg.out'):
    """Generate sum_crg.out file from corrected conformer probabilities"""

    print(f"\nGenerating {output_file}...")

    # Read original Protons and Electrons rows if file exists
    original_protons = None
    original_electrons = None

    import os
    if os.path.exists(original_sum_crg):
        with open(original_sum_crg, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('Protons'):
                    original_protons = line.strip()
                elif line.startswith('Electrons'):
                    original_electrons = line.strip()

    # Calculate net charges for each residue grouped by charge sign
    residue_charges = {}

    for res_id, res_conformers in residue_groups.items():
        res_type = res_conformers[0].res_type
        chain = res_conformers[0].chain
        res_num = res_conformers[0].res_num

        # Group conformers by charge sign (positive or negative only)
        charge_sign_groups = {'+': [], '-': []}

        for conf in res_conformers:
            charge_val = parse_charge_state(conf.charge_state)
            if charge_val > 0:
                charge_sign_groups['+'].append((conf, charge_val))
            elif charge_val < 0:
                charge_sign_groups['-'].append((conf, charge_val))
            # Skip neutral (0) charge states

        # For each charge sign group that has conformers, create a line
        for sign in ['+', '-']:
            conf_list = charge_sign_groups[sign]
            if not conf_list:
                continue

            # Create identifier: RESTYPE + SIGN + CHAIN + RESNUM + '_'
            identifier = f"{res_type}{sign}{chain}{res_num}_"

            # Calculate charge contributions at each pH
            charges = []
            for ph_idx in range(len(ph_values)):
                total_charge = sum(charge_val * conf.probabilities[ph_idx]
                                  for conf, charge_val in conf_list)
                charges.append(total_charge)

            residue_charges[identifier] = charges

    # Sort residue identifiers by chain and residue number
    def sort_key(identifier):
        chain = identifier[-6]
        res_num = identifier[-5:-1]
        return (chain, int(res_num))

    sorted_residues = sorted(residue_charges.keys(), key=sort_key)

    # Calculate net charges
    net_charges = [0.0] * len(ph_values)

    for res_id in sorted_residues:
        for ph_idx in range(len(ph_values)):
            net_charges[ph_idx] += residue_charges[res_id][ph_idx]

    # Write output file
    with open(output_file, 'w') as f:
        # Write header
        f.write(' ph            ')
        f.write(''.join(f'{ph:>5.1f} ' for ph in ph_values))
        f.write('\n')

        # Write residue data
        for res_id in sorted_residues:
            f.write(f'{res_id:<15}')
            for charge in residue_charges[res_id]:
                f.write(f'{charge:5.2f} ')
            f.write('\n')

        # Write separator
        f.write('-' * (15 + len(ph_values) * 6) + '\n')

        # Write Net_Charge row
        f.write(f'{"Net_Charge":<15}')
        for charge in net_charges:
            f.write(f'{charge:5.2f} ')
        f.write('\n')

        # Write original Protons and Electrons rows
        if original_protons:
            f.write(original_protons + '\n')
        else:
            print("Warning: Could not find Protons row in original sum_crg.out")

        if original_electrons:
            f.write(original_electrons + '\n')
        else:
            print("Warning: Could not find Electrons row in original sum_crg.out")

    print(f"Wrote sum_crg file to: {output_file}")


if __name__ == '__main__':
    import sys
    import os

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Apply Boltzmann entropy corrections to fort.38 file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
This script will:
  - Read the fort.38 file
  - Apply Boltzmann entropy corrections to non-amino acids (by default)
  - Penalize charge states with more conformers
  - Generate xts_fort.38 (corrected file)
  - Generate xts_sum_crg.out (corrected charge summary)
  - Generate entropy_correction.log (detailed log)

Examples:
  python xts_HETATMcorr.py                    # Process only non-amino acids
  python xts_HETATMcorr.py --doAAs          # Process all residues including amino acids
  python xts_HETATMcorr.py input.txt --doAAs
        '''
    )
    
    parser.add_argument(
        'input_file',
        nargs='?',
        default='fort.38',
        help='Path to fort.38 file (default: fort.38)'
    )
    
    parser.add_argument(
        '--doAAs',
        action='store_true',
        help='Also apply corrections to amino acids (default: only non-amino acids)'
    )
    
    args = parser.parse_args()
    input_file = args.input_file
    
    # Set AMINO_ACIDS based on flag
    if args.doAAs:
        AMINO_ACIDS = {}  # Empty set means no residues are excluded
        print("Mode: Correcting ALL residues (including amino acids)")
    else:
        AMINO_ACIDS = AMINO_ACIDS_LIST  # Use the full list to exclude amino acids
        print("Mode: Correcting non-amino acids only")

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found!")
        print("\nPlease ensure fort.38 exists in the current directory,")
        print("or specify the correct path as an argument.")
        sys.exit(1)

    output_file = 'xts_fort.38'
    sum_crg_file = 'xts_sum_crg.out'
    log_file = 'entropy_correction.log'

    print("=" * 60)
    print("BOLTZMANN ENTROPY CORRECTION FOR fort.38")
    print("=" * 60)
    print(f"Input file: {input_file}")
    print("\nMethod: Penalizing charge states with more conformers")
    print("Penalty: ΔG += RT × ln(N), where N = number of conformers")
    print()

    try:
        conformers, ph_values, residue_groups = process_fort38(input_file, output_file, log_file, AMINO_ACIDS)
        generate_sum_crg(conformers, ph_values, residue_groups, sum_crg_file)

        print("\n" + "=" * 60)
        print("COMPLETE!")
        print("=" * 60)
        print(f"\nGenerated files:")
        print(f"  - {output_file} (corrected probabilities)")
        print(f"  - {sum_crg_file} (corrected charge summary)")
        print(f"  - {log_file} (detailed correction log)")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
