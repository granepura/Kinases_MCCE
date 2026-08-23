# Kinases_MCCE

Data, simulation trees, and analysis scripts for the study
**"When do kinase inhibitors change protonation or tautomeric state upon binding?"**

Ranepura, Chowdhury, Rosenzweig, Rustenburg, López-Ríos de Castro, Mao, Chodera, Singh, and Gunner.

---

## Overview

The protonation and tautomeric states of a ligand and its target are first-order determinants of
binding affinity, yet most structure-based modeling workflows assign them heuristically and hold
them fixed. This project uses **MCCE4** (Multi-Conformer Continuum Electrostatics) to compute the
full Boltzmann distribution of protonation and tautomer states — for both partners simultaneously —
across **37 co-crystal structures spanning 18 FDA-approved inhibitors and 9 kinase domains**.

For every complex, three matched calculations are compared at pH 7.4:

| State | What is titrated |
|---|---|
| **Holo** | kinase + inhibitor (the parent calculation) |
| **Apo** | the same protein coordinates with the ligand removed |
| **Solution** | the inhibitor alone in continuum solvent |

The apo and solution runs are derived from the holo run *after* rotamer generation, so side-chain
conformers are held identical and the only difference between holo and apo is the presence of the
ligand. Inhibitor proton affinities and tautomer energies are taken from Schrödinger's Epik;
electrostatics are solved with DelPhi (ε = 4 protein / 80 solvent, 0.15 M salt, crystallographic
waters removed).

---

## Key results

**Inhibitor charge is dynamic and binding-site specific.**
Ensemble-average inhibitor charges in solution span 0 to +1.6 at pH 7.4; none are negative.
Inhibitors that are weakly charged in solution (q < 0.3) tend to *gain* positive charge on binding,
while those already substantially protonated (q ≥ 0.9) largely retain it. The correlation between
solution and bound charge is weak (r = 0.55, ρ = 0.54), reflecting the diversity of binding-site
electrostatic environments. Desolvation does not simply drive ligands toward neutrality — several
pockets stabilize the +1 form instead.

**The protein's net charge is buffered.**
Across all 37 complexes, the holo − apo protein charge difference is < |0.2| except for three EGFR
structures (two afatinib complexes lose 0.5–0.75 protons from Asp837; osimertinib-bound 4ZAU loses
1.05 protons across a network including Asp837, Glu711/736, and His805/893). This stability is not
the result of compensating shifts at different sites: the per-residue distributions themselves are
largely unchanged (N = 3508 titratable residues), with responses confined to a handful of residues.

**Responding residues are few, mostly Lys and His, and not always near the pocket.**
Strong responders (|Δq| ≥ 0.5) and affected residues (0.2 ≤ |Δq| < 0.5) sit between ~6 and ~16 Å
from the ligand center of mass. The largest single shift, Lys1205 in 4MKC (ceritinib), is ~10 Å
away — inhibitor binding can perturb protonation well beyond the first shell.

**Tautomer populations can invert on binding.**
Four inhibitors (axitinib, bosutinib, imatinib, ponatinib) have two tautomers within 2.8 kcal/mol in
solution, and this is where the largest redistributions occur: for imatinib bound to ABL and
ponatinib bound to DDR1, the lowest-energy *solution* tautomer is not the one populated when bound.
Binding shifts the relative probabilities of low-energy conformers (ΔG < 2.5 kcal/mol) without
recruiting higher-energy ones; typical ensemble free-energy changes are < 1 kcal/mol, but tautomers
that are degenerate in solution can differ by 3–6 kcal/mol when bound.

**DFG conformation is not predictive.**
Splitting the dataset by DFG-in vs. DFG-out label (bootstrapped ECDFs, 2000 resamples) shows no
strong conformational dependence for net-charge change in either partner. DFG-out complexes show
almost no inhibitor charge change; the broader spread among DFG-in complexes tracks the greater
structural diversity of that class. The DFG residues themselves never change protonation state
measurably in any structure studied.

**Practical implication.** Protonation and tautomer assignment should be treated as a
binding-site-specific property of each kinase–inhibitor complex, not inferred from solution-state
ligand energetics, kinase conformational label, or a single fixed-protonation model. Because the
protein response is local, dynamically titrating a limited binding-site shell may capture most of it
without titrating the whole protein.

---

## Repository layout

```
final_scripts/           curated analysis and run-setup scripts
  xts_corr.py            conformer-count entropy correction (see below)
  create_cof2.sh         builds the ligand-in-solution tree from the holo tree
  rm_cofs_step2_out.sh   builds the apo tree (ligand stripped after rotamer generation)
  rm_cofs.sh             older apo variant (ligand stripped before rotamer generation)
  plot_sumcrg_inhibitors_xts_Fig3.py    Figure 3: inhibitor charge, bound vs. solution
  plot_sumcrg_comparison_xts_Fig4.py    Figure 4: per-residue charge, holo vs. apo
  pdb_inhibitor.lst      PDB -> inhibitor -> ligand-code map (37 rows)

outlier-analysis/        Figure 4 structural overlay (PyMOL) and the outlier table
dfg-split-analysis/      Figure 5 DFG-in/DFG-out ECDFs with bootstrap CIs
test_xts_corr/           development lineage and test cases for the entropy correction
cof_tpl_mg/              snapshot of the MCCE run tree
  run_kin/               holo complexes, one directory per PDB ID
  run_prot2/             apo proteins (matched conformers; used for Fig. 4)
  run_cof2/              inhibitors alone in solution
  run_prot/              older apo variant, superseded
kinase_project-final-tables.xlsx   manuscript Table 1 and SI tables
```

---

## Reproducing the calculations

MCCE4 is run per structure from inside a run directory, via SLURM or directly:

```bash
cd cof_tpl_mg/run_prot2/2HYY
sbatch ../submit_mcce4_delphi.sh      # or: bash ../submit_mcce4_delphi.sh
```

Production settings:

```
step1.py -d 4 --noter --dry
step2.py -d 4 -l 1
step3.py -d 4 -s delphi -salt 0.15 --fly -p $CPUS -t $TMP
step4.py --xts --ms -i 7.4 -n 1
```

The derived trees (`run_cof2`, `run_prot2`, `run_prot`) inherit steps 1–2 from `run_kin` and are
submitted with `step1="f" step2="f"`. Re-running step 1 or 2 in a derived tree breaks the
shared-conformer premise that makes the holo/apo comparison meaningful.

### Entropy correction

`xts_corr.py` post-processes `fort.38` to remove the bias by which a charge state with more
conformers appears more probable purely from conformer count: it groups conformers by charge,
computes a Shannon entropy per group, adds it to each conformer's relative free energy, and
re-Boltzmanns.

```bash
cd <run_dir>/<PDB>
python3 xts_corr.py            # non-amino-acid residues only (as used in the paper)
```

It writes `xts_fort.38` and `xts_sum_crg.out` alongside the uncorrected files. **All published
figures read `xts_sum_crg.out`**; scripts without `_xts` in the name read the uncorrected output and
are retained only for reference.

### Analysis scripts

Run with `python3 <script>`. There are no command-line arguments — configuration is by editing the
`dir1` / `dir2` / `lst_file` paths at the top of each file, which must be pointed at your local copy
of the run tree before use. Requires `numpy`, `matplotlib`, `pandas`, and `seaborn`; the Figure 4
structural overlay additionally requires PyMOL.

Figure 4 classifies residues in two tiers by |q_holo − q_apo| (strong ≥ 0.5, affected 0.2–0.5);
Figure 3 flags outliers at 15% deviation from y = x; the ECDF bootstrap uses 2000 resamples.
The `_28pdb` script variants restrict the analysis to the 27-PDB subset matching the Excel tables
rather than the full 37 — check which N a reported number came from.

---

## Related resources

- **MCCE4-Alpha** (software and topology files): https://github.com/GunnerLab/MCCE4-Alpha
- **Epik inhibitor charges and parameters**: https://github.com/choderalab/mcce-charges/tree/master/epik_inhibitors

## Funding

NSF MCB-2141824 (M.R.G., G.R.); Damon Runyon Quantitative Biology Fellowship DRQ-14-22 and NCI K99
CA286801 (S.S.); NIH R35GM152017 and P30CA008748 (J.D.C.).

## Contact

Gehan A. Ranepura — granepura@gc.cuny.edu
