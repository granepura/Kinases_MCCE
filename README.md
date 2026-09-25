# Kinases_MCCE

Data, simulation trees, and analysis scripts for the study
**"When do kinase inhibitors change protonation or tautomeric state upon binding?"**

Ranepura, Chowdhury, Rosenzweig, Rustenburg, López-Ríos de Castro, Mao, Chodera, Singh, and Gunner.

---

## Overview

The protonation and tautomeric states of a ligand and its target are first-order determinants of
binding affinity, yet most structure-based modeling workflows assign them heuristically and hold
them fixed. This project uses **MCCE4** (Multi-Conformer Continuum Electrostatics) to compute the
full Boltzmann distribution of protonation and tautomer states, for both partners simultaneously,
across **37 co-crystal structures spanning 18 FDA-approved inhibitors and 9 kinase domains**.

For every complex, three matched calculations are compared at pH 7.4:

| State | Directory | What is titrated |
|---|---|---|
| **Holo** | `run_holo` | kinase + inhibitor (the parent calculation) |
| **Apo** | `run_apo` | the same protein coordinates with the ligand removed |
| **Solution** | `run_inhib` | the inhibitor alone in continuum solvent |

Apo is derived from holo *after* rotamer generation, so the side-chain conformers are identical and
the only difference between the two is the presence of the ligand. Inhibitor proton affinities and
tautomer energies come from Schrödinger's Epik; electrostatics are solved with DelPhi
(ε = 4 protein / 80 solvent, 0.15 M salt, crystallographic waters removed).

**Everything is run in triplicate.** `Trial01`, `Trial02` and `Trial03` are independent repeats of
the whole pipeline. This matters because two stages are stochastic: step 2's rotamer generation does
not reproduce between runs (no structure reproduces its conformer count across all three trials), and
step 4's Monte Carlo uses a per-trial seed. The trials therefore measure run-to-run reproducibility
of the entire calculation, not just the Monte Carlo.

---

## Key results

**Inhibitor charge is dynamic and binding-site specific.**
Ensemble-average inhibitor charges in solution span 0 to +1.6 at pH 7.4; none are negative.
Inhibitors that are weakly charged in solution (q < 0.3) tend to *gain* positive charge on binding,
while those already substantially protonated (q ≥ 0.9) largely retain it. The correlation between
solution and bound charge is weak, reflecting the diversity of binding-site electrostatic
environments. Desolvation does not simply drive ligands toward neutrality — several pockets
stabilize the +1 form instead.

**The protein's net charge is buffered.**
Across all 37 complexes the holo/apo charge difference is small, with the exception of a few EGFR
structures. This is not the result of compensating shifts at different sites: the per-residue
distributions themselves are largely unchanged (3668 titratable residues compared), with responses
confined to a handful of residues.

**Responding residues are few, mostly Lys and His, and not always near the pocket.**
The largest single shift, Lys1205 in 4MKC (ceritinib), is ~10 Å from the ligand — inhibitor binding
can perturb protonation well beyond the first shell. It reproduces to **SEM 0.000** across all three
trials.

**Tautomer populations can invert on binding.** Several inhibitors have two tautomers within a few
kcal/mol in solution, and this is where the largest redistributions occur.

**DFG conformation is not predictive.** Splitting by DFG-in vs DFG-out shows no strong
conformational dependence for net-charge change in either partner.

### What the replicates added

Reproducibility is predictable from the occupancies themselves, with no appeal to a pKa. The
quantity being averaged is a state occupancy, and for a site divided between two states with
occupancy `f` the variance is `f(1−f)`, largest when the two are evenly populated and vanishing
when either saturates:

| bound population at pH 7.4 | n | mean spread over 3 trials |
|---|---|---|
| most evenly divided, f(1−f) > 0.15 | 10 | 0.086 |
| intermediate | 4 | 0.040 |
| saturated in one state | 24 | **0.001** |

Two thirds of the ligands reproduce to ±0.01 across independent conformer sets. A large spread means
the population is genuinely divided, not that the calculation failed. 3CS9 nilotinib is the clearest
case, with a bound charge of 0.42, 0.63 and 0.19 in the three trials. Such sites should be reported
as a range rather than mean ± SEM.

Note that these runs use `TITR_STEPS = 1` at pH 7.4, so no titration is performed and no pKa is
determined by them; `pK.out` reports only `<7.4` or `>7.4` for every residue. An effective pKa can be
back-computed from the occupancy as `pKa = 7.4 - log10((1-f)/f)`, but that assumes the site has
exactly two charge states and is not a titration result, so label it as such wherever it is used.
See `CLAUDE.md` for the full treatment.

---

## Repository layout

```
kin-pdb/                  37 holo PDBs      | shared inputs, copied or symlinked
cof-pdb/                  37 ligand-only PDBs |   into every trial
pdb_inhibitor.lst         PDB -> inhibitor -> ligand code -> kinase (37 rows)

Trial01/  Trial02/  Trial03/          one independent repeat each
  run_holo/ run_apo/ run_inhib/       MCCE run trees, one directory per PDB
  0-prepare_run_apo.py                seed run_apo from run_holo
  1-run_xts_corr.py                   entropy correction, all three trees
  2-plot_sumcrg_inhibitors_xts_Fig3.py   Fig 3: inhibitor charge, bound vs solution
  3-plot_sumcrg_comparison_xts_Fig4A.py  Fig 4A: per-residue charge, holo vs apo
  RUNBOOK.md                          the order of operations + script fingerprints

scripts_Kinases_MCCE/                 canonical pipeline scripts
  setup_trial.sh  trial_config.sh     scaffold a trial; shared configuration
  submit_mcce4_template.sh            SBATCH/env template the submit scripts derive from
  make_holo_apo_step2_out.py          stepB in holo: split step2_out.pdb into holo/apo variants
  install_apo_step2_out.py            stepB in apo: verify the pair, link step2_out.pdb
  prune_kin-inhib_head3.py            stepC: prune inhibitor conformers, force ARG positive
  prepare_run_apo.py  run_xts_corr.py  plot_*.py    sources of the numbered trial scripts
  superseded/                         earlier versions and the published snapshot
    cof_tpl_mg/                       the published run tree (read-only history)
      kinase_project-final-tables_v0.xlsx   manuscript Table 1 and SI tables, published run

plot_trials_inhibitors_xts_Fig3.py    cross-trial Fig 3, mean ± SEM
plot_trials_comparison_xts_Fig4A.py   cross-trial Fig 4A, mean ± SEM
make_trials_tables_xlsx.py            builds tables_Trials/kinase_project-trials-tables.xlsx

plots_Trials_*/  tables_Trials/       cross-trial figures, CSVs and the workbook
outlier-analysis/  dfg-split-analysis/  test_xts_corr/    figure-specific analyses
```

---

## Reproducing the calculations

Requires MCCE4 (`mcce` on `PATH`), `pro_batch`, and Python with `numpy`, `matplotlib`, `openpyxl`.
Inhibitor topologies (`FMM.ftpl`, `IRE.ftpl`, …) and the conformer energies in `extra.tpl` ship with
MCCE4, so nothing needs staging per run.

### 1. Scaffold a trial

```bash
scripts_Kinases_MCCE/setup_trial.sh 1        # creates Trial01/ (--force to regenerate)
```

This copies the PDB folders, symlinks `pdb_inhibitor.lst`, writes the four numbered scripts, and
generates the submit scripts — which differ only in job name, step flags, `CPUS`, the hook scripts
and `MONTE_SEED` (Trial01 → 1001, Trial02 → 1002, Trial03 → 1003).

### 2. Run the three trees

Each tree is launched with `pro_batch` from inside its own directory. Clear a tree first if you are
restarting it — **this deletes all results in it**:

```bash
rm -rf 1* 2* 3* 4* 5* meta_bench pro_batch_* book.txt
```

**Solution** (shortest; independent of the other two, so a good first check):

```bash
cd Trial01/run_inhib
pro_batch cof-pdb -custom submit_mcce4.sh -job-name T01_inhib --skip-prerun
cat */mcce_timing.log | grep STEP4 | grep Success | wc -l      # 37 when done
```

**Holo, steps 1–2.** stepB here splits the finished `step2_out.pdb` into `holo_step2_out.pdb` and
`apo_step2_out.pdb` (the inhibitor deleted, matched on residue-name columns 18–20):

```bash
cd Trial01/run_holo
pro_batch kin-pdb -custom submit_mcce4_s1s2.sh -job-name T01_holo_s1s2 --skip-prerun
cat */mcce_timing.log | grep STEP2 | grep Success | wc -l      # 37
```

**Seed apo from holo.** Needs only steps 1–2, so it can run while holo's step 3 is still going:

```bash
cd Trial01 && ./0-prepare_run_apo.py           # --dry-run | --keep | 1XKK 2ITZ
```

It copies each `run_holo/<PDBID>` across — minus holo's own `step2_out.pdb` and any step 3/4
products — and links `step2_out.pdb` to `apo_step2_out.pdb`.

**Holo steps 3–4, and apo steps 3–4.** These are independent of each other and can run
concurrently. Apo's stepB re-checks that `apo_step2_out.pdb` is exactly `holo_step2_out.pdb` minus
the inhibitor before step 3 begins:

```bash
cd Trial01/run_holo
pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name T01_holo_s3s4 --skip-prerun
cd ../run_apo
pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name T01_apo_s3s4  --skip-prerun
cat */mcce_timing.log | grep STEP4 | grep Success | wc -l      # 37 each
```

Progress can also be checked with `pro_batch --check -job-name <name>` (`r` pending, `c` complete,
`e` error). Step 3 dominates at roughly 8 min/structure; step 4 takes ~15 s.

### 3. Entropy correction, then the per-trial figures

**step 4 does not apply the entropy correction** — `--xts` only switches on MCCE's internal entropy
term. The correction is a separate pass that removes the bias by which a charge state with more
conformers appears more probable purely from conformer count:

```bash
cd Trial01
./1-run_xts_corr.py                            # writes xts_sum_crg.out in all three trees
./2-plot_sumcrg_inhibitors_xts_Fig3.py         # --title to draw titles on the PNGs
./3-plot_sumcrg_comparison_xts_Fig4A.py
```

Every figure reads `xts_sum_crg.out`, so this must be done in all three trees or corrected numbers
would be compared against uncorrected ones.

### 4. Repeat for Trial02 and Trial03, then aggregate

```bash
./plot_trials_inhibitors_xts_Fig3.py           # --err sem|sd|range
./plot_trials_comparison_xts_Fig4A.py
./make_trials_tables_xlsx.py                   # -> tables_Trials/
```

The cross-trial scripts read the per-trial CSVs and the workbook reads both, so run them in that
order. All three discover `Trial*/` themselves — adding a Trial04 and re-running is all that is
needed to widen the averages, the error bars and the seed table.

### Things that silently corrupt the comparison

- **Never run step 1 or step 2 in `run_apo`.** It inherits holo's rotamers; regenerating them
  destroys the shared-conformer premise. Its submit script has `step1="f" step2="f"`.
- **Re-seed apo whenever holo is re-run.** step 2 is stochastic, so a new holo run has a different
  conformer set, and an apo tree left over from the previous one is silently mismatched.
- **stepC must run in all three trees.** `prune_kin-inhib_head3.py` does two unrelated edits —
  inhibitor conformer pruning *and* forcing ARG positive — so skipping it in apo would give holo and
  apo different ARG treatments.

`CLAUDE.md` documents these and the rest of the pipeline's invariants in detail.

---

## Related resources

- **MCCE4-Alpha** (software and topology files): https://github.com/GunnerLab/MCCE4-Alpha
- **Epik inhibitor charges and parameters**: https://github.com/choderalab/mcce-charges/tree/master/epik_inhibitors

## Funding

NSF MCB-2141824 (M.R.G., G.R.); Damon Runyon Quantitative Biology Fellowship DRQ-14-22 and NCI K99
CA286801 (S.S.); NIH R35GM152017 and P30CA008748 (J.D.C.).

## Contact

Gehan A. Ranepura — granepura@gc.cuny.edu
Sukrit Singh — sukrit.singh@choderalab.org
