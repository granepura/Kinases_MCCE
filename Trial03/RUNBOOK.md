# Trial03 runbook

Written: 2026-09-22 10:16:58
MONTE_SEED: 1003 (explicit)

Every trial is built and run the same way -- only the seed and the job names
differ.  The numbered scripts sit here at the trial root and work the rest out
for themselves; the pro_batch launches are run from inside each run directory.
Invariants and the reasoning behind them are in ../CLAUDE.md.

## Order of operations

Job names match the #SBATCH --job-name in each submit script so a trial's jobs are distinguishable in squeue.
--skip-prerun is used throughout: pro_batch's pre-run check is not needed here,
the PDBs are already curated.

To start a tree from scratch, clear it first (this deletes all results in it):
       rm -rf 1* 2* 3* 4* 5* meta_bench pro_batch_* book.txt

1. inhib, steps 1-4   (shortest; independent of holo/apo, so a good first check)
       cd run_inhib
       pro_batch cof-pdb -custom submit_mcce4.sh -job-name T03_inhib --skip-prerun
       cat */mcce_timing.log | grep STEP4 | wc -l        # 37 when done

2. holo, steps 1-2
       cd run_holo
       pro_batch kin-pdb -custom submit_mcce4_s1s2.sh -job-name T03_holo_s1s2 --skip-prerun
       cat */mcce_timing.log | grep STEP2 | wc -l        # 37 when done
   stepB = make_holo_apo_step2_out.py: splits the finished step2_out.pdb into
   holo_step2_out.pdb (exact copy) and apo_step2_out.pdb (inhibitor deleted).

3. seed run_apo from run_holo
       ./0-prepare_run_apo.py                 # --dry-run | --keep | 1XKK 2ITZ
   Copies each whole run_holo/<PDBID>, skipping holo's step2_out.pdb and any
   step3/4 products, then links step2_out.pdb -> apo_step2_out.pdb.
   Needs only step 2, so it can run while holo's step3/4 is still going.

4. holo, steps 3-4
       cd run_holo
       pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name T03_holo_s3s4 --skip-prerun
       cat */mcce_timing.log | grep STEP4 | wc -l        # 37 when done

5. apo, steps 3-4
       cd run_apo
       pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name T03_apo_s3s4 --skip-prerun
       cat */mcce_timing.log | grep STEP4 | wc -l        # 37 when done
   stepB = install_apo_step2_out.py: re-checks the step2 pair and the link.

   Steps 4 and 5 are independent of each other; run them concurrently.
   Status:  pro_batch --check -job-name <name>     (r pending, c done, e error)

6. entropy correction -- step4 does NOT do this
       ./1-run_xts_corr.py                    # --dry-run | -t run_inhib | --force
   Runs xts_corr.py in all three trees, producing xts_sum_crg.out (plus
   xts_fort.38, entropy_correction.log).  Every figure reads xts_sum_crg.out, so
   this must be done in all three trees or corrected numbers would be compared
   against uncorrected ones.  It prints one line per structure and lists the
   cause of any failure.

7. figures
       ./2-plot_sumcrg_inhibitors_xts_Fig3.py    # inhibitor: bound vs in solution
       ./3-plot_sumcrg_comparison_xts_Fig4A.py   # per residue: holo vs apo
   Add --title to draw titles on the PNGs.  Both write PNGs and a CSV carrying
   the trial name, into plots_Fig3_inhibitors_xts/ and
   plots_Fig4A_holo_vs_apo_xts/.  They only read the runs; nothing is modified.

## Checking one structure

       cat run_holo/1XKK/mcce_timing.log     per-step wall time, success/failure
       cat run_holo/1XKK/stepB.log           the holo/apo step2 split
       cat run_apo/1XKK/stepB.log            the pair check + the step2_out link
       cat run_holo/1XKK/stepC.log           the head3.lst edit
       ls  run_holo/1XKK/pK.out              exists => step4 finished

## Provenance -- sha256 of the scripts this trial was set up with
```
b0bd5af5096f42443e47d6a800cf1e4ca31eb8676e983f522284529bf20053f0  Trial03/0-prepare_run_apo.py
d6b540a22d945744fb738e81e9123a1eccde634e9117b16e9d5527a7cd5b21e9  Trial03/1-prepare_run_apo.py
f84e4ee956556962d3871dcc75fc4c2801fb49a5eff24c488c3529f8efec4211  Trial03/1-run_xts_corr.py
7f68dae8a59b8ad6d5f2bde0de08fb8579bff201583a5c76c07589820a6ec756  Trial03/2-plot_sumcrg_inhibitors_xts_Fig3.py
39de2252d15ceb296c82a74d262c38aa85c9deb3e8e38c0558428ed27bdb8bf5  Trial03/3-plot_sumcrg_comparison_xts_Fig4A.py
9c194507b76496d0ce1fb459c58c6a472eb57c77d91037bdece8e8c22a0a7221  scripts_Kinases_MCCE/make_holo_apo_step2_out.py
aaa5edb2e42dd35390d57a0abd7f3d1919a2e8d2c472e6164c574083a2461667  scripts_Kinases_MCCE/install_apo_step2_out.py
ea0b66ebc7e127eeca4c1ba2b26787d527e98e903cef1d4fd290dd8ebc1a9489  scripts_Kinases_MCCE/prune_kin-inhib_head3.py
2078a22d08fed503785d3e384a6981105ac4bb87cef83143ea473a2869eca4fa  scripts_Kinases_MCCE/trial_config.sh
```

Re-check these before comparing trials: the canonical scripts in
scripts_Kinases_MCCE/ change over time, and a trial's copies are the
record of what it actually ran.
