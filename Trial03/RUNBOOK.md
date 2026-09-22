# Trial03 runbook

Created: 2026-09-22 06:21:43
MONTE_SEED: 1003 (explicit)

## Shared scripts (sha256 at setup time -- re-check before comparing trials)
```
d6b540a22d945744fb738e81e9123a1eccde634e9117b16e9d5527a7cd5b21e9  /data/home/granepura/5-Kinases/Kinases_MCCE/Trial03/1-prepare_run_apo.py
a5375382a8c4a95b0ad22198ab5795d2a041ae524f6df9b0d67101feeacf5297  /data/home/granepura/5-Kinases/Kinases_MCCE/scripts_Kinases_MCCE/make_holo_apo_step2_out.py
0f8759a3df86fb8dab62fd89240706eb6d4ae70ea978555dd2c509ce1af3f684  /data/home/granepura/5-Kinases/Kinases_MCCE/scripts_Kinases_MCCE/install_apo_step2_out.py
3fe26e3b6ae1414d71702e2f4266fdca95285ac19fbdea76c018a56b4d91dfa1  /data/home/granepura/5-Kinases/Kinases_MCCE/scripts_Kinases_MCCE/prune_kin-inhib_head3.py
90bfc2c015ab74f4f51672ed0abe537c8ab016692565b1e2d06d9fdf655dc068  /home/granepura/5-Kinases/Kinases_MCCE/scripts_Kinases_MCCE/trial_config.sh
```

## Order of operations

1. holo, steps 1-2 -- builds the conformers and the coordinate frame everything else inherits
       cd run_holo
       pro_batch kin-pdb -custom submit_mcce4_s1s2.sh -job-name holo_s1s2 -j 15
   stepB here is make_holo_apo_step2_out.py: step3/step4 are off, so it runs last
   and splits the finished step2_out.pdb into holo_step2_out.pdb (exact copy) and
   apo_step2_out.pdb (inhibitor deleted).  Check each structure's stepB.log.

2. Seed apo from holo (needs only holo's step2_out.pdb, so it can run as soon as
   holo's steps 1-2 finish -- no need to wait for holo's step3/4).  Copies each
   whole run_holo/<PDBID> to run_apo/<PDBID> -- every file steps 1-2 left, plus
   both step2 variants -- replacing any that is already there, and skipping
   holo's step2_out.pdb and any step3/4 products.  It then links
   step2_out.pdb -> apo_step2_out.pdb, resetting the link if one exists.
       cd /data/home/granepura/5-Kinases/Kinases_MCCE/Trial03
       ./1-prepare_run_apo.py
       ./1-prepare_run_apo.py --dry-run   # inspect without writing
       ./1-prepare_run_apo.py 1XKK 2ITZ   # re-seed just these

3. holo, steps 3-4
       cd run_holo
       pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name holo_s3s4 -j 15

4. apo, steps 3-4
       cd run_apo
       pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name apo_s3s4 -j 15

5. inhib, steps 1-4 (independent: starts from cof-pdb, not carved from holo)
       cd run_inhib
       pro_batch cof-pdb -custom submit_mcce4.sh -job-name inhib -j 15

Steps 3 and 5 are independent of each other and of step 4; run them concurrently.
Check progress with:  pro_batch --check -job-name <name>

## The step2 chain

The ligand is deleted once, in the holo job, and both trees then share the same
two files.  Only step2_out.pdb differs between them:

    run_holo/<ID>/step2_out.pdb          what steps 1-2 produced = holo
         |  holo stepB: make_holo_apo_step2_out.py
         +-> holo_step2_out.pdb          exact copy of it
         +-> apo_step2_out.pdb           same file, inhibitor deleted

    1-prepare_run_apo.py copies the whole directory across, then:

    run_apo/<ID>/holo_step2_out.pdb      copied, for reference and the pair check
    run_apo/<ID>/apo_step2_out.pdb       copied -- the apo structure
    run_apo/<ID>/step2_out.pdb  ->  apo_step2_out.pdb   (relative symlink)

Deleting the ligand is the ONLY difference between the two structures -- step1
and step2 are off in run_apo, so the pocket is never repacked.  The code comes
from pdb_inhibitor.lst keyed on the directory name and is matched on the
residue-name columns 18-20, so other heteroatoms stay (2ITZ keeps its _CL).

apo's stepB is install_apo_step2_out.py.  It rewrites nothing in the normal
case: it checks that apo_step2_out.pdb holds none of this structure's inhibitor,
that it is exactly holo_step2_out.pdb minus those lines, and that step2_out.pdb
points at it -- resetting the link if it is missing, a plain file, or pointing
elsewhere.  That check runs immediately before the ~8 minutes of step3 that
depend on it.

driver_mcce4.sh logs a stepB failure without aborting the run, so when the pair
does not check out the script REMOVES step2_out.pdb.  step3 with step2="f" only
runs when that file exists, so the structure ends without a pK.out and
pro_batch --check flags it, instead of step3 computing something unvouched-for.
Check stepB.log and mcce_timing.log.

## Why stepC runs in all three

prune_kin-inhib_head3.py does two independent edits: the inhibitor conformer
pruning AND forcing ARG positive (neutral ARG -> FL=t).  apo has no inhibitor but
does have ARG, so skipping stepC there would give holo and apo different ARG
treatments and invalidate the comparison.  inhib has no ARG and gets only the
inhibitor edits.  All three therefore run stepC="t".

## Frames

run_holo and run_apo share one coordinate frame: apo reuses holo's step2_out.pdb
and its submit script has step1/step2 off, so nothing re-centers it.
run_inhib is built from cof-pdb through its own step1/step2 and therefore sits in
its own centered frame.  That is fine for pKa/charge -- it is an isolated-ligand
reference state -- but do not compare its coordinates to holo's.
