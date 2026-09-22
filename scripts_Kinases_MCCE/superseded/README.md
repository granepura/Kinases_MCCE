# Superseded run-setup scripts

Kept for reference only. They describe earlier layouts of the apo workflow and
will not work against the current `trial_config.sh` (the variables they source
are gone). Do not run them.

| file | replaced by | why |
|---|---|---|
| `derive_apo.sh` | `prepare_run_apo.py` (copied into each trial as `1-prepare_run_apo.py`) | same seeding in Python, no rsync dependency, lives inside the trial so the trial records the version it ran, and it links `step2_out.pdb` to `apo_step2_out.pdb` |
| `rm_inhib_step2_out.sh` | `make_holo_apo_step2_out.py` | the strip moved into holo's stepB, which writes both `holo_step2_out.pdb` and `apo_step2_out.pdb` once, instead of editing `step2_out.pdb` in place in run_apo |
| `rm_inhib_step2_out.py` | `make_holo_apo_step2_out.py` + `install_apo_step2_out.py` | same reason; run_apo no longer derives anything, it installs and checks what holo made |
| `backup_holo_step2_out.py` | `make_holo_apo_step2_out.py` | it only made the holo copy; the replacement makes both variants in one pass |

The one thing worth reading here is the column 18–20 residue-name matching in
`rm_inhib_step2_out.sh`, which the current scripts reproduce byte-for-byte.
