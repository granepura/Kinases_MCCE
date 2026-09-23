# Fig4A (per-residue charge, holo vs apo): how the statistics were computed

Written: 2026-09-23 07:06:58
Script : plot_trials_comparison_xts_Fig4A.py

## Replicates

n = 3 independent trials of the full MCCE pipeline (steps 1-4).

| Trial | MONTE_SEED |
|---|---|
| Trial01 | 1001 |
| Trial02 | 1002 |
| Trial03 | 1003 |

Both of the pipeline's stochastic stages differ between trials:

  * step2 rotamer generation does not reproduce between runs, so each
    trial has a different conformer set (0 of 37 holo structures match).
  * step4 Monte Carlo uses an explicit per-trial MONTE_SEED.

A trial is therefore the correct replicate unit; re-running step4 alone
in a finished tree would hold the conformers fixed and understate the
spread.  (The exception is a residue whose conformers happen to be
identical across trials -- there the split is pure Monte Carlo.)

## Statistics

| Quantity | Definition |
|---|---|
| mean | arithmetic mean over the n trials |
| sd | sample standard deviation, n-1 denominator |
| sem | sd / sqrt(n) -- **the plotted error bar** (--err sem) |
| min, max | the extreme trial values |
| delta | holo - apo computed WITHIN each trial, then averaged |

The difference is paired inside a trial because apo reuses holo's
step2 conformers there, so the pairing cancels the conformer set.
Combining two independent SEMs instead would overstate the error.

Error bars are +/- 1 SEM, NOT a confidence interval.  With n = 3 there
are 2 degrees of freedom, so a 95% CI would be mean +/- 4.30 x SEM.

## Is the mean a fair summary?

| Category | N | % |
|---|---|---|
| identical | 3179 | 86.7 |
| scattered | 488 | 13.3 |
| two-state | 1 | 0.0 |

Total points: 3668

  * identical -- every trial agreed exactly; SEM is 0 and no bar is
    drawn.  That is agreement, not a missing error bar.
  * scattered -- unimodal; mean +/- SEM is appropriate.
  * two-state -- the trials fall into two discrete states, so the mean
    is a value the simulation never produced.  Report the state
    populations instead; see the companion .tsv.
