multi.tempted.sim
================

Simulations comparing **multiTEMPTED** with **MEFISTO** (MOFA2) on
multi-omic longitudinal data. Each script is self-contained, runs
top-to-bottom, prints a summary, and writes its figures (and, where the
fits are slow, its results) into `output/`.

## What is being compared

Both methods decompose several subject × feature × time tensors into a
few interpretable components. They differ in which kind of variation
they allow, and the difference is a mirror image:

- **multiTEMPTED** gives every modality its own smooth temporal curve,
  but within a component all subjects share that curve, scaled by a
  subject loading. *Omics may follow different time courses; subjects
  share one.*
- **MEFISTO** gives every subject its own Gaussian-process trajectory,
  but a factor carries a single trajectory shared across modalities.
  *Subjects may follow different time courses; omics share one.*

The synthetic experiments target this asymmetry from both sides: first
where the truth is modality-specific (hard for MEFISTO), then where it
is shared across modalities (MEFISTO’s own ground). Real data comes
last.

## Scripts

| Script | Setting | What it shows | Figure |
|----|----|----|----|
| `01_validate.R` | multiTEMPTED only; 5 noise levels × 20 seeds | recovers planted components almost exactly at low noise, degrading smoothly as noise grows | [PDF](output/01_validate.pdf) |
| `02_compare_unaligned.R` | modality-specific curves, unaligned sampling times, 10 seeds | MEFISTO cannot express curves that differ across modalities, and needs far more compute | [PDF](output/02_compare_unaligned.pdf) |
| `03_compare_aligned.R` | shared curves, block-structured subjects, 100 seeds | even where MEFISTO can represent the truth, multiTEMPTED recovers groups and curves better, and far faster | [PDF](output/03_compare_aligned.pdf) |
| `04_compare_ipop.R` | real iPOP omics (cytokine, metabolome, lipid, protein) | multiTEMPTED’s subject embedding separates by sex more sharply | [PDF](output/04_compare_ipop.pdf) |

## Results

Recovery is scored as `|cor|` between an estimate and the planted truth;
misclassification is `1 -` group-assignment accuracy.

|  | multiTEMPTED | MEFISTO |
|----|----|----|
| **01** subject / feature / temporal recovery, noise 0.1 | 0.999 / 0.999 / 0.996 | n/a |
| **01** same, noise 4.0 | 0.281 / 0.262 / 0.431 | n/a |
| **02** temporal curve recovery | 0.990 | 0.671 |
| **03** subject misclassification | 0.000 | 0.472 |
| **03** temporal curve recovery | 0.971 | 0.738 |
| **04** max \|cor(score, sex)\| | 0.53 | 0.16 |

Across every comparison multiTEMPTED also fits in a small fraction of
the time MEFISTO takes: a gap of several orders of magnitude, widening
as the number of distinct sampling times grows. `03`’s figure shows this
directly as a per-seed compute-time panel alongside the two accuracy
panels.

## Run

``` r
source("01_validate.R")          # multiTEMPTED only; quick
source("02_compare_unaligned.R") # MEFISTO on unaligned data; slow
source("03_compare_aligned.R")   # long, resumable
source("04_compare_ipop.R")      # real data; MEFISTO is the slow step
```

Each script writes to an `output/` folder beside the script itself, so
results land in the project regardless of the session’s working
directory, and prints that destination on start-up.

## Reusing a finished run

The MEFISTO fits in `02`–`04` are slow, so each of those scripts saves
its results and reuses them on a re-run. Re-running then skips the fits
and goes straight to the figure, which makes re-plotting or restyling
free.

| Script | Saved to | When |
|----|----|----|
| `02_compare_unaligned.R` | `02_sims.csv`, `02_sims.rds` | once, after the last seed |
| `03_compare_aligned.R` | `03_sims.csv`, `03_sims.rds` | every `CKPT_EVERY` seeds |
| `04_compare_ipop.R` | `04_ipop.csv`, `04_ipop.rds` | once, after the run |

`03` is additionally **resumable mid-run**: it picks up at the first
missing seed, so an interruption costs at most `CKPT_EVERY` seeds.
Because every seed re-seeds the RNG from its own index, a resumed run
reproduces an uninterrupted one exactly. `02` and `04` save only on
completion, so an interruption there costs the whole run.

In each case the CSV is the readable artifact and the RDS carries the
exact values used for plotting. Delete them to force a fresh run; each
script refuses to reuse results saved under different settings,
reporting both configurations.

## Settings

`MEF_MAXITER`, `MEF_CONVERGENCE` and `MEF_OPTIMISE_GP` (top of 02–04)
trade speed for fairness: fewer iterations with `"fast"` and GP
optimisation off for quick checks; more iterations with `"slow"` and GP
optimisation on for the real comparison. `N_SEEDS` sets the scale.
multiTEMPTED’s cost is unaffected by any of these.

In `03`, `SAMPLING` chooses how sampling times are drawn: `"ALIGNED"`
(one shared grid), `"CLUSTERED"` (a cloud around nominal visit times,
with `CLUSTER_BY` deciding whether a subject’s modalities stay
co-measured), or `"UNALIGNED"` (independent random times). This
dominates MEFISTO’s runtime more than any other knob, because its GP
step additionally fits a group × group covariance whenever subjects
share sampling times; `"ALIGNED"` is therefore much the slowest and
impractical at full scale.

Plot aesthetics (theme, palette, point size) are set in a short block at
the top of each script rather than inside the package’s plotting
functions, so they can be changed without touching `multi.tempted`.

## Requirements

The `multi.tempted` and `MOFA2` packages, plus `ggplot2` and `patchwork`
for the figures. `BiocManager::install("MOFA2")` also pulls `basilisk`,
which supplies the Python backend; no manual setup needed.
