multi.tempted.sim
================

Reproducible simulation and analysis code for the manuscript comparing
**multiTEMPTED** (the
[`multi.tempted`](https://github.com/loulind/multi.tempted) package)
with **MEFISTO** (MOFA2) on multi-omic longitudinal data.

Every number and figure in the paper comes from a script here. The
manuscript carries the full specifications and results; this README is
just a map of the repo and how to run it.

## What is being compared

Both methods decompose several subject x feature x time tensors into a
few interpretable components, and they differ in which kind of variation
they allow:

- **multiTEMPTED** gives every modality its own smooth temporal curve,
  but within a component all subjects share that curve, scaled by a
  subject loading.
- **MEFISTO** gives every subject its own Gaussian-process trajectory,
  but a factor carries a single trajectory shared across modalities.

The experiments approach that asymmetry from both sides: first where the
truth is modality-specific (hard for MEFISTO), then where it is shared
across modalities (MEFISTO’s own ground), and finally on real data.

## Scripts

Each script is self-contained, runs top to bottom, prints a summary, and
writes to a folder beside the script rather than to `getwd()`.

| Script | Setting | Output |
|----|----|----|
| `01a_validate.R` | multiTEMPTED alone on data from the model it assumes, swept over noise | [PDF](output/01a_validate.pdf) |
| `01b_validate_larger.R` | the same validation at realistic feature dimensions and two cohort sizes | [PDF](output/01b_validate_larger.pdf) |
| `02_compare_unaligned.R` | modality-specific curves, unaligned sampling times | [PDF](output/02_compare_unaligned.pdf) |
| `03_compare_aligned.R` | shared curves, block-structured subjects | [PDF](output/03_compare_aligned.pdf) |
| `04_compare_ipop.R` | real iPOP omics (cytokine, metabolome, lipid, protein) | [PDF](output/04_compare_ipop.pdf) |
| `supp_sensitivity.R` | iPOP sensitivity to the modality weights and the smoothing penalty | [PDF](supp_output/supp_fig_weights.pdf), [table](supp_output/supp_tab_smoothing.tex) |

Scripts `01a` through `04` write to `output/`; `supp_sensitivity.R`
writes to `supp_output/`, including a LaTeX table ready to compile or to
drop into the supplement.

## Run

``` r
source("01a_validate.R")          # quick
source("01b_validate_larger.R")   # long, resumable
source("02_compare_unaligned.R")  # MEFISTO on unaligned data; slow
source("03_compare_aligned.R")    # long, resumable
source("04_compare_ipop.R")       # real data; MEFISTO is the slow step
source("supp_sensitivity.R")      # refits iPOP once per setting
```

Each script prints its output folder on start-up, so results land in the
project whatever the session’s working directory happens to be.

## Reusing a finished run

The slow scripts save their results and reuse them on a re-run, which
makes re-plotting and restyling free:

- `01b` and `03` checkpoint as they go and resume at the first missing
  run or seed, so an interruption costs at most `CKPT_EVERY` of them.
  Because every seed re-seeds the RNG from its own index, a resumed run
  reproduces an uninterrupted one exactly.
- `02` and `04` save once on completion; an interruption there costs the
  run.
- `01a` and `supp_sensitivity.R` refit from scratch every time.

In each case the CSV is the readable artifact and the RDS carries the
exact values used for plotting. Delete them to force a fresh run; a
script refuses to reuse results saved under different settings, and
reports both configurations.

## Settings

Knobs sit in a labelled block at the top of each script, above the
analysis. `MEF_MAXITER`, `MEF_CONVERGENCE` and `MEF_OPTIMISE_GP` trade
MEFISTO’s speed against fairness: loose settings for a quick check,
strict ones with GP optimisation on for the real comparison
(multiTEMPTED’s cost is unaffected by any of them). `N_SEEDS` sets the
scale, and in `03` the `SAMPLING` knob chooses how sampling times are
drawn, which dominates MEFISTO’s runtime more than anything else.

Plot aesthetics are set in a short block at the top of each script
rather than inside the package’s plotting functions, so they can be
changed without touching `multi.tempted`.

## Requirements

The `multi.tempted` and `MOFA2` packages, plus `ggplot2` and `patchwork`
for the figures. `BiocManager::install("MOFA2")` also pulls `basilisk`,
which supplies the Python backend; no manual setup needed.
