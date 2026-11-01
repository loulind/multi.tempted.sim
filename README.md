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
across modalities (MEFISTO’s own ground), and finally on real data. The
MOMS-PI pregnancy cohort is analysed with multiTEMPTED alone, jointly on
the vaginal microbiome and cervicovaginal cytokines.

## Scripts

Each script is self-contained, runs top to bottom, prints a summary, and
writes to a folder beside the script rather than to `getwd()`.

Scripts are numbered by the manuscript figure they produce; the `a`
scripts build that figure’s input data, and the `supp_` scripts produce
the supplementary figures and table.

| Script | Setting | Output |
|----|----|----|
| `01a_momspi_preprocess.R` | builds the MOMS-PI microbiome and cytokine data | `data/momspi/processed/` |
| `01_momspi.R` | MOMS-PI case study: joint vs microbiome-only fits (Figure 1) | [subjects](output/01_momspi_subject_loading.pdf), [features](output/01_momspi_feature_loading.pdf), [cytokines](output/01_momspi_cytokine_trajectory.pdf) |
| `02a_ipop_preprocess.R` | builds the iPOP data (`multi.tempted::ipop`) from the public download | `data/iPOP/processed/ipop.rds` |
| `02_ipop.R` | real iPOP omics (cytokine, metabolome, lipid, protein), multiTEMPTED vs MEFISTO (Figure 2) | [PDF](output/02_ipop.pdf) |
| `03_validate.R` | multiTEMPTED alone on data from the model it assumes, at realistic feature dimensions and two cohort sizes, swept over noise (Figure 3) | [PDF](output/03_validate.pdf) |
| `04_compare_unaligned.R` | modality-specific curves, unaligned sampling times (Figure 4) | [PDF](output/04_compare_unaligned.pdf) |
| `05_compare_aligned.R` | shared curves, block-structured subjects (Figure 5) | [PDF](output/05_compare_aligned.pdf) |
| `supp_momspi_stability.R` | MOMS-PI cytokine loading stability over 100 resamples of 80% of subjects (Figure S1) | [PDF](supp_output/supp_fig_momspi_stability.pdf) |
| `supp_sensitivity.R` | iPOP sensitivity to the modality weights (Figure S2) and the smoothing penalty (Table S1) | [PDF](supp_output/supp_fig_weights.pdf), [table](supp_output/supp_tab_smoothing.tex) |
| `supp_ablation.R` | multiTEMPTED vs its shared-trajectory restriction on Figure 5’s data over a noise sweep, no MEFISTO (Figure S3) | [PDF](supp_output/supp_fig_ablation.pdf) |

The numbered scripts write to `output/` (the `a` scripts write to
`data/`); the `supp_` scripts write to `supp_output/`, including a LaTeX
table ready to compile or to drop into the supplement.

## Data

The simulations generate their own data. The real-data inputs are
publicly available but are not tracked in this repository;
`data/momspi/` and `data/iPOP/` are ignored by git.

### MOMS-PI

Place these six files in `data/momspi/`:

| File | Contents | Source |
|----|----|----|
| `momspi16S_samp.rda` | 16S sample metadata (9,107 samples) | [HMP2Data](https://github.com/dozmorovlab/HMP2Data/) |
| `momspi16S_mtx.rda` | 16S OTU counts (7,665 taxa x 9,107 samples) | [HMP2Data](https://github.com/dozmorovlab/HMP2Data/) |
| `momspi16S_tax.rda` | 16S taxonomy (7,665 taxa x 7 ranks) | [HMP2Data](https://github.com/dozmorovlab/HMP2Data/) |
| `momspiCyto_samp.rda` | cytokine sample metadata (1,396 samples) | [HMP2Data](https://github.com/dozmorovlab/HMP2Data/) |
| `momspiCyto_mtx.rda` | cytokine concentrations (29 cytokines x 1,396 samples) | [HMP2Data](https://github.com/dozmorovlab/HMP2Data/) |
| `Vaginal_microbiome_log_ratios_per_sample_Pixu.csv` | pregnancy outcome (PTB / TB) per subject | [TEMPTED_paper](https://github.com/pixushi/TEMPTED_paper) |

The `.rda` objects are plain data frames and matrices, so the HMP2Data
package itself is not needed. In these objects the `subject_id` column
holds hashed identifiers, while the row names carry the original `EP...`
sample IDs that match the outcome file; the preprocessing script
extracts the subject from the row names.

`01a_momspi_preprocess.R` writes the analysis-ready data to
`data/momspi/processed/`, and `01` and `supp_momspi_stability.R` run it
automatically if that folder is missing.

### iPOP

The analysis-ready data ship with the `multi.tempted` package as `ipop`,
so `02` and `supp_sensitivity.R` need nothing else.
`02a_ipop_preprocess.R` documents how that object was built and rebuilds
it to `data/iPOP/processed/ipop.rds`. To rerun it, download these six
files from the iPOP data portal, <http://hmp2-data.stanford.edu/>
(Contrepois et al. 2020, *Cell*), into `data/iPOP/`:

| File | Contents |
|----|----|
| `NIHMS1588589-supplement-9.xlsx` | participant characteristics |
| `iPOP_visit_info.csv` | visit and sample information |
| `Exercise_Metabolomics.csv` | untargeted metabolomics (LC-MS) |
| `Exercise_Lipidomics.csv` | semi-targeted lipidomics (Lipidyzer) |
| `Exercise_targ_proteomics.csv` | targeted proteomics (immunoassays) |
| `Exercise_Proteomics.csv` | untargeted proteomics (SWATH-MS) |

## Run

``` r
source("01a_momspi_preprocess.R") # quick; needs the MOMS-PI inputs
source("01_momspi.R")             # quick
source("02a_ipop_preprocess.R")   # optional; needs the iPOP download
source("02_ipop.R")               # real data; MEFISTO is the slow step
source("03_validate.R")           # long, resumable
source("04_compare_unaligned.R")  # MEFISTO on unaligned data; slow
source("05_compare_aligned.R")    # long, resumable
source("supp_momspi_stability.R") # ~20 min for 100 refits; reuses saved results
source("supp_sensitivity.R")      # refits iPOP once per setting
source("supp_ablation.R")         # ~5 min; no MEFISTO
```

Each script prints its output folder on start-up, so results land in the
project whatever the session’s working directory happens to be.

## Reusing a finished run

The slow scripts save their results and reuse them on a re-run, which
makes re-plotting and restyling free:

- `03` and `05` checkpoint as they go and resume at the first missing
  run or seed, so an interruption costs at most `CKPT_EVERY` of them.
  Because every seed re-seeds the RNG from its own index, a resumed run
  reproduces an uninterrupted one exactly.
- `02`, `04` and `supp_momspi_stability.R` save once on completion; an
  interruption there costs the run.
- `01`, `supp_sensitivity.R` and `supp_ablation.R` refit from scratch
  every time.

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
scale, and in `05` the `SAMPLING` knob chooses how sampling times are
drawn, which dominates MEFISTO’s runtime more than anything else.

Plot aesthetics are set in a short block at the top of each script
rather than inside the package’s plotting functions, so they can be
changed without touching `multi.tempted`.

## Requirements

The `multi.tempted` and `MOFA2` packages, plus `ggplot2` (\>= 3.5.0) and
`patchwork` for the figures. `BiocManager::install("MOFA2")` also pulls
`basilisk`, which supplies the Python backend; no manual setup needed.
`02a_ipop_preprocess.R` also needs `PMCMRplus` and `readxl`.
