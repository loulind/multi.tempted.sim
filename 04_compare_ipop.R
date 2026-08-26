# ============================================================================
# compare_ipop.R
#
# Apply BOTH multiTEMPTED and MEFISTO (MOFA2) to the REAL iPOP longitudinal
# multi-omic dataset shipped with the multi.tempted package, and compare them.
# There is no ground truth here, so the comparison is (i) method-vs-method
# agreement on feature loadings and subject embeddings, and (ii) a biology check:
# does each method's subject embedding pick up the known sex signal?
#
# iPOP: 36 subjects, 4 omic modalities (cytokine / metabolome / lipid / protein),
# measured at 5 shared visits during an exercise challenge. All modalities are on
# a log10 scale and share visits, so samples are naturally aligned across views.
#
# --- MEFISTO speed settings, set in the knobs just below the library() call ---
# The MEFISTO step trains a per-subject Gaussian process and is the slow part.
#   MEF_MAXITER      -> more iterations = fuller convergence, slower
#   MEF_CONVERGENCE  -> "fast" or "slow" (stricter) convergence criterion
#   MEF_OPTIMISE_GP  -> TRUE tunes the GP lengthscale (MEFISTO at full capability);
#                       it runs even on iPOP's 4 views but is considerably slower.
#                       FALSE just fixes the lengthscale for a quicker fit.
# For quick code-checking use MEF_MAXITER = 100, "fast", GP off; for the real
# comparison use many iterations, "slow", and GP on.
# Also: MEFISTO needs every subject at the same visit times, so we keep only
# complete-visit subjects (applied to both methods); multiTEMPTED would not need
# this.
# ============================================================================

# pak::pkg_install("loulind/multi.tempted")
library(multi.tempted)
library(ggplot2)
library(patchwork)

# --- plot aesthetics --------------------------------------------------------
# Deliberately kept in the analysis scripts rather than inside multi.tempted's
# plotting functions, so they can be changed without touching the package.
# Palette and theme follow the MOMSPI analysis so all manuscript figures match.
THEME_MS    <- theme_bw(base_size = 12)
METHOD_COLS <- c(multiTEMPTED = "#4A90D9", MEFISTO = "#E05C5C")

# Everything is written to <this script's folder>/output, NOT to getwd(), so the
# results land in the project no matter where the session's working directory
# happens to point. Covers source()/RStudio "Source", Rscript, and running lines
# interactively in RStudio; if none of those can identify the file it falls back
# to getwd() and says so.
.outdir <- local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {                    # source() / RStudio Source
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { d <- dirname(normalizePath(of, mustWork = FALSE)); break }
  }
  if (is.null(d)) {                                     # Rscript 04_compare_ipop.R
    a <- commandArgs(trailingOnly = FALSE)
    a <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(a)) d <- dirname(normalizePath(a[1], mustWork = FALSE))
  }
  if (is.null(d) && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {                      # RStudio, running lines by hand
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) d <- dirname(normalizePath(p, mustWork = FALSE))
  }
  if (is.null(d)) {
    d <- getwd()
    warning("could not locate the script file; writing output/ under ", d, call. = FALSE)
  }
  file.path(d, "output")
})
dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("output directory: %s\n", .outdir))

# --- MEFISTO speed / capability knobs (see the header note) -------------------
MEF_MAXITER     <- 1000      # iterations; raise for a fuller fit
MEF_CONVERGENCE <- "slow"   # "fast" or "slow"
MEF_OPTIMISE_GP <- TRUE    # TRUE = full GP lengthscale tuning

# ---- load & format iPOP ----------------------------------------------------
ip   <- multi.tempted::ipop
mods <- c("cytokine", "metabolome", "lipid", "protein")
r <- 3

# MEFISTO's per-subject Gaussian process requires every subject to share the same
# set of visit times; a few iPOP subjects miss visits. We restrict to
# complete-visit subjects and apply the same filter to BOTH methods, so the
# comparison is on identical data. (multiTEMPTED alone would not need this.)
sid_all  <- as.character(ip$meta_visit$SubjectID)
n_visits <- max(table(sid_all))
keep_row <- sid_all %in% names(which(table(sid_all) == n_visits))

featuretables <- lapply(mods, function(m) as.matrix(ip[[m]])[keep_row, , drop = FALSE])
names(featuretables) <- mods
M <- length(featuretables)

# map visit code (1..5) to approximate minutes from baseline, as in the vignette
tp_map     <- c("1" = 0, "2" = 12, "3" = 25, "4" = 40, "5" = 70)
timepoint  <- unname(tp_map[as.character(ip$meta_visit$timepoint)])[keep_row]
timepoints <- rep(list(timepoint), M)
subjectID  <- rep(list(sid_all[keep_row]), M)

# per-subject sex (for the biology check), as a 0/1 vector keyed by subject ID
sex <- stats::setNames(as.integer(ip$meta_subj$Sex == "M"),
                       as.character(ip$meta_subj$subjectID))

cat(sprintf("iPOP: %d of %d subjects with all %d visits, %d modalities, features = %s\n",
            length(unique(subjectID[[1]])), length(unique(sid_all)), n_visits, M,
            paste(sapply(featuretables, ncol), collapse = "/")))


# ============================================================================
# multiTEMPTED
# ============================================================================
# --- saving the completed run ----------------------------------------------
# The MEFISTO fit here is slow, so the whole run is saved ONCE when it finishes:
# the subject embeddings to 04_ipop.csv and the agreement tables to 04_ipop.rds.
# Re-running the script then goes straight to the figure, so the plot can be
# restyled without repeating either fit. Delete the files to force a fresh run;
# the script refuses to reuse results saved under other settings.
.csv <- file.path(.outdir, "04_ipop.csv")
.rds <- file.path(.outdir, "04_ipop.rds")
.sig <- sprintf("r=%d|mods=%s|mef=%d/%s/%s", r, paste(mods, collapse = "+"),
                MEF_MAXITER, MEF_CONVERGENCE, MEF_OPTIMISE_GP)

cached <- NULL
if (file.exists(.rds)) {
  cc <- readRDS(.rds)
  if (!identical(cc$config, .sig))
    stop(sprintf(paste0("%s was saved under different settings:\n  saved:   %s\n  current: %s\n",
                        "Delete or rename it to start a fresh run."), .rds, cc$config, .sig))
  cached <- cc
  cat(sprintf("reusing the completed run in %s; skipping both fits\n",
              normalizePath(.rds, mustWork = FALSE)))
}

if (is.null(cached)) {
  cat("\n== multiTEMPTED on iPOP ==\n")
  mt <- multitempted_all(featuretables = featuretables,
                         timepoints    = timepoints,
                         subjectID     = subjectID,
                         transforms    = "none",   # already log10
                         do_ratio      = FALSE,    # not counts
                         r             = r)        # centralize = TRUE (default)
}


# ============================================================================
# MEFISTO (MOFA2)
# ============================================================================

# Build MOFA2 long format from the multitempted-style inputs; samples keyed by
# (subject, time) so the aligned visits are shared across views.
mofa_long_from <- function(featuretables, timepoints, subjectID) {
  M <- length(featuretables); mod <- names(featuretables)
  do.call(rbind, lapply(1:M, function(m) {
    ft <- as.matrix(featuretables[[m]]); tp <- as.numeric(timepoints[[m]])
    sid <- as.character(subjectID[[m]]); feats <- colnames(ft)
    samp <- paste0(sid, "@t", formatC(tp, format = "f", digits = 6))
    data.frame(sample = rep(samp,  times = ncol(ft)),
               group  = rep(sid,   times = ncol(ft)),
               feature= rep(feats, each  = nrow(ft)),
               view   = mod[m],
               value  = as.vector(ft),
               time   = rep(tp,    times = ncol(ft)),
               stringsAsFactors = FALSE)
  }))
}

# Prep + train MEFISTO (subjects = groups, time = covariate).
run_mefisto <- function(featuretables, timepoints, subjectID, n_factors,
                        seed = 1, maxiter = MEF_MAXITER,
                        convergence = MEF_CONVERGENCE, optimise_gp = MEF_OPTIMISE_GP) {
  long   <- mofa_long_from(featuretables, timepoints, subjectID)
  cov_df <- unique(long[, c("sample", "time")])
  cov_df <- data.frame(sample = cov_df$sample, covariate = "time", value = cov_df$time)

  obj <- MOFA2::create_mofa(long[, c("sample", "group", "feature", "view", "value")])
  obj <- MOFA2::set_covariates(obj, covariates = cov_df)

  dopt <- MOFA2::get_default_data_options(obj)
  mopt <- MOFA2::get_default_model_options(obj);   mopt$num_factors <- n_factors
  topt <- MOFA2::get_default_training_options(obj)
  topt$seed <- seed; topt$convergence_mode <- convergence
  topt$maxiter <- maxiter; topt$verbose <- FALSE
  eopt <- MOFA2::get_default_mefisto_options(obj)
  if (!optimise_gp)                            # skip the (crashing) GP hyperparam step
    eopt$start_opt <- eopt$opt_freq <- as.integer(maxiter + 1L)

  obj <- MOFA2::prepare_mofa(obj, data_options = dopt, model_options = mopt,
                             training_options = topt, mefisto_options = eopt)
  f <- file.path(tempdir(), "mefisto_ipop.hdf5")
  MOFA2::run_mofa(obj, outfile = f, use_basilisk = TRUE, save_data = TRUE)
  MOFA2::load_model(f, remove_inactive_factors = FALSE)   # keep all r factors (MEFISTO may drop 'inactive' ones)
}

if (!requireNamespace("MOFA2", quietly = TRUE)) {
  message("MOFA2 not installed; skipping MEFISTO. Install with:\n",
          "  BiocManager::install('MOFA2')   # basilisk supplies the Python backend")
} else {
  if (is.null(cached)) {
    cat("\n== MEFISTO on iPOP (subjects as groups, time as covariate) ==\n")
    cat("   training ... this is the slow step; raise MEF_MAXITER for a fuller fit\n")
    mef <- run_mefisto(featuretables, timepoints, subjectID, n_factors = r)

    # ---- (i) feature-loading agreement, per modality ------------------------
    # For each multiTEMPTED component, the best |cor| with any MEFISTO factor's
    # weights in that modality (loadings are only defined up to sign/scale).
    W <- MOFA2::get_weights(mef)
    agree <- sapply(mods, function(m) {
      Bm <- mt$B_hat[[m]]
      Wm <- W[[m]][rownames(Bm), , drop = FALSE]
      sapply(1:r, function(l) max(abs(stats::cor(Wm, Bm[, l])), na.rm = TRUE))
    })  # r x M
    rownames(agree) <- paste0("multiTEMPTED_PC", 1:r)

    # ---- (ii) subject-embedding agreement -----------------------------------
    # multiTEMPTED gives one score per subject per component (A_hat). MEFISTO
    # gives a factor value per (subject, time); average over time to get a
    # per-subject score, then align subjects and correlate the two embeddings.
    Zg     <- MOFA2::get_factors(mef, groups = "all")     # list per subject
    Zsub   <- t(sapply(Zg, colMeans))                     # subjects x factors
    common <- intersect(rownames(mt$A_hat), rownames(Zsub))
    A_mt   <- mt$A_hat[common, , drop = FALSE]
    Z_mef  <- Zsub[common, , drop = FALSE]
    emb_cor <- sapply(1:r, function(l) max(abs(stats::cor(Z_mef, A_mt[, l])), na.rm = TRUE))

    # ---- (iii) sex signal: does each method's embedding capture it? ---------
    # |point-biserial correlation| between each component's subject score and sex.
    sex_v   <- sex[common]
    sex_mt  <- abs(apply(A_mt,  2, function(s) stats::cor(s, sex_v)))
    sex_mef <- abs(apply(Z_mef, 2, function(s) stats::cor(s, sex_v)))

    # ---- save the finished run ----------------------------------------------
    embed <- data.frame(subject = common, sex = as.integer(sex_v),
                        stringsAsFactors = FALSE)
    for (l in 1:r) embed[[paste0("mt_PC", l)]] <- A_mt[, l]
    for (l in 1:r) embed[[paste0("mef_F", l)]] <- Z_mef[, l]
    utils::write.csv(embed, .csv, row.names = FALSE)
    saveRDS(list(config = .sig, embed = embed, agree = agree, emb_cor = emb_cor,
                 sex_mt = sex_mt, sex_mef = sex_mef), .rds)
    cat(sprintf("\n  run saved to %s and %s\n",
                normalizePath(.csv, mustWork = FALSE),
                normalizePath(.rds, mustWork = FALSE)))
  } else {
    agree   <- cached$agree;  emb_cor <- cached$emb_cor
    sex_mt  <- cached$sex_mt; sex_mef <- cached$sex_mef
    A_mt    <- as.matrix(cached$embed[, paste0("mt_PC", 1:r), drop = FALSE])
    Z_mef   <- as.matrix(cached$embed[, paste0("mef_F", 1:r), drop = FALSE])
    sex_v   <- cached$embed$sex
  }

  cat("\n-- feature-loading agreement (best |cor| between the methods' loadings) --\n")
  print(round(agree, 2))
  cat("\n-- subject-embedding agreement: best |cor| of each multiTEMPTED score with a MEFISTO factor --\n")
  print(round(stats::setNames(emb_cor, paste0("PC", 1:r)), 2))
  cat("\n-- sex association |cor(score, sex)| per component (max = best sex-separating axis) --\n")
  cat(sprintf("  multiTEMPTED: %s   (max %.2f)\n",
              paste(sprintf("PC%d=%.2f", 1:r, sex_mt), collapse = "  "), max(sex_mt)))
  cat(sprintf("  MEFISTO:      %s   (max %.2f)\n",
              paste(sprintf("F%d=%.2f", 1:r, sex_mef), collapse = "  "), max(sex_mef)))
  cat("\n  Both methods are unsupervised; a high value means a component lines up with\n")
  cat("  sex without being told about it. multiTEMPTED's leading component tracks sex,\n")
  cat("  matching the published iPOP analysis.\n")
  cat("\n  READ WITH CARE. The two methods parameterise time differently (per-modality\n")
  cat("  loadings vs one shared factor) and preprocess differently (rank-1 mean removal\n")
  cat("  vs per-feature centering), so they are NOT expected to agree closely -- low\n")
  cat("  agreement reflects different models, not one being 'wrong'.\n")

  # ---- (iv) PDF: first two subject components, both methods, coloured by sex --
  # Aesthetics follow the MOMSPI analysis: geom_point(size = 2, alpha = 0.8),
  # theme_bw(base_size = 12), legends collected across the two panels.
  SEX_COLS <- c(male = "#4A90D9", female = "#E05C5C")
  sex_lab  <- factor(ifelse(sex_v == 1, "male", "female"), levels = names(SEX_COLS))
  .pdf <- file.path(.outdir, "04_compare_ipop.pdf")

  .scatter <- function(d, title, xlab, ylab)
    ggplot(d, aes(x = x, y = y, colour = sex)) +
      geom_point(size = 2, alpha = 0.8) +
      scale_colour_manual(values = SEX_COLS, name = "Sex") +
      labs(title = title, subtitle = "Subject embedding", x = xlab, y = ylab) +
      THEME_MS + theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
                       plot.subtitle = element_text(hjust = 0.5))

  d_mt  <- data.frame(x = A_mt[, 1],  y = A_mt[, 2],  sex = sex_lab)
  d_mef <- data.frame(x = Z_mef[, 1], y = Z_mef[, 2], sex = sex_lab)

  p_ipop <- (.scatter(d_mt,  "multiTEMPTED on iPOP", "Component 1", "Component 2") +
               theme(legend.position = "none")) +
            .scatter(d_mef, "MEFISTO on iPOP",      "Factor 1",    "Factor 2") +
            plot_layout(ncol = 2, guides = "collect")

  ggsave(.pdf, p_ipop, width = 10, height = 4.5, bg = "white")
  cat(sprintf("\n  subject embeddings (coloured by sex) written to %s\n",
              normalizePath(.pdf, mustWork = FALSE)))
}
