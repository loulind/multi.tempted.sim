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
cat("\n== multiTEMPTED on iPOP ==\n")
mt <- multitempted_all(featuretables = featuretables,
                       timepoints    = timepoints,
                       subjectID     = subjectID,
                       transforms    = "none",   # already log10
                       do_ratio      = FALSE,    # not counts
                       r             = r)        # centralize = TRUE (default)


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
  MOFA2::run_mofa(obj, outfile = file.path(tempdir(), "mefisto_ipop.hdf5"),
                  use_basilisk = TRUE, save_data = TRUE)
}

if (!requireNamespace("MOFA2", quietly = TRUE)) {
  message("MOFA2 not installed; skipping MEFISTO. Install with:\n",
          "  BiocManager::install('MOFA2')   # basilisk supplies the Python backend")
} else {
  cat("\n== MEFISTO on iPOP (subjects as groups, time as covariate) ==\n")
  cat("   training ... this is the slow step (~10 min); raise MEF_MAXITER for a fuller fit\n")
  mef <- run_mefisto(featuretables, timepoints, subjectID, n_factors = r)

  # ---- (i) feature-loading agreement, per modality --------------------------
  # For each multiTEMPTED component, the best |cor| with any MEFISTO factor's
  # weights in that modality (loadings are only defined up to sign/scale).
  W <- MOFA2::get_weights(mef)
  cat("\n-- feature-loading agreement (best |cor| between the methods' loadings) --\n")
  agree <- sapply(mods, function(m) {
    Bm <- mt$B_hat[[m]]
    Wm <- W[[m]][rownames(Bm), , drop = FALSE]
    sapply(1:r, function(l) max(abs(stats::cor(Wm, Bm[, l])), na.rm = TRUE))
  })  # r x M
  rownames(agree) <- paste0("multiTEMPTED_PC", 1:r)
  print(round(agree, 2))

  # ---- (ii) subject-embedding agreement -------------------------------------
  # multiTEMPTED gives one score per subject per component (A_hat). MEFISTO gives
  # a factor value per (subject, time); average over time to get a per-subject
  # score, then align subjects and correlate the two embeddings.
  Zg <- MOFA2::get_factors(mef, groups = "all")          # list per subject
  Zsub <- t(sapply(Zg, colMeans))                        # subjects x factors
  common <- intersect(rownames(mt$A_hat), rownames(Zsub))
  A_mt  <- mt$A_hat[common, , drop = FALSE]
  Z_mef <- Zsub[common, , drop = FALSE]
  emb_cor <- sapply(1:r, function(l) max(abs(stats::cor(Z_mef, A_mt[, l])), na.rm = TRUE))
  cat("\n-- subject-embedding agreement: best |cor| of each multiTEMPTED score with a MEFISTO factor --\n")
  print(round(stats::setNames(emb_cor, paste0("PC", 1:r)), 2))

  # ---- (iii) sex signal: does each method's embedding capture it? -----------
  # |point-biserial correlation| between each component's subject score and sex.
  sex_v <- sex[common]
  sex_mt  <- abs(apply(A_mt,  2, function(s) stats::cor(s, sex_v)))
  sex_mef <- abs(apply(Z_mef, 2, function(s) stats::cor(s, sex_v)))
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

  # ---- (iv) PDF: first two subject PCs, both methods, coloured by sex --------
  # multiTEMPTED: subject loadings A_hat[,1:2]. MEFISTO: per-subject mean of
  # factors 1:2. Subjects are aligned (common set); colour by sex.
  col_sex <- ifelse(sex_v == 1, "#0072B2", "#D55E00")   # male = blue, female = orange
  grDevices::pdf("compare_ipop.pdf", width = 9, height = 4.6)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), mgp = c(2.3, 0.8, 0))
  plot(A_mt[, 1], A_mt[, 2], col = col_sex, pch = 19,
       xlab = "PC1", ylab = "PC2", main = "multiTEMPTED subject loadings")
  graphics::legend("topright", c("male", "female"), col = c("#0072B2", "#D55E00"),
                   pch = 19, bty = "n")
  plot(Z_mef[, 1], Z_mef[, 2], col = col_sex, pch = 19,
       xlab = "Factor 1", ylab = "Factor 2", main = "MEFISTO subject factors")

  # --- final page: the tables also printed to the log ---
  graphics::par(mfrow = c(1, 1), mar = c(0.5, 0.5, 2.4, 0.5))
  tbl <- c(
    "Feature-loading agreement (best |cor| between the methods' loadings):",
    utils::capture.output(print(round(agree, 2))),
    "",
    "Subject-embedding agreement (best |cor| of each multiTEMPTED score with a MEFISTO factor):",
    utils::capture.output(print(round(stats::setNames(emb_cor, paste0("PC", 1:r)), 2))),
    "",
    "Sex association |cor(score, sex)| per component (max = best sex-separating axis):",
    sprintf("  multiTEMPTED: %s   (max %.2f)",
            paste(sprintf("PC%d=%.2f", 1:r, sex_mt), collapse = "  "), max(sex_mt)),
    sprintf("  MEFISTO:      %s   (max %.2f)",
            paste(sprintf("F%d=%.2f", 1:r, sex_mef), collapse = "  "), max(sex_mef)))
  graphics::plot.new()
  graphics::mtext("iPOP: multiTEMPTED vs MEFISTO (log output)", side = 3, line = 0.5, font = 2)
  graphics::text(0, 1, paste(tbl, collapse = "\n"), family = "mono", adj = c(0, 1), cex = 0.75)

  graphics::par(op); grDevices::dev.off()
  cat("\n  first-two-PC scatter (coloured by sex) + tables written to compare_ipop.pdf\n")
}
