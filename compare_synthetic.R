# ============================================================================
# compare_synthetic.R
#
# Head-to-head of multiTEMPTED vs MEFISTO (MOFA2) on ONE fully synthetic dataset
# with a KNOWN block structure, formatted for both methods.
#
# Data design (block / biclustering structure):
#   * Subjects are split into r groups, one group per component.
#   * Each modality's features are split into r subsets, one per component.
#   * Component l is "on" only for group-l subjects and their subset-l features
#     (large loadings there, near-zero elsewhere).
#   * Every component uses the SAME temporal function across modalities, so both
#     models can represent the data (MEFISTO shares one factor trajectory across
#     views; multiTEMPTED shares the subject loading).
# So each component is a clean block: group-l subjects x subset-l features,
# evolving as xi_l(t). The question is whether each method recovers the temporal
# functions and separates the subject groups into the right components.
#
# The 3 subject GROUPS are LATENT: MEFISTO is given subjects-as-groups only (its
# per-subject GP), never told which subjects belong together. Both methods must
# discover the block structure.
#
# --- MEFISTO speed settings used here, and how to undo them -------------------
# The MEFISTO step is the slow part. Knobs are just below (edit them to taste):
#   MEF_MAXITER = 100      -> raise (e.g. 1000) for full convergence
#   MEF_CONVERGENCE "fast" -> set "slow" for a stricter convergence criterion
#   MEF_OPTIMISE_GP FALSE  -> the GP lengthscale optimisation. It DOES run with
#                             it TRUE (even with >=3 views), but it is much slower,
#                             especially with unaligned timepoints. Off here just
#                             for quick checking; turn it on for the real run.
#   SAMPLING "uniform"     -> how timepoints are placed (see below). Unaligned
#                             sampling is far slower for MEFISTO than "aligned".
# ============================================================================

library(multi.tempted)

# --- knobs --------------------------------------------------------------------
MEF_MAXITER     <- 1000        # MEFISTO iterations
MEF_CONVERGENCE <- "fast"     # "fast" or "slow"
MEF_OPTIMISE_GP <- TRUE      # TRUE = tune GP lengthscale (works but much slower)
SAMPLING        <- "aligned"  # "aligned"  : every subject & modality on one shared grid
                              # "clustered": small random jitter around a shared grid
                              # "uniform"  : fully random per subject & modality
                              #  ("aligned" gives MEFISTO co-measured samples across
                              #   modalities and is much faster; "uniform" is the
                              #   hardest, fully cross-modality-unaligned case.)


# ============================================================================
# GENERATOR
# ============================================================================
temporal_shape_library <- list(
  s_curve    = function(u) 1 / (1 + exp(-12 * (u - 0.5))),  # S-shaped
  increasing = function(u) exp(2 * u),                      # strictly increasing
  decreasing = function(u) exp(-2 * u)                      # strictly decreasing
)
.trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

# Block-structured data from the multiTEMPTED model with shared temporal dynamics.
generate_block_data <- function(n = 15, M = 3, p = 24, r = 3, n_timepoints = 8,
                                shapes = temporal_shape_library, lambda = c(8, 6, 4),
                                noise_sd = 0.1, offgroup = 0.05, sampling = "uniform",
                                time_range = c(0, 1), seed = 1) {
  set.seed(seed)
  mod <- paste0("mod", 1:M); subj <- sprintf("subj%02d", 1:n); PC <- paste0("PC", 1:r)
  a_end <- time_range[1]; b_end <- time_range[2]

  # subjects and features assigned to r blocks (one per component)
  subj_group <- rep(1:r, length.out = n)
  feat_group <- rep(1:r, length.out = p)

  # block subject loadings: large in-group, near-zero elsewhere
  A <- matrix(offgroup * abs(stats::rnorm(n * r)), n, r)
  for (l in 1:r) { g <- subj_group == l; A[g, l] <- 1 + 0.2 * abs(stats::rnorm(sum(g))) }
  A <- apply(A, 2, function(x) x / sqrt(sum(x^2))); dimnames(A) <- list(subj, PC)

  # block feature loadings per modality
  B <- lapply(1:M, function(m) {
    Bm <- matrix(offgroup * abs(stats::rnorm(p * r)), p, r)
    for (l in 1:r) { g <- feat_group == l; Bm[g, l] <- 1 + 0.2 * abs(stats::rnorm(sum(g))) }
    Bm <- apply(Bm, 2, function(x) x / sqrt(sum(x^2)))
    dimnames(Bm) <- list(sprintf("%s_feat%02d", mod[m], 1:p), PC); Bm
  }); names(B) <- mod

  Lambda <- matrix(lambda[1:r], M, r, byrow = TRUE, dimnames = list(mod, PC))

  # shared, L2-normalised temporal functions (same curve per component in every modality)
  fine <- seq(a_end, b_end, length.out = 2001)
  xi <- lapply(1:r, function(l) {
    g <- function(t) shapes[[l]]((t - a_end) / (b_end - a_end))
    nrm <- sqrt(.trapz(fine, g(fine)^2)); function(t) g(t) / nrm })

  grid <- seq(a_end, b_end, length.out = n_timepoints)
  cw   <- 0.3 * (b_end - a_end) / (n_timepoints - 1)   # cluster half-width
  draw_times <- function() {
    if (sampling == "aligned")   grid
    else if (sampling == "clustered")
      sort(pmin(pmax(grid + stats::runif(n_timepoints, -cw, cw), a_end), b_end))
    else sort(stats::runif(n_timepoints, a_end, b_end))
  }

  featuretables <- timepoints <- subjectID <- vector("list", M)
  for (m in 1:M) {
    per <- vector("list", n); tvec <- numeric(0); svec <- character(0)
    for (i in 1:n) {
      ts <- draw_times()
      Xi <- vapply(1:r, function(l) xi[[l]](ts), numeric(length(ts)))
      signal <- B[[m]] %*% (t(Xi) * (Lambda[m, ] * A[i, ]))
      per[[i]] <- t(signal + matrix(stats::rnorm(p * length(ts), sd = noise_sd), p))
      tvec <- c(tvec, ts); svec <- c(svec, rep(subj[i], length(ts)))
    }
    ft <- do.call(rbind, per); colnames(ft) <- rownames(B[[m]])
    featuretables[[m]] <- ft; timepoints[[m]] <- tvec; subjectID[[m]] <- svec
  }
  names(featuretables) <- names(timepoints) <- names(subjectID) <- mod

  list(featuretables = featuretables, timepoints = timepoints, subjectID = subjectID,
       truth = list(A = A, B = B, Lambda = Lambda, xi = xi,
                    subj_group = subj_group, feat_group = feat_group),
       params = list(n = n, M = M, p = p, r = r, n_timepoints = n_timepoints,
                     noise_sd = noise_sd, sampling = sampling, time_range = time_range,
                     modality_names = mod))
}


# ============================================================================
# MEFISTO helpers
# ============================================================================

# Long-format MOFA input; samples keyed by (subject, time) -- so modalities
# sampled at the SAME time share a sample (aligned => co-measured), while
# modalities at different times give separate, single-view samples.
mofa_long_from <- function(featuretables, timepoints, subjectID) {
  M <- length(featuretables); mod <- names(featuretables)
  do.call(rbind, lapply(1:M, function(m) {
    ft <- as.matrix(featuretables[[m]]); tp <- as.numeric(timepoints[[m]]); sid <- as.character(subjectID[[m]])
    samp <- paste0(sid, "@t", formatC(tp, format = "f", digits = 6))
    data.frame(sample = rep(samp, times = ncol(ft)), group = rep(sid, times = ncol(ft)),
               feature = rep(colnames(ft), each = nrow(ft)), view = mod[m],
               value = as.vector(ft), time = rep(tp, times = ncol(ft)), stringsAsFactors = FALSE)
  }))
}

run_mefisto <- function(featuretables, timepoints, subjectID, n_factors,
                        maxiter = MEF_MAXITER, convergence = MEF_CONVERGENCE,
                        optimise_gp = MEF_OPTIMISE_GP, seed = 1) {
  long   <- mofa_long_from(featuretables, timepoints, subjectID)
  cov_df <- unique(long[, c("sample", "time")])
  cov_df <- data.frame(sample = cov_df$sample, covariate = "time", value = cov_df$time)
  obj  <- MOFA2::create_mofa(long[, c("sample", "group", "feature", "view", "value")])
  obj  <- MOFA2::set_covariates(obj, covariates = cov_df)      # groups = subjects, covariate = time
  dopt <- MOFA2::get_default_data_options(obj)
  mopt <- MOFA2::get_default_model_options(obj);  mopt$num_factors <- n_factors
  topt <- MOFA2::get_default_training_options(obj)
  topt$seed <- seed; topt$convergence_mode <- convergence; topt$maxiter <- maxiter; topt$verbose <- FALSE
  eopt <- MOFA2::get_default_mefisto_options(obj)
  if (!optimise_gp) eopt$start_opt <- eopt$opt_freq <- as.integer(maxiter + 1L)
  obj <- MOFA2::prepare_mofa(obj, data_options = dopt, model_options = mopt,
                             training_options = topt, mefisto_options = eopt)
  MOFA2::run_mofa(obj, outfile = file.path(tempdir(), "mefisto_compare.hdf5"),
                  use_basilisk = TRUE, save_data = TRUE)
}

# Split each MEFISTO factor into subject scale (u) and temporal shape (v) via an
# SVD of its (subject x time-bin) matrix, mirroring multiTEMPTED's a * xi.
mefisto_decompose <- function(model, r, n_bins = 12) {
  Z  <- do.call(rbind, MOFA2::get_factors(model, groups = "all"))    # samples x factors
  sm <- MOFA2::samples_metadata(model); ord <- match(rownames(Z), sm$sample)
  tvec <- sm$time[ord]; gvec <- sm$group[ord]; subs <- sort(unique(gvec))
  br  <- seq(min(tvec), max(tvec), length.out = n_bins + 1)
  bin <- cut(tvec, br, include.lowest = TRUE, labels = FALSE)
  centers <- (utils::head(br, -1) + utils::tail(br, -1)) / 2
  shapes <- matrix(0, n_bins, r); uscale <- matrix(0, length(subs), r, dimnames = list(subs, NULL))
  for (k in 1:r) {
    Mmat <- matrix(NA_real_, length(subs), n_bins)
    for (si in seq_along(subs)) { sel <- gvec == subs[si]
      for (b in unique(bin[sel])) Mmat[si, b] <- mean(Z[sel & bin == b, k]) }
    keep <- colSums(!is.na(Mmat)) > 0
    Mk <- Mmat[, keep, drop = FALSE]
    for (b in seq_len(ncol(Mk))) { na <- is.na(Mk[, b]); if (any(na)) Mk[na, b] <- mean(Mk[!na, b]) }
    sv <- svd(Mk, nu = 1, nv = 1); shapes[keep, k] <- sv$v[, 1]; uscale[, k] <- sv$u[, 1]
  }
  list(centers = centers, shapes = shapes, uscale = uscale)
}


# ============================================================================
# 1. GENERATE + SUMMARISE
# ============================================================================
cat("== Generating one block-structured dataset ==\n")
sim <- generate_block_data(n = 15, M = 3, p = 24, r = 3, n_timepoints = 8,
                           sampling = SAMPLING, seed = 1)
pr <- sim$params
gsz <- table(sim$truth$subj_group); fsz <- table(sim$truth$feat_group)
cat(sprintf("  %d subjects in %d groups (sizes %s); %d modalities; %d features/modality\n",
            pr$n, pr$r, paste(gsz, collapse = "/"), pr$M, pr$p))
cat(sprintf("  each component = one subject group x one feature subset (sizes %s per modality)\n",
            paste(fsz, collapse = "/")))
cat(sprintf("  temporal functions shared across modalities: %s\n",
            paste(names(temporal_shape_library)[1:pr$r], collapse = ", ")))
cat(sprintf("  sampling = '%s', %d timepoints/subject, noise_sd = %.2f\n",
            pr$sampling, pr$n_timepoints, pr$noise_sd))


# ============================================================================
# 2. FIT BOTH METHODS
# ============================================================================
cat("\n== Fitting multiTEMPTED ==\n")
mt <- multitempted_all(sim$featuretables, sim$timepoints, sim$subjectID,
                       transforms = "none", do_ratio = FALSE, centralize = FALSE,
                       smooth = 1e-4, r = pr$r)

if (!requireNamespace("MOFA2", quietly = TRUE)) stop("Install MOFA2: BiocManager::install('MOFA2')")
cat(sprintf("\n== Fitting MEFISTO (maxiter=%d, %s, GP-opt=%s) ==\n",
            MEF_MAXITER, MEF_CONVERGENCE, MEF_OPTIMISE_GP))
mef  <- run_mefisto(sim$featuretables, sim$timepoints, sim$subjectID, n_factors = pr$r)
mdec <- mefisto_decompose(mef, pr$r)


# ============================================================================
# 3. COMPARE: functional estimation + group separation
# ============================================================================
r <- pr$r; M <- pr$M; grp <- sim$truth$subj_group; A_true <- sim$truth$A

# per-subject scores for each method
mt_load <- mt$A_hat                                          # subjects x r
Zsub    <- t(sapply(MOFA2::get_factors(mef, groups = "all"), colMeans))  # subjects x factors
Zsub    <- Zsub[rownames(mt_load), , drop = FALSE]

# ONE matching drives everything: assign each subject to its argmax component,
# then choose the component->group labelling that maximises agreement. This same
# component-per-group mapping is used for both group separation and functional
# recovery, so the two never disagree.
.all_perms <- function(v) if (length(v) == 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i)
    lapply(.all_perms(v[-i]), function(p) c(v[i], p))))
assign_map <- function(load) {
  pred  <- apply(abs(load), 1, which.max)                       # subject -> component
  perms <- .all_perms(1:r)
  bp    <- perms[[which.max(vapply(perms, function(pm) mean(pm[pred] == grp), numeric(1)))]]
  list(pred = pred,
       comp_of_group = vapply(1:r, function(l) which(bp == l), integer(1)),  # group -> component
       accuracy = mean(bp[pred] == grp))
}
mt_map  <- assign_map(mt_load); mt_ci  <- mt_map$comp_of_group
mef_map <- assign_map(Zsub);    mef_ci <- mef_map$comp_of_group

# --- group separation: per-group recall (fraction of a group's subjects whose
#     top component is the one matched to that group) ---
mt_recall  <- vapply(1:r, function(l) mean(mt_map$pred[grp == l]  == mt_ci[l]),  numeric(1))
mef_recall <- vapply(1:r, function(l) mean(mef_map$pred[grp == l] == mef_ci[l]), numeric(1))

# --- functional estimation: temporal recovery, same matching ---
grid_mt     <- mt$time_Zeta[[1]]
xi_true_mt  <- sapply(1:r, function(l) sim$truth$xi[[l]](grid_mt))
xi_true_mef <- sapply(1:r, function(l) sim$truth$xi[[l]](mdec$centers))
mt_temporal  <- sapply(1:r, function(l) mean(sapply(1:M, function(m)
                  abs(stats::cor(mt$Zeta_hat[[m]][, mt_ci[l]], xi_true_mt[, l])))))
mef_temporal <- sapply(1:r, function(l) abs(stats::cor(mdec$shapes[, mef_ci[l]], xi_true_mef[, l])))

cat("\n== Recovery vs ground truth (temporal = |cor| with true curve;",
    "recall = fraction of the group's subjects put on its component) ==\n")
tab <- data.frame(
  group        = paste0("group", 1:r),
  shape        = names(temporal_shape_library)[1:r],
  mT_temporal  = round(mt_temporal, 3),  MEF_temporal = round(mef_temporal, 3),
  mT_recall    = round(mt_recall, 2),    MEF_recall   = round(mef_recall, 2),
  row.names = NULL)
print(tab, row.names = FALSE)
cat(sprintf("\n  overall group-assignment accuracy:  multiTEMPTED %.2f | MEFISTO %.2f\n",
            mt_map$accuracy, mef_map$accuracy))


# ============================================================================
# 4. PDF: temporal-function recovery of each method
# ============================================================================
unit  <- function(v) { v <- v - mean(v); if (sqrt(sum(v^2)) > 0) v / sqrt(sum(v^2)) else v }
align <- function(est, ref) if (stats::cor(est, ref) < 0) -est else est
fine  <- seq(pr$time_range[1], pr$time_range[2], length.out = 200)

grDevices::pdf("compare_synthetic.pdf", width = 3.2 * r, height = 3.4)
op <- graphics::par(mfrow = c(1, r), mar = c(3.2, 3.2, 2.4, 1), mgp = c(1.9, 0.6, 0))
for (l in 1:r) {
  tru  <- unit(sim$truth$xi[[l]](fine))
  mtv  <- align(unit(mt$Zeta_hat[[1]][, mt_ci[l]]), unit(sim$truth$xi[[l]](grid_mt)))
  mefv <- align(unit(mdec$shapes[, mef_ci[l]]), unit(xi_true_mef[, l]))
  yl <- range(c(tru, mtv, mefv))
  plot(fine, tru, type = "l", lwd = 2.5, col = "black", ylim = yl,
       xlab = "time", ylab = "temporal loading (scaled)",
       main = sprintf("PC%d: %s", l, tab$shape[l]))
  graphics::lines(grid_mt, mtv, lwd = 2, lty = 2, col = "#D55E00")
  graphics::lines(mdec$centers, mefv, lwd = 2, lty = 3, col = "#0072B2")
  if (l == 1) graphics::legend("topleft", c("true", "multiTEMPTED", "MEFISTO"),
                               lwd = 2, lty = c(1, 2, 3), col = c("black", "#D55E00", "#0072B2"),
                               bty = "n", cex = 0.85)
}

# --- final page: the recovery table also printed to the log ---
graphics::par(mfrow = c(1, 1), mar = c(0.5, 0.5, 2.4, 0.5))
tbl <- c(
  sprintf("Block data: %d subjects in %d groups, %d modalities, sampling = '%s'",
          pr$n, pr$r, pr$M, pr$sampling),
  "",
  "Recovery (temporal = |cor| with true curve; recall = fraction of a group's",
  "subjects put on its component):",
  utils::capture.output(print(tab, row.names = FALSE)),
  "",
  sprintf("overall group-assignment accuracy:  multiTEMPTED %.2f | MEFISTO %.2f",
          mt_map$accuracy, mef_map$accuracy))
graphics::plot.new()
graphics::mtext("multiTEMPTED vs MEFISTO (log output)", side = 3, line = 0.5, font = 2)
graphics::text(0, 1, paste(tbl, collapse = "\n"), family = "mono", adj = c(0, 1), cex = 0.8)

graphics::par(op); grDevices::dev.off()
cat("\n  temporal-recovery overlay + table written to compare_synthetic.pdf\n")
