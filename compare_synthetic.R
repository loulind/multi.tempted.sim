# ============================================================================
# compare_synthetic.R
#
# Head-to-head of multiTEMPTED vs MEFISTO (MOFA2) on ONE fully synthetic dataset,
# formatted two ways -- once for each method -- so the comparison is as direct as
# possible (identical underlying values and ground truth; only the container and
# the way each method ingests time differ).
#
# Design choices that make the comparison FAIR:
#   * SHARED temporal dynamics across modalities. multiTEMPTED allows a different
#     time curve per modality; MEFISTO uses one factor trajectory shared across
#     views. Making the true curves shared means BOTH models can represent the
#     data, so neither is handicapped by construction.
#   * Both formats carry the SAME time UNALIGNMENT across subjects AND modalities
#     (each subject, and each modality, is sampled at its own random times).
#
# The truth is generated from the multiTEMPTED model (shared subject loading a,
# per-modality feature loadings b, shared temporal functions xi). Because the
# curves are shared, this is exactly representable by MEFISTO too: its weights W
# play the role of b, and its factor Z_i(t) = a_i * xi(t) (a subject-scaled copy
# of the shared curve).
#
# --- MEFISTO time-saving settings used here, and how to undo them -------------
#   maxiter = 100        -> raise (e.g. 1000) for full convergence
#   convergence "fast"   -> set MEF_CONVERGENCE <- "slow" for a stricter fit
#   n = 12 subjects      -> MEFISTO's per-subject GP scales poorly; raising n
#                           makes it much slower (multiTEMPTED does not care)
# ============================================================================

library(multi.tempted)

MEF_MAXITER <- 100
MEF_CONVERGENCE  <- "fast"
MEF_OPTIMISE_GP  <- TRUE
SEED <- 1


# ============================================================================
# GENERATOR  (copy of the generator from validate_synthetic.R)
# ============================================================================

temporal_shape_library <- list(
  s_curve = function(u) 1 / (1 + exp(-12 * (u - 0.5))),  # S-shaped
  increasing = function(u) exp(2 * u), # strictly increasing
  decreasing = function(u) exp(-2 * u) # strictly decreasing
)

.trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

.make_loadings <- function(nr, r) qr.Q(qr(matrix(stats::rnorm(nr * r), nr, r)))[, 1:r, drop = FALSE]

# Generate synthetic data from the multiTEMPTED model with dynamics SHARED across
# modalities and uniformly-random (unaligned) sampling times.
generate_shared_data <- function(n = 12, M = 2, p = 20, r = 3, n_timepoints = 10,
                                 shapes = temporal_shape_library,
                                 lambda = c(8, 5, 3), noise_sd = 0.1,
                                 time_range = c(0, 1), seed = SEED) {
  set.seed(seed)
  if (length(shapes) < r) stop("need >= r shapes")
  mod  <- paste0("mod", 1:M); subj <- sprintf("subj%02d", 1:n); PC <- paste0("PC", 1:r)

  Lambda <- matrix(lambda[1:r], M, r, byrow = TRUE, dimnames = list(mod, PC))
  A <- .make_loadings(n, r); dimnames(A) <- list(subj, PC)
  B <- lapply(1:M, function(m) {
    Bm <- .make_loadings(p, r)
    dimnames(Bm) <- list(sprintf("%s_feat%02d", mod[m], 1:p), PC); Bm
  }); names(B) <- mod

  # L2-normalise the r shared shapes on the interval T
  a_end <- time_range[1]; b_end <- time_range[2]; fine <- seq(a_end, b_end, length.out = 2001)
  xi <- lapply(1:r, function(l) {
    g <- function(t) shapes[[l]]((t - a_end) / (b_end - a_end))
    nrm <- sqrt(.trapz(fine, g(fine)^2)); function(t) g(t) / nrm
  })

  featuretables <- timepoints <- subjectID <- vector("list", M)
  for (m in 1:M) {
    per <- vector("list", n); tvec <- numeric(0); svec <- character(0)
    for (i in 1:n) {
      ts <- sort(stats::runif(n_timepoints, a_end, b_end))      # unaligned per subject & modality
      Xi <- vapply(1:r, function(l) xi[[l]](ts), numeric(length(ts)))  # q x r
      coef   <- Lambda[m, ] * A[i, ]
      signal <- B[[m]] %*% (t(Xi) * coef)                       # p x q
      per[[i]] <- t(signal + matrix(stats::rnorm(p * length(ts), sd = noise_sd), p))
      tvec <- c(tvec, ts); svec <- c(svec, rep(subj[i], length(ts)))
    }
    ft <- do.call(rbind, per); colnames(ft) <- rownames(B[[m]])
    featuretables[[m]] <- ft; timepoints[[m]] <- tvec; subjectID[[m]] <- svec
  }
  names(featuretables) <- names(timepoints) <- names(subjectID) <- mod

  list(featuretables = featuretables, timepoints = timepoints, subjectID = subjectID,
       truth = list(A = A, B = B, Lambda = Lambda, xi = xi),
       params = list(n = n, M = M, p = p, r = r, n_timepoints = n_timepoints,
                     noise_sd = noise_sd, time_range = time_range, modality_names = mod))
}


# ============================================================================
# MEFISTO helpers
# ============================================================================

# Long-format MOFA input; samples keyed by (subject, modality, time). Because the
# modalities are sampled at different times, each sample carries ONE modality --
# MEFISTO must link modalities through its temporal GP alone.
mofa_long_from <- function(featuretables, timepoints, subjectID) {
  M <- length(featuretables); mod <- names(featuretables)
  do.call(rbind, lapply(1:M, function(m) {
    ft <- as.matrix(featuretables[[m]]); tp <- as.numeric(timepoints[[m]]); sid <- as.character(subjectID[[m]])
    samp <- paste0(sid, "@", mod[m], "t", formatC(tp, format = "f", digits = 5))
    data.frame(sample = rep(samp, times = ncol(ft)), group = rep(sid, times = ncol(ft)),
               feature = rep(colnames(ft), each = nrow(ft)), view = mod[m],
               value = as.vector(ft), time = rep(tp, times = ncol(ft)), stringsAsFactors = FALSE)
  }))
}

run_mefisto <- function(featuretables, timepoints, subjectID, n_factors,
                        maxiter = MEF_MAXITER, convergence = MEF_CONVERGENCE,
                        optimise_gp = MEF_OPTIMISE_GP, seed = SEED) {
  long   <- mofa_long_from(featuretables, timepoints, subjectID)
  cov_df <- unique(long[, c("sample", "time")])
  cov_df <- data.frame(sample = cov_df$sample, covariate = "time", value = cov_df$time)

  obj  <- MOFA2::create_mofa(long[, c("sample", "group", "feature", "view", "value")])
  obj  <- MOFA2::set_covariates(obj, covariates = cov_df)
  dopt <- MOFA2::get_default_data_options(obj)
  mopt <- MOFA2::get_default_model_options(obj);  mopt$num_factors <- n_factors
  topt <- MOFA2::get_default_training_options(obj)
  topt$seed <- seed; topt$convergence_mode <- convergence; topt$maxiter <- maxiter; topt$verbose <- FALSE
  eopt <- MOFA2::get_default_mefisto_options(obj)
  eopt$model_groups   <- FALSE   # default TRUE
  eopt$frac_inducing  <- 0.5     # default 0.75
  if (!optimise_gp) eopt$start_opt <- eopt$opt_freq <- as.integer(maxiter + 1L)
  obj <- MOFA2::prepare_mofa(obj, data_options = dopt, model_options = mopt,
                             training_options = topt, mefisto_options = eopt)
  obj@mefisto_options$sparseGP <- TRUE
  MOFA2::run_mofa(obj, outfile = file.path(tempdir(), "mefisto_compare.hdf5"),
                  use_basilisk = TRUE, save_data = TRUE)
}

# Split each MEFISTO factor into a subject scale (u) and a temporal shape (v) via
# an SVD of its (subject x time-bin) matrix, mirroring multiTEMPTED's a * xi.
mefisto_decompose <- function(model, r, n_bins = 15) {
  Z  <- do.call(rbind, MOFA2::get_factors(model, groups = "all")) # samples x factors
  sm <- MOFA2::samples_metadata(model) 
  ord  <- match(rownames(Z), sm$sample)
  tvec <- sm$time[ord]; gvec <- sm$group[ord]
  subs <- sort(unique(gvec))
  br   <- seq(min(tvec), max(tvec), length.out = n_bins + 1)
  bin  <- cut(tvec, br, include.lowest = TRUE, labels = FALSE)
  centers <- (utils::head(br, -1) + utils::tail(br, -1)) / 2

  shapes <- matrix(0, n_bins, r); uscale <- matrix(0, length(subs), r, dimnames = list(subs, NULL))
  for (k in 1:r) {
    Mmat <- matrix(NA_real_, length(subs), n_bins)
    for (si in seq_along(subs)) {
      sel <- gvec == subs[si]
      for (b in unique(bin[sel])) Mmat[si, b] <- mean(Z[sel & bin == b, k])
    }
    keep <- colSums(!is.na(Mmat)) > 0                    # drop empty bins
    Mk <- Mmat[, keep, drop = FALSE]
    for (b in seq_len(ncol(Mk))) { na <- is.na(Mk[, b]); if (any(na)) Mk[na, b] <- mean(Mk[!na, b]) }
    sv <- svd(Mk, nu = 1, nv = 1)
    shapes[keep, k] <- sv$v[, 1]; uscale[, k] <- sv$u[, 1]
  }
  list(centers = centers, shapes = shapes, uscale = uscale)
}


# ============================================================================
# 1. GENERATE ONE DATASET, FORMAT FOR BOTH METHODS
# ============================================================================
cat("== Generating one synthetic dataset (shared dynamics, unaligned times) ==\n")
sim <- generate_shared_data(n = 12, M = 2, p = 20, r = 3, n_timepoints = 10,
                            shapes = temporal_shape_library, lambda = c(8, 5, 3),
                            noise_sd = 0.1, seed = SEED)
pr  <- sim$params
long <- mofa_long_from(sim$featuretables, sim$timepoints, sim$subjectID)

cat("\n-- the two formats of the SAME data --\n")
cat(sprintf("  SHARED: %d subjects, %d modalities, %d features/modality, %d timepoints/subject,\n",
            pr$n, pr$M, pr$p, pr$n_timepoints))
cat(sprintf("          r=%d components, noise_sd=%.2f; times drawn uniformly at random per subject\n",
            pr$r, pr$noise_sd))
cat(sprintf("          AND per modality (so times are unaligned across both).\n"))
cat(sprintf("  multiTEMPTED format: %d feature tables, each %d samples x %d features, with a\n",
            pr$M, pr$n * pr$n_timepoints, pr$p))
cat("          per-modality time vector (time handled natively).\n")
cat(sprintf("  MEFISTO format:      one long table, %d rows, %d samples keyed by (subject,modality,\n",
            nrow(long), length(unique(long$sample))))
cat("          time); each sample has ONE modality (others missing), so MEFISTO links\n")
cat("          modalities only through its temporal GP.\n")
cat("  DIFFERENCE: representation only -- identical values, subjects, times, and ground truth.\n")


# ============================================================================
# 2. FIT BOTH METHODS
# ============================================================================
cat("\n== Fitting multiTEMPTED ==\n")
mt <- multitempted_all(sim$featuretables, sim$timepoints, sim$subjectID,
                       transforms = "none", do_ratio = FALSE, centralize = FALSE,
                       smooth = 1e-4, r = pr$r)

if (!requireNamespace("MOFA2", quietly = TRUE)) {
  stop("MOFA2 not installed. Install with: BiocManager::install('MOFA2')")
}
cat("\n== Fitting MEFISTO (GP hyperparam optimisation; a few minutes) ==\n")
mef <- run_mefisto(sim$featuretables, sim$timepoints, sim$subjectID, n_factors = pr$r)
mdec <- mefisto_decompose(mef, pr$r)
Wm <- MOFA2::get_weights(mef)


# ============================================================================
# 3. RECOVERY METRICS (both methods vs the SAME ground truth)
# ============================================================================
# All loadings are identified up to sign/scale, so we use absolute correlation.
r <- pr$r; M <- pr$M; mods <- pr$modality_names
A_true <- sim$truth$A

# match estimated components to true components greedily by subject loadings
greedy_match <- function(est_subj) {         # est_subj: n x r
  cm <- abs(stats::cor(est_subj, A_true)); mt <- integer(r); used <- logical(r)
  for (k in 1:r) { c <- cm[k, ]; c[used] <- -Inf; j <- which.max(c); mt[k] <- j; used[j] <- TRUE }
  mt   # mt[k] = true component matched to estimated column k
}
inv <- function(perm) { o <- integer(length(perm)); o[perm] <- seq_along(perm); o }  # true -> est col

mt_perm <- inv(greedy_match(mt$A_hat))       # mt_perm[l] = multiTEMPTED column for true comp l

# MEFISTO: for each true component, the factor whose subject-scale best matches.
# On fully cross-modality-unaligned data MEFISTO's factors are largely single-
# modality (see the integration note below), so we match by best fit rather than
# assume one integrated factor spans both modalities.
mef_k <- sapply(1:r, function(l) which.max(abs(stats::cor(mdec$uscale, A_true[, l]))))

# best |cor| over ALL factors' weights in a modality (robust to zero-weight cols)
safe_cor_max <- function(Wmat, b) {
  cs <- suppressWarnings(abs(stats::cor(Wmat, b)))
  if (all(is.na(cs))) 0 else max(cs, na.rm = TRUE)
}

# true temporal curves on grids
grid_mt  <- mt$time_Zeta[[1]]
xi_true_mt  <- sapply(1:r, function(l) sim$truth$xi[[l]](grid_mt))
xi_true_mef <- sapply(1:r, function(l) sim$truth$xi[[l]](mdec$centers))

metrics <- data.frame(
  component = paste0("PC", 1:r),
  shape     = names(temporal_shape_library)[1:r],
  # feature loadings (avg over modalities); MEFISTO: best-matching factor per modality
  mT_feature  = sapply(1:r, function(l) mean(sapply(1:M, function(m)
                  abs(stats::cor(mt$B_hat[[m]][, mt_perm[l]], sim$truth$B[[m]][, l]))))),
  MEF_feature = sapply(1:r, function(l) mean(sapply(1:M, function(m)
                  safe_cor_max(Wm[[m]][rownames(sim$truth$B[[m]]), , drop = FALSE], sim$truth$B[[m]][, l])))),
  # subject loadings
  mT_subject  = sapply(1:r, function(l) abs(stats::cor(mt$A_hat[, mt_perm[l]], A_true[, l]))),
  MEF_subject = sapply(1:r, function(l) abs(stats::cor(mdec$uscale[, mef_k[l]], A_true[, l]))),
  # temporal function
  mT_temporal  = sapply(1:r, function(l) mean(sapply(1:M, function(m)
                   abs(stats::cor(mt$Zeta_hat[[m]][, mt_perm[l]], xi_true_mt[, l]))))),
  MEF_temporal = sapply(1:r, function(l) abs(stats::cor(mdec$shapes[, mef_k[l]], xi_true_mef[, l]))),
  row.names = NULL)

cat("\n== Recovery vs ground truth (absolute correlation; 1.000 = perfect) ==\n")
print(cbind(metrics[, 1:2], round(metrics[, -(1:2)], 3)), row.names = FALSE)

# integration: share of each MEFISTO factor's weight energy in its top modality
energy   <- sapply(mods, function(v) colSums(Wm[[v]]^2)) # factors x M
frac_top <- apply(energy, 1, function(e) max(e) / sum(e))
cat(sprintf("\n-- MEFISTO integration: mean weight-energy share in each factor's top modality: %.2f (1/M = %.2f) --\n",
            mean(frac_top), 1 / M))
cat("   Near 1 means factors are single-modality: with no samples co-measured across\n")
cat("   modalities (fully unaligned), MEFISTO cannot fuse them, whereas multiTEMPTED\n")
cat("   links modalities through the shared subject loading regardless of timing.\n")


# ============================================================================
# 4. PDF: temporal-function recovery of each method
# ============================================================================
unit <- function(v) { v <- v - mean(v); if (sqrt(sum(v^2)) > 0) v / sqrt(sum(v^2)) else v }
align <- function(est, ref) if (stats::cor(est, ref) < 0) -est else est
fine  <- seq(pr$time_range[1], pr$time_range[2], length.out = 200)

grDevices::pdf("compare_synthetic.pdf", width = 3.2 * r, height = 3.4)
op <- graphics::par(mfrow = c(1, r), mar = c(3.2, 3.2, 2.4, 1), mgp = c(1.9, 0.6, 0))
for (l in 1:r) {
  tru  <- unit(sim$truth$xi[[l]](fine))
  mtv  <- align(unit(mt$Zeta_hat[[1]][, mt_perm[l]]), unit(sim$truth$xi[[l]](grid_mt)))
  mefv <- align(unit(mdec$shapes[, mef_k[l]]), unit(xi_true_mef[, l]))
  yl <- range(c(tru, mtv, mefv))
  plot(fine, tru, type = "l", lwd = 2.5, col = "black", ylim = yl,
       xlab = "time", ylab = "temporal loading (scaled)",
       main = sprintf("PC%d: %s", l, metrics$shape[l]))
  graphics::lines(grid_mt, mtv, lwd = 2, lty = 2, col = "#D55E00")
  graphics::lines(mdec$centers, mefv, lwd = 2, lty = 3, col = "#0072B2")
  if (l == 1) graphics::legend("topleft", c("true", "multiTEMPTED", "MEFISTO"),
                               lwd = 2, lty = c(1, 2, 3), col = c("black", "#D55E00", "#0072B2"),
                               bty = "n", cex = 0.85)
}
graphics::par(op); grDevices::dev.off()
cat("\n  temporal-recovery overlay written to compare_synthetic.pdf\n")
