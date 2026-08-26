# ============================================================================
# exp1_unaligned_mxf.R   (storyline point 2)
#
# Setting that is HARD for MEFISTO and natural for multiTEMPTED:
#   * m x f temporal functions: every (modality, component) pair has its OWN
#     curve, so the temporal dynamics differ ACROSS modalities.
#   * unaligned timepoints: each subject and each modality is sampled at its own
#     random times (no samples co-measured across modalities).
#
# Over a handful of random seeds we compare, for each method:
#   (a) computation time, and
#   (b) function estimation -- how well the estimated temporal curves match truth.
#
# Expected story: multiTEMPTED is fast and recovers every curve; MEFISTO is slow
# on unaligned data and, because a MEFISTO factor carries ONE trajectory shared
# across views, it cannot reproduce curves that differ across modalities.
#
# --- MEFISTO / run settings (efficient here; make them fair for the real run) --
#   N_SEEDS         -> number of random datasets (<= 10 in the paper)
#   MEF_MAXITER 100 -> raise (e.g. 1000) for full convergence
#   MEF_CONVERGENCE -> "fast" (checking) or "slow" (fair)
#   MEF_OPTIMISE_GP -> FALSE (checking; much faster). TRUE tunes the GP length-
#                      scale -- MEFISTO's fair setting, but far slower, ESPECIALLY
#                      on the unaligned data used here.
# ============================================================================

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
  if (is.null(d)) {                                     # Rscript 02_compare_unaligned.R
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

N_SEEDS         <- 10
MEF_MAXITER     <- 1000
MEF_CONVERGENCE <- "slow"
MEF_OPTIMISE_GP <- TRUE

# fixed data dimensions (change to taste; see README)
N <- 12; M <- 3; P <- 24; R <- 3; NT <- 8; NOISE <- 0.1; LAMBDA <- c(8, 6, 4)


# ============================================================================
# helpers
# ============================================================================
# nine distinct shapes on u in [0,1]; modality m, component l uses shape (m-1)*R+l
shape_library <- list(
  s_curve    = function(u) 1 / (1 + exp(-12 * (u - 0.5))),
  concave    = function(u) 1 - 4 * (u - 0.5)^2,
  decreasing = function(u) exp(-2 * u),
  increasing = function(u) exp(2 * u),
  m_shaped   = function(u) abs(sin(2 * pi * u)),
  convex     = function(u) 4 * (u - 0.5)^2,
  plateau    = function(u) 1 - exp(-8 * u),
  w_shaped   = function(u) 1 - abs(sin(2 * pi * u)),
  sine       = function(u) sin(2 * pi * u))
.trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
.ortho <- function(nr, r) qr.Q(qr(matrix(stats::rnorm(nr * r), nr, r)))[, 1:r, drop = FALSE]

# m x f generator with unaligned (uniform) sampling
gen_mxf <- function(seed) {
  set.seed(seed)
  mod <- paste0("mod", 1:M); subj <- sprintf("subj%02d", 1:N); PC <- paste0("PC", 1:R)
  A <- .ortho(N, R); dimnames(A) <- list(subj, PC)
  B <- lapply(1:M, function(m) { Bm <- .ortho(P, R)
    dimnames(Bm) <- list(sprintf("%s_feat%02d", mod[m], 1:P), PC); Bm }); names(B) <- mod
  Lambda <- matrix(LAMBDA[1:R], M, R, byrow = TRUE, dimnames = list(mod, PC))
  fine <- seq(0, 1, length.out = 2001)
  xi <- lapply(1:M, function(m) lapply(1:R, function(l) {
    f <- shape_library[[(m - 1) * R + l]]; nrm <- sqrt(.trapz(fine, f(fine)^2))
    function(t) f(t) / nrm }))                                  # xi[[m]][[l]]
  ft <- tp <- sid <- vector("list", M)
  for (m in 1:M) {
    per <- vector("list", N); tv <- numeric(0); sv <- character(0)
    for (i in 1:N) {
      ts <- sort(stats::runif(NT, 0, 1))                        # unaligned
      Xi <- vapply(1:R, function(l) xi[[m]][[l]](ts), numeric(length(ts)))
      per[[i]] <- t(B[[m]] %*% (t(Xi) * (Lambda[m, ] * A[i, ])) +
                    matrix(stats::rnorm(P * length(ts), sd = NOISE), P))
      tv <- c(tv, ts); sv <- c(sv, rep(subj[i], length(ts)))
    }
    x <- do.call(rbind, per); colnames(x) <- rownames(B[[m]])
    ft[[m]] <- x; tp[[m]] <- tv; sid[[m]] <- sv
  }
  names(ft) <- names(tp) <- names(sid) <- mod
  list(featuretables = ft, timepoints = tp, subjectID = sid,
       truth = list(A = A, B = B, xi = xi), modality_names = mod)
}

# MEFISTO: long format keyed by (subject,time); subjects are groups, time is covariate
mofa_long_from <- function(ft, tp, sid) {
  do.call(rbind, lapply(seq_along(ft), function(m) {
    x <- as.matrix(ft[[m]]); t <- as.numeric(tp[[m]]); s <- as.character(sid[[m]])
    data.frame(sample = rep(paste0(s, "@t", formatC(t, format = "f", digits = 6)), times = ncol(x)),
               group = rep(s, times = ncol(x)), feature = rep(colnames(x), each = nrow(x)),
               view = names(ft)[m], value = as.vector(x), time = rep(t, times = ncol(x)),
               stringsAsFactors = FALSE) }))
}
run_mefisto <- function(ft, tp, sid, n_factors) {
  long <- mofa_long_from(ft, tp, sid); cd <- unique(long[, c("sample", "time")])
  obj <- MOFA2::create_mofa(long[, c("sample", "group", "feature", "view", "value")])
  obj <- MOFA2::set_covariates(obj, data.frame(sample = cd$sample, covariate = "time", value = cd$time))
  mo <- MOFA2::get_default_model_options(obj); mo$num_factors <- n_factors
  to <- MOFA2::get_default_training_options(obj)
  to$seed <- 1; to$convergence_mode <- MEF_CONVERGENCE; to$maxiter <- MEF_MAXITER; to$verbose <- FALSE
  eo <- MOFA2::get_default_mefisto_options(obj)
  if (!MEF_OPTIMISE_GP) eo$start_opt <- eo$opt_freq <- as.integer(MEF_MAXITER + 1L)
  obj <- MOFA2::prepare_mofa(obj, data_options = MOFA2::get_default_data_options(obj),
                             model_options = mo, training_options = to, mefisto_options = eo)
  f <- file.path(tempdir(), "mef.hdf5")
  MOFA2::run_mofa(obj, outfile = f, use_basilisk = TRUE, save_data = TRUE)
  MOFA2::load_model(f, remove_inactive_factors = FALSE)   # keep all r factors (MEFISTO may drop 'inactive' ones)
}
# per-factor (subject-scale, temporal-shape) via SVD of the subject x time-bin matrix
mefisto_decompose <- function(model, r, n_bins = 12) {
  Z <- do.call(rbind, MOFA2::get_factors(model, groups = "all")); sm <- MOFA2::samples_metadata(model)
  ord <- match(rownames(Z), sm$sample); tv <- sm$time[ord]; gv <- sm$group[ord]; subs <- sort(unique(gv))
  br <- seq(min(tv), max(tv), length.out = n_bins + 1); bin <- cut(tv, br, include.lowest = TRUE, labels = FALSE)
  centers <- (utils::head(br, -1) + utils::tail(br, -1)) / 2
  shapes <- matrix(0, n_bins, r); uscale <- matrix(0, length(subs), r, dimnames = list(subs, NULL))
  for (k in 1:r) {
    Mm <- matrix(NA_real_, length(subs), n_bins)
    for (si in seq_along(subs)) { sel <- gv == subs[si]; for (b in unique(bin[sel])) Mm[si, b] <- mean(Z[sel & bin == b, k]) }
    keep <- colSums(!is.na(Mm)) > 0; Mk <- Mm[, keep, drop = FALSE]
    for (b in seq_len(ncol(Mk))) { na <- is.na(Mk[, b]); if (any(na)) Mk[na, b] <- mean(Mk[!na, b]) }
    if (ncol(Mk) == 0 || stats::sd(as.vector(Mk)) < 1e-9) next   # inactive factor: leave zeros
    sv <- svd(Mk, nu = 1, nv = 1); shapes[keep, k] <- sv$v[, 1]; uscale[, k] <- sv$u[, 1]
  }
  list(centers = centers, shapes = shapes, uscale = uscale)
}
.match <- function(load, A) sapply(1:ncol(A), function(l) which.max(abs(stats::cor(load, A[, l]))))


# ============================================================================
# run over seeds
# ============================================================================
if (!requireNamespace("MOFA2", quietly = TRUE)) stop("Install MOFA2: BiocManager::install('MOFA2')")
cat(sprintf("== exp1: m x f (modality-specific) curves, unaligned; %d seeds ==\n", N_SEEDS))

# --- saving the completed run ----------------------------------------------
# The MEFISTO fits here take hours, so the whole run is saved ONCE when the last
# seed finishes: metrics to 02_sims.csv and everything the figure needs to
# 02_sims.rds. Re-running the script then goes straight to the figure, so the
# plot can be restyled without repeating the fits. Delete the files to force a
# fresh run; the script refuses to reuse results saved under other settings.
.csv <- file.path(.outdir, "02_sims.csv")
.rds <- file.path(.outdir, "02_sims.rds")
.sig <- sprintf("N=%d|M=%d|P=%d|R=%d|NT=%d|noise=%g|lam=%s|seeds=%d|mef=%d/%s/%s",
                N, M, P, R, NT, NOISE, paste(LAMBDA[1:R], collapse = "+"),
                N_SEEDS, MEF_MAXITER, MEF_CONVERGENCE, MEF_OPTIMISE_GP)

cached <- NULL
if (file.exists(.rds)) {
  cc <- readRDS(.rds)
  if (!identical(cc$config, .sig))
    stop(sprintf(paste0("%s was saved under different settings:\n  saved:   %s\n  current: %s\n",
                        "Delete or rename it to start a fresh run."), .rds, cc$config, .sig))
  cached <- cc
  cat(sprintf("reusing the completed run in %s (%d seeds); skipping the fits\n",
              normalizePath(.rds, mustWork = FALSE), nrow(cc$R_all)))
}

res <- vector("list", N_SEEDS); keep_seed1 <- NULL
for (s in if (is.null(cached)) seq_len(N_SEEDS) else integer(0)) {
  sim <- gen_mxf(s)
  t_mt <- system.time(mt <- multitempted_all(sim$featuretables, sim$timepoints, sim$subjectID,
             transforms = "none", do_ratio = FALSE, centralize = FALSE, smooth = 1e-4, r = R))[["elapsed"]]
  t_mef <- system.time(mef <- run_mefisto(sim$featuretables, sim$timepoints, sim$subjectID, n_factors = R))[["elapsed"]]
  mdec <- mefisto_decompose(mef, R)

  # function estimation per (modality, component): |cor(estimate, true xi_ml)|
  mt_k  <- .match(mt$A_hat, sim$truth$A)
  mef_k <- .match(mdec$uscale, sim$truth$A)
  gmt <- mt$time_Zeta[[1]]; gmf <- mdec$centers
  mt_fe <- mef_fe <- matrix(NA_real_, M, R)
  for (m in 1:M) for (l in 1:R) {
    mt_fe[m, l]  <- abs(stats::cor(mt$Zeta_hat[[m]][, mt_k[l]], sim$truth$xi[[m]][[l]](gmt)))
    mef_fe[m, l] <- abs(stats::cor(mdec$shapes[, mef_k[l]],     sim$truth$xi[[m]][[l]](gmf)))  # ONE shape per component
  }
  res[[s]] <- data.frame(seed = s, mt_time = t_mt, mef_time = t_mef,
                         mt_fe = mean(mt_fe), mef_fe = mean(mef_fe))
  if (s == 1) keep_seed1 <- list(sim = sim, mt = mt, mdec = mdec, mt_k = mt_k, mef_k = mef_k,
                                 gmt = gmt, gmf = gmf, mt_fe = mt_fe, mef_fe = mef_fe)
  cat(sprintf("  seed %d: time mT=%.1fs MEF=%.1fs | function-est mT=%.3f MEF=%.3f\n",
              s, t_mt, t_mef, mean(mt_fe), mean(mef_fe)))
}
R_all <- if (is.null(cached)) do.call(rbind, res) else cached$R_all

summ <- data.frame(
  metric = c("compute time (s)", "function estimation (|cor|)"),
  multiTEMPTED = c(sprintf("%.1f (median)", median(R_all$mt_time)), sprintf("%.3f", mean(R_all$mt_fe))),
  MEFISTO      = c(sprintf("%.1f (median)", median(R_all$mef_time)), sprintf("%.3f", mean(R_all$mef_fe))))
cat("\n== summary over seeds ==\n"); print(summ, row.names = FALSE)


# ============================================================================
# figure: estimated vs true curves (seed 1) + summary table page
# ============================================================================
.unit  <- function(v) { v <- v - mean(v); if (sqrt(sum(v^2)) > 0) v / sqrt(sum(v^2)) else v }
.align <- function(e, ref) if (stats::cor(e, ref) < 0) -e else e
fine   <- seq(0, 1, length.out = 200)

.pdf <- file.path(.outdir, "02_compare_unaligned.pdf")

if (is.null(cached)) {
  k1 <- keep_seed1
  curve_df <- do.call(rbind, lapply(1:M, function(m) do.call(rbind, lapply(1:R, function(l) {
    tru  <- .unit(k1$sim$truth$xi[[m]][[l]](fine))
    mtv  <- .align(.unit(k1$mt$Zeta_hat[[m]][, k1$mt_k[l]]),
                   .unit(k1$sim$truth$xi[[m]][[l]](k1$gmt)))
    mefv <- .align(.unit(k1$mdec$shapes[, k1$mef_k[l]]),
                   .unit(k1$sim$truth$xi[[m]][[l]](k1$gmf)))
    pan  <- sprintf("Modality %d - PC %d", m, l)
    rbind(data.frame(time = fine,   value = tru,  series = "true",         panel = pan),
          data.frame(time = k1$gmt, value = mtv,  series = "multiTEMPTED", panel = pan),
          data.frame(time = k1$gmf, value = mefv, series = "MEFISTO",      panel = pan))
  }))))
  curve_df$series <- factor(curve_df$series, levels = c("true", names(METHOD_COLS)))
  # save the finished run: metrics as CSV, plus everything the figure needs
  utils::write.csv(R_all, .csv, row.names = FALSE)
  saveRDS(list(config = .sig, R_all = R_all, curve_df = curve_df), .rds)
  cat(sprintf("  run saved to %s and %s\n",
              normalizePath(.csv, mustWork = FALSE),
              normalizePath(.rds, mustWork = FALSE)))
} else {
  curve_df <- cached$curve_df
}

p_curves <- ggplot(curve_df, aes(x = time, y = value, colour = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ panel, ncol = R, scales = "free_y") +
  scale_colour_manual(values = c(true = "black", METHOD_COLS), name = NULL) +
  scale_linetype_manual(values = c(true = "solid", multiTEMPTED = "22", MEFISTO = "42"),
                        name = NULL) +
  labs(title    = "Estimated vs true temporal loadings",
       subtitle = "modality-specific curves, unaligned sampling times (seed 1)",
       x = "Time", y = "Loading") +
  THEME_MS + theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle   = element_text(hjust = 0.5),
                   legend.position = "bottom")

ggsave(.pdf, p_curves, width = 3.1 * R, height = 2.6 * M + 1.2, bg = "white")
cat(sprintf("\n  curves written to %s\n", normalizePath(.pdf, mustWork = FALSE)))
