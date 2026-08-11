# ============================================================================
# exp2_aligned_1xf.R   (storyline point 3)
#
# Setting that FAVORS MEFISTO -- yet multiTEMPTED still wins:
#   * 1 x f temporal functions: each component uses the SAME curve in every
#     modality (dynamics identical across modalities -> MEFISTO can represent it).
#   * sampling times set by SAMPLING (see below); the default "ALIGNED" gives
#     every subject and modality one shared time grid, so all modalities are
#     co-measured (MEFISTO's natural case).
#   * block structure: subjects split into r groups (one per component) with a
#     matching feature subset per modality; component l is "on" only for group-l
#     subjects and subset-l features.
#
# MEFISTO is run with group = subject (its per-subject GP); it is NOT told which
# subjects share a group. Over many random seeds we compare both methods on:
#   (a) group loading -- misclassification error (1 - group-assignment accuracy),
#   (b) function estimation -- |cor| of the estimated temporal curve with truth.
#
# Expected story: even on MEFISTO's home turf, multiTEMPTED classifies subjects
# and estimates the curves at least as well, usually better.
#
# --- sampling design -------------------------------------------------------
#   SAMPLING "ALIGNED"   -> every subject and modality sampled at exactly the
#                           same NT grid times (the original behaviour).
#            "CLUSTERED" -> samples land in a small cloud around the nominal grid
#                           times (sd = TIME_JITTER), drawn independently for each
#                           subject AND modality, so nothing is exactly co-measured
#                           but the visit schedule is still recognisable.
#            "UNALIGNED" -> every subject x modality draws its own random times,
#                           exactly as in 02_compare_unaligned.R.
#   NOTE: this drives MEFISTO's cost far more than any other knob. Its GP step
#   optimises a group x group covariance only when groups share covariate values,
#   so "ALIGNED" is by far the most expensive of the three (measured: a single GP
#   hyperparameter update takes ~1 h aligned vs ~10 min for 02's unaligned data).
#
# --- checkpointing ---------------------------------------------------------
#   Long runs are resumable. Finished seeds are appended to output/03_sims.csv
#   every CKPT_EVERY seeds; re-running the script reads that file and picks up
#   at the first missing seed (break between 90 and 100 -> only 10 seeds left).
#   Each seed re-seeds the RNG from its own index, so the recovered run gives
#   bit-identical metrics to an uninterrupted one. Delete the CSV to start over;
#   the script refuses to resume onto a checkpoint written under other settings.
#
# --- run / MEFISTO settings (efficient here; make them fair for the real run) --
#   N_SEEDS 100     -> paper scale; drop for a quick check
#   MEF_MAXITER 100 -> raise (e.g. 1000) for full convergence
#   MEF_CONVERGENCE -> "fast" (checking) or "slow" (fair)
#   MEF_OPTIMISE_GP -> FALSE (checking; faster). TRUE tunes the GP lengthscale
#                      (MEFISTO's fair setting; see the cost note above).
# ============================================================================

library(multi.tempted)

SAMPLING        <- "CLUSTERED"          # "ALIGNED" | "CLUSTERED" | "UNALIGNED"
TIME_JITTER     <- 0.02               # sd of the visit-time noise when CLUSTERED

N_SEEDS         <- 100
CKPT_EVERY      <- 10                 # flush results to output/03_sims.csv every n seeds
MEF_MAXITER     <- 1000
MEF_CONVERGENCE <- "slow"
MEF_OPTIMISE_GP <- TRUE

N <- 15; M <- 3; P <- 24; R <- 3; NT <- 8; NOISE <- 1.0; LAMBDA <- c(8, 6, 4); OFFGROUP <- 0.05

stopifnot(SAMPLING %in% c("ALIGNED", "CLUSTERED", "UNALIGNED"))


# ============================================================================
# helpers
# ============================================================================
shapes <- list(s_curve    = function(u) 1 / (1 + exp(-12 * (u - 0.5))),
               increasing = function(u) exp(2 * u),
               decreasing = function(u) exp(-2 * u))               # shared across modalities
.trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

# block generator, shared (1 x f) dynamics, ALIGNED sampling
gen_block <- function(seed) {
  set.seed(seed)
  mod <- paste0("mod", 1:M); subj <- sprintf("subj%02d", 1:N); PC <- paste0("PC", 1:R)
  subj_group <- rep(1:R, length.out = N); feat_group <- rep(1:R, length.out = P)
  A <- matrix(OFFGROUP * abs(stats::rnorm(N * R)), N, R)
  for (l in 1:R) { g <- subj_group == l; A[g, l] <- 1 + 0.2 * abs(stats::rnorm(sum(g))) }
  A <- apply(A, 2, function(x) x / sqrt(sum(x^2))); dimnames(A) <- list(subj, PC)
  B <- lapply(1:M, function(m) {
    Bm <- matrix(OFFGROUP * abs(stats::rnorm(P * R)), P, R)
    for (l in 1:R) { g <- feat_group == l; Bm[g, l] <- 1 + 0.2 * abs(stats::rnorm(sum(g))) }
    Bm <- apply(Bm, 2, function(x) x / sqrt(sum(x^2)))
    dimnames(Bm) <- list(sprintf("%s_feat%02d", mod[m], 1:P), PC); Bm }); names(B) <- mod
  Lambda <- matrix(LAMBDA[1:R], M, R, byrow = TRUE)
  fine <- seq(0, 1, length.out = 2001)
  xi <- lapply(1:R, function(l) { nrm <- sqrt(.trapz(fine, shapes[[l]](fine)^2)); function(t) shapes[[l]](t) / nrm })
  grid <- seq(0, 1, length.out = NT)                                # nominal visit schedule
  # Sampling times, drawn per (modality, subject). CLUSTERED jitter is reflected
  # at 0 and 1 rather than clamped, to avoid piling samples up on the two
  # endpoints (which would re-align them across subjects and modalities).
  .reflect01 <- function(u) { u <- abs(u); 1 - abs(1 - u) }
  ft <- tp <- sid <- vector("list", M)
  for (m in 1:M) {
    per <- vector("list", N); tv <- numeric(0); sv <- character(0)
    for (i in 1:N) {
      ts <- switch(SAMPLING,
                   ALIGNED   = grid,
                   CLUSTERED = sort(.reflect01(grid + stats::rnorm(NT, sd = TIME_JITTER))),
                   UNALIGNED = sort(stats::runif(NT, 0, 1)))
      Xi <- vapply(1:R, function(l) xi[[l]](ts), numeric(NT))       # truth at the ACTUAL times
      per[[i]] <- t(B[[m]] %*% (t(Xi) * (Lambda[m, ] * A[i, ])) +
                      matrix(stats::rnorm(P * NT, sd = NOISE), P))
      tv <- c(tv, ts); sv <- c(sv, rep(subj[i], NT))
    }
    x <- do.call(rbind, per); colnames(x) <- rownames(B[[m]]); ft[[m]] <- x; tp[[m]] <- tv; sid[[m]] <- sv
  }
  names(ft) <- names(tp) <- names(sid) <- mod
  list(featuretables = ft, timepoints = tp, subjectID = sid,
       truth = list(A = A, xi = xi, subj_group = subj_group))
}

mofa_long_from <- function(ft, tp, sid) {
  do.call(rbind, lapply(seq_along(ft), function(m) {
    x <- as.matrix(ft[[m]]); t <- as.numeric(tp[[m]]); s <- as.character(sid[[m]])
    data.frame(sample = rep(paste0(s, "@t", formatC(t, format = "f", digits = 6)), times = ncol(x)),
               group = rep(s, times = ncol(x)), feature = rep(colnames(x), each = nrow(x)),
               view = names(ft)[m], value = as.vector(x), time = rep(t, times = ncol(x)), stringsAsFactors = FALSE) }))
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
mefisto_decompose <- function(model, r, n_bins = 12) {
  Z <- do.call(rbind, MOFA2::get_factors(model, groups = "all")); sm <- MOFA2::samples_metadata(model)
  ord <- match(rownames(Z), sm$sample); tv <- sm$time[ord]; gv <- sm$group[ord]; subs <- sort(unique(gv))
  br <- seq(min(tv), max(tv), length.out = n_bins + 1); bin <- cut(tv, br, include.lowest = TRUE, labels = FALSE)
  centers <- (utils::head(br, -1) + utils::tail(br, -1)) / 2
  shp <- matrix(0, n_bins, r); usc <- matrix(0, length(subs), r, dimnames = list(subs, NULL))
  for (k in 1:r) {
    Mm <- matrix(NA_real_, length(subs), n_bins)
    for (si in seq_along(subs)) { sel <- gv == subs[si]; for (b in unique(bin[sel])) Mm[si, b] <- mean(Z[sel & bin == b, k]) }
    keep <- colSums(!is.na(Mm)) > 0; Mk <- Mm[, keep, drop = FALSE]
    for (b in seq_len(ncol(Mk))) { na <- is.na(Mk[, b]); if (any(na)) Mk[na, b] <- mean(Mk[!na, b]) }
    if (ncol(Mk) == 0 || stats::sd(as.vector(Mk)) < 1e-9) next   # inactive factor: leave zeros
    sv <- svd(Mk, nu = 1, nv = 1); shp[keep, k] <- sv$v[, 1]; usc[, k] <- sv$u[, 1] }
  list(centers = centers, shapes = shp, uscale = usc)
}
.perms <- function(v) if (length(v) == 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(.perms(v[-i]), function(p) c(v[i], p))))
# misclassification error via best-matching argmax assignment; also returns component-per-group
misclass <- function(load, grp, r) {
  pred <- apply(abs(load), 1, which.max); pm <- .perms(1:r)
  bp <- pm[[which.max(vapply(pm, function(q) mean(q[pred] == grp), numeric(1)))]]
  list(err = 1 - mean(bp[pred] == grp), comp_of_group = vapply(1:r, function(l) which(bp == l), integer(1)))
}


# ============================================================================
# run over seeds
# ============================================================================
if (!requireNamespace("MOFA2", quietly = TRUE)) stop("Install MOFA2: BiocManager::install('MOFA2')")
cat(sprintf("== exp2: 1 x f (shared) curves, %s sampling, block structure; %d seeds ==\n",
            tolower(SAMPLING), N_SEEDS))

# --- checkpointing ----------------------------------------------------------
# Every seed is self-contained: gen_block(s) calls set.seed(s) before drawing
# anything, so seed s always produces the same data and the same fit no matter
# what ran before it. That is what makes a resumed run agree with a start-to-
# finish run. Finished seeds are flushed to output/03_sims.csv every CKPT_EVERY
# seeds; on start-up we read that file and only run the seeds still missing.
dir.create("output", showWarnings = FALSE)
.csv <- file.path("output", "03_sims.csv")
.rds <- file.path("output", "03_sims.rds")   # same rows, exact doubles (see .write_ckpt)

# Config stamp: refuse to resume onto a checkpoint written under different
# settings, which would silently mix incompatible simulations.
.sig <- sprintf("%s|jit=%g|N=%d|M=%d|P=%d|R=%d|NT=%d|noise=%g|off=%g|lam=%s|mef=%d/%s/%s",
                SAMPLING, TIME_JITTER, N, M, P, R, NT, NOISE, OFFGROUP,
                paste(LAMBDA[1:R], collapse = "+"),   # no commas: keeps the CSV one-field-per-column
                MEF_MAXITER, MEF_CONVERGENCE, MEF_OPTIMISE_GP)

# The CSV is the readable artifact, but no decimal format R writes round-trips a
# double exactly (as.numeric mis-rounds long digit strings), so a CSV-only resume
# would shift recovered metrics by ~1e-16. The .rds alongside it stores the same
# rows as exact doubles and is what we actually resume from; the CSV still gates
# whether a resume happens at all, so deleting it starts a clean run either way.
.write_ckpt <- function(tab) {
  tab <- tab[order(tab$seed), , drop = FALSE]
  row.names(tab) <- NULL
  utils::write.csv(tab, .csv, row.names = FALSE)
  saveRDS(tab, .rds)
  invisible(tab)
}

done <- NULL
if (file.exists(.csv)) {
  prev <- utils::read.csv(.csv, stringsAsFactors = FALSE)
  if (file.exists(.rds)) {
    exact <- readRDS(.rds)
    if (identical(sort(exact$seed), sort(prev$seed))) prev <- exact
  }
  if (nrow(prev) > 0) {
    if (any(prev$config != .sig))
      stop(sprintf(paste0("%s was written under different settings:\n  checkpoint: %s\n  current:    %s\n",
                          "Delete or rename it to start a fresh run."),
                   .csv, prev$config[1], .sig))
    done <- prev[prev$seed >= 1 & prev$seed <= N_SEEDS, , drop = FALSE]
    done <- done[!duplicated(done$seed), , drop = FALSE]
  }
}
todo <- setdiff(seq_len(N_SEEDS), done$seed)

if (!is.null(done) && nrow(done) > 0)
  cat(sprintf("resuming from %s: %d/%d seeds already done, %d to run\n",
              .csv, nrow(done), N_SEEDS, length(todo)))
if (length(todo) == 0) cat("all seeds already present in the checkpoint; nothing to run\n")

res <- if (is.null(done)) list() else split(done, seq_len(nrow(done)))
for (s in todo) {
  cat(sprintf("[sim %d/%d] seed %d\n", s, N_SEEDS, s)); utils::flush.console()
  sim <- gen_block(s); grp <- sim$truth$subj_group
  t_mt <- system.time(mt <- multitempted_all(sim$featuretables, sim$timepoints, sim$subjectID,
                                             transforms = "none", do_ratio = FALSE, centralize = FALSE, smooth = 1e-4, r = R))[["elapsed"]]
  t_mef <- system.time(mef <- run_mefisto(sim$featuretables, sim$timepoints, sim$subjectID, n_factors = R))[["elapsed"]]
  mdec <- mefisto_decompose(mef, R)
  Zsub <- t(sapply(MOFA2::get_factors(mef, groups = "all"), colMeans))[rownames(mt$A_hat), , drop = FALSE]
  
  mmt <- misclass(mt$A_hat, grp, R); mmf <- misclass(Zsub, grp, R)
  # function estimation (shared curve per component)
  gmt <- mt$time_Zeta[[1]]; gmf <- mdec$centers
  .acor <- function(x, y) { v <- suppressWarnings(abs(stats::cor(x, y))); if (is.na(v)) 0 else v }
  mt_fe <- mean(vapply(1:R, function(l) mean(vapply(1:M, function(m)
    .acor(mt$Zeta_hat[[m]][, mmt$comp_of_group[l]], sim$truth$xi[[l]](gmt)), numeric(1))), numeric(1)))
  mef_fe <- mean(vapply(1:R, function(l)
    .acor(mdec$shapes[, mmf$comp_of_group[l]], sim$truth$xi[[l]](gmf)), numeric(1)))
  
  res[[as.character(s)]] <- data.frame(seed = s, mt_time = t_mt, mef_time = t_mef,
                                       mt_miscl = mmt$err, mef_miscl = mmf$err,
                                       mt_fe = mt_fe, mef_fe = mef_fe, config = .sig,
                                       stringsAsFactors = FALSE)
  cat(sprintf("  -> miscl mT=%.2f MEF=%.2f | function-est mT=%.3f MEF=%.3f | time mT=%.1fs MEF=%.1fs\n",
              mmt$err, mmf$err, mt_fe, mef_fe, t_mt, t_mef))
  
  if (s %% CKPT_EVERY == 0 || s == N_SEEDS) {
    .write_ckpt(do.call(rbind, res))
    cat(sprintf("  -- checkpoint: %d seeds saved to %s\n",
                length(res), normalizePath(.csv, mustWork = FALSE)))
  }
}
R_all <- do.call(rbind, res)
R_all <- R_all[order(R_all$seed), , drop = FALSE]   # seed order, so a resumed run
row.names(R_all) <- NULL                            # summarises exactly like a full one
.write_ckpt(R_all)

summ <- data.frame(
  metric = c("misclassification error", "function estimation (|cor|)", "compute time (s)"),
  multiTEMPTED = c(sprintf("%.3f (+/-%.3f)", mean(R_all$mt_miscl), sd(R_all$mt_miscl)),
                   sprintf("%.3f (+/-%.3f)", mean(R_all$mt_fe),    sd(R_all$mt_fe)),
                   sprintf("%.1f", median(R_all$mt_time))),
  MEFISTO      = c(sprintf("%.3f (+/-%.3f)", mean(R_all$mef_miscl), sd(R_all$mef_miscl)),
                   sprintf("%.3f (+/-%.3f)", mean(R_all$mef_fe),    sd(R_all$mef_fe)),
                   sprintf("%.1f", median(R_all$mef_time))))
cat("\n== summary over seeds (mean +/- sd) ==\n"); print(summ, row.names = FALSE)


# ============================================================================
# figure: boxplots across seeds + summary table page
# ============================================================================
dir.create("output", showWarnings = FALSE)
.pdf <- file.path("output", "03_compare_aligned.pdf")
.dl <- grDevices::dev.list(); if (!is.null(.dl)) for (.d in .dl[names(.dl) == "pdf"]) grDevices::dev.off(.d)
grDevices::pdf(.pdf, width = 9, height = 4.6)
graphics::par(mfrow = c(1, 2), mar = c(3, 4, 3, 1), mgp = c(2.3, 0.8, 0))
cols <- c("#D55E00", "#0072B2")
boxplot(list(multiTEMPTED = R_all$mt_miscl, MEFISTO = R_all$mef_miscl), col = cols,
        ylab = "misclassification error", main = "group loading", ylim = c(0, max(0.05, R_all$mef_miscl)))
boxplot(list(multiTEMPTED = R_all$mt_fe, MEFISTO = R_all$mef_fe), col = cols,
        ylab = "function estimation |cor|", main = "temporal curve recovery", ylim = c(min(R_all$mef_fe, 0.8), 1))
# table page
graphics::par(mfrow = c(1, 1), mar = c(0.5, 0.5, 2.4, 0.5))
tbl <- c(sprintf("exp2: 1 x f shared curves, block structure; %d seeds (N=%d, M=%d, p=%d, r=%d)",
                 N_SEEDS, N, M, P, R),
         sprintf("sampling: %s%s", SAMPLING,
                 if (SAMPLING == "CLUSTERED") sprintf(" (time jitter sd=%g)", TIME_JITTER) else ""),
         sprintf("MEFISTO: maxiter=%d, %s, GP-opt=%s", MEF_MAXITER, MEF_CONVERGENCE, MEF_OPTIMISE_GP), "",
         "Summary over seeds:", utils::capture.output(print(summ, row.names = FALSE)))
graphics::plot.new(); graphics::mtext("exp2: misclassification & function estimation (log output)", side = 3, line = 0.5, font = 2)
graphics::text(0, 1, paste(tbl, collapse = "\n"), family = "mono", adj = c(0, 1), cex = 0.85)
grDevices::dev.off()
cat(sprintf("\n  boxplots + table written to %s\n", normalizePath(.pdf, mustWork = FALSE)))
