# ============================================================================
# 01b_validate_larger.R
#
# The validation of 01a_validate.R repeated at feature dimensions closer to real
# omics data: each modality keeps its 20 signal features and gains features with
# zero loadings in every component, giving p = 500, 100 and 30. Only a small
# fraction of features therefore carries signal, as in practice.
#
# Two cohort sizes are run over the same noise sweep, so the effect of sample
# size at these dimensions can be read off directly.
#
# Reported per noise level and cohort size:
#   * recovery of the subject, feature and temporal loadings (|cor| with truth),
#   * the proportion of the signal features ranked among the largest absolute
#     estimated loadings,
#   * runtime per run.
#
# Writes output/01b_validate_larger.pdf. Runs are flushed to
# output/01b_validate_larger_runs.csv every CKPT_EVERY runs and the script
# resumes at the first missing run, so a long sweep can be interrupted safely.
#
# The generator and the recovery metric come from 01a_validate.R rather than
# being copied, so the two scripts cannot drift apart.
# ============================================================================

library(multi.tempted)
library(ggplot2)
library(patchwork)

.here <- local({
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of, mustWork = FALSE)))
  }
  a <- commandArgs(trailingOnly = FALSE)
  a <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(a)) return(dirname(normalizePath(a[1], mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  warning("could not locate 01b_validate_larger.R; using the working directory",
          call. = FALSE)
  getwd()
})
.outdir <- file.path(.here, "output")
dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("output directory: %s\n", .outdir))

THEME_MS     <- theme_bw(base_size = 12)
LOADING_COLS <- c(subject = "#4A90D9", feature = "#E05C5C", temporal = "#5CB85C")

# ---- knobs -----------------------------------------------------------------
N_SUBS       <- c(30, 200)                       # cohort sizes, small first
P_TARGET     <- c(500, 100, 30)                  # features per modality
SIGNAL_P     <- 20                               # signal features per modality
NOISE_LEVELS <- c(0.1, 0.5, 1.0, 2.0, 4.0)
N_SEEDS      <- 20
SMOOTH       <- 1e-4
MAXITER      <- 100     # iterations per component, passed to multitempted_all
EPSILON      <- 1e-4    # convergence threshold on the loading change
CKPT_EVERY   <- 10
TIMING_PASS  <- TRUE    # re-time one seed per (n, noise) sequentially

# ---- generator and recovery metric from 01a_validate.R ---------------------
# Evaluated into a private environment so 01a's own output-directory and theme
# definitions cannot overwrite the ones set above.
.a <- local({
  src <- readLines(file.path(.here, "01a_validate.R"))
  cut <- grep("^# ------ RUN: noise study", src)
  if (!length(cut)) cut <- grep("^NOISE_LEVELS", src)[1]
  env <- new.env(parent = globalenv())
  eval(parse(text = src[1:(cut[1] - 1)]), envir = env)
  env
})
generate_multitempted_data <- .a$generate_multitempted_data
evaluate_recovery          <- .a$evaluate_recovery

# ---- append zero-loading (noise-only) features ------------------------------
pad_features <- function(sim, p_target, noise_sd, seed) {
  set.seed(seed * 1000L + 7L)
  for (m in seq_len(sim$params$M)) {
    extra <- p_target[m] - ncol(sim$featuretables[[m]])
    stopifnot(extra >= 0)
    if (extra == 0) next
    nr <- nrow(sim$featuretables[[m]])
    noise <- matrix(stats::rnorm(nr * extra, sd = noise_sd), nr, extra)
    colnames(noise) <- sprintf("%s_null%04d", names(sim$featuretables)[m], seq_len(extra))
    sim$featuretables[[m]] <- cbind(sim$featuretables[[m]], noise)
    Bm <- sim$truth$B[[m]]
    sim$truth$B[[m]] <- rbind(Bm, matrix(0, extra, ncol(Bm),
                                         dimnames = list(colnames(noise), colnames(Bm))))
  }
  sim
}

# proportion of the SIGNAL_P true signal features ranked in the top SIGNAL_P by
# |estimated loading|, averaged over components, per modality
signal_detection <- function(sim, fit, match) {
  M <- sim$params$M; r <- sim$params$r
  vapply(seq_len(M), function(m) {
    Bhat <- fit$B_hat[[m]]; Btrue <- sim$truth$B[[m]]
    mean(vapply(seq_len(r), function(lh) {
      truth_idx <- which(abs(Btrue[, match[lh]]) > 0)
      top_idx   <- order(abs(Bhat[, lh]), decreasing = TRUE)[seq_along(truth_idx)]
      length(intersect(truth_idx, top_idx)) / length(truth_idx)
    }, numeric(1)))
  }, numeric(1))
}

# ---- checkpointing ----------------------------------------------------------
.csv <- file.path(.outdir, "01b_validate_larger_runs.csv")
.sig <- sprintf("n=%s|p=%s|signal=%d|seeds=%d|smooth=%g|iter=%d/%g|a=sqrt(n)|mt-only",
                paste(N_SUBS, collapse = "+"), paste(P_TARGET, collapse = "+"),
                SIGNAL_P, N_SEEDS, SMOOTH, MAXITER, EPSILON)

done <- NULL
if (file.exists(.csv)) {
  prev <- utils::read.csv(.csv, stringsAsFactors = FALSE)
  if (nrow(prev) > 0) {
    if (any(prev$config != .sig))
      stop(sprintf(paste0("%s was written under different settings:\n  saved:   %s\n",
                          "  current: %s\nDelete it to start a fresh run."),
                   .csv, prev$config[1], .sig))
    done <- prev
    if (!"parallel" %in% names(done)) done$parallel <- FALSE
    cat(sprintf("resuming: %d runs already done\n", nrow(done)))
  }
}

grid <- expand.grid(seed = seq_len(N_SEEDS), noise = NOISE_LEVELS, n = N_SUBS)
grid$key <- sprintf("n%d_%g_%d", grid$n, grid$noise, grid$seed)
todo <- if (is.null(done)) grid else grid[!(grid$key %in% done$key), , drop = FALSE]
cat(sprintf("== %s subjects x %d noise levels x %d seeds; %d runs to go ==\n",
            paste(N_SUBS, collapse = "/"), length(NOISE_LEVELS), N_SEEDS, nrow(todo)))

# Runs are independent and each is seeded from its own (noise, seed) pair, so
# evaluating them in parallel gives results identical to running in sequence.
NCORES <- max(1L, min(6L, parallel::detectCores() - 2L))

fit_one <- function(n_sub, nl, s) {
  sim <- generate_multitempted_data(n = n_sub, M = 3, p = SIGNAL_P, r = 3,
                                    noise_sd = nl, seed = s)
  sim <- pad_features(sim, P_TARGET, nl, s)

  # The package reports "Converged: dif = X after N iterations" even when the
  # loop exits because it reached maxiter, so capture the messages and judge
  # convergence from dif against EPSILON rather than from the wording.
  msgs <- character(0)
  el <- system.time(
    fit <- withCallingHandlers(
      multitempted_all(sim$featuretables, sim$timepoints, sim$subjectID,
                       transforms = "none", do_ratio = FALSE, centralize = FALSE,
                       smooth = SMOOTH, maxiter = MAXITER, epsilon = EPSILON,
                       r = sim$params$r),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
      }))[["elapsed"]]

  mm    <- regmatches(msgs, regexec("dif = ([0-9.e+-]+) after ([0-9]+) iterations", msgs))
  mm    <- mm[lengths(mm) == 3]
  iters <- as.integer(vapply(mm, `[`, character(1), 3))
  difs  <- as.numeric(vapply(mm, `[`, character(1), 2))
  list(sim = sim, fit = fit, secs = el,
       max_iters = if (length(iters)) max(iters) else NA_integer_,
       n_at_maxiter = sum(iters >= MAXITER & difs > EPSILON))
}

run_one <- function(k) {
  n_sub <- todo$n[k]; nl <- todo$noise[k]; s <- todo$seed[k]
  o   <- fit_one(n_sub, nl, s)
  rp  <- evaluate_recovery(o$sim, o$fit)
  det <- signal_detection(o$sim, o$fit, rp$match)
  data.frame(key = todo$key[k], n = n_sub, noise = nl, seed = s, secs = o$secs,
             max_iters = o$max_iters, n_at_maxiter = o$n_at_maxiter,
             subject  = mean(rp$per_component$cor_A),
             feature  = mean(rp$per_component$cor_B_mean),
             temporal = mean(rp$per_component$cor_Z_mean),
             det_mod1 = det[1], det_mod2 = det[2], det_mod3 = det[3],
             det_mean = mean(det), config = .sig,
             parallel = NCORES > 1L, stringsAsFactors = FALSE)
}

res <- if (is.null(done)) list() else split(done, seq_len(nrow(done)))
if (nrow(todo) > 0) {
  cat(sprintf("running %d jobs on %d cores\n", nrow(todo), NCORES))
  chunks <- split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / CKPT_EVERY))
  for (ci in seq_along(chunks)) {
    ch  <- chunks[[ci]]
    out <- parallel::mclapply(ch, function(k) tryCatch(run_one(k), error = function(e) e),
                              mc.cores = NCORES)
    for (j in seq_along(ch)) {
      r1 <- out[[j]]
      if (inherits(r1, "error"))
        stop(sprintf("run %s failed: %s", todo$key[ch[j]], conditionMessage(r1)))
      res[[todo$key[ch[j]]]] <- r1
      cat(sprintf("  [%s] subject=%.3f feature=%.3f temporal=%.3f | signal=%.2f | %.0fs\n",
                  r1$key, r1$subject, r1$feature, r1$temporal, r1$det_mean, r1$secs))
    }
    utils::write.csv(do.call(rbind, res), .csv, row.names = FALSE)
    cat(sprintf("  -- checkpoint: %d runs saved (chunk %d/%d)\n",
                length(res), ci, length(chunks)))
    utils::flush.console()
  }
}

R_all <- do.call(rbind, res)
R_all <- R_all[order(R_all$n, R_all$noise, R_all$seed), , drop = FALSE]
row.names(R_all) <- NULL
utils::write.csv(R_all, .csv, row.names = FALSE)

# ---- clean timing pass ------------------------------------------------------
# Runs above may have been evaluated in parallel, which inflates each run's wall
# clock through contention. Re-time one seed per (n, noise) with nothing else
# running, and quote those as the runtime.
.tcsv <- file.path(.outdir, "01b_validate_larger_timing.csv")
timing <- if (file.exists(.tcsv)) utils::read.csv(.tcsv, stringsAsFactors = FALSE) else NULL
if (TIMING_PASS) {
  want <- expand.grid(noise = NOISE_LEVELS, n = N_SUBS)[, c("n", "noise")]
  have <- if (is.null(timing)) character(0) else sprintf("%d_%g", timing$n, timing$noise)
  need <- want[!(sprintf("%d_%g", want$n, want$noise) %in% have), , drop = FALSE]
  if (nrow(need) > 0) {
    cat(sprintf("\n== clean timing pass: %d sequential run(s) ==\n", nrow(need)))
    add <- do.call(rbind, lapply(seq_len(nrow(need)), function(i) {
      n_sub <- need$n[i]; nl <- need$noise[i]
      el <- fit_one(n_sub, nl, 1)$secs
      cat(sprintf("  n=%d noise=%g  %.1f s\n", n_sub, nl, el)); utils::flush.console()
      data.frame(n = n_sub, noise = nl, seed = 1, secs = el)
    }))
    timing <- if (is.null(timing)) add else rbind(timing, add)
    timing <- timing[order(timing$n, timing$noise), , drop = FALSE]
    utils::write.csv(timing, .tcsv, row.names = FALSE)
  }
}

# ---- summary ---------------------------------------------------------------
agg <- function(d, col) tapply(d[[col]], d$noise, mean)
sdv <- function(d, col) tapply(d[[col]], d$noise, stats::sd)

for (n_sub in N_SUBS) {
  d <- R_all[R_all$n == n_sub, , drop = FALSE]
  tsub <- if (is.null(timing)) NULL else timing[timing$n == n_sub, , drop = FALSE]
  summ <- data.frame(
    noise_sd = NOISE_LEVELS,
    subject  = sprintf("%.3f (+/-%.3f)", agg(d, "subject"),  sdv(d, "subject")),
    feature  = sprintf("%.3f (+/-%.3f)", agg(d, "feature"),  sdv(d, "feature")),
    temporal = sprintf("%.3f (+/-%.3f)", agg(d, "temporal"), sdv(d, "temporal")),
    signal_detected = sprintf("%.3f (+/-%.3f)", agg(d, "det_mean"), sdv(d, "det_mean")),
    secs_per_run = if (is.null(tsub)) sprintf("%.1f", agg(d, "secs"))
                   else sprintf("%.1f", tsub$secs[match(NOISE_LEVELS, tsub$noise)]),
    max_iters    = sprintf("%.0f", agg(d, "max_iters")),
    runs_hitting_maxiter = sprintf("%d/%d",
      tapply(d$n_at_maxiter > 0, d$noise, sum), tapply(d$seed, d$noise, length)))
  cat(sprintf("\n== n = %d, p = %s: recovery by noise level (mean +/- sd over seeds) ==\n",
              n_sub, paste(P_TARGET, collapse = "/")))
  print(summ, row.names = FALSE)
}
cat(sprintf("\nas-run total: %.1f core-hours over %d runs (%d evaluated in parallel)\n",
            sum(R_all$secs) / 3600, nrow(R_all), sum(R_all$parallel %in% TRUE)))
cat(sprintf("iteration cap (maxiter = %d): %d of %d runs had at least one component\n",
            MAXITER, sum(R_all$n_at_maxiter > 0), nrow(R_all)))
cat(sprintf("  stop on the cap rather than on the convergence threshold (epsilon = %g).\n",
            EPSILON))

# ---- figure: estimated curves (left) and combined recovery (right) ---------
# Left: the temporal loadings recovered at the smallest cohort and the lowest
# noise level. Right: recovery accuracy for both cohort sizes in one plot,
# colour for the loading type and line type for the sample size.
ref_n  <- min(N_SUBS)
ref_nl <- min(NOISE_LEVELS)
cat(sprintf("\nrefitting n = %d at noise %g for the curve panel\n", ref_n, ref_nl))
ref <- fit_one(ref_n, ref_nl, 1)
ref_rp <- evaluate_recovery(ref$sim, ref$fit)

.unit <- function(v) v / sqrt(sum(v^2))
Mv <- ref$sim$params$M; rv <- ref$sim$params$r; mtch <- ref_rp$match
curve_df <- do.call(rbind, lapply(1:Mv, function(m) do.call(rbind, lapply(1:rv, function(lh) {
  tg  <- ref$fit$time_Zeta[[m]]
  est <- .unit(ref$fit$Zeta_hat[[m]][, lh])
  tru <- .unit(ref$sim$truth$xi[[m]][[mtch[lh]]](tg))
  if (stats::cor(est, tru) < 0) est <- -est
  pan <- sprintf("Modality %d - PC %d (%s)", m, lh,
                 ref$sim$truth$shape_names[m, mtch[lh]])
  rbind(data.frame(time = tg, value = tru, series = "true",         panel = pan),
        data.frame(time = tg, value = est, series = "multiTEMPTED", panel = pan))
}))))
curve_df$series <- factor(curve_df$series, levels = c("true", "multiTEMPTED"))

p_curves <- ggplot(curve_df, aes(x = time, y = value, colour = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ panel, ncol = rv, scales = "free_y") +
  scale_colour_manual(values = c(true = "black", multiTEMPTED = "#4A90D9"), name = NULL) +
  scale_linetype_manual(values = c(true = "solid", multiTEMPTED = "22"), name = NULL) +
  labs(title = sprintf("Estimated vs true temporal loadings (n = %d, noise %g)",
                       ref_n, ref_nl),
       x = "Time", y = "Loading") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom")

rec_df <- do.call(rbind, lapply(N_SUBS, function(n_sub) {
  d <- R_all[R_all$n == n_sub, , drop = FALSE]
  do.call(rbind, lapply(names(LOADING_COLS), function(nm)
    data.frame(noise = NOISE_LEVELS, loading = nm, n = n_sub,
               mean = as.numeric(agg(d, nm)), sd = as.numeric(sdv(d, nm)))))
}))
rec_df$loading <- factor(rec_df$loading, levels = names(LOADING_COLS))
rec_df$cohort  <- factor(sprintf("n = %d", rec_df$n),
                         levels = sprintf("n = %d", N_SUBS))

p_rec <- ggplot(rec_df, aes(x = noise, y = mean, colour = loading, linetype = cohort)) +
  geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
                width = 0.04, linewidth = 0.3, alpha = 0.55, linetype = "solid") +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = cohort), size = 2, alpha = 0.85) +
  scale_x_log10(breaks = NOISE_LEVELS, labels = NOISE_LEVELS) +
  scale_colour_manual(values = LOADING_COLS, name = "Loading") +
  scale_linetype_manual(values = c("solid", "22"), name = "Cohort") +
  scale_shape_manual(values = c(16, 17), name = "Cohort") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = sprintf("Recovery accuracy (p = %s)", paste(P_TARGET, collapse = "/")),
       x = "Noise SD (log scale)", y = "|cor| with truth") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom", legend.box = "vertical")

p_fig <- (p_curves | p_rec) + plot_layout(widths = c(1.6, 1)) +
  plot_annotation(tag_levels = "A")

.pdf <- file.path(.outdir, "01b_validate_larger.pdf")
ggsave(.pdf, p_fig, width = 14, height = 6.4, bg = "white")
cat(sprintf("\n  figure written to %s\n  runs written to   %s\n",
            normalizePath(.pdf, mustWork = FALSE), normalizePath(.csv, mustWork = FALSE)))
