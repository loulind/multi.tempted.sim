# ============================================================================
# supp_ablation.R
#
# PURPOSE (reviewer 1, comment 5)
# Ablation of the modality-specific temporal loadings. Data are generated
# exactly as in 05_compare_aligned.R (shared curves across modalities,
# co-measured visits, block-structured subjects), but over a range of noise
# levels, and only two models are fitted:
#   * multiTEMPTED, with modality-specific temporal loadings;
#   * multiTEMPTED (shared trajectory), with ONE temporal loading per
#     component shared by all modalities. With co-measured visits, equal weights and a common
#     smoothing penalty, this is exactly multiTEMPTED applied to the
#     concatenated features of all modalities.
# MEFISTO is not run. At NOISE = 1 and seeds 1..100 the multiTEMPTED results
# reproduce those of 05_compare_aligned.R exactly. (Supplementary Figure S3)
#
# OUTPUT (supp_output/)
#   supp_ablation.csv        one row per (noise level, seed), both models
#   supp_fig_ablation.pdf    misclassification, temporal and subject loading
#                            recovery against noise
# Runs in about 5 minutes.
# ============================================================================

library(multi.tempted)
library(ggplot2)
library(patchwork)

NOISES  <- c(1, 1.5, 2, 2.5, 3, 4)
N_SEEDS <- 100
# Output goes to supp_output/ beside this script, NOT getwd().
.outdir <- local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {                    # source() / RStudio Source
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { d <- dirname(normalizePath(of, mustWork = FALSE)); break }
  }
  if (is.null(d)) {                                     # Rscript supp_ablation.R
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
  file.path(d, "supp_output")
})

dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)

# --- copied verbatim from 05_compare_aligned.R (data generator + scoring);
# NOISE is set per level in the loop below ---
SAMPLING <- "CLUSTERED"; TIME_JITTER <- 0.02; CLUSTER_BY <- "visit"
N <- 15; M <- 3; P <- 24; R <- 3; NT <- 8; LAMBDA <- c(8, 6, 4); OFFGROUP <- 0.05
shapes <- list(s_curve    = function(u) 1 / (1 + exp(-12 * (u - 0.5))),
               increasing = function(u) exp(2 * u),
               decreasing = function(u) exp(-2 * u))
.trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
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
  grid <- seq(0, 1, length.out = NT)
  .reflect01 <- function(u) { u <- abs(u); 1 - abs(1 - u) }
  .jitter <- function() sort(.reflect01(grid + stats::rnorm(NT, sd = TIME_JITTER)))
  visit <- if (SAMPLING == "CLUSTERED" && CLUSTER_BY == "visit") lapply(1:N, function(i) .jitter()) else NULL
  ft <- tp <- sid <- vector("list", M)
  for (m in 1:M) {
    per <- vector("list", N); tv <- numeric(0); sv <- character(0)
    for (i in 1:N) {
      ts <- switch(SAMPLING, ALIGNED = grid,
                   CLUSTERED = if (CLUSTER_BY == "visit") visit[[i]] else .jitter(),
                   UNALIGNED = sort(stats::runif(NT, 0, 1)))
      Xi <- vapply(1:R, function(l) xi[[l]](ts), numeric(NT))
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
.perms <- function(v) if (length(v) == 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(.perms(v[-i]), function(p) c(v[i], p))))
misclass <- function(load, grp, r) {
  pred <- apply(abs(load), 1, which.max); pm <- .perms(1:r)
  bp <- pm[[which.max(vapply(pm, function(q) mean(q[pred] == grp), numeric(1)))]]
  list(err = 1 - mean(bp[pred] == grp), comp_of_group = vapply(1:r, function(l) which(bp == l), integer(1)))
}
.acor <- function(x, y) { v <- suppressWarnings(abs(stats::cor(x, y))); if (is.na(v)) 0 else v }

# --- new: fit full and restricted (shared-trajectory) models on one dataset ---
fit_pair <- function(seed) {
  sim <- gen_block(seed); grp <- sim$truth$subj_group
  q <- function(...) suppressMessages(multitempted_all(...))
  mt <- q(sim$featuretables, sim$timepoints, sim$subjectID,
          transforms = "none", do_ratio = FALSE, centralize = FALSE, smooth = 1e-4, r = R)
  # multiTEMPTED (shared trajectory): one temporal loading per component shared by all modalities.
  # With co-measured visits, equal weights and a common penalty, this is exactly
  # multiTEMPTED applied to the concatenated features of all modalities.
  stopifnot(all(vapply(2:M, function(m) identical(sim$timepoints[[m]], sim$timepoints[[1]]) &&
                                         identical(sim$subjectID[[m]], sim$subjectID[[1]]), logical(1))))
  rt <- q(list(concat = do.call(cbind, sim$featuretables)), list(concat = sim$timepoints[[1]]),
          list(concat = sim$subjectID[[1]]),
          transforms = "none", do_ratio = FALSE, centralize = FALSE, smooth = 1e-4, r = R)
  mmt <- misclass(mt$A_hat, grp, R); mrt <- misclass(rt$A_hat, grp, R)
  gmt <- mt$time_Zeta[[1]]; grt <- rt$time_Zeta[[1]]
  mt_fe <- mean(vapply(1:R, function(l) mean(vapply(1:M, function(m)
    .acor(mt$Zeta_hat[[m]][, mmt$comp_of_group[l]], sim$truth$xi[[l]](gmt)), numeric(1))), numeric(1)))
  rt_fe <- mean(vapply(1:R, function(l)
    .acor(rt$Zeta_hat[[1]][, mrt$comp_of_group[l]], sim$truth$xi[[l]](grt)), numeric(1)))
  # subject loadings, continuous: |cor| with the true loading column, after matching
  mt_sub <- mean(vapply(1:R, function(l) .acor(mt$A_hat[, mmt$comp_of_group[l]], sim$truth$A[, l]), numeric(1)))
  rt_sub <- mean(vapply(1:R, function(l) .acor(rt$A_hat[, mrt$comp_of_group[l]], sim$truth$A[, l]), numeric(1)))
  data.frame(seed = seed, noise = NOISE, mt_miscl = mmt$err, rt_miscl = mrt$err,
             mt_fe = mt_fe, rt_fe = rt_fe, mt_sub = mt_sub, rt_sub = rt_sub)
}

# ============================================================================
# run
# ============================================================================
out <- list()
for (NOISE in NOISES) {
  cat(sprintf("noise %g: %d seeds\n", NOISE, N_SEEDS)); utils::flush.console()
  out[[as.character(NOISE)]] <- do.call(rbind, lapply(seq_len(N_SEEDS), fit_pair))
}
res <- do.call(rbind, out)
utils::write.csv(res, file.path(.outdir, "supp_ablation.csv"), row.names = FALSE)

r <- res
# ============================================================================
# figure: line summaries over noise; panel titles and y-axis metrics follow
# 05_compare_aligned.pdf (Figure 5), plus subject loading recovery
# ============================================================================
THEME_MS    <- theme_bw(base_size = 12)
METHOD_COLS <- c(multiTEMPTED = "#4A90D9", "multiTEMPTED (shared trajectory)" = "#D98E2B")
long <- rbind(
  data.frame(noise = r$noise, method = "multiTEMPTED",      miscl = r$mt_miscl, fe = r$mt_fe, sub = r$mt_sub),
  data.frame(noise = r$noise, method = "multiTEMPTED (shared trajectory)", miscl = r$rt_miscl, fe = r$rt_fe, sub = r$rt_sub))
long$method <- factor(long$method, levels = names(METHOD_COLS))
.by <- split(long, list(long$noise, long$method), drop = TRUE)
# misclassification: mean +/- SE (median is 0 for most settings);
# |cor| metrics: median with interquartile range
mis <- do.call(rbind, lapply(.by, function(d) { m <- mean(d$miscl); se <- sd(d$miscl) / sqrt(nrow(d))
  data.frame(noise = d$noise[1], method = d$method[1], y = m, lo = m - se, hi = m + se) }))
.med <- function(v) do.call(rbind, lapply(.by, function(d)
  data.frame(noise = d$noise[1], method = d$method[1], y = median(d[[v]]),
             lo = unname(quantile(d[[v]], .25)), hi = unname(quantile(d[[v]], .75)))))
pd <- position_dodge(width = 0.12)
.linepanel <- function(df, ylab, title, ylim) {
  ggplot(df, aes(noise, y, colour = method, shape = method, linetype = method)) +
    geom_line(linewidth = 0.6, position = pd) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5, linetype = "solid", position = pd) +
    geom_point(size = 2.4, position = pd) +
    scale_colour_manual(values = METHOD_COLS, name = NULL) +
    scale_shape_manual(values = c(16, 17), name = NULL) +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    scale_x_continuous(breaks = sort(unique(df$noise))) + coord_cartesian(ylim = ylim) +
    labs(title = title, x = "Noise SD", y = ylab) +
    THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
p_box <- (.linepanel(mis,          "Misclassification error",   "Group loading",           c(0, 0.5)) |
          .linepanel(.med("fe"),   "Function estimation |cor|", "Temporal curve recovery", c(0.4, 1)) |
          .linepanel(.med("sub"),  "Subject loading |cor|",     "Subject loading recovery", c(0.2, 1))) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(.outdir, "supp_fig_ablation.pdf"), p_box, width = 12, height = 4.4, bg = "white")
