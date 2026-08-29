# ============================================================================
# 01_validate.R   (storyline point 1)
#
# Validate multiTEMPTED: data are generated from the *exact* model the method
# assumes, and we check that it recovers the planted components -- and how that
# recovery degrades as noise grows.
#
# Model: For modality m = 1..M, subject i = 1..n, feature j = 1..p_m, and time
# t in T_mi (a subset of the continuous interval T):
#
#     Y_mijt = sum_{l=1}^r  lambda_ml * a_il * b_mjl * xi_ml(t)  +  Z_mijt
#
# ...with identifiability constraints (sum_i a_il^2 = sum_j b_mjl^2 =
# integral_T xi_ml(t)^2 dt = 1) and error Z_mijt ~ N(0, noise_sd^2). The subject
# loadings a_l are shared across modalities; b_ml, xi_ml, lambda_ml are
# modality-specific.
#
# The run at the bottom is a NOISE STUDY: several noise levels x many seeds
# (default 5 x 20 = 100 runs). It reports recovery accuracy (|cor| of the
# estimated subject / feature / temporal loadings with truth) per noise level
# and writes output/01_validate.pdf: one figure with the estimated-vs-true curves
# at the lowest noise level (A) beside the recovery-vs-noise plot (B).
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
  if (is.null(d)) {                                     # Rscript 01_validate.R
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

# ---- temporal loading function shapes -------
# Each shape is a function of u in [0, 1] (the generator maps the real time
# interval T onto [0, 1] before evaluating, so shapes are independent of the
# chosen time_range).
#
# By default the generator walks this library in order and hands each modality
# the next r shapes, so every (modality, component) pair gets a UNIQUE curve:
# with M = r = 3, mod1 = shapes 1-3, mod2 = shapes 4-6, mod3 = shapes 7-9.
#
# The order is deliberate. Each modality gets one monotone, one single-bend and
# one multi-bend curve. And the two affine pairs (concave/convex and
# m_shaped/w_shaped, where one is 1 - the other) are split across different
# modalities: two affinely dependent temporal loadings inside the SAME modality
# would be degenerate. Swap the whole thing via `temporal_funcs`.
temporal_shape_library <- list(
  # -- modality 1 --
  s_curve    = function(u) 1 / (1 + exp(-12 * (u - 0.5))), # S-shaped
  concave    = function(u) 1 - 4 * (u - 0.5)^2,            # single peak (cap)
  decreasing = function(u) exp(-2 * u),                    # strictly decreasing
  # -- modality 2 --
  increasing = function(u) exp(2 * u),                     # strictly increasing
  m_shaped   = function(u) abs(sin(2 * pi * u)),           # two peaks
  convex     = function(u) 4 * (u - 0.5)^2,                # single trough (cup)
  # -- modality 3 --
  plateau    = function(u) 1 - exp(-8 * u),                # fast rise, then flat
  w_shaped   = function(u) 1 - abs(sin(2 * pi * u)),       # two troughs
  sine       = function(u) sin(2 * pi * u)                 # one full oscillation
)


# ---- small helpers (preceding dot hides function from global env) -------

# Trapezoidal integral of y over the (sorted) grid x.
.trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

# Random loadings with the required normalization. `orthonormal = TRUE` gives
# orthonormal columns (cleanest identifiability); otherwise just unit-norm
# columns. Needs nr >= r.
.make_loadings <- function(nr, r, orthonormal) {
  X <- matrix(stats::rnorm(nr * r), nr, r)
  if (orthonormal) {
    qr.Q(qr(X))[, 1:r, drop = FALSE]
  } else {
    apply(X, 2, function(z) z / sqrt(sum(z^2)))
  }
}


# ------- GENERATE DATA -------

#' Generate synthetic data from the multiTEMPTED model.
#'
#' @param n Number of subjects (shared across modalities).
#' @param M Number of modalities.
#' @param p Features per modality: scalar or length-M vector.
#' @param r Number of latent components.
#' @param n_timepoints Sampled timepoints per subject: scalar or length-M vector.
#' @param time_range Length-2 numeric giving the continuous interval T.
#' @param temporal_funcs NULL to give every (modality, component) pair a unique
#'   shape from `temporal_shape_library`; or a list of r functions of u in [0,1]
#'   reused for every modality; or a length-M list of such lists to set each
#'   modality's shapes explicitly.
#' @param lambda Component scalings: NULL (auto), a length-r vector (same
#'   across modalities), or an M x r matrix.
#' @param lambda_max,lambda_decay Used when lambda = NULL:
#'   lambda_l = lambda_max * lambda_decay^(l - 1) so decreases wrt l=1..r.
#' @param noise_sd Gaussian noise SD: scalar or length-M vector.
#' @param sampling "uniform" (default) draws each subject's timepoints uniformly
#'   at random from time_range, so subjects (and modalities) are sampled at
#'   unaligned times. "grid" puts every subject on the same evenly spaced grid.
#' @param orthonormal  Orthonormalise the A and B loadings (TRUE) or only
#'   unit-normalise them (FALSE).
#' @param modality_names Optional length-M character vector.
#' @param seed Integer seed for full reproducibility.
#' @return A list with `featuretables`, `timepoints`, `subjectID` (prep for
#'   `multitempted_all()`), the original components under `truth`, and the
#'   settings under `params`.
generate_multitempted_data <- function(
    n = 30,
    M = 3,
    p = 20,
    r = 3,
    n_timepoints = 8,
    time_range = c(0, 1),
    temporal_funcs = NULL,
    lambda = NULL,
    lambda_max = 8,
    lambda_decay = 0.6,
    noise_sd = 0.1,
    sampling = c("uniform", "grid"),
    orthonormal = TRUE,
    modality_names = NULL,
    seed = 1) {
  
  sampling <- match.arg(sampling)
  set.seed(seed)
  
  # -- recycle per-modality scalars and validate --
  if (length(p) == 1)            p <- rep(p, M)
  if (length(n_timepoints) == 1) n_timepoints <- rep(n_timepoints, M)
  if (length(noise_sd) == 1)     noise_sd <- rep(noise_sd, M)
  if (length(p) != M)            stop("'p' must have length 1 or M.")
  if (length(n_timepoints) != M) stop("'n_timepoints' must have length 1 or M.")
  if (length(noise_sd) != M)     stop("'noise_sd' must have length 1 or M.")
  if (n < r)                     stop("Need n >= r subjects.")
  if (any(p < r))                stop("Each modality needs p >= r features.")
  if (any(n_timepoints < 2))     stop("Each subject needs >= 2 timepoints.")
  
  if (is.null(modality_names)) modality_names <- paste0("mod", 1:M)
  subj_names <- sprintf("subj%02d", 1:n)
  PCnames <- paste0("PC", 1:r)
  
  # -- Modality scaling (lambda) (M x r matrix) --
  if (is.null(lambda)) {
    lam_vec <- lambda_max * lambda_decay^(1:r - 1)
    Lambda_true <- matrix(lam_vec, M, r, byrow = TRUE)
  } else if (is.matrix(lambda)) {
    if (!all(dim(lambda) == c(M, r))) stop("'lambda' matrix must be M x r.")
    Lambda_true <- lambda
  } else {
    if (length(lambda) != r) stop("'lambda' vector must have length r.")
    Lambda_true <- matrix(lambda, M, r, byrow = TRUE)
  }
  dimnames(Lambda_true) <- list(modality_names, PCnames)
  
  # -- Shared subject loadings A (n x r) and feature loadings B_m (p_m x r) --
  A_true <- .make_loadings(n, r, orthonormal)
  dimnames(A_true) <- list(subj_names, PCnames)
  
  B_true <- lapply(1:M, function(m) {
    Bm <- .make_loadings(p[m], r, orthonormal)
    dimnames(Bm) <- list(sprintf("%s_feat%02d", modality_names[m], seq_len(p[m])),
                         PCnames)
    Bm
  })
  names(B_true) <- modality_names
  
  # -- Temporal loadings (Xi): an M-list of r-lists of L2-normalized functions --
  if (is.null(temporal_funcs)) {
    # Unique shape per (modality, component): modality m takes the next r shapes
    # from the library. Only repeats if M*r outruns the library.
    L <- length(temporal_shape_library)
    if (M * r > L)
      warning(sprintf("M*r = %d exceeds the %d shapes in the library, so shapes repeat.",
                      M * r, L))
    idx <- (1:(M * r) - 1) %% L + 1
    base_shapes <- lapply(1:M, function(m)
      temporal_shape_library[idx[((m - 1) * r + 1):(m * r)]])
  } else if (is.function(temporal_funcs[[1]])) {
    if (length(temporal_funcs) < r) stop("'temporal_funcs' needs at least r functions.")
    base_shapes <- rep(list(temporal_funcs[1:r]), M)
  } else {
    if (length(temporal_funcs) != M) stop("'temporal_funcs' must be length M (per modality).")
    base_shapes <- lapply(temporal_funcs, function(fl) {
      if (length(fl) < r) stop("Each modality's temporal_funcs needs at least r functions.")
      fl[1:r]
    })
  }

  # Record which shape landed in each (modality, component) slot, for labelling.
  shape_names <- matrix(NA_character_, M, r, dimnames = list(modality_names, PCnames))
  for (m in 1:M) {
    nm <- names(base_shapes[[m]])
    if (!is.null(nm)) shape_names[m, ] <- nm
  }

  a_end <- time_range[1]; b_end <- time_range[2]
  to_unit <- function(t) (t - a_end) / (b_end - a_end)
  fine_t  <- seq(a_end, b_end, length.out = 2001)
  # Turn a shape on [0,1] into a function on T with integral_T xi^2 dt = 1.
  l2_normalise <- function(shape) {
    g <- function(t) shape(to_unit(t))
    nrm <- sqrt(.trapz(fine_t, g(fine_t)^2))
    function(t) g(t) / nrm
  }
  xi_true <- lapply(base_shapes, function(fl) lapply(fl, l2_normalise))
  names(xi_true) <- modality_names
  
  # -- Sampling times (t in T_mi): list over modalities of list over subjects --
  draw_times <- function(m) {
    q <- n_timepoints[m]
    if (sampling == "uniform") {
      # q times drawn uniformly from T, independently per subject, so subjects
      # (and modalities) land on unaligned timepoints.
      lapply(1:n, function(i) sort(stats::runif(q, a_end, b_end)))
    } else {
      g <- seq(a_end, b_end, length.out = q)
      rep(list(g), n) # every subject on the same even grid
    }
  }
  times_by_mod <- lapply(1:M, draw_times)
  
  # -- assemble featuretables (samples x features) + metadata per modality --
  featuretables <- vector("list", M)
  timepoints <- vector("list", M)
  subjectID <- vector("list", M)
  
  for (m in 1:M) {
    per_subj_vals <- vector("list", n)
    tvec <- numeric(0); svec <- character(0)
    for (i in 1:n) {
      ts <- times_by_mod[[m]][[i]]
      q  <- length(ts)
      # xi_ml(t) on this subject's times: q x r
      Xi <- vapply(1:r, function(l) xi_true[[m]][[l]](ts), numeric(q))
      if (q == 1) Xi <- matrix(Xi, nrow = 1)
      coef <- Lambda_true[m, ] * A_true[i, ] # length r: lambda_ml * a_il
      signal <- B_true[[m]] %*% (t(Xi) * coef) # p_m x q rank-r signal
      noise  <- matrix(stats::rnorm(p[m] * q, sd = noise_sd[m]), nrow = p[m])
      per_subj_vals[[i]] <- t(signal + noise) # q x p_m (samples x features)
      tvec <- c(tvec, ts)
      svec <- c(svec, rep(subj_names[i], q))
    }
    ft <- do.call(rbind, per_subj_vals)
    colnames(ft) <- rownames(B_true[[m]])
    featuretables[[m]] <- ft
    timepoints[[m]] <- tvec
    subjectID[[m]] <- svec
  }
  names(featuretables) <- modality_names
  names(timepoints) <- modality_names
  names(subjectID) <- modality_names
  
  list(
    featuretables = featuretables,
    timepoints    = timepoints,
    subjectID     = subjectID,
    truth = list(A = A_true, B = B_true, Lambda = Lambda_true,
                 xi = xi_true, shape_names = shape_names,
                 times = times_by_mod, time_range = time_range),
    params = list(n = n, M = M, p = p, r = r, n_timepoints = n_timepoints,
                  noise_sd = noise_sd, sampling = sampling,
                  orthonormal = orthonormal, seed = seed,
                  modality_names = modality_names)
  )
}


# ------- RECOVERY -------

#' Compare a fitted decomposition against the planted ground truth.
#'
#' Loadings are only identified up to sign and (for temporal loadings) scale, so
#' all comparisons use absolute correlation. Recovered components are matched to
#' true components greedily by the shared subject loadings A.
#'
#' @param sim Output of `generate_multitempted_data()`.
#' @param fit Output of `multitempted_all()` / `multi_tempted_decomp()`.
#' @return A list with `per_component` (data frame), `cor_B`/`cor_Zeta`
#'   (component x modality matrices), the `match` permutation, and the final
#'   accumulated R^2 per modality.
evaluate_recovery <- function(sim, fit) {
  r <- sim$params$r
  M <- sim$params$M
  
  # greedy match: recovered component lh -> true component match[lh], by |cor(A)|
  cor_A_mat <- abs(stats::cor(fit$A_hat, sim$truth$A)) # r_hat x r_true
  match <- integer(r); used <- logical(r)
  for (lh in 1:r) {
    cand <- cor_A_mat[lh, ]; cand[used] <- -Inf
    j <- which.max(cand); match[lh] <- j; used[j] <- TRUE
  }
  
  cor_A <- vapply(1:r,
                  function(lh) abs(stats::cor(fit$A_hat[, lh], sim$truth$A[, match[lh]])),
                  numeric(1))
  
  cor_B <- vapply(1:M, function(m)
    vapply(1:r, function(lh)
      abs(stats::cor(fit$B_hat[[m]][, lh], sim$truth$B[[m]][, match[lh]])),
      numeric(1)),
    numeric(r)) # r x M
  
  cor_Zeta <- vapply(1:M, function(m) {
    tgrid <- fit$time_Zeta[[m]]
    vapply(1:r, function(lh) {
      xi_vals <- sim$truth$xi[[m]][[match[lh]]](tgrid)
      abs(stats::cor(fit$Zeta_hat[[m]][, lh], xi_vals))
    }, numeric(1))
  }, numeric(r)) # r x M
  if (r == 1) { cor_B <- matrix(cor_B, 1); cor_Zeta <- matrix(cor_Zeta, 1) }
  
  per_component <- data.frame(
    recovered   = paste0("PC", 1:r),
    true        = paste0("PC", match),
    cor_A       = round(cor_A, 4),
    cor_B_mean  = round(rowMeans(cor_B), 4),
    cor_B_min   = round(apply(cor_B, 1, min), 4),
    cor_Z_mean  = round(rowMeans(cor_Zeta), 4),
    cor_Z_min   = round(apply(cor_Zeta, 1, min), 4),
    row.names   = NULL
  )
  
  list(per_component = per_component,
       cor_B = cor_B, cor_Zeta = cor_Zeta, match = match,
       final_r2 = fit$accum_r_square[, r])
}


# ------ PLOT ESTIMATED OVER TRUE ------

#' Plot estimated vs true temporal loadings, one panel per modality x component.
#' @param file NULL to draw on the current device, or a path to save a PDF.
plot_temporal_recovery <- function(sim, fit, report = NULL, file = NULL) {
  if (is.null(report)) report <- evaluate_recovery(sim, fit)
  r <- sim$params$r; M <- sim$params$M; match <- report$match
  if (!is.null(file)) { grDevices::pdf(file, width = 2.6 * r, height = 2.4 * M); on.exit(grDevices::dev.off()) }
  op <- graphics::par(mfrow = c(M, r), mar = c(3, 3, 2, 1), mgp = c(1.6, 0.5, 0))
  on.exit(graphics::par(op), add = TRUE)
  unit <- function(v) v / sqrt(sum(v^2))
  for (m in 1:M) for (lh in 1:r) {
    tgrid <- fit$time_Zeta[[m]]
    est   <- unit(fit$Zeta_hat[[m]][, lh])
    tru   <- unit(sim$truth$xi[[m]][[match[lh]]](tgrid))
    if (stats::cor(est, tru) < 0) est <- -est # align sign for display
    yl <- range(c(est, tru))
    shape <- sim$truth$shape_names[m, match[lh]] # true curve in this slot
    plot(tgrid, tru, type = "l", lwd = 2, col = "black", ylim = yl,
         xlab = "time", ylab = "loading",
         main = sprintf("%s  PC%d (%s)", sim$params$modality_names[m], lh, shape))
    graphics::lines(tgrid, est, lwd = 2, lty = 2, col = "red")
    if (m == 1 && lh == 1)
      graphics::legend("topleft", c("true", "estimated"), lwd = 2,
                       lty = c(1, 2), col = c("black", "red"), bty = "n", cex = 0.8)
  }
}


# ============================================================================
# NOISE STUDY: how well multiTEMPTED recovers the truth as noise increases.
# For each of several noise levels we run many random seeds and record the
# recovery accuracy (|cor| of the estimated subject / feature / temporal
# loadings with the truth). Fixing the seed and only scaling noise_sd keeps the
# planted structure identical across noise levels -- a controlled experiment.
# ============================================================================

NOISE_LEVELS   <- c(0.1, 0.5, 1.0, 2.0, 4.0)   # 5 increasing noise SDs
N_SEEDS        <- 20                           # seeds per level (5 x 20 = 100 runs)
smooth_penalty <- 1e-4                         # RKHS penalty suited to unaligned sampling

cat(sprintf("== multiTEMPTED validation: %d noise levels x %d seeds = %d runs ==\n",
            length(NOISE_LEVELS), N_SEEDS, length(NOISE_LEVELS) * N_SEEDS))

res <- vector("list", length(NOISE_LEVELS) * N_SEEDS); i <- 0
low <- NULL                                    # keep the cleanest run for the overlay figure
for (nl in NOISE_LEVELS) {
  for (s in 1:N_SEEDS) {
    sim <- generate_multitempted_data(r = 3, M = 3, noise_sd = nl, seed = s)
    fit <- suppressMessages(multitempted_all(
      sim$featuretables, sim$timepoints, sim$subjectID,
      transforms = "none", do_ratio = FALSE, centralize = FALSE,
      smooth = smooth_penalty, r = sim$params$r))
    rp <- evaluate_recovery(sim, fit)
    i <- i + 1
    res[[i]] <- data.frame(noise = nl, seed = s,
                           subject  = mean(rp$per_component$cor_A),
                           feature  = mean(rp$per_component$cor_B_mean),
                           temporal = mean(rp$per_component$cor_Z_mean))
    if (nl == NOISE_LEVELS[1] && s == 1) low <- list(sim = sim, fit = fit, report = rp)
  }
  cat(sprintf("  noise_sd=%.2f done (%d seeds)\n", nl, N_SEEDS))
}
R_all <- do.call(rbind, res)

# summarise per noise level (mean over seeds, with sd)
mns <- function(col) tapply(R_all[[col]], R_all$noise, mean)
sds <- function(col) tapply(R_all[[col]], R_all$noise, stats::sd)
summ <- data.frame(
  noise_sd = NOISE_LEVELS,
  subject  = sprintf("%.3f (+/-%.3f)", mns("subject"),  sds("subject")),
  feature  = sprintf("%.3f (+/-%.3f)", mns("feature"),  sds("feature")),
  temporal = sprintf("%.3f (+/-%.3f)", mns("temporal"), sds("temporal")),
  row.names = NULL)
cat("\n== recovery accuracy (|cor| with truth) by noise level, mean +/- sd over seeds ==\n")
print(summ, row.names = FALSE)


# ------ PDF figures (ggplot2, white background, manuscript theme) ------
LOADING_COLS <- c(subject = "#4A90D9", feature = "#E05C5C", temporal = "#5CB85C")

rec_df <- do.call(rbind, lapply(names(LOADING_COLS), function(nm)
  data.frame(noise = NOISE_LEVELS, loading = nm,
             mean = as.numeric(mns(nm)), sd = as.numeric(sds(nm)))))
rec_df$loading <- factor(rec_df$loading, levels = names(LOADING_COLS))

p_recovery <- ggplot(rec_df, aes(x = noise, y = mean, colour = loading)) +
  geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
                width = 0.05, linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2, alpha = 0.8) +
  scale_x_log10(breaks = NOISE_LEVELS, labels = NOISE_LEVELS) +
  scale_colour_manual(values = LOADING_COLS, name = "Loading") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title    = "multiTEMPTED recovery decreases with noise",
       subtitle = sprintf("%d seeds per noise level; mean +/- sd", N_SEEDS),
       x = "Noise SD (log scale)", y = "Recovery accuracy (|cor| with truth)") +
  THEME_MS + theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = element_text(hjust = 0.5))

# left panel: estimated vs true curves at the lowest noise level
.unit <- function(v) v / sqrt(sum(v^2))
Mv <- low$sim$params$M; rv <- low$sim$params$r; mtch <- low$report$match
curve_df <- do.call(rbind, lapply(1:Mv, function(m) do.call(rbind, lapply(1:rv, function(lh) {
  tg  <- low$fit$time_Zeta[[m]]
  est <- .unit(low$fit$Zeta_hat[[m]][, lh])
  tru <- .unit(low$sim$truth$xi[[m]][[mtch[lh]]](tg))
  if (stats::cor(est, tru) < 0) est <- -est
  pan <- sprintf("Modality %d - PC %d (%s)", m, lh, low$sim$truth$shape_names[m, mtch[lh]])
  rbind(data.frame(time = tg, value = tru, series = "true",         panel = pan),
        data.frame(time = tg, value = est, series = "multiTEMPTED", panel = pan))
}))))
curve_df$series <- factor(curve_df$series, levels = c("true", "multiTEMPTED"))

p_curves <- ggplot(curve_df, aes(x = time, y = value, colour = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ panel, ncol = rv, scales = "free_y") +
  scale_colour_manual(values = c(true = "black", multiTEMPTED = "#4A90D9"), name = NULL) +
  scale_linetype_manual(values = c(true = "solid", multiTEMPTED = "22"), name = NULL) +
  labs(title    = "Estimated vs true temporal loadings",
       subtitle = sprintf("lowest noise level (SD = %g)", NOISE_LEVELS[1]),
       x = "Time", y = "Loading") +
  THEME_MS + theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = element_text(hjust = 0.5),
                   legend.position = "bottom")

# ---- one figure: the two panels side by side, each with its own bottom legend --
# patchwork rather than gridExtra: grid.arrange does not align panel regions
# across subplots, so the two plotting areas would not line up.
out_pdf <- file.path(.outdir, "01_validate.pdf")
.leg_bottom <- theme(legend.position = "bottom", legend.justification = "center")
p_fig <- ((p_curves + .leg_bottom) | (p_recovery + .leg_bottom)) +
  plot_layout(widths = c(1.6, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(out_pdf, p_fig, width = 14, height = 6, bg = "white")
cat(sprintf("\n  validation figure written to %s\n",
            normalizePath(out_pdf, mustWork = FALSE)))
