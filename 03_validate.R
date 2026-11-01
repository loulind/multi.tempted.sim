# ============================================================================
# 03_validate.R
#
# Validation of multiTEMPTED on data generated from the model it assumes (each
# (modality, component) pair has its own temporal shape, and each subject and
# modality its own sampling times), at feature dimensions closer to real omics
# data: each modality has 20 signal features plus features with zero loadings in
# every component, giving p = 500, 100 and 30. Only a small
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
# Writes output/03_validate.pdf. Runs are flushed to
# output/03_validate_runs.csv every CKPT_EVERY runs and the script
# resumes at the first missing run, so a long sweep can be interrupted safely.
#
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
  warning("could not locate 03_validate.R; using the working directory",
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

# ---- data generator and recovery metric -------------------------------------
# (moved here verbatim from the former 01a_validate.R, which this script
# used to source at run time)

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
# Cohort size at which this parameterisation matches the original one, where A
# had unit norm and lambda was undivided.
A_NORM_REF <- 30

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
    # Divided by sqrt(A_NORM_REF) to pair with subject loadings of L2 norm
    # sqrt(n) below: together these keep the per-subject signal independent of
    # n, and reproduce the original parameterisation exactly at n = A_NORM_REF.
    lam_vec <- lambda_max * lambda_decay^(1:r - 1) / sqrt(A_NORM_REF)
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
  # A is scaled to L2 norm sqrt(n) rather than 1, so a single subject's loading
  # does not shrink as the cohort grows. B keeps unit norm.
  A_true <- .make_loadings(n, r, orthonormal) * sqrt(n)
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
.csv <- file.path(.outdir, "03_validate_runs.csv")
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
.tcsv <- file.path(.outdir, "03_validate_timing.csv")
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
  mod  <- sprintf("Modality %d (p = %d)", m, P_TARGET[m])
  comp <- sprintf("Component %d", lh)
  rbind(data.frame(time = tg, value = tru, series = "true",
                   modality = mod, component = comp),
        data.frame(time = tg, value = est, series = "multiTEMPTED",
                   modality = mod, component = comp))
}))))
curve_df$series    <- factor(curve_df$series, levels = c("true", "multiTEMPTED"))
curve_df$modality  <- factor(curve_df$modality,
                             levels = sprintf("Modality %d (p = %d)", 1:Mv, P_TARGET[1:Mv]))
curve_df$component <- factor(curve_df$component, levels = sprintf("Component %d", 1:rv))

p_curves <- ggplot(curve_df, aes(x = time, y = value, colour = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  facet_grid(modality ~ component) +
  scale_colour_manual(values = c(true = "black", multiTEMPTED = "#4A90D9"), name = NULL) +
  scale_linetype_manual(values = c(true = "solid", multiTEMPTED = "22"), name = NULL) +
  labs(title = sprintf("Estimated vs true temporal loadings (n = %d, noise %g)",
                       ref_n, ref_nl),
       x = "Time", y = "Loading") +
  # Shared y range across rows, clipped rather than filtered so the lines stay
  # continuous at the boundary. Strip text is shrunk to fit the row labels,
  # which are the longest, inside their strips.
  coord_cartesian(ylim = c(-0.23, 0.23)) +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom",
                   strip.text = element_text(size = 9))

rec_df <- do.call(rbind, lapply(N_SUBS, function(n_sub) {
  d <- R_all[R_all$n == n_sub, , drop = FALSE]
  do.call(rbind, lapply(names(LOADING_COLS), function(nm)
    data.frame(noise = NOISE_LEVELS, loading = nm, n = n_sub,
               mean = as.numeric(agg(d, nm)), sd = as.numeric(sdv(d, nm)))))
}))
rec_df$loading <- factor(rec_df$loading, levels = names(LOADING_COLS))
rec_df$cohort  <- factor(sprintf("n = %d", rec_df$n),
                         levels = sprintf("n = %d", N_SUBS))

# A small horizontal dodge separates the six series; the width is in log10
# units because the x scale is logged, and the noise levels sit ~0.4 apart.
.dodge <- position_dodge(width = 0.22)
p_rec <- ggplot(rec_df, aes(x = noise, y = mean, colour = loading, linetype = cohort,
                            group = interaction(loading, cohort))) +
  geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
                width = 0.03, linewidth = 0.3, alpha = 0.55, linetype = "solid",
                position = .dodge) +
  geom_line(linewidth = 0.7, position = .dodge) +
  geom_point(aes(shape = cohort), size = 2, alpha = 0.85, position = .dodge) +
  scale_x_log10(breaks = NOISE_LEVELS, labels = NOISE_LEVELS) +
  scale_colour_manual(values = LOADING_COLS, name = "Loading") +
  scale_linetype_manual(values = c("solid", "22"), name = "Sample size") +
  scale_shape_manual(values = c(16, 17), name = "Sample size") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = sprintf("Recovery accuracy (p = %s)", paste(P_TARGET, collapse = "/")),
       x = "Noise SD (log scale)", y = "|cor| with truth") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom", legend.box = "vertical")

p_fig <- (p_curves | p_rec) + plot_layout(widths = c(1.6, 1)) +
  plot_annotation(tag_levels = "A")

.pdf <- file.path(.outdir, "03_validate.pdf")
ggsave(.pdf, p_fig, width = 14, height = 6.4, bg = "white")
cat(sprintf("\n  figure written to %s\n  runs written to   %s\n",
            normalizePath(.pdf, mustWork = FALSE), normalizePath(.csv, mustWork = FALSE)))
