# ============================================================================
# simulation.R
#
# A fully synthetic simulation analysis for multiTEMPTED. 
# Data are gen'd from the *exact* model assumed in the method, to see if 
# the method can recover the original components.
#
# Model: For modality m = 1..M, subject i = 1..n, feature j = 1..p_m, and time
# t in T_mi (a subset of the continuous interval T):
#
#     Y_mijt = sum_{l=1}^r  lambda_ml * a_il * b_mjl * xi_ml(t)  +  Z_mijt
#
# ...with identifiability constraints, for every l and m,
#
#     sum_i a_il^2 = 1,   sum_j b_mjl^2 = 1,   integral_T xi_ml(t)^2 dt = 1,
#
# ....and error Z_mijt ~ N(0, noise_sd^2). The subject loadings a_l are
# shared across modalities; b_ml, xi_ml, and lambda_ml are modality specific.
#
# The generator returns `featuretables`, `timepoints`, and `subjectID` in
# exactly the shape `multitempted_all()` expects, plus the original
# components under `$truth` so recovery can be checked.
#
# Everything is driven by generate_multitempted_data(); see its arguments to
# modify components.
# ============================================================================

# pak::pkg_install("loulind/multi.tempted")
library(multi.tempted)

# ---- default temporal loading shapes -------
# Each shape is a function of u in [0, 1] (the generator maps the real time
# interval T onto [0, 1] before evaluating, so shapes are independent of the
# chosen time_range). Component 1 is S-shaped, 2 strictly increasing, 3
# strictly decreasing. Swap these via `temporal_funcs`.
temporal_shapes_default <- list(
  s_curve = function(u) 1 / (1 + exp(-12 * (u - 0.5))), # S-shaped on 0 to 1
  increasing = function(u) exp(2 * u), # strictly incr on 0 to 1
  decreasing = function(u) exp(-2 * u) # strictly decr on 0 to 1
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


# ============================================================================
# GENERATE DATA
# ============================================================================

#' Generate synthetic data from the multiTEMPTED model.
#'
#' @param n Number of subjects (shared across modalities).
#' @param M Number of modalities.
#' @param p Features per modality: scalar or length-M vector.
#' @param r Number of latent components.
#' @param n_timepoints Sampled timepoints per subject: scalar or length-M vector.
#' @param time_range Length-2 numeric giving the continuous interval T.
#' @param temporal_funcs NULL for the defaults (S / increasing / decreasing);
#'   or a list of r functions of u in [0,1] used for every modality; or a
#'   length-M list of such lists to give each modality its own shapes.
#' @param lambda Component scalings: NULL (auto), a length-r vector (same
#'   across modalities), or an M x r matrix.
#' @param lambda_max,lambda_decay Used when lambda = NULL:
#'   lambda_l = lambda_max * lambda_decay^(l - 1) so decreases wrt l=1..r.
#' @param noise_sd Gaussian noise SD: scalar or length-M vector.
#' @param sampling "grid" (evenly spaced times, shared by all subjects) or
#'   "random" (irregular times drawn per subject).
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
    sampling = c("grid", "random"),
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
    if (length(temporal_shapes_default) < r)
      stop("Provide 'temporal_funcs': fewer than r default shapes are defined.")
    base_shapes <- rep(list(temporal_shapes_default[1:r]), M)
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
    if (sampling == "grid") {
      g <- seq(a_end, b_end, length.out = q)
      rep(list(g), n) # every subject same grid
    } else {
      lapply(1:n, function(i)
        sort(unique(stats::runif(q, a_end, b_end)))) # irregular per subject
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
                 xi = xi_true, times = times_by_mod, time_range = time_range),
    params = list(n = n, M = M, p = p, r = r, n_timepoints = n_timepoints,
                  noise_sd = noise_sd, sampling = sampling,
                  orthonormal = orthonormal, seed = seed,
                  modality_names = modality_names)
  )
}


# ============================================================================
# RECOVERY ANALYSIS
# ============================================================================

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


# ============================================================================
# OPTIONAL: overlay true vs estimated temporal loadings
# ============================================================================

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
    plot(tgrid, tru, type = "l", lwd = 2, col = "black", ylim = yl,
         xlab = "time", ylab = "loading",
         main = sprintf("%s  PC%d", sim$params$modality_names[m], lh))
    graphics::lines(tgrid, est, lwd = 2, lty = 2, col = "red")
    if (m == 1 && lh == 1)
      graphics::legend("topleft", c("true", "estimated"), lwd = 2,
                       lty = c(1, 2), col = c("black", "red"), bty = "n", cex = 0.8)
  }
}


# ============================================================================
# DEMO (runs only when the file is executed directly, e.g. `Rscript
# synthetic_analysis.R`; skipped when the file is source()d to load functions)
# ============================================================================

cat("== Generating synthetic multiTEMPTED data ==\n")
sim <- generate_multitempted_data(r = 3, M=4, seed = 1)
pr <- sim$params
cat(sprintf("  M=%d modalities, n=%d subjects, p=%s features, %s timepoints/subj\n",
            pr$M, pr$n, paste(pr$p, collapse = "/"),
            paste(pr$n_timepoints, collapse = "/")))
cat("  temporal shapes: PC1 S-curve, PC2 increasing, PC3 decreasing; ")
cat(sprintf("noise_sd=%s, sampling=%s\n",
            paste(pr$noise_sd, collapse = "/"), pr$sampling))

cat("\n-- featuretables format (what multitempted_all consumes) --\n")
ft1 <- sim$featuretables[[1]]
cat(sprintf("  featuretables[[1]] '%s': %d samples x %d features\n",
            pr$modality_names[1], nrow(ft1), ncol(ft1)))
cat(sprintf("  timepoints[[1]] length %d, subjectID[[1]] length %d (%d unique)\n",
            length(sim$timepoints[[1]]), length(sim$subjectID[[1]]),
            length(unique(sim$subjectID[[1]]))))

cat("\n== Fitting multitempted_all (centralize=FALSE: data ARE the model) ==\n")
fit <- multitempted_all(
  featuretables = sim$featuretables,
  timepoints    = sim$timepoints,
  subjectID     = sim$subjectID,
  transforms    = "none", # already-scaled real values
  do_ratio      = FALSE, # not counts
  centralize    = FALSE, # no separate mean term in the model
  r             = pr$r)

cat("\n== Recovery (absolute correlations; 1.000 = perfect) ==\n")
report <- evaluate_recovery(sim, fit)
print(report$per_component, row.names = FALSE)
cat(sprintf("\n  final accumulated R^2 per modality: %s\n",
            paste(sprintf("%.3f", report$final_r2), collapse = ", ")))
cat(sprintf("  min subject-loading correlation across components: %.4f\n",
            min(report$per_component$cor_A)))

# save a visual overlay next to the script
out_pdf <- "recovery.pdf"
plot_temporal_recovery(sim, fit, report, file = out_pdf)
cat(sprintf("  temporal-loading overlay written to %s\n", out_pdf))
