# ============================================================================
# supp_tab_01.R
#
# Sensitivity of the iPOP analysis to two tuning choices: the modality weights
# w_m and the kernel ridge penalty C_K^(m).
#
# Everything else matches 04_compare_ipop.R: same subject filter, same
# preprocessing, r = 3. The reference fit is equal weights with C_K = 1e-4.
#
# Weighting schemes, both computed on the tensors that actually enter the
# decomposition, i.e. AFTER rank-1 mean removal by svd_centralize():
#   equal       w_m = 1
#   frobenius   w_m proportional to 1 / ||Y_m||_F^2       (equalises total sum of squares)
#   mfa         w_m proportional to 1 / sigma_1(Y_m)^2    (equalises the leading direction,
#                                                          as in multiple factor analysis)
# Both schemes are rescaled to average 1 across modalities, because the penalty
# term is NOT multiplied by the weights: leaving them unscaled would change the
# effective smoothing and confound the two effects being tested.
#
# Writes supp_output/supp_tab_01.csv and supp_tab_01_weights.csv.
# ============================================================================

library(multi.tempted)

# Everything is written to <this script's folder>/supp_output, not to getwd().
.outdir <- local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { d <- dirname(normalizePath(of, mustWork = FALSE)); break }
  }
  if (is.null(d)) {
    a <- commandArgs(trailingOnly = FALSE)
    a <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(a)) d <- dirname(normalizePath(a[1], mustWork = FALSE))
  }
  if (is.null(d) && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) d <- dirname(normalizePath(p, mustWork = FALSE))
  }
  if (is.null(d)) {
    d <- getwd()
    warning("could not locate the script file; writing supp_output/ under ", d, call. = FALSE)
  }
  file.path(d, "supp_output")
})
dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("output directory: %s\n", .outdir))

r        <- 3
R_SVD    <- 1                                   # rank of the removed mean structure
SMOOTH0  <- 1e-4                                # reference penalty
SMOOTHS  <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)     # sweep, equal weights
set.seed(1)


# ---- iPOP data and preprocessing (identical to 04_compare_ipop.R) ----------
ip   <- multi.tempted::ipop
mods <- c("cytokine", "metabolome", "lipid", "protein")

sid_all  <- as.character(ip$meta_visit$SubjectID)
n_visits <- max(table(sid_all))
keep_row <- sid_all %in% names(which(table(sid_all) == n_visits))

featuretables <- lapply(mods, function(m) as.matrix(ip[[m]])[keep_row, , drop = FALSE])
names(featuretables) <- mods
M <- length(featuretables)

tp_map     <- c("1" = 0, "2" = 12, "3" = 25, "4" = 40, "5" = 70)
timepoint  <- unname(tp_map[as.character(ip$meta_visit$timepoint)])[keep_row]
timepoints <- rep(list(timepoint), M)
subjectID  <- rep(list(sid_all[keep_row]), M)

sex <- stats::setNames(as.integer(ip$meta_subj$Sex == "M"),
                       as.character(ip$meta_subj$subjectID))

cat(sprintf("iPOP: %d subjects, %d modalities, features = %s\n",
            length(unique(subjectID[[1]])), M,
            paste(sapply(featuretables, ncol), collapse = "/")))


# ---- format + centralise, mirroring multitempted_all() --------------------
datlists <- lapply(1:M, function(m)
  format_tempted(featuretable = featuretables[[m]], timepoint = timepoints[[m]],
                 subjectID = subjectID[[m]], threshold = 0.95, pseudo = NULL,
                 transform = "none"))          # iPOP is already log10
names(datlists) <- mods
ref_subjs <- names(datlists[[1]])
for (m in 1:M) datlists[[m]] <- datlists[[m]][ref_subjs]

datlists_decomp <- svd_centralize(datlists, r = R_SVD)$datlists

# mode-2 unfolding: p_m by (sum_i q_im), features as rows (row 1 holds the times)
unfold <- function(dl_m) {
  p_m <- nrow(dl_m[[1]]) - 1L
  do.call(cbind, lapply(dl_m, function(x) x[2:(p_m + 1), , drop = FALSE]))
}
unfolded <- lapply(datlists_decomp, unfold)

frob2 <- vapply(unfolded, function(X) sum(X^2), numeric(1))
sig12 <- vapply(unfolded, function(X) max(svd(X, nu = 0, nv = 0)$d)^2, numeric(1))

rescale <- function(w) w / mean(w)              # average 1 across modalities
WEIGHTS <- list(equal     = rescale(rep(1, M)),
                frobenius = rescale(1 / frob2),
                mfa       = rescale(1 / sig12))
for (nm in names(WEIGHTS)) names(WEIGHTS[[nm]]) <- mods

cat("\n-- modality weights (rescaled to average 1) --\n")
wtab <- data.frame(modality = mods,
                   p = sapply(featuretables, ncol),
                   frob_norm_sq = round(frob2, 1),
                   leading_sv_sq = round(sig12, 1),
                   w_equal = round(WEIGHTS$equal, 4),
                   w_frobenius = round(WEIGHTS$frobenius, 4),
                   w_mfa = round(WEIGHTS$mfa, 4), row.names = NULL)
print(wtab, row.names = FALSE)


# ---- fits ------------------------------------------------------------------
fit_one <- function(w, smooth) {
  # maxiter/epsilon must match multitempted_all()'s defaults, NOT
  # multi_tempted_decomp()'s own (maxiter = 20), or the reference fit stops
  # early: the iPOP components need 24-46 iterations to converge.
  res <- multi_tempted_decomp(datlists = datlists_decomp, r = r, smooth = smooth,
                              maxiter = 100, epsilon = 1e-4,
                              weights = unname(w))
  A <- res$A_hat
  rownames(A) <- ref_subjs
  A
}

cat(sprintf("\nreference fit: equal weights, C_K = %g\n", SMOOTH0))
A0 <- fit_one(WEIGHTS$equal, SMOOTH0)
sex_v <- sex[rownames(A0)]
stopifnot(!any(is.na(sex_v)))

settings <- c(
  lapply(names(WEIGHTS), function(nm) list(kind = "weights", label = nm,
                                           w = WEIGHTS[[nm]], smooth = SMOOTH0)),
  lapply(SMOOTHS, function(sm) list(kind = "smoothing", label = format(sm, scientific = TRUE),
                                    w = WEIGHTS$equal, smooth = sm)))

rows <- lapply(settings, function(st) {
  A <- fit_one(st$w, st$smooth)
  cors <- vapply(1:r, function(l) abs(stats::cor(A[, l], A0[, l])), numeric(1))
  data.frame(varied = st$kind, setting = st$label,
             weights = paste(sprintf("%.3f", st$w), collapse = "/"),
             C_K = st$smooth,
             cor_PC1 = cors[1], cor_PC2 = cors[2], cor_PC3 = cors[3],
             sex_cor_PC1 = abs(stats::cor(A[, 1], sex_v)),
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)

cat("\n-- Supplementary Table S1: sensitivity of the iPOP fit --\n")
cat("   cor_PCk = |cor| of component k's subject loadings with the reference fit\n")
cat("   sex_cor_PC1 = |cor| between the leading subject component and sex\n\n")
print(transform(tab, cor_PC1 = round(cor_PC1, 4), cor_PC2 = round(cor_PC2, 4),
                cor_PC3 = round(cor_PC3, 4), sex_cor_PC1 = round(sex_cor_PC1, 4)),
      row.names = FALSE)

cat(sprintf("\nminimum |cor| with the reference fit, leading component: %.4f\n", min(tab$cor_PC1)))
cat(sprintf("minimum |cor| with the reference fit, any component:      %.4f\n",
            min(c(tab$cor_PC1, tab$cor_PC2, tab$cor_PC3))))
cat(sprintf("sex association of the leading component: %.4f to %.4f\n",
            min(tab$sex_cor_PC1), max(tab$sex_cor_PC1)))

utils::write.csv(tab,  file.path(.outdir, "supp_tab_01.csv"), row.names = FALSE)
utils::write.csv(wtab, file.path(.outdir, "supp_tab_01_weights.csv"), row.names = FALSE)
cat(sprintf("\n  written to %s\n  and        %s\n",
            normalizePath(file.path(.outdir, "supp_tab_01.csv"), mustWork = FALSE),
            normalizePath(file.path(.outdir, "supp_tab_01_weights.csv"), mustWork = FALSE)))
