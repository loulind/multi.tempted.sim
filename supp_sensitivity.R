# ============================================================================
# supp_sensitivity.R
#
# Sensitivity of the iPOP analysis to two tuning choices: the modality weights
# w_m and the kernel ridge penalty C_K^(m).
#
# Everything else matches 02_ipop.R: same subject filter, same
# preprocessing, r = 3. The reference fit is equal weights with C_K = 1e-4.
#
# Weights are perturbed one modality at a time by each of WEIGHT_MULT, then
# renormalised to sum to the number of modalities, so that the effective
# smoothing is unchanged.
#
# Each perturbed fit is compared with the reference on the subject, feature and
# temporal loadings, and on how strongly the subject loadings track sex. The sex
# measure takes the largest absolute correlation across components, since a
# perturbation can move the sex signal from one component to another.
#
# Before comparing, the perturbed components are matched to the reference ones,
# so that a reordering is not mistaken for a loss of agreement. See compare().
#
# The weight sweep is reported as a figure, the smoothing sweep as a LaTeX
# table. Writes supp_output/supp_fig_weights.pdf, supp_tab_smoothing.tex and the
# underlying CSVs.
# ============================================================================

library(multi.tempted)
library(ggplot2)
library(patchwork)

# Output goes to supp_output/ beside this script, not getwd().
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


# ---- iPOP data and preprocessing (identical to 02_ipop.R) ----------
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

# ---- weight settings: local perturbation around equal weights --------------
# For modality k, multiply its weight by c and renormalise so the weights sum to
# M. c = 1 is the reference (equal weights).
WEIGHT_MULT <- c(2, 1.5, 0.75, 0.5)

perturb <- function(k, mult) {
  w <- rep(1, M); w[k] <- mult
  w * M / sum(w)
}

weight_settings <- c(
  list(list(label = "equal", modality = "equal", mult = 1, w = rep(1, M))),
  unlist(lapply(seq_len(M), function(k) lapply(WEIGHT_MULT, function(cm)
    list(label = sprintf("%s x%g", mods[k], cm), modality = mods[k],
         mult = cm, w = perturb(k, cm)))), recursive = FALSE))

wtab <- do.call(rbind, lapply(weight_settings, function(st)
  data.frame(setting = st$label, modality = st$modality, mult = st$mult,
             t(setNames(round(st$w, 4), mods)), check.names = FALSE,
             stringsAsFactors = FALSE)))
cat("\n-- modality weights (each row sums to M) --\n")
print(wtab, row.names = FALSE)

# ---- fits ------------------------------------------------------------------
# maxiter and epsilon as in the multitempted_all() defaults
fit_one <- function(w, smooth) {
  res <- multi_tempted_decomp(datlists = datlists_decomp, r = r, smooth = smooth,
                              maxiter = 100, epsilon = 1e-4, weights = unname(w))
  rownames(res$A_hat) <- ref_subjs
  res
}

cat(sprintf("\nreference fit: equal weights, C_K = %g\n", SMOOTH0))
ref   <- fit_one(rep(1, M), SMOOTH0)
sex_v <- sex[rownames(ref$A_hat)]
stopifnot(!any(is.na(sex_v)))

# All r! orderings of the components, so they can be matched without pulling in
# an assignment-problem dependency (r = 3 here, so six candidates).
.perms <- function(n) {
  if (n == 1L) return(list(1L))
  do.call(c, lapply(seq_len(n), function(i)
    lapply(.perms(n - 1L), function(p) c(i, p + (p >= i)))))
}

# |cor| with the reference fit, component by component, for each kind of
# loading; feature and temporal are averaged over modalities.
#
# Perturbed components are first matched to the reference by the permutation
# that maximises agreement of the subject loadings, and the same permutation is
# applied to the feature and temporal loadings. Unmatched values are also returned.
compare <- function(res) {
  # C[i, j] = |cor(reference component i, perturbed component j)|
  C     <- abs(stats::cor(ref$A_hat, res$A_hat))
  cands <- .perms(r)
  ord   <- cands[[which.max(vapply(cands, function(q)
    sum(diag(C[, q, drop = FALSE])), numeric(1)))]]

  A1 <- res$A_hat[, ord, drop = FALSE]
  B1 <- lapply(res$B_hat, function(b) b[, ord, drop = FALSE])
  Z1 <- lapply(res$Zeta_hat, function(z) z[, ord, drop = FALSE])

  per_comp <- function(X0, X1) abs(diag(stats::cor(X0, X1)))
  avg_mod  <- function(L0, L1) rowMeans(vapply(1:M, function(m)
    per_comp(L0[[m]], L1[[m]]), numeric(r)))

  # a perturbation can move the sex signal between components, so take the best
  sex_max <- max(vapply(1:r, function(l)
    abs(stats::cor(res$A_hat[, l], sex_v)), numeric(1)))

  list(order = ord, swapped = !all(ord == seq_len(r)), sex_max = sex_max,
       subject  = per_comp(ref$A_hat, A1),
       feature  = avg_mod(ref$B_hat, B1),
       temporal = avg_mod(ref$Zeta_hat, Z1),
       subject0  = per_comp(ref$A_hat, res$A_hat),
       feature0  = avg_mod(ref$B_hat, res$B_hat),
       temporal0 = avg_mod(ref$Zeta_hat, res$Zeta_hat))
}

summarise_fit <- function(res, extra) {
  cp <- compare(res)
  cbind(extra, data.frame(
    subject = mean(cp$subject), feature = mean(cp$feature), temporal = mean(cp$temporal),
    subject_PC1 = cp$subject[1], subject_PC2 = cp$subject[2], subject_PC3 = cp$subject[3],
    feature_PC1 = cp$feature[1], temporal_PC1 = cp$temporal[1],
    sex_max = cp$sex_max,
    order = paste(cp$order, collapse = "-"), swapped = cp$swapped,
    subject_unmatched = mean(cp$subject0), feature_unmatched = mean(cp$feature0),
    temporal_unmatched = mean(cp$temporal0), stringsAsFactors = FALSE))
}

# ---- weight perturbation ----------------------------------------------------
wrows <- lapply(weight_settings, function(st)
  summarise_fit(fit_one(st$w, SMOOTH0),
                data.frame(setting = st$label, modality = st$modality,
                           mult = st$mult, stringsAsFactors = FALSE)))
wres <- do.call(rbind, wrows)

# ---- smoothing sweep --------------------------------------------------------
srows <- lapply(SMOOTHS, function(sm)
  summarise_fit(fit_one(rep(1, M), sm),
                data.frame(C_K = sm, stringsAsFactors = FALSE)))
sres <- do.call(rbind, srows)

sex_ref <- wres$sex_max[wres$setting == "equal"]

cat("\n-- weight perturbation: |cor| with the equal-weight fit --\n")
print(transform(wres[, c("setting", "mult", "subject", "feature", "temporal", "sex_max")],
                subject = round(subject, 4), feature = round(feature, 4),
                temporal = round(temporal, 4), sex_max = round(sex_max, 4)),
      row.names = FALSE)
cat("\n-- smoothing sweep: |cor| with the C_K = 1e-4 fit --\n")
print(transform(sres[, c("C_K", "subject", "feature", "temporal", "sex_max")],
                subject = round(subject, 4), feature = round(feature, 4),
                temporal = round(temporal, 4), sex_max = round(sex_max, 4)),
      row.names = FALSE)
cat("\n-- component matching: ordering chosen on the subject loadings --\n")
print(wres[, c("setting", "order", "swapped", "subject_unmatched", "subject")],
      row.names = FALSE, digits = 4)
cat(sprintf("fits whose components were reordered: %d of %d\n",
            sum(wres$swapped), nrow(wres)))
cat(sprintf("\nsex association at equal weights: %.4f\n", sex_ref))
cat(sprintf("weights:   loading agreement >= %.4f ; sex %.4f to %.4f\n",
            min(c(wres$subject, wres$feature, wres$temporal)),
            min(wres$sex_max), max(wres$sex_max)))
cat(sprintf("smoothing: loading agreement >= %.4f ; sex %.4f to %.4f\n",
            min(c(sres$subject, sres$feature, sres$temporal)),
            min(sres$sex_max), max(sres$sex_max)))

utils::write.csv(wres, file.path(.outdir, "supp_sensitivity_weights.csv"), row.names = FALSE)
utils::write.csv(wtab, file.path(.outdir, "supp_sensitivity_weight_values.csv"), row.names = FALSE)

# ---- smoothing sweep: supplementary table ----------------------------------
# written as a standalone LaTeX document
.tex <- file.path(.outdir, "supp_tab_smoothing.tex")
writeLines(c(
  "\\documentclass{article}",
  "\\usepackage{array}",
  "",
  "\\begin{document}",
  "",
  "\\begin{table}[ht]",
  "\\centering",
  paste0("\\caption{Sensitivity of the iPOP fit to the temporal smoothing ",
         "penalty $C_K$.}"),
  "\\label{tab:supp-smoothing}",
  "\\begin{tabular}{lcccc}",
  "\\hline",
  "$C_K$ & Subject & Feature & Temporal & Sex association \\\\",
  "\\hline",
  sprintf("$10^{%d}$%s & %.4f & %.4f & %.4f & %.4f \\\\",
          round(log10(sres$C_K)),
          ifelse(sres$C_K == SMOOTH0, " (reference)", ""),
          sres$subject, sres$feature, sres$temporal, sres$sex_max),
  "\\hline",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\end{document}"), .tex)


# ---- figures ----------------------------------------------------------------
THEME_MS <- theme_bw(base_size = 12)
LOAD_COLS <- c(subject = "#4A90D9", feature = "#E05C5C", temporal = "#5CB85C")
MULTS <- sort(unique(wres$mult))

long_loadings <- function(d, xcol) {
  do.call(rbind, lapply(names(LOAD_COLS), function(nm)
    data.frame(x = d[[xcol]], loading = nm, value = d[[nm]], stringsAsFactors = FALSE)))
}

# weights: one facet per perturbed modality, reference (mult = 1) shared
wl <- do.call(rbind, lapply(mods, function(md) {
  d <- wres[wres$modality %in% c(md, "equal"), , drop = FALSE]
  d$mult[d$modality == "equal"] <- 1
  out <- long_loadings(d, "mult"); out$modality <- md; out
}))
wl$loading  <- factor(wl$loading, levels = names(LOAD_COLS))
wl$modality <- factor(wl$modality, levels = mods)

# horizontal dodge (in log2 units) to separate overlapping series
.wdodge <- position_dodge(width = 0.16)

p_w1 <- ggplot(wl, aes(x = x, y = value, colour = loading)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7, position = .wdodge) +
  geom_point(size = 2, alpha = 0.85, position = .wdodge) +
  facet_wrap(~ modality, nrow = 1) +
  scale_x_continuous(transform = "log2", breaks = MULTS, labels = MULTS) +
  scale_colour_manual(values = LOAD_COLS, name = "Loading") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Correlation with the equal-weight fit loadings (iPOP)",
       x = "Weight multiplier (log2 scale)", y = "|cor|") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom")

# faceted like panel A, with the equal-weight fit shown at mult = 1 in every facet
ws <- wres[wres$modality != "equal", c("modality", "mult", "sex_max")]
ws <- rbind(ws, data.frame(modality = mods, mult = 1, sex_max = sex_ref))
ws$modality <- factor(ws$modality, levels = mods)
ws <- ws[order(ws$modality, ws$mult), ]

p_w2 <- ggplot(ws, aes(x = mult, y = sex_max)) +
  geom_hline(yintercept = sex_ref, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7, colour = "#4A90D9") +
  geom_point(size = 2, alpha = 0.85, colour = "#4A90D9") +
  facet_wrap(~ modality, nrow = 1) +
  scale_x_continuous(transform = "log2", breaks = MULTS, labels = MULTS) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Association of the subject loadings with sex (iPOP)",
       subtitle = sprintf("dashed line: equal weights (%.2f)", sex_ref),
       x = "Weight multiplier (log2 scale)", y = "max |cor(component, sex)|") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = element_text(hjust = 0.5),
                   legend.position = "none")

ggsave(file.path(.outdir, "supp_fig_weights.pdf"),
       (p_w1 / p_w2) + plot_layout(heights = c(1.1, 1)) +
         plot_annotation(tag_levels = "A"),
       width = 12, height = 8.5, bg = "white")

cat(sprintf("\n  figure written to %s\n  table written to  %s\n",
            normalizePath(file.path(.outdir, "supp_fig_weights.pdf"), mustWork = FALSE),
            normalizePath(.tex, mustWork = FALSE)))
