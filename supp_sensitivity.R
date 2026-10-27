# ============================================================================
# supp_sensitivity.R
#
# Sensitivity of the iPOP analysis to two tuning choices: the modality weights
# w_m and the kernel ridge penalty C_K^(m).
#
# Everything else matches 04_compare_ipop.R: same subject filter, same
# preprocessing, r = 3. The reference fit is equal weights with C_K = 1e-4.
#
# Weights are perturbed locally around equal weighting: one modality at a time
# is multiplied by each of WEIGHT_MULT, then all weights are renormalised to sum
# to the number of modalities. Renormalising matters because the penalty term is
# not multiplied by the weights, so leaving them unscaled would also change the
# effective smoothing and confound the two effects.
#
# Each perturbed fit is compared with the reference on the subject, feature and
# temporal loadings, and on how strongly the subject loadings track sex. The sex
# measure takes the largest absolute correlation across components, since a
# perturbation can move the sex signal from one component to another.
#
# Writes supp_output/supp_fig_weights.pdf, supp_fig_smoothing.pdf and the
# underlying CSVs.
# ============================================================================

library(multi.tempted)
library(ggplot2)
library(patchwork)

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
# maxiter/epsilon must match multitempted_all()'s defaults, NOT
# multi_tempted_decomp()'s own (maxiter = 20), or the reference fit stops early:
# the iPOP components need 24-46 iterations to converge.
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

# |cor| with the reference, component by component, for each kind of loading.
# Feature and temporal loadings are averaged over modalities.
compare <- function(res) {
  subj <- vapply(1:r, function(l) abs(stats::cor(res$A_hat[, l], ref$A_hat[, l])), numeric(1))
  feat <- vapply(1:r, function(l) mean(vapply(1:M, function(m)
    abs(stats::cor(res$B_hat[[m]][, l], ref$B_hat[[m]][, l])), numeric(1))), numeric(1))
  temp <- vapply(1:r, function(l) mean(vapply(1:M, function(m)
    abs(stats::cor(res$Zeta_hat[[m]][, l], ref$Zeta_hat[[m]][, l])), numeric(1))), numeric(1))
  # a perturbation can move the sex signal between components, so take the best
  sex_max <- max(vapply(1:r, function(l)
    abs(stats::cor(res$A_hat[, l], sex_v)), numeric(1)))
  list(subject = subj, feature = feat, temporal = temp, sex_max = sex_max)
}

summarise_fit <- function(res, extra) {
  cp <- compare(res)
  cbind(extra, data.frame(
    subject = mean(cp$subject), feature = mean(cp$feature), temporal = mean(cp$temporal),
    subject_PC1 = cp$subject[1], subject_PC2 = cp$subject[2], subject_PC3 = cp$subject[3],
    feature_PC1 = cp$feature[1], temporal_PC1 = cp$temporal[1],
    sex_max = cp$sex_max, stringsAsFactors = FALSE))
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
cat(sprintf("\nsex association at equal weights: %.4f\n", sex_ref))
cat(sprintf("weights:   loading agreement >= %.4f ; sex %.4f to %.4f\n",
            min(c(wres$subject, wres$feature, wres$temporal)),
            min(wres$sex_max), max(wres$sex_max)))
cat(sprintf("smoothing: loading agreement >= %.4f ; sex %.4f to %.4f\n",
            min(c(sres$subject, sres$feature, sres$temporal)),
            min(sres$sex_max), max(sres$sex_max)))

utils::write.csv(wres, file.path(.outdir, "supp_sensitivity_weights.csv"), row.names = FALSE)
utils::write.csv(sres, file.path(.outdir, "supp_sensitivity_smoothing.csv"), row.names = FALSE)
utils::write.csv(wtab, file.path(.outdir, "supp_sensitivity_weight_values.csv"), row.names = FALSE)


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

p_w1 <- ggplot(wl, aes(x = x, y = value, colour = loading)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7) + geom_point(size = 2, alpha = 0.85) +
  facet_wrap(~ modality, nrow = 1) +
  scale_x_continuous(breaks = MULTS, labels = MULTS) +
  scale_colour_manual(values = LOAD_COLS, name = "Loading") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Correlation with the equal-weight fit loadings (iPOP)",
       x = "Weight multiplier", y = "|cor|") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom")

# Faceted like panel A, so the two rows share x breaks and panel widths. The
# equal-weight fit is carried into every facet at mult = 1, as it is in panel
# A, so each line runs through the reference rather than jumping across it.
ws <- wres[wres$modality != "equal", c("modality", "mult", "sex_max")]
ws <- rbind(ws, data.frame(modality = mods, mult = 1, sex_max = sex_ref))
ws$modality <- factor(ws$modality, levels = mods)
ws <- ws[order(ws$modality, ws$mult), ]

p_w2 <- ggplot(ws, aes(x = mult, y = sex_max)) +
  geom_hline(yintercept = sex_ref, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7, colour = "#4A90D9") +
  geom_point(size = 2, alpha = 0.85, colour = "#4A90D9") +
  facet_wrap(~ modality, nrow = 1) +
  scale_x_continuous(breaks = MULTS, labels = MULTS) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Association of the subject loadings with sex (iPOP)",
       subtitle = sprintf("dashed line: equal weights (%.2f)", sex_ref),
       x = "Weight multiplier", y = "max |cor(component, sex)|") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = element_text(hjust = 0.5),
                   legend.position = "none")

ggsave(file.path(.outdir, "supp_fig_weights.pdf"),
       (p_w1 / p_w2) + plot_layout(heights = c(1.1, 1)) +
         plot_annotation(tag_levels = "A"),
       width = 12, height = 8.5, bg = "white")

# smoothing
sl <- long_loadings(sres, "C_K")
sl$loading <- factor(sl$loading, levels = names(LOAD_COLS))
p_s1 <- ggplot(sl, aes(x = x, y = value, colour = loading)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.7) + geom_point(size = 2, alpha = 0.85) +
  scale_x_log10(breaks = SMOOTHS, labels = format(SMOOTHS, scientific = TRUE)) +
  scale_colour_manual(values = LOAD_COLS, name = "Loading") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Agreement with the reference fit",
       x = expression(C[K]~"(log scale)"), y = "|cor| with reference") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom")

p_s2 <- ggplot(sres, aes(x = C_K, y = sex_max)) +
  geom_hline(yintercept = sex_ref, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.7, colour = "#4A90D9") +
  geom_point(size = 2, alpha = 0.85, colour = "#4A90D9") +
  scale_x_log10(breaks = SMOOTHS, labels = format(SMOOTHS, scientific = TRUE)) +
  expand_limits(y = 0) +
  labs(title = "Sex association of the subject loadings",
       subtitle = sprintf("dashed line: reference C_K (%.2f)", sex_ref),
       x = expression(C[K]~"(log scale)"), y = "max |cor(component, sex)|") +
  THEME_MS + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(.outdir, "supp_fig_smoothing.pdf"),
       (p_s1 | p_s2) + plot_annotation(tag_levels = "A"),
       width = 11, height = 4.8, bg = "white")

cat(sprintf("\n  figures written to %s\n  and                %s\n",
            normalizePath(file.path(.outdir, "supp_fig_weights.pdf"), mustWork = FALSE),
            normalizePath(file.path(.outdir, "supp_fig_smoothing.pdf"), mustWork = FALSE)))
