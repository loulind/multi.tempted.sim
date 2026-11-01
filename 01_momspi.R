# ============================================================================
# 01_momspi.R
#
# PURPOSE
# MOMS-PI case study (manuscript Figure 1). Fit multiTEMPTED jointly to the
# vaginal microbiome (CLR) and cervicovaginal cytokines (log2), and, for
# comparison, to the microbiome alone, each with r = 3 components.
#
# OUTPUT (output/)
#   01_momspi_subject_loading.pdf      Figure 1a: first two subject components,
#       joint vs microbiome-only, coloured by pregnancy outcome
#   01_momspi_feature_loading.pdf      Figure 1b: top 5% microbiome feature
#       loadings for components 1 and 2 of the joint fit
#   01_momspi_cytokine_trajectory.pdf  Figure 1c: smoothed IL-1ra and GM-CSF
#       trajectories by pregnancy outcome
#   01_momspi.rds                      both fits, for downstream use
#
# Requires the processed data from 01a_momspi_preprocess.R, which is run
# automatically if it is missing.
# ============================================================================

library(multi.tempted)
library(ggplot2)
library(patchwork)

# --- plot aesthetics (same palette and theme as the other manuscript figures)
THEME_MS     <- theme_bw(base_size = 12)
OUTCOME_COLS <- c(PTB = "#E05C5C", TB = "#4A90D9")

# Everything is read from and written to folders beside this script, NOT
# getwd(), so it runs the same from source()/RStudio "Source" or Rscript.
.scriptdir <- local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {                    # source() / RStudio Source
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { d <- dirname(normalizePath(of, mustWork = FALSE)); break }
  }
  if (is.null(d)) {                                     # Rscript 01_momspi.R
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
    warning("could not locate the script file; using ", d, call. = FALSE)
  }
  d
})
.datadir <- file.path(.scriptdir, "data", "momspi", "processed")
.outdir  <- file.path(.scriptdir, "output")
dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("output directory: %s\n", .outdir))

if (!file.exists(file.path(.datadir, "MOMSPI_microbiome_clr.csv"))) {
  cat("processed MOMS-PI data not found; running 01a_momspi_preprocess.R\n")
  source(file.path(.scriptdir, "01a_momspi_preprocess.R"), local = new.env())
}


# ============================================================================
# load data and build inputs
# ============================================================================
.rd <- function(f) utils::read.csv(file.path(.datadir, f), row.names = 1)
clr_16s   <- .rd("MOMSPI_microbiome_clr.csv")
cyto_log2 <- .rd("MOMSPI_cytokine_log2.csv")
meta_16s  <- .rd("MOMSPI_microbiome_meta.csv")
meta_cyto <- .rd("MOMSPI_cytokine_meta.csv")
taxa      <- .rd("MOMSPI_microbiome_taxa.csv")

# Readable taxon labels: "Genus species", else genus, else [family] or
# [order]; labels shared by several OTUs get the OTU ID appended. read.csv
# prefixes numeric OTU IDs with "X", which is stripped first.
otu_ids <- sub("^X", "", setdiff(colnames(clr_16s), "Other"))
.label <- function(otu) {
  if (!otu %in% rownames(taxa)) return(otu)
  v <- vapply(c("Genus", "Species", "Family", "Order"), function(k) {
    x <- trimws(as.character(taxa[otu, k])); if (is.na(x)) "" else x }, character(1))
  if (v["Genus"] != "" && v["Species"] != "") return(paste(v["Genus"], v["Species"]))
  if (v["Genus"]  != "") return(unname(v["Genus"]))
  if (v["Family"] != "") return(paste0("[", v["Family"], "]"))
  if (v["Order"]  != "") return(paste0("[", v["Order"], "]"))
  otu
}
raw_labels <- vapply(otu_ids, .label, character(1))
dup <- duplicated(raw_labels) | duplicated(raw_labels, fromLast = TRUE)
raw_labels[dup] <- paste0(raw_labels[dup], "_", otu_ids[dup])
colnames(clr_16s) <- c(raw_labels, "Other")

tp_16s  <- as.numeric(meta_16s$visit_number);  sid_16s  <- as.character(meta_16s$subject_ep)
tp_cyto <- as.numeric(meta_cyto$visit_number); sid_cyto <- as.character(meta_cyto$subject_ep)
group_subID <- meta_16s[!duplicated(meta_16s$subject_ep), c("subject_ep", "Group")]
colnames(group_subID) <- c("subjectID", "Group")


# ============================================================================
# fits
# ============================================================================
cat("== joint fit: microbiome + cytokine ==\n")
out_joint <- multitempted_all(
  featuretables = list(microbiome = as.matrix(clr_16s), cytokine = as.matrix(cyto_log2)),
  timepoints    = list(microbiome = tp_16s,  cytokine = tp_cyto),
  subjectID     = list(microbiome = sid_16s, cytokine = sid_cyto),
  transforms = "none", do_ratio = FALSE, smooth = 1e-4, r = 3)

cat("== microbiome-only fit ==\n")
out_micro <- multitempted_all(
  featuretables = list(microbiome = as.matrix(clr_16s)),
  timepoints    = list(microbiome = tp_16s),
  subjectID     = list(microbiome = sid_16s),
  transforms = "none", do_ratio = FALSE, smooth = 1e-4, r = 3)

saveRDS(list(joint = out_joint[c("A_hat", "B_hat", "Zeta_hat", "time_Zeta", "Lambda")],
             micro = out_micro[c("A_hat", "B_hat", "Zeta_hat", "time_Zeta", "Lambda")],
             group = group_subID),
        file.path(.outdir, "01_momspi.rds"))


# ============================================================================
# Figure 1a: subject loadings, joint vs microbiome-only
# ============================================================================
.subject_plot <- function(A, title) {
  df <- data.frame(subjectID = rownames(A), PC1 = A[, 1], PC2 = A[, 2])
  df <- merge(df, group_subID, by = "subjectID")
  ggplot(df, aes(PC1, PC2, colour = Group)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_colour_manual(values = OUTCOME_COLS, name = "Outcome") +
    labs(title = title, subtitle = "Subject Loadings", x = "Component 1", y = "Component 2") +
    THEME_MS
}
p_subj <- (.subject_plot(out_joint$A_hat, "multiTEMPTED on Microbiome + Cytokine") +
             theme(legend.position = "none")) +
  .subject_plot(out_micro$A_hat, "TEMPTED on Microbiome Only") +
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(.outdir, "01_momspi_subject_loading.pdf"), p_subj, width = 10, height = 4.5)


# ============================================================================
# Figure 1b: top 5% microbiome feature loadings, components 1 and 2 (joint)
# ============================================================================
.feature_plot <- function(out, component, pct = 0.05, modality = "microbiome") {
  b    <- out$B_hat[[modality]]
  vals <- b[, component]
  idx  <- order(abs(vals), decreasing = TRUE)[seq_len(max(1, round(pct * length(vals))))]
  df   <- data.frame(feature = rownames(b)[idx], loading = vals[idx])
  df$feature   <- factor(df$feature, levels = df$feature[order(df$loading)])
  df$direction <- ifelse(df$loading > 0, "positive", "negative")
  ggplot(df, aes(x = loading, y = feature, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = c(positive = "#4A90D9", negative = "#E05C5C"), guide = "none") +
    labs(title = paste0("Microbiome loadings - Component ", component), x = "Loading", y = NULL) +
    theme_bw(base_size = 11)
}
p_feat <- .feature_plot(out_joint, 1) | .feature_plot(out_joint, 2)
ggsave(file.path(.outdir, "01_momspi_feature_loading.pdf"), p_feat, width = 10, height = 3)


# ============================================================================
# Figure 1c: IL-1ra and GM-CSF trajectories by pregnancy outcome
# ============================================================================
# read.csv turns "IL-1ra" / "GM-CSF" into "IL.1ra" / "GM.CSF"
il1ra_col <- grep("IL.1ra", colnames(cyto_log2), value = TRUE)[1]
gmcsf_col <- grep("GM.CSF", colnames(cyto_log2), value = TRUE)[1]
cyto_long <- rbind(
  data.frame(cytokine = "IL-1ra", log2_conc = cyto_log2[[il1ra_col]],
             visit_number = meta_cyto$visit_number, Group = meta_cyto$Group),
  data.frame(cytokine = "GM-CSF", log2_conc = cyto_log2[[gmcsf_col]],
             visit_number = meta_cyto$visit_number, Group = meta_cyto$Group))
cyto_long$cytokine <- factor(cyto_long$cytokine, levels = c("IL-1ra", "GM-CSF"))
p_cyto <- ggplot(cyto_long, aes(visit_number, log2_conc, colour = Group, fill = Group)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, alpha = 0.2) +
  facet_wrap(~cytokine, scales = "free_y") +
  scale_colour_manual(values = OUTCOME_COLS, name = "Outcome") +
  scale_fill_manual(values = OUTCOME_COLS, name = "Outcome") +
  scale_x_continuous(breaks = 1:6) +
  labs(x = "Visit number", y = expression(log[2] ~ "concentration"),
       title = "IL-1ra and GM-CSF over time by pregnancy outcome") +
  THEME_MS + theme(legend.position = "bottom")
ggsave(file.path(.outdir, "01_momspi_cytokine_trajectory.pdf"), p_cyto, width = 10, height = 4)

# top cytokines per component of the joint fit (the basis for Figure 1c)
b_cyto <- out_joint$B_hat$cytokine
cat("\n== top two cytokines by |loading|, joint fit ==\n")
for (k in seq_len(ncol(b_cyto))) {
  top <- order(abs(b_cyto[, k]), decreasing = TRUE)[1:2]
  cat(sprintf("  component %d: %s (%.3f), %s (%.3f)\n", k, rownames(b_cyto)[top[1]],
              b_cyto[top[1], k], rownames(b_cyto)[top[2]], b_cyto[top[2], k]))
}
cat(sprintf("\nfigures written to %s\n", normalizePath(.outdir, mustWork = FALSE)))
