# ============================================================================
# supp_momspi_stability.R
#
# PURPOSE (reviewer 1, comment 7)
# Stability of the MOMS-PI cytokine feature loadings (Supplementary Figure S1).
# In the joint microbiome + cytokine fit of 01_momspi.R, IL-1ra and GM-CSF are
# the two cytokines with the largest absolute loadings, with opposite signs, in
# each of the two leading components. This script refits the joint model on
# 100 resamples, each drawing 80% of the subjects without replacement
# stratified by birth outcome (PTB / TB), and records how often that ranking
# holds. "Top two with opposite signs" is invariant to the sign indeterminacy
# of the decomposition, so no sign alignment across resamples is needed.
#
# OUTPUT (supp_output/)
#   supp_momspi_stability.csv       one row per (resample, component)
#   supp_fig_momspi_stability.pdf   IL-1ra vs GM-CSF loadings, components 1-2
#
# The 100 refits take roughly 15-20 minutes. If supp_momspi_stability.csv
# already exists it is reused and only the figure is redrawn; delete it to
# recompute. Requires the processed data from 01a_momspi_preprocess.R, which is
# run automatically if it is missing.
# ============================================================================

library(multi.tempted)
library(ggplot2)
library(patchwork)

N_RESAMPLE <- 100
FRAC       <- 0.8
SEED       <- 2024
R          <- 3

# Everything is read from and written to folders beside this script, NOT
# getwd(), so it runs the same from source()/RStudio "Source" or Rscript.
.scriptdir <- local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {                    # source() / RStudio Source
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { d <- dirname(normalizePath(of, mustWork = FALSE)); break }
  }
  if (is.null(d)) {                                     # Rscript supp_momspi_stability.R
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
.outdir  <- file.path(.scriptdir, "supp_output")
dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("output directory: %s\n", .outdir))
.csv <- file.path(.outdir, "supp_momspi_stability.csv")


# ============================================================================
# resampling (skipped if results are already saved)
# ============================================================================
if (file.exists(.csv)) {
  cat(sprintf("reusing %s (delete it to recompute)\n", .csv))
  res <- utils::read.csv(.csv, stringsAsFactors = FALSE)
} else {
  if (!file.exists(file.path(.datadir, "MOMSPI_microbiome_clr.csv"))) {
    cat("processed MOMS-PI data not found; running 01a_momspi_preprocess.R\n")
    source(file.path(.scriptdir, "01a_momspi_preprocess.R"), local = new.env())
  }
  .rd <- function(f) utils::read.csv(file.path(.datadir, f), row.names = 1)
  clr_16s   <- .rd("MOMSPI_microbiome_clr.csv")
  cyto_log2 <- .rd("MOMSPI_cytokine_log2.csv")
  meta_16s  <- .rd("MOMSPI_microbiome_meta.csv")
  meta_cyto <- .rd("MOMSPI_cytokine_meta.csv")
  tp_16s  <- as.numeric(meta_16s$visit_number);  sid_16s  <- as.character(meta_16s$subject_ep)
  tp_cyto <- as.numeric(meta_cyto$visit_number); sid_cyto <- as.character(meta_cyto$subject_ep)

  group_subID <- meta_16s[!duplicated(meta_16s$subject_ep), c("subject_ep", "Group")]
  ptb_subs <- group_subID$subject_ep[group_subID$Group == "PTB"]
  tb_subs  <- group_subID$subject_ep[group_subID$Group == "TB"]
  n_ptb <- round(FRAC * length(ptb_subs)); n_tb <- round(FRAC * length(tb_subs))
  cat(sprintf("resampling %d PTB of %d and %d TB of %d subjects, %d times\n",
              n_ptb, length(ptb_subs), n_tb, length(tb_subs), N_RESAMPLE))

  # read.csv turns "IL-1ra" / "GM-CSF" into "IL.1ra" / "GM.CSF"
  il1ra_col <- grep("IL.1ra", colnames(cyto_log2), value = TRUE)[1]
  gmcsf_col <- grep("GM.CSF", colnames(cyto_log2), value = TRUE)[1]

  set.seed(SEED)
  rows <- vector("list", N_RESAMPLE)
  for (i in seq_len(N_RESAMPLE)) {
    sel <- c(sample(ptb_subs, n_ptb, replace = FALSE), sample(tb_subs, n_tb, replace = FALSE))
    i16 <- which(sid_16s %in% sel); icy <- which(sid_cyto %in% sel)
    fit <- suppressMessages(multitempted_all(
      featuretables = list(microbiome = as.matrix(clr_16s[i16, ]), cytokine = as.matrix(cyto_log2[icy, ])),
      timepoints    = list(microbiome = tp_16s[i16],  cytokine = tp_cyto[icy]),
      subjectID     = list(microbiome = sid_16s[i16], cytokine = sid_cyto[icy]),
      transforms = "none", do_ratio = FALSE, smooth = 1e-4, r = R))
    b <- fit$B_hat$cytokine
    rows[[i]] <- do.call(rbind, lapply(seq_len(R), function(k) {
      top2 <- rownames(b)[order(abs(b[, k]), decreasing = TRUE)[1:2]]
      data.frame(resample = i, component = k,
                 top2     = il1ra_col %in% top2 && gmcsf_col %in% top2,
                 opp_sign = sign(b[il1ra_col, k]) != sign(b[gmcsf_col, k]),
                 il1ra    = b[il1ra_col, k], gmcsf = b[gmcsf_col, k])
    }))
    if (i %% 10 == 0) { cat(sprintf("  %d/%d resamples done\n", i, N_RESAMPLE)); utils::flush.console() }
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, .csv, row.names = FALSE)
}

res$stable <- res$top2 & res$opp_sign
cat("\n== IL-1ra and GM-CSF top two with opposite signs (% of resamples) ==\n")
print(stats::aggregate(cbind(pct = 100 * stable) ~ component, data = res, FUN = mean), row.names = FALSE)


# ============================================================================
# figure: components 1 and 2
# ============================================================================
STATUS_COLS <- c("Top two, opposite signs" = "#4A90D9", "Otherwise" = "#BBBBBB")
res$status  <- factor(ifelse(res$stable, names(STATUS_COLS)[1], names(STATUS_COLS)[2]),
                      levels = names(STATUS_COLS))
.panel <- function(k) {
  d <- res[res$component == k, ]
  ggplot(d, aes(x = il1ra, y = gmcsf, colour = status)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_colour_manual(values = STATUS_COLS, name = NULL, drop = FALSE) +
    labs(title = paste0("Component ", k),
         subtitle = sprintf("Top two with opposite signs in %.0f%% of resamples", 100 * mean(d$stable)),
         x = "IL-1ra loading", y = "GM-CSF loading") +
    theme_bw(base_size = 11)
}
p <- (.panel(1) | .panel(2)) + plot_layout(guides = "collect") +
  plot_annotation(title = sprintf("IL-1ra vs GM-CSF loadings across %d resamples", N_RESAMPLE)) &
  theme(legend.position = "bottom")
ggsave(file.path(.outdir, "supp_fig_momspi_stability.pdf"), p, width = 10, height = 4.4, bg = "white")
cat(sprintf("\nfigure written to %s\n", normalizePath(file.path(.outdir, "supp_fig_momspi_stability.pdf"), mustWork = FALSE)))
