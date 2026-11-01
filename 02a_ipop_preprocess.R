# ============================================================================
# 02a_ipop_preprocess.R
#
# PURPOSE
# Build the iPOP exercise data analysed in 02_ipop.R and
# supp_sensitivity.R. The result is the `ipop` object shipped with the
# multi.tempted package (multi.tempted::ipop); this script documents how it
# was made and rebuilds it from the public download.
#
# INPUTS (data/iPOP/; NOT included in this repository)
# Download from the iPOP data portal, http://hmp2-data.stanford.edu/
# (Contrepois et al. 2020, Cell, "Molecular choreography of acute exercise"):
#   NIHMS1588589-supplement-9.xlsx   participant characteristics
#   iPOP_visit_info.csv              visit / sample information
#   Exercise_Metabolomics.csv        untargeted metabolomics (LC-MS)
#   Exercise_Lipidomics.csv          semi-targeted lipidomics (Lipidyzer)
#   Exercise_targ_proteomics.csv     targeted proteomics (immunoassays)
#   Exercise_Proteomics.csv          untargeted proteomics (SWATH-MS)
#
# STEPS
#   * visits: keep the exercise session and code its five time points 1-5
#   * metabolome: drop metabolites with extreme outliers (generalized ESD test,
#     up to 20 outliers, minimum p-value < 1e-19), then log10
#   * lipidome: group lipid species into 13 classes, with TAG further split by
#     fatty-acid chain length (MCFA <= 12, LCFA 13-21, VLCFA >= 22 carbons) and,
#     for LCFA, by saturation; within each group keep the principal component
#     scores of the scaled log10 data that explain < 95% cumulative variance
#   * targeted proteomics ("cytokine"): drop one sample with many missing
#     values, then log10
#   * untargeted proteomics ("protein"): log10
#   * keep samples present in all four modalities
#
# OUTPUT
#   data/iPOP/processed/ipop.rds   list(meta_subj, meta_visit, cytokine,
#                                       metabolome, lipid, protein)
#
# REQUIRES: PMCMRplus (generalized ESD test) and readxl.
# ============================================================================

library(PMCMRplus)
library(readxl)

# Everything is read from and written to folders beside this script, NOT
# getwd(), so it runs the same from source()/RStudio "Source" or Rscript.
.scriptdir <- local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {                    # source() / RStudio Source
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { d <- dirname(normalizePath(of, mustWork = FALSE)); break }
  }
  if (is.null(d)) {                                     # Rscript 02a_ipop_preprocess.R
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
.indir  <- file.path(.scriptdir, "data", "iPOP")
.outdir <- file.path(.indir, "processed")
.need <- c("NIHMS1588589-supplement-9.xlsx", "iPOP_visit_info.csv", "Exercise_Metabolomics.csv",
           "Exercise_Lipidomics.csv", "Exercise_targ_proteomics.csv", "Exercise_Proteomics.csv")
.missing <- .need[!file.exists(file.path(.indir, .need))]
if (length(.missing))
  stop("download the iPOP files from http://hmp2-data.stanford.edu/ into ", .indir,
       "\n  missing: ", paste(.missing, collapse = ", "), call. = FALSE)
dir.create(.outdir, showWarnings = FALSE, recursive = TRUE)
.rd <- function(f) utils::read.csv(file.path(.indir, f), header = TRUE, row.names = 1)


# ============================================================================
# subjects and visits
# ============================================================================
# participant table: header on row 3, one row per participant (36)
meta_subj <- as.data.frame(read_excel(file.path(.indir, "NIHMS1588589-supplement-9.xlsx"),
                                      sheet = 1, skip = 2, n_max = 36, trim_ws = FALSE))
colnames(meta_subj) <- make.names(colnames(meta_subj))
rownames(meta_subj) <- meta_subj$Participant.ID
colnames(meta_subj)[1] <- "subjectID"
meta_subj <- meta_subj[sort(meta_subj$subjectID), ]

visits <- .rd("iPOP_visit_info.csv")
visits$visitID <- rownames(visits)
visits <- visits[visits$Event == "Exercise", ]
old_val <- unique(visits$Event_Note1)                   # baseline, 2, 15, 30 min, 1 h
visits$timepoint <- match(visits$Event_Note1, old_val)
cat("exercise visits:", nrow(visits), "\n")


# ============================================================================
# metabolome: remove features with extreme outliers, log10
# ============================================================================
dat_metab <- .rd("Exercise_Metabolomics.csv")
min_gesd_p <- apply(dat_metab, 2, function(x) min(gesdTest(x, 20)$p.value))
log_metab  <- log10(dat_metab[, min_gesd_p > 10e-20])
cat(sprintf("metabolites: %d of %d retained\n", ncol(log_metab), ncol(dat_metab)))


# ============================================================================
# lipidome: class-wise principal components
# ============================================================================
dat_lipid <- .rd("Exercise_Lipidomics.csv")
types <- c("TAG", "CE", "CER", "DAG", "DCER", "FFA", "HCER", "LCER", "LPC", "LPE", "PC", "PE", "SM")
list_lipid <- lapply(seq_along(types), function(jj) {
  pat <- paste0("(?<![A-Za-z])", types[jj])
  if (jj > 1) pat <- paste0(pat, "[.]")
  dat_lipid[, grep(pat, colnames(dat_lipid), perl = TRUE)]
})
names(list_lipid) <- types

# split TAG by fatty-acid chain length and, within LCFA, by saturation
TAG_names     <- colnames(list_lipid$TAG)
fa_carbons    <- as.numeric(sub(".*FA(\\d+)\\.\\d+", "\\1", TAG_names))
double_bonds  <- as.numeric(sub("TAG\\d+\\.(\\d+)\\.FA\\d+\\.\\d+", "\\1", TAG_names))
carbon_subgrp <- ifelse(fa_carbons <= 12, "MCFA", ifelse(fa_carbons < 22, "LCFA", "VLCFA"))
satur_subgrp  <- ifelse(double_bonds == 0, "unsaturated", ifelse(double_bonds == 1, "mono", "poly"))
TAG <- list_lipid$TAG
list_lipid$TAG.MCFA               <- TAG[, carbon_subgrp == "MCFA"]
list_lipid$TAG.VLCFA              <- TAG[, carbon_subgrp == "VLCFA"]
list_lipid$TAG.LCFA.unsaturated   <- TAG[, carbon_subgrp == "LCFA" & satur_subgrp == "unsaturated"]
list_lipid$TAG.LCFA.monosaturated <- TAG[, carbon_subgrp == "LCFA" & satur_subgrp == "mono"]
list_lipid$TAG.LCFA.polysaturated <- TAG[, carbon_subgrp == "LCFA" & satur_subgrp == "poly"]
list_lipid$TAG <- NULL

# PC scores of the scaled log10 data explaining < 95% cumulative variance,
# signed so that each component's loadings sum to a positive value
get_pc <- function(x) {
  s  <- svd(scale(log10(x)))
  kp <- cumsum(s$d^2) / sum(s$d^2) < 0.95
  score <- sweep(s$u[, kp, drop = FALSE], 2, sign(colSums(s$v[, kp, drop = FALSE])), FUN = "*")
  dimnames(score) <- list(rownames(x), paste0("pc", seq_len(ncol(score))))
  score
}
log_lipid_pc <- lapply(list_lipid, get_pc)
for (g in names(log_lipid_pc)) colnames(log_lipid_pc[[g]]) <- paste0(g, ".", colnames(log_lipid_pc[[g]]))
log_lipid_pc <- do.call(cbind, log_lipid_pc)
cat(sprintf("lipids: %d species summarized by %d class-wise PCs\n", ncol(dat_lipid), ncol(log_lipid_pc)))


# ============================================================================
# targeted proteomics ("cytokine") and untargeted proteomics ("protein")
# ============================================================================
dat_cytokine <- .rd("Exercise_targ_proteomics.csv")
dat_cytokine <- dat_cytokine[rownames(dat_cytokine) != "ZOZOW1T-E12", ]   # mostly missing
log_cytokine <- log10(dat_cytokine)
log_protein  <- log10(.rd("Exercise_Proteomics.csv"))


# ============================================================================
# samples present in all four modalities
# ============================================================================
sID <- sort(Reduce(intersect, list(rownames(log_cytokine), rownames(log_metab),
                                   rownames(log_lipid_pc), rownames(log_protein))))
ipop <- list(meta_subj  = meta_subj,
             meta_visit = visits[sID, ],
             cytokine   = log_cytokine[sID, ],
             metabolome = log_metab[sID, ],
             lipid      = log_lipid_pc[sID, ],
             protein    = log_protein[sID, ])
print(sapply(ipop, dim))
saveRDS(ipop, file.path(.outdir, "ipop.rds"))
cat(sprintf("written to %s\n", normalizePath(file.path(.outdir, "ipop.rds"), mustWork = FALSE)))
