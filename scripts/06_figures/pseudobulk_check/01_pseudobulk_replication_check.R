# =========================================================
# Section 8: Pseudoreplication diagnostics
# Environment: SEURAT session
# =========================================================
# Addresses the criticism that FindMarkers on nuclei from one animal versus
# nuclei from another treats nuclei as replicates.
#
# The criticism is correct for the reported p-values. What it does NOT
# automatically invalidate is the RANKING used for GSEA, which is built
# from avg_log2FC rather than from p-values. This script tests whether the
# ranking itself is driven by infection or by between-animal variation.
#
# NOTE ON A BAD ARGUMENT: "the signal appears in many clusters, so it
# cannot be a single-animal artefact" is NOT valid. Every cluster comes
# from the same animals, so a whole-animal difference appears in all of
# them. Do not use that argument.
#
# BOTH REPLICATED TIMEPOINTS ARE NOW RUN. 24-hpiA is back in the analysis,
# so 6 and 24 hpi each have two animals and each supports a 2-vs-2
# pseudobulk comparison.
#
# READ TEST B WITH THE qPCR TABLE OPEN. Host-normalised qPCR shows that at
# both early timepoints one animal of the pair carried little or no virus
# (6-hpiD undetectable; 24-hpiA at 9.7e-7 relative to host, against 0.313
# in 24-hpiJ). So B is NOT a pure between-animal measure at either
# timepoint: it mixes individual variation with a real difference in
# infection. Test A, the two uninfected controls, is the only comparison
# in the dataset where no infection difference is possible, and it remains
# the decisive negative control.
#
# FOUR TESTS
#
#   A. CONTROL vs CONTROL (Homogenate vs Uninfected)
#      Two animals, neither infected. If the translational suppression
#      signature appears here, it is between-animal variation. Decisive.
#
#   B. ANIMAL vs ANIMAL within a timepoint, for 6 and 24 hpi
#      An upper bound on between-animal variation, inflated by the
#      infection difference described above.
#
#   C. EACH ANIMAL SEPARATELY vs pooled control
#      If the suppression reproduces in both animals of a pair
#      independently, it is reproducible across individuals.
#
#   D. PSEUDOBULK, 2 vs 2, at 6 and at 24 hpi
#      The properly replicated comparisons: two control animals against
#      two infected animals, aggregated to one profile per animal per
#      cluster. Underpowered, but the fold change is animal-level rather
#      than nucleus-level. Note that at each timepoint one of the two
#      infected animals carries little virus, so D is conservative: a
#      signature that survives it is not an artefact of aggregation.
#
# Outputs are written in the same format as 04b, so 05c can be pointed at
# this directory and run unchanged.

library(Seurat)
library(Matrix)
library(dplyr)
library(readr)
library(tibble)

# =========================================================
# CONFIG
# =========================================================


SEU_RDS <- "/home/pdewari/Documents/parse_2025/seurat_2025/seu_obj_umap_18d_6r_3kRes.rds"
CLUSTER_IDS <- c(
  "Cluster 0", "Cluster 1",
  "Gill ciliary cells", "Hepatopancreas cells",
  "Gill neuroepithelial cells", "Gill cell type 1",
  "Hyalinocytes", "Haemocyte cell type 1",
  "Mantle cell type 1", "Cluster 9",
  "Vesicular haemocytes", "Immature haemocytes",
  "Macrophage-like cells", "Adductor muscle cells",
  "Mantle cell type 2", "Mantle epithelial cells",
  "Gill cell type 2", "Small granule cells")

SEURAT_OBJ <- readRDS(SEU_RDS)

names(CLUSTER_IDS) <- levels(SEURAT_OBJ)
SEURAT_OBJ <- RenameIdents(SEURAT_OBJ, CLUSTER_IDS)
SEURAT_OBJ$cluster_annotation <- as.character(Idents(SEURAT_OBJ))

CLUSTER_COL   <- "cluster_annotation"
CONDITION_COL <- "condition_new"
SAMPLE_COL    <- "sample"

CONTROL_A <- "Uninfected"     # unchallenged
CONTROL_B <- "Homogenate"     # mock-challenged

## Every timepoint with two animals. Add a pair here and all of B, C and D
## pick it up; nothing else needs changing.
TIMEPOINTS <- list(
  "6hpi"  = c("6-hpiA",  "6-hpiD"),
  "24hpi" = c("24-hpiA", "24-hpiJ")
)

MIN_PCT         <- 0.1        # matches the pathway-analysis DE run
MIN_CELLS_GROUP <- 20

# restrict to the clusters that carry the finding; NULL = all
CLUSTERS_TO_RUN <- c("Cluster 1", "Hepatopancreas cells", "Gill ciliary cells",
                     "Immature haemocytes", "Mantle cell type 1",
                     "Small granule cells", "Adductor muscle cells",
                     "Haemocyte cell type 1")

ASSAY  <- "RNA"
OUTDIR <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/new_figures_8samples/pseudobulk_checks"

# =========================================================
dir.create(file.path(OUTDIR, "full"), showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
clean_id  <- function(x) sub("^(\\S+)\\s.*$", "\\1", x)

get_counts <- function(obj, assay) {
  tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
           error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = "counts"))
}

DefaultAssay(SEURAT_OBJ) <- ASSAY
if (inherits(SEURAT_OBJ[[ASSAY]], "Assay5")) {
  if (length(SeuratObject::Layers(SEURAT_OBJ[[ASSAY]], search = "data")) > 1)
    SEURAT_OBJ[[ASSAY]] <- SeuratObject::JoinLayers(SEURAT_OBJ[[ASSAY]])
}

meta <- SEURAT_OBJ@meta.data

## Fail early and by name rather than producing empty comparisons later.
wanted <- c(CONTROL_A, CONTROL_B, unlist(TIMEPOINTS))
missing <- setdiff(wanted, unique(as.character(meta[[SAMPLE_COL]])))
if (length(missing))
  stop("samples not found in ", SAMPLE_COL, ": ", paste(missing, collapse = ", "))

clusters <- if (is.null(CLUSTERS_TO_RUN)) {
  sort(unique(as.character(meta[[CLUSTER_COL]])))
} else {
  intersect(CLUSTERS_TO_RUN, unique(as.character(meta[[CLUSTER_COL]])))
}

cat("Clusters:", length(clusters), "\n")
cat("Replicated timepoints:", paste(names(TIMEPOINTS), collapse = ", "), "\n\n")

log_rows <- list(); i_log <- 0L

# ---------------------------------------------------------
# helper: one FindMarkers comparison, written in 04b format
# ---------------------------------------------------------
run_pair <- function(cl, cells_1, cells_2, label_1, label_2, tag) {

  if (length(cells_1) < MIN_CELLS_GROUP || length(cells_2) < MIN_CELLS_GROUP) {
    cat(sprintf("    %-26s SKIPPED (%d vs %d cells)\n", tag,
                length(cells_1), length(cells_2)))
    return(NULL)
  }

  sub <- subset(SEURAT_OBJ, cells = c(cells_1, cells_2))
  grp <- ifelse(colnames(sub) %in% cells_1, label_1, label_2)
  names(grp) <- colnames(sub)
  Idents(sub) <- factor(grp)

  # ident.1 = the "control-like" side, matching the 04b convention:
  # positive avg_log2FC = higher in ident.1
  res <- tryCatch(
    FindMarkers(sub, ident.1 = label_1, ident.2 = label_2,
                logfc.threshold = 0, min.pct = MIN_PCT),
    error = function(e) { cat("    FAILED:", conditionMessage(e), "\n"); NULL })

  if (is.null(res) || nrow(res) == 0) return(NULL)

  res <- res %>%
    rownames_to_column("gene_full") %>%
    mutate(gene_id = clean_id(gene_full),
           cluster = safe_name(cl),
           stage   = tag)

  write_tsv(res, file.path(OUTDIR, "full",
                           paste0(safe_name(cl), "_", tag, "_full.tsv")))

  n_sig <- sum(res$p_val_adj < 0.05 & abs(res$avg_log2FC) >= 1)
  cat(sprintf("    %-26s tested: %5d | |lfc|>=1 & padj<0.05: %4d\n",
              tag, nrow(res), n_sig))

  data.frame(cluster = cl, comparison = tag,
             n_1 = length(cells_1), n_2 = length(cells_2),
             n_tested = nrow(res), n_sig = n_sig, stringsAsFactors = FALSE)
}

# =========================================================
# TESTS A, B, C
# =========================================================
for (cl in clusters) {

  cat("==================================================\n", cl, "\n")
  cells_cl <- rownames(meta)[as.character(meta[[CLUSTER_COL]]) == cl]
  smp <- function(s) intersect(cells_cl, rownames(meta)[as.character(meta[[SAMPLE_COL]]) == s])

  ctrlA <- smp(CONTROL_A); ctrlB <- smp(CONTROL_B)
  ctrl_all <- c(ctrlA, ctrlB)

  # A — control vs control: the decisive negative control, run once
  r <- run_pair(cl, ctrlA, ctrlB, "ctrlA", "ctrlB", "A_ctrl_vs_ctrl")
  if (!is.null(r)) { i_log <- i_log + 1L; log_rows[[i_log]] <- r }

  for (tp in names(TIMEPOINTS)) {
    a <- TIMEPOINTS[[tp]][1]; b <- TIMEPOINTS[[tp]][2]
    cellsA <- smp(a); cellsB <- smp(b)

    # B — animal vs animal within the timepoint
    r <- run_pair(cl, cellsA, cellsB, "animalA", "animalB",
                  paste0("B_", tp, "_between_animals"))
    if (!is.null(r)) { i_log <- i_log + 1L; log_rows[[i_log]] <- r }

    # C — each animal separately against the pooled controls
    for (s in TIMEPOINTS[[tp]]) {
      r <- run_pair(cl, ctrl_all, smp(s), "control", "infected",
                    paste0("C_ctrl_vs_", safe_name(s)))
      if (!is.null(r)) { i_log <- i_log + 1L; log_rows[[i_log]] <- r }
    }
  }
}

# =========================================================
# TEST D — PSEUDOBULK, 2 controls vs 2 infected, per timepoint
# =========================================================
# Counts are summed per animal per cluster, converted to CPM, and the
# animal-level log2 fold change computed between group means. This is a
# biological-replicate effect size: nuclei contribute to their animal's
# profile but are not themselves treated as replicates.
#
# With n = 2 per group no per-gene significance is claimed; the output is
# a ranking for GSEA, which is exactly what is needed to ask whether the
# translational signature survives animal-level aggregation.

counts_all <- get_counts(SEURAT_OBJ, ASSAY)
gene_ids   <- clean_id(rownames(counts_all))

for (tp in names(TIMEPOINTS)) {

  cat("\n==================================================\n")
  cat("TEST D - pseudobulk 2 vs 2 at", tp, "\n")
  cat("==================================================\n")

  pb_animals <- c(CONTROL_A, CONTROL_B, TIMEPOINTS[[tp]])
  pb_group   <- c("control", "control", "infected", "infected")
  tag        <- paste0("D_pseudobulk_", tp)

  for (cl in clusters) {

    cells_cl <- rownames(meta)[as.character(meta[[CLUSTER_COL]]) == cl]

    prof <- sapply(pb_animals, function(s) {
      cells <- intersect(cells_cl, rownames(meta)[as.character(meta[[SAMPLE_COL]]) == s])
      if (length(cells) < MIN_CELLS_GROUP) return(rep(NA_real_, nrow(counts_all)))
      v <- Matrix::rowSums(counts_all[, cells, drop = FALSE])
      1e6 * v / sum(v)                                   # CPM
    })

    if (any(is.na(prof[1, ]))) {
      cat(sprintf("  %-24s skipped (an animal has too few nuclei)\n", cl)); next
    }

    # expression filter: detected in this cluster at all
    keep <- rowMeans(prof) > 1
    prof <- prof[keep, , drop = FALSE]
    gid  <- gene_ids[keep]

    mean_ctrl <- rowMeans(prof[, pb_group == "control",  drop = FALSE])
    mean_inf  <- rowMeans(prof[, pb_group == "infected", drop = FALSE])

    # sign convention matches 04b: positive = higher in CONTROL
    lfc <- log2((mean_ctrl + 1) / (mean_inf + 1))

    res <- data.frame(gene_full = rownames(prof), gene_id = gid,
                      avg_log2FC = lfc,
                      ctrl_cpm = mean_ctrl, inf_cpm = mean_inf,
                      p_val = NA_real_, p_val_adj = NA_real_,
                      cluster = safe_name(cl), stage = tag,
                      stringsAsFactors = FALSE)

    write_tsv(res, file.path(OUTDIR, "full",
                             paste0(safe_name(cl), "_", tag, "_full.tsv")))

    cat(sprintf("  %-24s %5d genes | median |lfc| %.2f\n",
                cl, nrow(res), median(abs(res$avg_log2FC))))
  }
}

# =========================================================
write_tsv(bind_rows(log_rows), file.path(OUTDIR, "pseudorep_check_log.tsv"))

cat("\nOutputs:", file.path(OUTDIR, "full"), "\n")
cat("\nNext, in the go-enrich session, point 05c at this directory:\n")
cat('  de_dir <- "', OUTDIR, '"\n', sep = "")
cat("and run it. Then read GSEA_tidy.tsv by the `stage` column:\n\n")
cat("  A_ctrl_vs_ctrl            -> translational terms here mean the\n")
cat("                               signature is between-animal variation.\n")
cat("                               The test that matters most.\n")
cat("  B_6hpi_between_animals    -> upper bound on between-animal variation;\n")
cat("  B_24hpi_between_animals      inflated, because one animal of each pair\n")
cat("                               carried little or no virus (see qPCR).\n")
cat("  C_ctrl_vs_<sample>        -> does the signature appear in each animal\n")
cat("                               independently?\n")
cat("  D_pseudobulk_6hpi         -> does it survive animal-level aggregation,\n")
cat("  D_pseudobulk_24hpi           at each replicated timepoint?\n\n")
cat("05d filters GSEA_tidy.tsv by these stage labels - update its list,\n")
cat("since B and D are no longer single comparisons.\n")

