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
# FOUR TESTS
#
#   A. CONTROL vs CONTROL (Homogenate vs Uninfected)
#      Two animals, neither infected. If the translational suppression
#      signature appears here, it is between-animal variation. This is the
#      decisive negative control.
#
#   B. INFECTED vs INFECTED at one timepoint (6-hpiA vs 6-hpiD)
#      Two animals, same treatment. Measures the magnitude of
#      between-animal variation directly.
#
#   C. EACH 6 hpi ANIMAL SEPARATELY vs pooled control
#      If the suppression reproduces in both animals independently, it is
#      reproducible across individuals at that timepoint.
#
#   D. PSEUDOBULK, 2 vs 2, at 6 hpi
#      The only properly replicated comparison in the dataset: two control
#      animals against two 6 hpi animals, aggregated to one profile per
#      animal per cluster. Underpowered, but the fold change is
#      animal-level, not nucleus-level. If the signature survives here it
#      answers the criticism directly.
#
# Outputs are written in the same format as 04b, so 05c can be pointed at
# this directory and run unchanged.

library(Seurat)
library(Matrix)
library(dplyr)
library(readr)
library(tibble)

seurat_obj <- readRDS("/home/pdewari/Documents/parse_2025/seurat_2025/seu_obj_umap_18d_6r_3kRes.rds")

# Rename clusters with final cell-type annotations
new.cluster.ids <- c(
  "Cluster 0", "Cluster 1",
  "Gill ciliary cells", "Hepatopancreas cells",
  "Gill neuroepithelial cells", "Gill cell type 1",
  "Hyalinocytes", "Haemocyte cell type 1",
  "Mantle cell type 1", "Cluster 9",
  "Vesicular haemocytes", "Immature haemocytes",
  "Macrophage like cells", "Adductor muscle cells",
  "Mantle cell type 2", "Mantle epithelial cells",
  "Gill cell type 2", "Small granule cells"
)

names(new.cluster.ids) <- levels(seurat_obj)

seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)


###############
# Add simplified condition column
seurat_obj$condition_new <- case_when(
  seurat_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seurat_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seurat_obj$sample == "24-hpiA" ~ "mid?",
  seurat_obj$sample == "24-hpiJ" ~ "24hpi",
  seurat_obj$sample == "72-hpiJ" ~ "72hpi",
  seurat_obj$sample == "96-hpiE" ~ "96hpi"
)

# Set factor order
seurat_obj$condition_new <- factor(
  seurat_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove 24-hpiA ("mid?")
seurat_obj_clean <- subset(seurat_obj, subset = condition_new != "mid?")
rm(seurat_obj)

# Drop unused factor levels
seurat_obj_clean$condition_new <- droplevels(seurat_obj_clean$condition_new)

seurat_obj_clean$cluster_annotation <- as.character(Idents(seurat_obj_clean))

head(seurat_obj_clean)


SEURAT_OBJ <- seurat_obj_clean

# =========================================================
# CONFIG
# =========================================================
#SEURAT_OBJ <- seurat_obj_clean

if (!"cluster_annotation" %in% colnames(SEURAT_OBJ@meta.data)) {
  SEURAT_OBJ$cluster_annotation <- as.character(Idents(SEURAT_OBJ))
}

CLUSTER_COL   <- "cluster_annotation"
CONDITION_COL <- "condition_new"
SAMPLE_COL    <- "sample"

CONTROL_A <- "Uninfected"     # unchallenged
CONTROL_B <- "Homogenate"     # mock-challenged
INF_A     <- "6-hpiA"
INF_B     <- "6-hpiD"

MIN_PCT         <- 0.1        # matches the pathway-analysis DE run
MIN_CELLS_GROUP <- 20

# restrict to the clusters that carry the finding; NULL = all
CLUSTERS_TO_RUN <- c("Cluster 1", "Hepatopancreas cells", "Gill ciliary cells",
                     "Immature haemocytes", "Mantle cell type 1",
                     "Small granule cells", "Adductor muscle cells",
                     "Haemocyte cell type 1")

ASSAY  <- "RNA"
OUTDIR <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_pseudorep_checks"

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

clusters <- if (is.null(CLUSTERS_TO_RUN)) {
  sort(unique(as.character(meta[[CLUSTER_COL]])))
} else {
  intersect(CLUSTERS_TO_RUN, unique(as.character(meta[[CLUSTER_COL]])))
}

cat("Clusters:", length(clusters), "\n\n")

log_rows <- list(); i_log <- 0L

# ---------------------------------------------------------
# helper: one FindMarkers comparison, written in 04b format
# ---------------------------------------------------------
run_pair <- function(cl, cells_1, cells_2, label_1, label_2, tag) {

  if (length(cells_1) < MIN_CELLS_GROUP || length(cells_2) < MIN_CELLS_GROUP) {
    cat(sprintf("    %-22s SKIPPED (%d vs %d cells)\n", tag,
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
  cat(sprintf("    %-22s tested: %5d | |lfc|>=1 & padj<0.05: %4d\n",
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
  infA  <- smp(INF_A);     infB  <- smp(INF_B)
  ctrl_all <- c(ctrlA, ctrlB)

  # A — control vs control: the decisive negative control
  r <- run_pair(cl, ctrlA, ctrlB, "ctrlA", "ctrlB", "A_ctrl_vs_ctrl")
  if (!is.null(r)) { i_log <- i_log + 1L; log_rows[[i_log]] <- r }

  # B — infected vs infected at the same timepoint
  r <- run_pair(cl, infA, infB, "infA", "infB", "B_6hpiA_vs_6hpiD")
  if (!is.null(r)) { i_log <- i_log + 1L; log_rows[[i_log]] <- r }

  # C — each 6 hpi animal separately against pooled control
  r <- run_pair(cl, ctrl_all, infA, "control", "infA", "C_ctrl_vs_6hpiA")
  if (!is.null(r)) { i_log <- i_log + 1L; log_rows[[i_log]] <- r }
  r <- run_pair(cl, ctrl_all, infB, "control", "infB", "C_ctrl_vs_6hpiD")
  if (!is.null(r)) { i_log <- i_log + 1L; log_rows[[i_log]] <- r }
}

# =========================================================
# TEST D — PSEUDOBULK, 2 controls vs 2 infected at 6 hpi
# =========================================================
# Counts are summed per animal per cluster, converted to CPM, and the
# animal-level log2 fold change computed between group means. This is a
# biological-replicate effect size: nuclei contribute to their animal's
# profile but are not themselves treated as replicates.
#
# With n = 2 per group no per-gene significance is claimed; the output is
# a ranking for GSEA, which is exactly what is needed to ask whether the
# translational signature survives animal-level aggregation.

cat("\n==================================================\n")
cat("TEST D — pseudobulk 2 vs 2 at 6 hpi\n")
cat("==================================================\n")

counts_all <- get_counts(SEURAT_OBJ, ASSAY)
gene_ids   <- clean_id(rownames(counts_all))

pb_animals <- c(CONTROL_A, CONTROL_B, INF_A, INF_B)
pb_group   <- c("control", "control", "infected", "infected")

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
                    cluster = safe_name(cl), stage = "D_pseudobulk_6hpi",
                    stringsAsFactors = FALSE)

  write_tsv(res, file.path(OUTDIR, "full",
                           paste0(safe_name(cl), "_D_pseudobulk_6hpi_full.tsv")))

  cat(sprintf("  %-24s %5d genes | median |lfc| %.2f\n",
              cl, nrow(res), median(abs(res$avg_log2FC))))
}

# =========================================================
write_tsv(bind_rows(log_rows), file.path(OUTDIR, "pseudorep_check_log.tsv"))

cat("\nOutputs:", file.path(OUTDIR, "full"), "\n")
cat("\nNext, in the go-enrich session, point 05c at this directory:\n")
cat('  de_dir <- "', OUTDIR, '"\n', sep = "")
cat("and run it. Then read GSEA_tidy.tsv by the `stage` column:\n\n")
cat("  A_ctrl_vs_ctrl    -> translational terms here mean the signature is\n")
cat("                       between-animal variation. This is the test that\n")
cat("                       matters most.\n")
cat("  B_6hpiA_vs_6hpiD  -> magnitude of between-animal variation.\n")
cat("  C_ctrl_vs_6hpiA   -> does the signature appear in each infected\n")
cat("  C_ctrl_vs_6hpiD      animal independently?\n")
cat("  D_pseudobulk_6hpi -> does it survive animal-level aggregation?\n")
