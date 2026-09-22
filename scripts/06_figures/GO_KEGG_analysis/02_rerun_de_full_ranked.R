# =========================================================
# Section 4b: Re-run control-vs-stage DE, full ranked output
# Environment: SEURAT session
# =========================================================
# ONE PASS THAT PRODUCES THREE THINGS
#
#   1. Full ranked statistics per comparison  -> input for GSEA (05c)
#   2. The tested gene set per cluster        -> the ORA universe (05b)
#   3. DEG lists at |log2FC| >= 1             -> identical to your current
#                                                manuscript lists
#
# The only change from your original call is logfc.threshold = 1 -> 0:
#
#   ORIGINAL: FindMarkers(seu_obj_clean, logfc.threshold = 1, min.pct = 0.25,
#                         ident.1 = control_label, ident.2 = target_label)
#   HERE:     FindMarkers(sub,           logfc.threshold = 0, min.pct = 0.25,
#                         ident.1 = CONTROL_LABEL, ident.2 = stage)
#
# YOUR EXISTING RESULTS ARE NOT INVALIDATED.
# Seurat computes p_val_adj as Bonferroni over ALL features in the assay,
# not over the genes that passed the fold-change prefilter. Genes at
# |log2FC| >= 1 therefore get identical p_val, p_val_adj, pct.1 and pct.2
# whether the threshold was 1 or 0. This re-run reproduces your current
# DEG lists exactly and adds everything below the cut.
#
# WHY BOTHER
# logfc.threshold = 1 is a 2-fold cut. It is a sound choice for naming
# individual genes in the text, but it discards most of the pathway
# signal: a coordinated 1.5-fold shift across thirty genes in one pathway
# is a real biological response and is invisible to an analysis that only
# sees 2-fold changes. GSEA on the full ranked list recovers exactly that,
# and needs no background at all. You keep the 2-fold lists for the
# gene-level narrative and use the ranked lists for the pathway figure.
#
# DIRECTION CONVENTION — PRESERVED FROM YOUR ORIGINAL SCRIPT
# ident.1 = control, ident.2 = target. So:
#   positive avg_log2FC = higher in CONTROL  = "up_in_control"
#   negative avg_log2FC = higher in INFECTED = "up_in_target"
# Every downstream script assumes this. Do not swap the idents.
#
# RUNTIME: roughly 72 comparisons (18 clusters x 4 stages) with no
# fold-change prefilter. Budget 40-90 minutes. Set CLUSTERS_TO_RUN to the
# three responsive clusters if you only need those.

library(Seurat)
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
SEURAT_OBJ <- seurat_obj_clean     # same object as the original DE run

# RenameIdents() changes Idents(), NOT a metadata column. Create it if absent.
if (!"cluster_annotation" %in% colnames(SEURAT_OBJ@meta.data)) {
  SEURAT_OBJ$cluster_annotation <- as.character(Idents(SEURAT_OBJ))
  message("Created cluster_annotation from Idents()")
}

CLUSTER_COL   <- "cluster_annotation"
CONDITION_COL <- "condition_new"   # your simplified column
CONTROL_LABEL <- "control"         # pools Homogenate + Uninfected
STAGES        <- c("6hpi", "24hpi", "72hpi", "96hpi")

# NULL = every cluster. Or restrict — LIVE names, with spaces:
#   c("Cluster 1", "Hepatopancreas cells", "Gill ciliary cells")
CLUSTERS_TO_RUN <- NULL

# >>> PATHWAY-ANALYSIS RUN — min.pct RELAXED FROM 0.25 TO 0.1 <<<
#
# At min.pct = 0.25 only 500-1,200 genes were testable per comparison, and
# after restricting to GO-annotated genes the ranked lists fell to 450-750.
# GSEA requires at least MIN_GS_SIZE of a term's genes to be present in the
# ranked list, so most GO terms were never tested at all.
#
# 0.1 gives ranked lists several times larger. This run is FOR PATHWAY
# ANALYSIS ONLY and writes to its own directory: the 0.25 results, which
# reproduce the DEG counts reported in the manuscript, are untouched.
#
# The DEG lists this run writes will therefore NOT match the manuscript's
# counts — do not use them for Table 1 or for gene-level statements.
MIN_PCT          <- 0.1

LOGFC_FOR_LISTS  <- 1      # matches your original DE call (2-fold)
PADJ_CUTOFF      <- 0.05

# Wilcoxon on ~10 nuclei per side gives unstable p-values. Your smallest
# clusters (Small granule cells: 367 cells across 5 conditions, Cluster 9:
# 714) will have thin per-condition groups. 20 is a safer floor; the log
# records every skipped comparison so nothing disappears silently.
MIN_CELLS_GROUP  <- 20

ASSAY  <- "RNA"
# separate from the min.pct = 0.25 run, which stays as it is
OUTDIR <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_full_ranked_minpct01"

# =========================================================
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "universes"),  showWarnings = FALSE)
dir.create(file.path(OUTDIR, "full"),       showWarnings = FALSE)
dir.create(file.path(OUTDIR, "deg_lists"),  showWarnings = FALSE)

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
clean_id  <- function(x) sub("^(\\S+)\\s.*$", "\\1", x)

DefaultAssay(SEURAT_OBJ) <- ASSAY

# Seurat v5: FindMarkers errors on split layers. Join them first.
if (inherits(SEURAT_OBJ[[ASSAY]], "Assay5")) {
  lyrs <- SeuratObject::Layers(SEURAT_OBJ[[ASSAY]], search = "data")
  if (length(lyrs) > 1) {
    cat("Seurat v5 object with split layers — running JoinLayers()\n")
    SEURAT_OBJ[[ASSAY]] <- SeuratObject::JoinLayers(SEURAT_OBJ[[ASSAY]])
  }
}

meta <- SEURAT_OBJ@meta.data
stopifnot(CLUSTER_COL %in% colnames(meta), CONDITION_COL %in% colnames(meta))

all_clusters <- sort(unique(as.character(meta[[CLUSTER_COL]])))
clusters_to_run <- if (is.null(CLUSTERS_TO_RUN)) all_clusters else
  intersect(CLUSTERS_TO_RUN, all_clusters)

cat("Clusters to run:", length(clusters_to_run), "\n")
cat("Stages:", paste(STAGES, collapse = ", "), "\n\n")

run_log <- list(); i_log <- 0L
t_start <- Sys.time()

for (cl in clusters_to_run) {

  cells_cl <- rownames(meta)[as.character(meta[[CLUSTER_COL]]) == cl]
  if (length(cells_cl) == 0) next

  cat("\n==================================================\n")
  cat("CLUSTER:", cl, "(", length(cells_cl), "cells )\n")
  cat("==================================================\n")

  sub <- subset(SEURAT_OBJ, cells = cells_cl)
  Idents(sub) <- sub@meta.data[[CONDITION_COL]]

  tested_union <- character(0)

  for (stage in STAGES) {

    n_ctrl <- sum(as.character(Idents(sub)) == CONTROL_LABEL)
    n_tgt  <- sum(as.character(Idents(sub)) == stage)

    if (n_ctrl < MIN_CELLS_GROUP || n_tgt < MIN_CELLS_GROUP) {
      cat(sprintf("  %-6s SKIPPED (control: %d, %s: %d)\n", stage, n_ctrl, stage, n_tgt))
      i_log <- i_log + 1L
      run_log[[i_log]] <- data.frame(
        cluster = cl, stage = stage, n_control = n_ctrl, n_target = n_tgt,
        n_tested = NA_integer_, n_up_control = NA_integer_, n_up_target = NA_integer_,
        skipped = TRUE, stringsAsFactors = FALSE)
      next
    }

    t0 <- Sys.time()

    res <- tryCatch(
      FindMarkers(sub,
                  ident.1         = CONTROL_LABEL,
                  ident.2         = stage,
                  logfc.threshold = 0,          # <-- the only change
                  min.pct         = MIN_PCT),
      error = function(e) { cat("    FAILED:", conditionMessage(e), "\n"); NULL }
    )

    if (is.null(res) || nrow(res) == 0) next

    # write the SANITISED cluster name so ORA (which derives cluster names
    # from universe filenames) and GSEA (which reads this column) agree
    res <- res %>%
      rownames_to_column(var = "gene_full") %>%
      mutate(gene_id = clean_id(gene_full),
             cluster = safe_name(cl),
             stage   = stage)

    tag <- paste0(safe_name(cl), "_control_vs_", stage)

    # --- 1. full ranked table (GSEA input) ---
    write_tsv(res, file.path(OUTDIR, "full", paste0(tag, "_full.tsv")))

    # --- 2. tested genes contribute to this cluster's universe ---
    tested_union <- union(tested_union, res$gene_id)

    # --- 3. DEG lists at the manuscript threshold ---
    up_control <- res %>%
      filter(p_val_adj < PADJ_CUTOFF, avg_log2FC >=  LOGFC_FOR_LISTS) %>% pull(gene_id) %>% unique()
    up_target  <- res %>%
      filter(p_val_adj < PADJ_CUTOFF, avg_log2FC <= -LOGFC_FOR_LISTS) %>% pull(gene_id) %>% unique()

    writeLines(up_control, file.path(OUTDIR, "deg_lists", paste0(tag, "_de_genes_up_in_control.txt")))
    writeLines(up_target,  file.path(OUTDIR, "deg_lists", paste0(tag, "_de_genes_up_in_target.txt")))

    el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
    cat(sprintf("  %-6s tested: %5d | up_control: %4d | up_target: %4d  (%ds)\n",
                stage, nrow(res), length(up_control), length(up_target), el))

    i_log <- i_log + 1L
    run_log[[i_log]] <- data.frame(
      cluster = cl, stage = stage, n_control = n_ctrl, n_target = n_tgt,
      n_tested = nrow(res), n_up_control = length(up_control),
      n_up_target = length(up_target), skipped = FALSE, stringsAsFactors = FALSE)
  }

  # --- per-cluster ORA universe: union of tested genes across the four
  # --- comparisons, fixed so the timepoints stay mutually comparable ---
  if (length(tested_union)) {
    write_csv(data.frame(gene_id = sort(tested_union)),
              file.path(OUTDIR, "universes", paste0("universe_", safe_name(cl), ".csv")))
    cat("  universe:", length(tested_union), "genes\n")
  }

  # write the log incrementally — this run is long enough to be interrupted
  write_tsv(bind_rows(run_log), file.path(OUTDIR, "de_rerun_log.tsv"))
}

log_df <- bind_rows(run_log)
write_tsv(log_df, file.path(OUTDIR, "de_rerun_log.tsv"))

cat("\n==================================================\n")
cat("DONE in", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "minutes\n")
cat("==================================================\n\n")
print(as.data.frame(log_df))

cat("\nOutputs:\n")
cat("  full ranked tables :", file.path(OUTDIR, "full"), "\n")
cat("  DEG lists          :", file.path(OUTDIR, "deg_lists"), "\n")
cat("  per-cluster universes:", file.path(OUTDIR, "universes"), "\n")

# ---------------------------------------------------------
# VERIFICATION — run this before trusting anything downstream
# ---------------------------------------------------------
# Compare a re-run DEG list against the equivalent original file. They
# should be identical. If they are not, the ident order or the condition
# labels differ between the two scripts.
#
#   old <- readLines("<de_plots_full_ann_30032026>/72hpi/Cluster_1_control_vs_72hpi/Cluster_1_control_vs_72hpi_de_genes_up_in_target.txt")
#   old <- unique(sub("^(\\S+)\\s.*$", "\\1", old))
#   new <- readLines(file.path(OUTDIR, "deg_lists", "Cluster_1_control_vs_72hpi_de_genes_up_in_target.txt"))
#   length(old); length(new); length(intersect(old, new))
