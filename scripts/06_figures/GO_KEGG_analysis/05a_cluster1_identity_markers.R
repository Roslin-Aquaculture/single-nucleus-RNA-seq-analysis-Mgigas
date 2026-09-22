# =========================================================
# Section 6a: Cluster 1 identity — marker tables + universes
# Environment: SEURAT session (not go-enrich)
# =========================================================
# Answers "what biology defines Cluster 1?" — a different question from
# "what changes in Cluster 1 during infection", and one the DE-stage
# enrichment (04/05/05b) cannot answer at all.
#
# THREE PROBLEMS THIS FIXES IN THE 02/03 MARKER PATH
#
#   1. INFECTION-STATE CONTAMINATION.
#      df_all_markers was computed across ALL cells — control and infected
#      pooled. So Cluster 1's "identity" signature contains the very
#      infection-response genes (Ced-3, TRIMs, E3 ligases, Malt1) that the
#      DE analysis reports separately. For a cluster you are still trying
#      to annotate, that is circular: the cluster looks immune-activated
#      partly because infected cells are in it. Identity markers are
#      computed here on CONTROL CELLS ONLY. The infection response is then
#      a genuinely independent result.
#
#   2. ONE-VS-ALL DILUTION.
#      FindAllMarkers compares each cluster against every other cell. If
#      Cluster 1 is haemocyte-adjacent, every gene it shares with the
#      haemocyte and macrophage-like populations is suppressed — those are
#      exactly the genes that would tell you what it is. Targeted pairwise
#      comparisons against the candidate identities are far sharper.
#
#   3. WRONG UNIVERSE, AGAIN.
#      A pairwise Cluster_1-vs-Haemocytes comparison has its own eligible
#      gene set: the union of what is testable in those two clusters. Not
#      the whole object.
#
# Outputs feed 07b_cluster1_identity_enrichment.R.

library(Seurat)
library(Matrix)
library(dplyr)
library(readr)

# =========================================================
# CONFIG
# =========================================================
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
#####################################


table(seurat_obj_clean$condition_new, useNA = "ifany")
table(seurat_obj_clean$sample, seurat_obj_clean$condition_new)
table(seurat_obj_clean$cluster_annotation)
ncol(seurat_obj_clean)
##############################################

SEURAT_OBJ <- seurat_obj_clean



# RenameIdents() changes Idents(), NOT a metadata column. Create the column
# from the active identities if it is not already there.
if (!"cluster_annotation" %in% colnames(SEURAT_OBJ@meta.data)) {
  SEURAT_OBJ$cluster_annotation <- as.character(Idents(SEURAT_OBJ))
  message("Created cluster_annotation from Idents()")
}

CLUSTER_COL   <- "cluster_annotation"
CONDITION_COL <- "condition_new"      # your simplified column
CONTROL_LABEL <- "control"            # pools Homogenate + Uninfected

# >>> LIVE object names — spaces and a hyphen, NOT the underscored form
# >>> used in df_all_markers and the DE output directories. Filenames are
# >>> sanitised automatically by safe_name().
TARGET_CLUSTER <- "Cluster 1"

# Candidate identities, taken from the actual cluster labels in
# df_all_markers_20Jan26.xlsx. Your dataset has SIX haemocyte-lineage
# populations, which is exactly why one-vs-all dilutes Cluster 1: every
# gene it shares with any of them is suppressed by the "all others"
# comparison. The pairwise contrasts below are what resolve it.
CANDIDATE_CLUSTERS <- c("Haemocyte cell type 1",
                        "Immature haemocytes",
                        "Vesicular haemocytes",
                        "Hyalinocytes",
                        "Macrophage like cells",
                        "Small granule cells",
                        "Hepatopancreas cells",
                        "Gill ciliary cells")

# Matches your marker call: FindAllMarkers(min.pct = 0.20, ...)
MIN_PCT       <- 0.20

# CONFIRMED FROM df_all_markers_20Jan26.xlsx: min |avg_log2FC| = 0.1004,
# so the applied threshold was 0.1, NOT the 0.25 that was intended. The
# original call passed `log2fc.threshold = 0.25`, which is not a Seurat
# argument (it is `logfc.threshold`, no "2"); R could not partial-match it,
# so it was absorbed into `...` and Seurat's v5 default of 0.1 applied.
#
# Set to 0.1 here so the control-only markers are directly comparable with
# the existing all-cells marker table. See the note in the reply about what
# this means for the Methods text.
LOGFC_THRESH  <- 0.1
MIN_CELLS     <- 20      # skip a cluster with fewer control cells than this

ASSAY  <- "RNA"
OUTDIR <- "/home/pdewari/Documents/parse_2025/seurat_2025/cluster1_identity"

# =========================================================
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
clean_id  <- function(x) sub("^(\\S+)\\s.*$", "\\1", x)

get_counts <- function(obj, assay) {
  tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
           error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = "counts"))
}

meta <- SEURAT_OBJ@meta.data
stopifnot(CLUSTER_COL %in% colnames(meta), CONDITION_COL %in% colnames(meta))

cat("Clusters in object:\n")
print(table(meta[[CLUSTER_COL]]))

# ---------------------------------------------------------
# 6a.1 Control-only subset
# ---------------------------------------------------------
control_cells <- rownames(meta)[as.character(meta[[CONDITION_COL]]) == CONTROL_LABEL]
cat("\nControl cells:", length(control_cells), "\n")
if (length(control_cells) == 0) stop("No cells matched CONTROL_LABEL = ", CONTROL_LABEL)

ctrl <- subset(SEURAT_OBJ, cells = control_cells)
Idents(ctrl) <- ctrl@meta.data[[CLUSTER_COL]]

ctrl_counts <- get_counts(ctrl, ASSAY)
ctrl_meta   <- ctrl@meta.data
gene_ids    <- clean_id(rownames(ctrl_counts))

cat("Control cells per cluster:\n")
print(table(Idents(ctrl)))

keep_clusters <- names(which(table(Idents(ctrl)) >= MIN_CELLS))
cat("\nClusters with >=", MIN_CELLS, "control cells:", length(keep_clusters), "\n")
if (!TARGET_CLUSTER %in% keep_clusters) {
  stop(TARGET_CLUSTER, " has fewer than ", MIN_CELLS,
       " control cells — identity markers cannot be computed on controls alone.")
}

ctrl <- subset(ctrl, idents = keep_clusters)

# ---------------------------------------------------------
# 6a.2 One-vs-all identity markers, CONTROL CELLS ONLY
# ---------------------------------------------------------
cat("\nRunning FindAllMarkers on control cells...\n")
markers_ctrl <- FindAllMarkers(
  ctrl,
  only.pos        = TRUE,
  min.pct         = MIN_PCT,
  logfc.threshold = LOGFC_THRESH
)

# normalise cluster names to the underscored form on write, so everything
# downstream (07b, figures, the DE pipeline) shares one naming convention
markers_ctrl <- markers_ctrl %>%
  mutate(gene_id = clean_id(gene),
         cluster = safe_name(as.character(cluster)))
write_tsv(markers_ctrl, file.path(OUTDIR, "markers_control_only.tsv"))
cat("Control-only markers written:", nrow(markers_ctrl), "rows\n")

# universe for one-vs-all markers within control cells:
# genes reaching min.pct in at least one retained cluster
pct_by_cluster <- sapply(keep_clusters, function(cl) {
  cells <- rownames(ctrl_meta)[as.character(ctrl_meta[[CLUSTER_COL]]) == cl]
  cells <- intersect(cells, colnames(ctrl_counts))
  if (length(cells) == 0) return(rep(0, nrow(ctrl_counts)))
  Matrix::rowSums(ctrl_counts[, cells, drop = FALSE] > 0) / length(cells)
})
universe_markers <- unique(gene_ids[apply(pct_by_cluster, 1, max) >= MIN_PCT])
write_csv(data.frame(gene_id = universe_markers),
          file.path(OUTDIR, "universe_markers_control.csv"))
cat("Marker universe (control cells):", length(universe_markers), "genes\n")

# ---------------------------------------------------------
# 6a.3 STABILITY CHECK — does pooling infected cells change the signature?
# ---------------------------------------------------------
# If Cluster 1's top markers differ substantially between control-only and
# all-cells, its apparent identity is partly infection state. Report this
# number; it is a direct answer to a reviewer asking whether Cluster 1 is
# a cell type or an activation state.
# df_all_markers uses underscored cluster names; markers_ctrl is already
# normalised above. Compare on the sanitised form so both sides match.
TARGET_SAFE <- safe_name(TARGET_CLUSTER)
markers <- readxl::read_excel("/home/pdewari/Documents/parse_2025/seurat_2025/df_all_markers_20Jan26.xlsx") %>%
  mutate(gene_id = sub("^(\\S+)\\s.*$", "\\1", Gene.ID))

if (exists("markers")) {
  top_ctrl <- markers_ctrl %>%
    filter(cluster == TARGET_SAFE, p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>% head(100) %>% pull(gene_id)

  top_all <- markers %>%
    mutate(cluster = safe_name(as.character(cluster))) %>%
    filter(cluster == TARGET_SAFE, p_val_adj < 0.05, avg_log2FC > LOGFC_THRESH) %>%
    arrange(desc(avg_log2FC)) %>% head(100) %>% pull(gene_id)

  ov <- length(intersect(top_ctrl, top_all))
  cat("\n--- marker stability for", TARGET_CLUSTER, "---\n")
  cat("Top-100 overlap, control-only vs all-cells:", ov, "/ 100\n")
  cat("Control-only unique:", paste(head(setdiff(top_ctrl, top_all), 15), collapse = ", "), "\n")
  cat("All-cells unique   :", paste(head(setdiff(top_all, top_ctrl), 15), collapse = ", "), "\n")

  write_tsv(data.frame(
    metric = c("top100_overlap", "control_only_unique", "all_cells_unique"),
    value  = c(ov, length(setdiff(top_ctrl, top_all)), length(setdiff(top_all, top_ctrl)))
  ), file.path(OUTDIR, "marker_stability.tsv"))
} else {
  cat("\n(Load your df_all_markers table as `markers` to get the stability check.)\n")
}

# ---------------------------------------------------------
# 6a.4 PER-ANIMAL CONSISTENCY — is Cluster 1 a cell type or a batch?
# ---------------------------------------------------------
# You already noted that sub-clustering Cluster 1 partitions almost
# entirely by individual animal. Before calling it a cell type, check
# whether its markers reproduce in each control animal separately.
# Set SAMPLE_COL to your per-animal metadata column to enable this.
# Your per-animal column is `sample`. Controls are two animals (Homogenate
# = mock, Uninfected = unchallenged), so this check is meaningful.
SAMPLE_COL <- "sample"

if (!is.null(SAMPLE_COL) && SAMPLE_COL %in% colnames(ctrl@meta.data)) {
  animals <- unique(as.character(ctrl@meta.data[[SAMPLE_COL]]))
  per_animal <- list()
  for (a in animals) {
    cells_a <- rownames(ctrl@meta.data)[as.character(ctrl@meta.data[[SAMPLE_COL]]) == a]
    sub_a <- subset(ctrl, cells = cells_a)
    if (sum(Idents(sub_a) == TARGET_CLUSTER) < MIN_CELLS) next
    m <- tryCatch(
      FindMarkers(sub_a, ident.1 = TARGET_CLUSTER, only.pos = TRUE,
                  min.pct = MIN_PCT, logfc.threshold = LOGFC_THRESH),
      error = function(e) NULL
    )
    if (!is.null(m) && nrow(m) > 0) {
      per_animal[[a]] <- clean_id(rownames(m)[m$p_val_adj < 0.05])
    }
  }
  if (length(per_animal) > 1) {
    core <- Reduce(intersect, per_animal)
    cat("\n--- per-animal marker consistency ---\n")
    for (a in names(per_animal)) cat(sprintf("  %-16s %5d markers\n", a, length(per_animal[[a]])))
    cat("  Shared by ALL control animals:", length(core), "\n")
    write_csv(data.frame(gene_id = core), file.path(OUTDIR, "cluster1_core_markers_all_animals.csv"))
    cat("  -> use this core set for the identity claim; it cannot be a single-animal artefact.\n")
  }
}

# ---------------------------------------------------------
# 6a.5 Targeted pairwise comparisons + matched universes
# ---------------------------------------------------------
present_candidates <- intersect(CANDIDATE_CLUSTERS, keep_clusters)
missing <- setdiff(CANDIDATE_CLUSTERS, keep_clusters)
if (length(missing)) cat("\nCandidates not usable (absent or too few control cells):",
                         paste(missing, collapse = ", "), "\n")

pair_summary <- data.frame()

for (cand in present_candidates) {

  cat("\n", TARGET_CLUSTER, "vs", cand, "\n")

  res <- FindMarkers(ctrl, ident.1 = TARGET_CLUSTER, ident.2 = cand,
                     min.pct = MIN_PCT, logfc.threshold = LOGFC_THRESH)
  res$gene_id <- clean_id(rownames(res))

  f <- file.path(OUTDIR, paste0("pairwise_", safe_name(TARGET_CLUSTER),
                                "_vs_", safe_name(cand), ".tsv"))
  write_tsv(res, f)

  # matched universe: genes testable in EITHER of the two clusters
  cells_t <- rownames(ctrl_meta)[as.character(ctrl_meta[[CLUSTER_COL]]) == TARGET_CLUSTER]
  cells_c <- rownames(ctrl_meta)[as.character(ctrl_meta[[CLUSTER_COL]]) == cand]
  cells_t <- intersect(cells_t, colnames(ctrl_counts))
  cells_c <- intersect(cells_c, colnames(ctrl_counts))

  pct_t <- Matrix::rowSums(ctrl_counts[, cells_t, drop = FALSE] > 0) / length(cells_t)
  pct_c <- Matrix::rowSums(ctrl_counts[, cells_c, drop = FALSE] > 0) / length(cells_c)
  uni   <- unique(gene_ids[pmax(pct_t, pct_c) >= MIN_PCT])

  write_csv(data.frame(gene_id = uni),
            file.path(OUTDIR, paste0("universe_", safe_name(TARGET_CLUSTER),
                                     "_vs_", safe_name(cand), ".csv")))

  n_up   <- sum(res$p_val_adj < 0.05 & res$avg_log2FC >  LOGFC_THRESH)
  n_down <- sum(res$p_val_adj < 0.05 & res$avg_log2FC < -LOGFC_THRESH)

  cat(sprintf("  up in %s: %d | up in %s: %d | universe: %d\n",
              TARGET_CLUSTER, n_up, cand, n_down, length(uni)))

  pair_summary <- rbind(pair_summary, data.frame(
    target = TARGET_CLUSTER, candidate = cand,
    n_cells_target = length(cells_t), n_cells_candidate = length(cells_c),
    n_up_target = n_up, n_up_candidate = n_down,
    n_universe = length(uni), stringsAsFactors = FALSE
  ))
}

write_tsv(pair_summary, file.path(OUTDIR, "pairwise_summary.tsv"))

cat("\n--- pairwise summary ---\n")
print(pair_summary)
cat("\nREAD THIS: the candidate with the FEWEST differentially expressed\n")
cat("genes against", TARGET_CLUSTER, "is its closest transcriptional neighbour.\n")
cat("That number is the single most informative output for annotation —\n")
cat("more than any enrichment term.\n")

cat("\nAll outputs in:", OUTDIR, "\n")
cat("Next: 07b_cluster1_identity_enrichment.R (go-enrich session)\n")


######
