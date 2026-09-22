# =========================================================
# Section 6b: Cluster 1 identity — enrichment
# Environment: go-enrich (local machine)
# =========================================================
# Requires:
#   - 01b_build_reference_tables_v2.R run in this session
#   - 07a_cluster1_identity_markers.R already run in the Seurat session
#
# This is the replacement for the 02/03 marker-enrichment path for the
# specific purpose of annotating Cluster 1. It differs from 02/03 in
# three ways:
#   - uses the v2 reference tables (propagated GO, split ontologies,
#     deduplicated KEGG)
#   - uses control-only markers, so infection response does not leak into
#     the identity signature
#   - adds targeted pairwise enrichment against candidate identities,
#     each against its own matched universe
#
# WHAT ENRICHMENT CAN AND CANNOT TELL YOU HERE
# Enrichment is supporting evidence for a cell-type call, not the call
# itself. "Immune system process is enriched" narrows Cluster 1 to the
# immune compartment; it does not distinguish a haemocyte subtype from a
# macrophage-like cell. The pairwise DE gene counts from 07a, canonical
# marker expression, and cross-referencing to published C. gigas
# single-cell atlases are what actually make the annotation. Treat this
# script's output as the paragraph of supporting biology, not the answer.

library(readr)
library(dplyr)
library(tidyr)
library(clusterProfiler)

REF_DIR <- "/home/pdewari/eggnog/results/full_proteome_20260820_124651/reference_tables_v2"

t2g_BP         <- readRDS(file.path(REF_DIR, "t2g_BP.rds"))
t2g_MF         <- readRDS(file.path(REF_DIR, "t2g_MF.rds"))
t2g_CC         <- readRDS(file.path(REF_DIR, "t2g_CC.rds"))
term2name_go   <- readRDS(file.path(REF_DIR, "term2name_go.rds"))
kegg_term2gene <- readRDS(file.path(REF_DIR, "kegg_term2gene.rds"))
kegg_names     <- readRDS(file.path(REF_DIR, "kegg_names.rds"))
background_global <- readRDS(file.path(REF_DIR, "background_global.rds"))



# =========================================================
# CONFIG
# =========================================================
identity_dir <- "/home/pdewari/Documents/parse_2025/seurat_2025/cluster1_identity"
outdir       <- "/home/pdewari/eggnog/results/plots/cluster1_identity"

TARGET_CLUSTER <- "Cluster_1"
PADJ_CUTOFF    <- 0.05
MIN_GENES      <- 15
MIN_GS_SIZE    <- 10
MAX_GS_SIZE    <- 500
LOGFC_THRESH   <- 0.25

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

ontologies <- list(
  GO_BP = list(t2g = t2g_BP,         t2n = term2name_go),
  GO_MF = list(t2g = t2g_MF,         t2n = term2name_go),
  GO_CC = list(t2g = t2g_CC,         t2n = term2name_go),
  KEGG  = list(t2g = kegg_term2gene, t2n = kegg_names)
)

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE),
         function(y) as.numeric(y[1]) / as.numeric(y[2]), numeric(1))
}

run_one <- function(genes, t2g, t2n, universe) {
  uni <- intersect(universe, unique(t2g[[2]]))
  gl  <- intersect(genes, uni)
  if (length(gl) < MIN_GENES) return(NULL)

  res <- tryCatch(
    enricher(gene = gl, universe = uni, TERM2GENE = t2g, TERM2NAME = t2n,
             pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
             minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE),
    error = function(e) { cat("    failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(res)) return(NULL)

  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)

  df %>%
    mutate(fold_enrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio),
           n_input = length(gl)) %>%
    filter(p.adjust < PADJ_CUTOFF)
}

results <- list(); i <- 0L

# =========================================================
# 6b.1 One-vs-all identity enrichment (control cells only)
# =========================================================
markers_ctrl <- read_tsv(file.path(identity_dir, "markers_control_only.tsv"),
                         show_col_types = FALSE)
universe_markers <- read_csv(file.path(identity_dir, "universe_markers_control.csv"),
                             show_col_types = FALSE)$gene_id

target_genes <- markers_ctrl %>%
  filter(cluster == TARGET_CLUSTER, p_val_adj < 0.05, avg_log2FC > LOGFC_THRESH) %>%
  pull(gene_id) %>% unique()

cat("=== One-vs-all identity:", TARGET_CLUSTER, "===\n")
cat("Marker genes (control only):", length(target_genes), "\n")
cat("Universe:", length(universe_markers), "\n")

for (onto in names(ontologies)) {
  r <- run_one(target_genes, ontologies[[onto]]$t2g, ontologies[[onto]]$t2n, universe_markers)
  n <- if (is.null(r)) 0 else nrow(r)
  cat(sprintf("  %-6s -> %3d terms\n", onto, n))
  if (!is.null(r)) {
    i <- i + 1L
    results[[i]] <- r %>% mutate(analysis = "one_vs_all", cluster = TARGET_CLUSTER,
                                 contrast = "vs_all_other_clusters",
                                 ontology = onto, .before = 1)
  }
}

# --- context: the same for every other cluster, so Cluster 1's profile
# --- can be read against its neighbours rather than in isolation
RUN_ALL_CLUSTERS <- TRUE
if (RUN_ALL_CLUSTERS) {
  other_clusters <- setdiff(unique(markers_ctrl$cluster), TARGET_CLUSTER)
  cat("\n=== Context: one-vs-all for", length(other_clusters), "other clusters ===\n")
  for (cl in other_clusters) {
    g <- markers_ctrl %>%
      filter(cluster == cl, p_val_adj < 0.05, avg_log2FC > LOGFC_THRESH) %>%
      pull(gene_id) %>% unique()
    if (length(g) < MIN_GENES) { cat(sprintf("  %-24s skipped (%d genes)\n", cl, length(g))); next }
    n_tot <- 0
    for (onto in names(ontologies)) {
      r <- run_one(g, ontologies[[onto]]$t2g, ontologies[[onto]]$t2n, universe_markers)
      if (!is.null(r)) {
        i <- i + 1L
        results[[i]] <- r %>% mutate(analysis = "one_vs_all", cluster = cl,
                                     contrast = "vs_all_other_clusters",
                                     ontology = onto, .before = 1)
        n_tot <- n_tot + nrow(r)
      }
    }
    cat(sprintf("  %-24s %4d genes -> %3d terms\n", cl, length(g), n_tot))
  }
}

# =========================================================
# 6b.2 Targeted pairwise enrichment vs candidate identities
# =========================================================
pair_files <- list.files(identity_dir,
                         pattern = paste0("^pairwise_", safe_name(TARGET_CLUSTER), "_vs_.*\\.tsv$"),
                         full.names = TRUE)

cat("\n=== Pairwise contrasts:", length(pair_files), "===\n")

for (pf in pair_files) {
  cand <- sub("\\.tsv$", "", sub(paste0("^pairwise_", safe_name(TARGET_CLUSTER), "_vs_"), "",
                                 basename(pf)))
  uf <- file.path(identity_dir, paste0("universe_", safe_name(TARGET_CLUSTER), "_vs_", cand, ".csv"))
  if (!file.exists(uf)) { cat("  missing universe for", cand, "— skipping\n"); next }

  res_pair <- read_tsv(pf, show_col_types = FALSE)
  uni_pair <- read_csv(uf, show_col_types = FALSE)$gene_id

  up_target <- res_pair %>% filter(p_val_adj < 0.05, avg_log2FC >  LOGFC_THRESH) %>% pull(gene_id)
  up_cand   <- res_pair %>% filter(p_val_adj < 0.05, avg_log2FC < -LOGFC_THRESH) %>% pull(gene_id)

  cat(sprintf("\n  %s vs %s  (up in target: %d | up in candidate: %d)\n",
              TARGET_CLUSTER, cand, length(up_target), length(up_cand)))

  for (side in c("up_in_target", "up_in_candidate")) {
    g <- if (side == "up_in_target") up_target else up_cand
    for (onto in names(ontologies)) {
      r <- run_one(g, ontologies[[onto]]$t2g, ontologies[[onto]]$t2n, uni_pair)
      if (!is.null(r)) {
        i <- i + 1L
        results[[i]] <- r %>% mutate(analysis = "pairwise", cluster = TARGET_CLUSTER,
                                     contrast = paste0(side, "_vs_", cand),
                                     ontology = onto, .before = 1)
        cat(sprintf("    %-16s %-6s -> %3d terms\n", side, onto, nrow(r)))
      }
    }
  }
}

# =========================================================
# 6b.3 Outputs
# =========================================================
tidy_identity <- if (length(results)) bind_rows(results) else data.frame()
write_tsv(tidy_identity, file.path(outdir, "cluster1_identity_tidy.tsv"))

if (nrow(tidy_identity)) {

  # What is enriched in Cluster 1 that is NOT enriched in the candidates?
  # This is the discriminating signal — terms Cluster 1 shares with a
  # candidate say nothing about which it is.
  ova <- tidy_identity %>% filter(analysis == "one_vs_all")

  target_terms <- ova %>% filter(cluster == TARGET_CLUSTER) %>% pull(ID) %>% unique()

  shared <- ova %>%
    filter(cluster != TARGET_CLUSTER, ID %in% target_terms) %>%
    count(ID, name = "n_other_clusters_sharing")

  discriminating <- ova %>%
    filter(cluster == TARGET_CLUSTER) %>%
    dplyr::select(ontology, ID, Description, fold_enrichment, p.adjust, Count) %>%
    left_join(shared, by = "ID") %>%
    mutate(n_other_clusters_sharing = coalesce(n_other_clusters_sharing, 0L)) %>%
    arrange(n_other_clusters_sharing, p.adjust)

  write_tsv(discriminating, file.path(outdir, "cluster1_discriminating_terms.tsv"))

  cat("\n=== Most discriminating terms for", TARGET_CLUSTER, "===\n")
  cat("(enriched in Cluster 1, shared with fewest other clusters)\n\n")
  print(as.data.frame(head(discriminating, 20)))
}

cat("\nTidy table:", file.path(outdir, "cluster1_identity_tidy.tsv"), "\n")
cat("Discriminating terms:", file.path(outdir, "cluster1_discriminating_terms.tsv"), "\n")

#############

eggnog_full <- readr::read_tsv(
  "/home/pdewari/eggnog/results/full_proteome_20260820_124651/full_proteome.emapper.annotations",
  comment = "##", show_col_types = FALSE) %>%
  dplyr::mutate(gene_id = sub("\\..*$", "", sub("^transcript:", "", `#query`)))

ann <- eggnog_full %>%
  dplyr::select(gene_id, Preferred_name, Description, PFAMs, COG_category) %>%
  dplyr::distinct()

tidy_identity %>%
  filter(analysis == "one_vs_all", cluster == "Cluster_1", ID %in% key) %>%
  dplyr::select(ID, Description_term = Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::distinct(geneID) %>%
  left_join(ann, by = c("geneID" = "gene_id")) %>%
  print(n = 30)

#######
birc_all <- eggnog_full %>%
  filter(Preferred_name == "birc2" | grepl("\\bBIR\\b", PFAMs)) %>%
  pull(gene_id) %>% unique()
length(birc_all)

# which are Cluster 1 control markers
c1_markers <- markers_ctrl %>% filter(cluster == "Cluster_1", p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene_id)
intersect(birc_all, c1_markers)

de_dir <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_plots_full_ann_30032026"

inf <- unique(unlist(lapply(c("6hpi","24hpi","72hpi","96hpi"), function(s) {
  d <- list.dirs(file.path(de_dir, s), full.names = TRUE, recursive = FALSE)
  d <- d[grepl("^Cluster.?1_control_vs_", basename(d))]
  if (!length(d)) return(character(0))
  f <- list.files(d, pattern = "de_genes_up_in_target\\.txt$", full.names = TRUE)
  if (!length(f)) return(character(0))
  sub("^(\\S+)\\s.*$", "\\1", readLines(f[1]))
})))

length(inf)   # expect 308
# which are infection-responsive
intersect(birc_all, inf)

##############
markers_ctrl %>%
  filter(gene_id %in% birc_all, p_val_adj < 0.05, avg_log2FC > 0) %>%
  dplyr::select(cluster, gene_id) %>%
  arrange(cluster) %>%
  print(n = 60)
###
for (cl in c("Hepatopancreas_cells", "Gill_ciliary_cells")) {
  g <- unique(unlist(lapply(c("6hpi","24hpi","72hpi","96hpi"), function(s) {
    d <- list.dirs(file.path(de_dir, s), full.names = TRUE, recursive = FALSE)
    d <- d[grepl(paste0("^", cl, "_control_vs_"), basename(d))]
    if (!length(d)) return(character(0))
    f <- list.files(d, pattern = "de_genes_up_in_target\\.txt$", full.names = TRUE)
    if (!length(f)) return(character(0))
    sub("^(\\S+)\\s.*$", "\\1", readLines(f[1]))
  })))
  cat(cl, ":", paste(intersect(birc_all, g), collapse = ", "), "\n")
}

#
sum(birc_all %in% universe_markers)   # how many of the 47 were testable at all
#####
eggnog_full %>%
  filter(gene_id %in% birc_all) %>%
  count(PFAMs, sort = TRUE) %>%
  print(n = 20)

##############
markers <- readr::read_tsv(
  "/home/pdewari/Documents/parse_2025/seurat_2025/df_all_markers_20Jan26_new.tsv",
  show_col_types = FALSE) %>%
  dplyr::mutate(gene_id = sub("^(\\S+)\\s.*$", "\\1", Gene.ID),
                across(any_of(c("p_val", "p_val_adj")), as.numeric))

nrow(markers)                 # expect 8960
class(markers$p_val_adj)      # want "numeric"


combined_c1 <- markers %>%
  mutate(cluster = safe_name(as.character(cluster))) %>%
  filter(cluster == "Cluster_1", avg_log2FC > 0) %>%
  pull(gene_id) %>% unique()

intersect(c("G19035","G20059","G25799"), combined_c1)   # inducible IAPs as "markers"?
intersect(birc_all, combined_c1)                        # which IAPs appear

