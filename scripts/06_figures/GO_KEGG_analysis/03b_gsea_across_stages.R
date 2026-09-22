# =========================================================
# Session reload — go-enrich environment
# =========================================================
# Run this after restarting R in the go-enrich conda env, before
# 05b / 05c / 06 (or any 07b follow-up work).
#
#   source("/path/to/00_reload_go_enrich_session.R")
#
# Loads the reference tables built by 01b plus, optionally, the objects
# used in the Cluster 1 identity work. Prints a summary so you can see at
# a glance that everything arrived.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(ggplot2)
  library(stringr)
})

# =========================================================
# PATHS — edit if yours differ
# =========================================================
REF_DIR      <- "/home/pdewari/eggnog/results/full_proteome_20260820_124651/reference_tables_v2"
EGGNOG_FILE  <- "/home/pdewari/eggnog/results/full_proteome_20260820_124651/full_proteome.emapper.annotations"
# 04b output used for PATHWAY analysis (min.pct = 0.1)
DE_DIR       <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_full_ranked_minpct01"
# the min.pct = 0.25 run, kept for gene-level results that match the manuscript
DE_DIR_025   <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_full_ranked"
OLD_DE_DIR   <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_plots_full_ann_30032026"
IDENTITY_DIR <- "/home/pdewari/Documents/parse_2025/seurat_2025/cluster1_identity"     # 07a output
MARKERS_TSV  <- "/home/pdewari/Documents/parse_2025/seurat_2025/df_all_markers_20Jan26_new.tsv"

# Set FALSE to skip the Cluster 1 identity objects (not needed for 05b/05c/06)
LOAD_IDENTITY <- TRUE

# =========================================================
# SHARED HELPERS  (same definitions used across the pipeline)
# =========================================================
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
clean_id  <- function(x) sub("^(\\S+)\\s.*$", "\\1", x)

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE),
         function(y) as.numeric(y[1]) / as.numeric(y[2]), numeric(1))
}

# =========================================================
# 1. REFERENCE TABLES (required for 05b / 05c / 07b)
# =========================================================
stopifnot(dir.exists(REF_DIR))

t2g_BP            <- readRDS(file.path(REF_DIR, "t2g_BP.rds"))
t2g_MF            <- readRDS(file.path(REF_DIR, "t2g_MF.rds"))
t2g_CC            <- readRDS(file.path(REF_DIR, "t2g_CC.rds"))
term2name_go      <- readRDS(file.path(REF_DIR, "term2name_go.rds"))
kegg_term2gene    <- readRDS(file.path(REF_DIR, "kegg_term2gene.rds"))
kegg_names        <- readRDS(file.path(REF_DIR, "kegg_names.rds"))
background_global <- readRDS(file.path(REF_DIR, "background_global.rds"))

cat("=== reference tables ===\n")
cat(sprintf("  GO BP : %6d terms | %5d genes\n",
            dplyr::n_distinct(t2g_BP$GO_ID), dplyr::n_distinct(t2g_BP$gene_id)))
cat(sprintf("  GO MF : %6d terms | %5d genes\n",
            dplyr::n_distinct(t2g_MF$GO_ID), dplyr::n_distinct(t2g_MF$gene_id)))
cat(sprintf("  GO CC : %6d terms | %5d genes\n",
            dplyr::n_distinct(t2g_CC$GO_ID), dplyr::n_distinct(t2g_CC$gene_id)))
cat(sprintf("  KEGG  : %6d paths | %5d genes\n",
            dplyr::n_distinct(kegg_term2gene$pathway_id),
            dplyr::n_distinct(kegg_term2gene$gene_id)))
cat(sprintf("  global background: %d genes\n", length(background_global)))

# =========================================================
# 2. 04b OUTPUT — check it is where 05b / 05c expect it
# =========================================================
cat("\n=== 04b output (min.pct = 0.1, for GSEA) ===\n")
if (dir.exists(DE_DIR)) {
  n_full <- length(list.files(file.path(DE_DIR, "full"), pattern = "_full\\.tsv$"))
  cat(sprintf("  ranked tables: %3d\n", n_full))
  if (n_full == 0) cat("  WARNING: no ranked tables — 05c (GSEA) will stop.\n")
} else {
  cat("  NOT FOUND:", DE_DIR, "\n  Run 04b in the Seurat session first.\n")
}

# ORA was dropped: with ~25 annotated genes per thresholded list against
# 5,234 testable BP terms it could not clear multiple-testing correction.
# Make sure no stale ORA output is left where 06 would pick it up and plot
# it alongside GSEA results from a different DE run.
stale <- "/home/pdewari/eggnog/results/enrichment/ora_across_stages/ORA_tidy.tsv"
if (file.exists(stale)) {
  cat("\n  WARNING: stale ORA output present at\n    ", stale, "\n")
  cat("  06 will plot it alongside GSEA. Move or delete it first.\n")
}

# =========================================================
# 3. OPTIONAL — Cluster 1 identity objects
# =========================================================
if (LOAD_IDENTITY) {
  
  cat("\n=== identity objects ===\n")
  
  # eggNOG annotation: gene symbols, descriptions, domains
  if (file.exists(EGGNOG_FILE)) {
    eggnog_full <- read_tsv(EGGNOG_FILE, comment = "##", show_col_types = FALSE) %>%
      mutate(gene_id = sub("\\..*$", "", sub("^transcript:", "", `#query`)))
    
    ann <- eggnog_full %>%
      dplyr::select(gene_id, Preferred_name, Description, PFAMs, COG_category) %>%
      distinct()
    
    birc_all <- eggnog_full %>%
      filter(Preferred_name == "birc2" | grepl("\\bBIR\\b", PFAMs)) %>%
      pull(gene_id) %>% unique()
    
    cat(sprintf("  eggnog_full: %d rows | ann: %d genes | birc_all: %d IAPs\n",
                nrow(eggnog_full), nrow(ann), length(birc_all)))
  } else {
    cat("  eggNOG annotation not found at", EGGNOG_FILE, "\n")
  }
  
  # control-only markers from 07a
  f <- file.path(IDENTITY_DIR, "markers_control_only.tsv")
  if (file.exists(f)) {
    markers_ctrl     <- read_tsv(f, show_col_types = FALSE)
    universe_markers <- read_csv(file.path(IDENTITY_DIR, "universe_markers_control.csv"),
                                 show_col_types = FALSE)$gene_id
    cat(sprintf("  markers_ctrl: %d rows | universe_markers: %d genes\n",
                nrow(markers_ctrl), length(universe_markers)))
  } else {
    cat("  07a output not found at", IDENTITY_DIR, "\n")
  }
  
  # combined-object markers — force numeric p-values.
  # Reading this table from the .xlsx returns p_val_adj as CHARACTER, which
  # makes `p_val_adj < 0.05` a string comparison that silently keeps only
  # the LEAST significant rows. The coercion below prevents that.
  if (file.exists(MARKERS_TSV)) {
    markers <- read_tsv(MARKERS_TSV, show_col_types = FALSE) %>%
      mutate(gene_id = clean_id(Gene.ID),
             across(any_of(c("p_val", "p_val_adj")), as.numeric))
    cat(sprintf("  markers: %d rows | p_val_adj is %s\n",
                nrow(markers), class(markers$p_val_adj)[1]))
  }
  
  # Cluster 1 infection-induced genes, from the ORIGINAL DE output
  if (dir.exists(OLD_DE_DIR)) {
    inf <- unique(unlist(lapply(c("6hpi","24hpi","72hpi","96hpi"), function(s) {
      d <- list.dirs(file.path(OLD_DE_DIR, s), full.names = TRUE, recursive = FALSE)
      d <- d[grepl("^Cluster.?1_control_vs_", basename(d))]
      if (!length(d)) return(character(0))
      fl <- list.files(d, pattern = "de_genes_up_in_target\\.txt$", full.names = TRUE)
      if (!length(fl)) return(character(0))
      clean_id(readLines(fl[1]))
    })))
    cat(sprintf("  inf (Cluster 1 up in infected): %d genes\n", length(inf)))
  }
  
  # enrichment results from 07b, if already written
  f <- "/home/pdewari/eggnog/results/plots/cluster1_identity/cluster1_identity_tidy.tsv"
  if (file.exists(f)) {
    tidy_identity <- read_tsv(f, show_col_types = FALSE)
    cat(sprintf("  tidy_identity: %d rows\n", nrow(tidy_identity)))
  }
}

cat("\nReady. Next: 05b -> 05c -> 06\n")



# =========================================================
# Section 5c: GSEA across stages
# Environment: go-enrich (local machine)
# =========================================================
# Requires in this session:
#   - 01b_build_reference_tables_v2.R  (t2g_BP/MF/CC, term2name_go,
#                                       kegg_term2gene, kegg_names)
# Requires on disk:
#   - 04b_rerun_de_full_ranked.R output (full/ ranked tables)
#
# WHY GSEA IS THE BETTER ANALYSIS FOR "WHAT PATHWAYS CHANGE OVER TIME"
#
#   - No arbitrary threshold. Your DE used logfc.threshold = 1, so ORA
#     only ever sees 2-fold changes. A coordinated 1.5-fold shift across
#     thirty genes of one pathway is a real response and is invisible to
#     ORA. GSEA ranks every expressed gene and finds exactly that.
#   - No background question. The ranked list IS the universe, so the
#     per-cluster-vs-global problem does not arise at all. Whatever a
#     reviewer thinks about background choice, it cannot touch this
#     result.
#   - Better suited to snRNA-seq effect sizes, which are modest and
#     spread across many genes rather than concentrated in a few.
#
# >>> DIRECTION CONVENTION — READ THIS <<<
# 04b keeps your original ident order (ident.1 = control), so a POSITIVE
# avg_log2FC means higher in CONTROL. That is the opposite of how a
# volcano plot is normally read, so this script NEGATES the ranking
# statistic. In every output below:
#
#      NES > 0  =  enriched in INFECTED  =  "up_in_target"
#      NES < 0  =  enriched in CONTROL   =  "up_in_control"
#
# The `direction` column states this explicitly for every row. Do not
# re-derive it from avg_log2FC without accounting for the negation.

library(readr)
library(dplyr)
library(tidyr)
library(clusterProfiler)

# =========================================================
# CONFIG
# =========================================================
de_dir <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_full_ranked_minpct01"
outdir <- "/home/pdewari/eggnog/results/enrichment/gsea_across_stages"

stages <- c("6hpi", "24hpi", "72hpi", "96hpi")

PADJ_CUTOFF  <- 0.05

# Lowered from 10. With the relaxed min.pct the ranked lists are larger,
# but many GO terms still contribute only a handful of genes to them. 5 is
# a common floor for modestly sized ranked lists; state it in Methods.
MIN_GS_SIZE  <- 5

MAX_GS_SIZE  <- 500
MIN_RANKED   <- 200     # skip a comparison with too few ranked genes
SEED         <- 42

# Ranking metric:
#   "logfc"  -> avg_log2FC (simple, interpretable, the usual choice)
#   "signed_p" -> sign(log2FC) * -log10(p_val); more sensitive but
#                 p-values from per-nucleus Wilcoxon are extreme and tied,
#                 so use it only as a sensitivity check
RANK_METRIC <- "logfc"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

set.seed(SEED)

ontologies <- list(
  GO_BP = list(t2g = t2g_BP,         t2n = term2name_go),
  GO_MF = list(t2g = t2g_MF,         t2n = term2name_go),
  GO_CC = list(t2g = t2g_CC,         t2n = term2name_go),
  KEGG  = list(t2g = kegg_term2gene, t2n = kegg_names)
)

# =========================================================
# BUILD RANKED LIST
# =========================================================
make_ranked <- function(res) {

  stat <- if (RANK_METRIC == "signed_p") {
    sign(res$avg_log2FC) * -log10(pmax(res$p_val, .Machine$double.xmin))
  } else {
    res$avg_log2FC
  }

  d <- data.frame(gene_id = res$gene_id, stat = stat, stringsAsFactors = FALSE) %>%
    filter(!is.na(stat), is.finite(stat))

  # one value per gene (cleaned IDs can collide); keep the strongest
  d <- d %>%
    group_by(gene_id) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()

  # NEGATE so positive = up in infected (see the header note)
  v <- setNames(-d$stat, d$gene_id)
  sort(v, decreasing = TRUE)
}

# =========================================================
# MAIN LOOP
# =========================================================
full_files <- list.files(file.path(de_dir, "full"), pattern = "_full\\.tsv$", full.names = TRUE)
if (length(full_files) == 0) stop("No ranked tables found — run 04b first.")

cat("Ranked tables found:", length(full_files), "\n")
cat("Ranking metric:", RANK_METRIC, "(negated: NES > 0 = up in infected)\n\n")

results <- list(); logs <- list(); i_r <- 0L; i_l <- 0L

for (f in full_files) {

  res <- read_tsv(f, show_col_types = FALSE)
  cl    <- unique(res$cluster)[1]
  stage <- unique(res$stage)[1]

  ranked <- make_ranked(res)

  cat(sprintf("%-26s %-6s  ranked genes: %5d\n", cl, stage, length(ranked)))

  if (length(ranked) < MIN_RANKED) {
    cat("    too few ranked genes — skipping\n")
    next
  }

  for (onto in names(ontologies)) {

    gres <- tryCatch(
      GSEA(geneList      = ranked,
           TERM2GENE     = ontologies[[onto]]$t2g,
           TERM2NAME     = ontologies[[onto]]$t2n,
           minGSSize     = MIN_GS_SIZE,
           maxGSSize     = MAX_GS_SIZE,
           pvalueCutoff  = 1,          # filter on p.adjust below
           pAdjustMethod = "BH",
           eps           = 0,          # accurate small p-values
           seed          = TRUE,
           verbose       = FALSE),
      error = function(e) { cat("    ", onto, "failed:", conditionMessage(e), "\n"); NULL }
    )

    n_sig <- 0L
    if (!is.null(gres)) {
      df <- as.data.frame(gres)
      if (nrow(df) > 0) {
        df <- df %>%
          filter(p.adjust < PADJ_CUTOFF) %>%
          mutate(direction = ifelse(NES > 0, "up_in_target", "up_in_control"))
        n_sig <- nrow(df)
        if (n_sig > 0) {
          i_r <- i_r + 1L
          results[[i_r]] <- df %>%
            mutate(method = "GSEA", cluster = cl, stage = stage,
                   ontology = onto, .before = 1)
        }
      }
    }

    cat(sprintf("    %-6s -> %3d terms\n", onto, n_sig))

    i_l <- i_l + 1L
    logs[[i_l]] <- data.frame(
      cluster = cl, stage = stage, ontology = onto,
      n_ranked = length(ranked), n_sig_terms = n_sig,
      rank_metric = RANK_METRIC, stringsAsFactors = FALSE)
  }

  if (length(results)) write_tsv(bind_rows(results), file.path(outdir, "GSEA_tidy.tsv"))
  if (length(logs))    write_tsv(bind_rows(logs),    file.path(outdir, "gsea_run_log.tsv"))
}

tidy_gsea <- if (length(results)) bind_rows(results) else data.frame()
write_tsv(tidy_gsea, file.path(outdir, "GSEA_tidy.tsv"))
write_tsv(bind_rows(logs), file.path(outdir, "gsea_run_log.tsv"))

# =========================================================
# SUMMARIES
# =========================================================
if (nrow(tidy_gsea)) {

  cat("\n--- terms per cluster x direction ---\n")
  print(as.data.frame(
    tidy_gsea %>% count(cluster, direction, ontology) %>%
      pivot_wider(names_from = direction, values_from = n, values_fill = 0)
  ))

  # Which pathways are SHARED across the responsive clusters, and which
  # are cluster-specific? This is the table that decides whether the
  # merged results sections are organised by shared core vs specific.
  shared <- tidy_gsea %>%
    group_by(ontology, ID, Description, direction) %>%
    summarise(n_clusters = n_distinct(cluster),
              clusters   = paste(sort(unique(cluster)), collapse = "; "),
              n_stages   = n_distinct(stage),
              stages     = paste(sort(unique(stage)), collapse = "; "),
              median_NES = median(NES),
              best_padj  = min(p.adjust), .groups = "drop") %>%
    arrange(desc(n_clusters), best_padj)

  write_tsv(shared, file.path(outdir, "GSEA_shared_vs_specific.tsv"))

  cat("\n--- most broadly shared terms (the conserved core) ---\n")
  print(as.data.frame(head(shared %>% filter(ontology == "GO_BP"), 15)))

  # supplementary table with leading-edge genes intact
  supp <- tidy_gsea %>%
    arrange(cluster, direction, ontology, stage, p.adjust) %>%
    transmute(Cluster = cluster, Direction = direction, Ontology = ontology,
              Stage = stage, Term_ID = ID, Term = Description,
              Set_size = setSize, NES = round(NES, 3),
              P_value = pvalue, P_adjusted = p.adjust,
              Leading_edge_genes = core_enrichment)
  write_tsv(supp, file.path(outdir, "AdditionalFile_GSEA_full.tsv"))
}

cat("\nGSEA tidy table:", file.path(outdir, "GSEA_tidy.tsv"), "\n")
cat("Shared vs specific:", file.path(outdir, "GSEA_shared_vs_specific.tsv"), "\n")
