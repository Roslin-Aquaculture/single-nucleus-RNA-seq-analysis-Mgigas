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

library(readr); library(dplyr); library(tidyr); library(clusterProfiler)

REF_DIR <- "/home/pdewari/eggnog/results/full_proteome_20260820_124651/reference_tables_v2"
t2g_BP         <- readRDS(file.path(REF_DIR, "t2g_BP.rds"))
t2g_MF         <- readRDS(file.path(REF_DIR, "t2g_MF.rds"))
t2g_CC         <- readRDS(file.path(REF_DIR, "t2g_CC.rds"))
term2name_go   <- readRDS(file.path(REF_DIR, "term2name_go.rds"))
kegg_term2gene <- readRDS(file.path(REF_DIR, "kegg_term2gene.rds"))
kegg_names     <- readRDS(file.path(REF_DIR, "kegg_names.rds"))

de_dir <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_pseudorep_checks"
outdir <- "/home/pdewari/eggnog/results/enrichment/gsea_pseudorep_checks"

cat("ranked tables found:",
    length(list.files(file.path(de_dir, "full"), pattern = "_full\\.tsv$")), "\n")



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


##########################
tg <- read_tsv("/home/pdewari/eggnog/results/enrichment/gsea_pseudorep_checks/GSEA_tidy.tsv",
               show_col_types = FALSE)

TRANS <- "ribosom|translat|peptide biosynth|peptide metab|amide biosynth|tRNA|aminoacyl"

cat("\n=== 1. terms per comparison ===\n")
print(as.data.frame(count(tg, stage, direction)))

cat("\n=== 2. TRANSLATIONAL terms per comparison ===\n")
print(as.data.frame(
  tg %>% filter(grepl(TRANS, Description, ignore.case = TRUE)) %>%
    count(stage, direction, name = "n_translational")))

cat("\n=== 3. A_ctrl_vs_ctrl — translational terms (the decisive test) ===\n")
a <- tg %>% filter(stage == "A_ctrl_vs_ctrl", grepl(TRANS, Description, ignore.case = TRUE))
if (nrow(a) == 0) cat("NONE — signature absent between two control animals\n") else
  print(as.data.frame(a[, c("cluster","direction","Description","NES","p.adjust")]))

cat("\n=== 4. D_pseudobulk_6hpi — translational terms (animal-level) ===\n")
d <- tg %>% filter(stage == "D_pseudobulk_6hpi", grepl(TRANS, Description, ignore.case = TRUE))
if (nrow(d) == 0) cat("NONE\n") else
  print(as.data.frame(d[, c("cluster","direction","Description","NES","p.adjust")]))

cat("\n=== 5. C — does it appear in each 6 hpi animal separately? ===\n")
print(as.data.frame(
  tg %>% filter(stage %in% c("C_ctrl_vs_6hpiA","C_ctrl_vs_6hpiD"),
                grepl(TRANS, Description, ignore.case = TRUE)) %>%
    count(stage, cluster, direction)))

######

core <- read_csv("/home/pdewari/Documents/parse_2025/seurat_2025/cluster1_identity/cluster1_core_markers_all_animals.csv",
                 show_col_types = FALSE)$gene_id
length(core)
intersect(c("G19384","G19385","G19421","G19423","G19424"), core)

