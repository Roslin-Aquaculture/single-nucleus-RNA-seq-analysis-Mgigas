# =========================================================
# Section 5b: Over-representation analysis across stages (ORA)
# Environment: go-enrich (local machine)
# =========================================================
# Replaces 05_de_all_stages_batch.R.
#
# Requires in this session:
#   - 01b_build_reference_tables_v2.R  (t2g_BP/MF/CC, term2name_go,
#                                       kegg_term2gene, kegg_names,
#                                       background_global)
# Requires on disk:
#   - 04b_rerun_de_full_ranked.R output (DEG lists + per-cluster universes)
#
# CHANGES FROM v1
#   1. Per-cluster universe instead of the single global background.
#   2. Four ontologies (GO BP / MF / CC, KEGG) instead of one mixed GO
#      table plus a doubled KEGG table.
#   3. Background sensitivity: every test is run twice, once per-cluster
#      and once global, so you can see which terms are real and which were
#      artefacts of the inflated background.
#   4. One tidy master table for the figure and for drafting.
#   5. Explicit p.adjust filtering rather than relying on clusterProfiler's
#      cutoff semantics, which have shifted between versions.
#   6. Mapped-gene reporting, so an underpowered comparison is visible
#      rather than silently empty.
#
# ORA vs GSEA: this script tests the 2-fold DEG lists, matching the gene
# names in your manuscript text. 05c runs GSEA on the full ranked lists
# and will be more sensitive. Run both; use ORA where you name genes and
# GSEA for the pathway figure if it tells the cleaner story.

library(readr)
library(dplyr)
library(tidyr)
library(clusterProfiler)

# =========================================================
# CONFIG
# =========================================================
de_dir <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_full_ranked_minpct01"
outdir <- "/home/pdewari/eggnog/results/enrichment/ora_across_stages"

stages     <- c("6hpi", "24hpi", "72hpi", "96hpi")
directions <- c("up_in_control", "up_in_target")

MIN_DE_GENES <- 15    # floor on ANNOTATED genes, not raw list length
PADJ_CUTOFF  <- 0.05
MIN_GS_SIZE  <- 10
MAX_GS_SIZE  <- 500

RUN_BACKGROUND_SENSITIVITY <- TRUE

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

read_deg <- function(cl, stage, direction) {
  f <- file.path(de_dir, "deg_lists",
                 paste0(safe_name(cl), "_control_vs_", stage,
                        "_de_genes_", direction, ".txt"))
  if (!file.exists(f) || file.info(f)$size == 0) return(character(0))
  g <- readLines(f)
  unique(sub("^(\\S+)\\s.*$", "\\1", g[g != "" & !is.na(g)]))
}

universe_cache <- new.env(parent = emptyenv())
get_universe <- function(cl) {
  key <- safe_name(cl)
  if (!is.null(universe_cache[[key]])) return(universe_cache[[key]])
  f <- file.path(de_dir, "universes", paste0("universe_", key, ".csv"))
  if (!file.exists(f)) stop("No universe for cluster '", cl, "' — run 04b first.")
  u <- read_csv(f, show_col_types = FALSE)$gene_id
  universe_cache[[key]] <- u
  u
}

run_one <- function(genes, t2g, t2n, universe) {
  uni <- intersect(universe, unique(t2g[[2]]))
  gl  <- intersect(genes, uni)
  if (length(gl) < MIN_DE_GENES) {
    return(list(res = NULL, n_mapped = length(gl), n_universe = length(uni)))
  }
  res <- tryCatch(
    enricher(gene = gl, universe = uni, TERM2GENE = t2g, TERM2NAME = t2n,
             pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
             minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE),
    error = function(e) { cat("      failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(res)) return(list(res = NULL, n_mapped = length(gl), n_universe = length(uni)))
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(list(res = NULL, n_mapped = length(gl), n_universe = length(uni)))
  df <- df %>%
    mutate(fold_enrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio)) %>%
    filter(p.adjust < PADJ_CUTOFF)
  if (nrow(df) == 0) return(list(res = NULL, n_mapped = length(gl), n_universe = length(uni)))
  list(res = df, n_mapped = length(gl), n_universe = length(uni))
}

# =========================================================
# MAIN LOOP
# =========================================================
universe_files <- list.files(file.path(de_dir, "universes"), pattern = "^universe_.*\\.csv$")
clusters <- sub("\\.csv$", "", sub("^universe_", "", universe_files))
cat("Clusters with a universe:", length(clusters), "\n\n")

all_results <- list(); log_entries <- list(); i_res <- 0L; i_log <- 0L
bg_modes <- if (RUN_BACKGROUND_SENSITIVITY) c("per_cluster", "global") else "per_cluster"

for (cl in clusters) {

  cat("==================================================\n")
  cat("CLUSTER:", cl, "\n")
  cat("==================================================\n")

  for (stage in stages) {
    for (direction in directions) {

      de_genes <- read_deg(cl, stage, direction)
      cat(sprintf("  %-6s %-14s  DE genes: %5d\n", stage, direction, length(de_genes)))

      for (bg_mode in bg_modes) {
        universe <- if (bg_mode == "per_cluster") get_universe(cl) else background_global

        for (onto in names(ontologies)) {
          out <- run_one(de_genes, ontologies[[onto]]$t2g, ontologies[[onto]]$t2n, universe)
          n_terms <- if (is.null(out$res)) 0L else nrow(out$res)

          if (bg_mode == "per_cluster") {
            cat(sprintf("        %-6s mapped: %4d / uni: %5d -> %3d terms\n",
                        onto, out$n_mapped, out$n_universe, n_terms))
          }

          if (!is.null(out$res)) {
            i_res <- i_res + 1L
            all_results[[i_res]] <- out$res %>%
              mutate(method = "ORA", cluster = cl, stage = stage,
                     direction = direction, ontology = onto,
                     background = bg_mode, .before = 1)
          }

          i_log <- i_log + 1L
          log_entries[[i_log]] <- data.frame(
            cluster = cl, stage = stage, direction = direction, ontology = onto,
            background = bg_mode, n_de_genes = length(de_genes),
            n_de_mapped = out$n_mapped, n_universe = out$n_universe,
            n_sig_terms = n_terms, below_min_genes = out$n_mapped < MIN_DE_GENES,
            stringsAsFactors = FALSE)
        }
      }
    }
  }

  if (length(log_entries)) write_tsv(bind_rows(log_entries), file.path(outdir, "ora_run_log.tsv"))
  if (length(all_results)) write_tsv(bind_rows(all_results), file.path(outdir, "ORA_tidy.tsv"))
}

run_log <- bind_rows(log_entries)
tidy_all <- if (length(all_results)) bind_rows(all_results) else data.frame()
write_tsv(run_log,  file.path(outdir, "ora_run_log.tsv"))
write_tsv(tidy_all, file.path(outdir, "ORA_tidy.tsv"))

# ---- how many comparisons were actually testable ----
cat("\n--- power check ---\n")
print(as.data.frame(
  run_log %>% filter(background == "per_cluster") %>%
    group_by(ontology) %>%
    summarise(comparisons = n(),
              testable    = sum(!below_min_genes),
              with_terms  = sum(n_sig_terms > 0), .groups = "drop")
))
cat("\nIf `testable` is low, the 2-fold DEG lists are too short for ORA.\n")
cat("That is what 05c (GSEA) is for — it uses every gene and needs no cut.\n")

# ---- background sensitivity: read this before writing ----
if (RUN_BACKGROUND_SENSITIVITY && nrow(tidy_all)) {
  sens <- tidy_all %>%
    dplyr::select(cluster, stage, direction, ontology, ID, Description, background, p.adjust) %>%
    pivot_wider(names_from = background, values_from = p.adjust, names_prefix = "padj_") %>%
    mutate(verdict = case_when(
      !is.na(padj_per_cluster) & !is.na(padj_global) ~ "both_robust",
       is.na(padj_per_cluster) & !is.na(padj_global) ~ "global_only_LIKELY_ARTEFACT",
      !is.na(padj_per_cluster) &  is.na(padj_global) ~ "per_cluster_only_newly_revealed",
      TRUE ~ NA_character_))
  write_tsv(sens, file.path(outdir, "background_sensitivity.tsv"))
  cat("\n--- background sensitivity ---\n")
  print(as.data.frame(count(sens, verdict)))
}

cat("\nORA tidy table:", file.path(outdir, "ORA_tidy.tsv"), "\n")

