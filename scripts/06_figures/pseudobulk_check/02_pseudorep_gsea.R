# =========================================================
# Section 5c (pseudoreplication copy): GSEA engine
# Environment: go-enrich (local machine)
# =========================================================
# A DEDICATED COPY of 05c_gsea_across_stages.R, reserved for the
# pseudoreplication diagnostics. The analysis is byte-for-byte the same;
# only the input and output directories differ. Keeping a separate copy
# means the two pipelines cannot overwrite each other's paths - the one
# failure mode that produces a plausible-looking wrong answer rather than
# an error.
#
# Run it through 05d_gsea_pseudorep_checks.R, not directly.
#
# Requires in this session:
#   - reference tables from 01b - see 00_reload_go_enrich_session.R
# Requires on disk:
#   - 08_pseudoreplication_checks.R output
#
# >>> DIRECTION CONVENTION - READ THIS <<<
# 04b and 08 keep the original ident order (ident.1 = control), so a
# POSITIVE avg_log2FC means higher in CONTROL. This script NEGATES the
# ranking statistic, so in every output below:
#
#      NES > 0  =  enriched in INFECTED  =  "up_in_target"
#      NES < 0  =  enriched in CONTROL   =  "up_in_control"
#
# The `direction` column states this for every row. Do not re-derive it
# from avg_log2FC without accounting for the negation.

library(readr)
library(dplyr)
library(tidyr)
library(clusterProfiler)

# =========================================================
# CONFIG - pseudoreplication paths, fixed
# =========================================================
de_dir <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/new_figures_8samples/pseudobulk_checks"
outdir <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/new_figures_8samples/gsea_pseudobulk_checks"

REF_DIR <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/new_figures_8samples/reference_tables_3_Oct26"

PADJ_CUTOFF  <- 0.05
MIN_GS_SIZE  <- 5       # matches the main run; state it in Methods
MAX_GS_SIZE  <- 500
MIN_RANKED   <- 200     # skip a comparison with too few ranked genes
SEED         <- 42

## Must stay "logfc" here. The test D pseudobulk tables carry p_val = NA by
## design - there are no nucleus-level p-values to report at animal level -
## so "signed_p" would drop every D comparison silently.
RANK_METRIC <- "logfc"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

t2g_BP            <- readRDS(file.path(REF_DIR, "t2g_BP.rds"))
t2g_MF            <- readRDS(file.path(REF_DIR, "t2g_MF.rds"))
t2g_CC            <- readRDS(file.path(REF_DIR, "t2g_CC.rds"))
term2name_go      <- readRDS(file.path(REF_DIR, "term2name_go.rds"))
kegg_term2gene    <- readRDS(file.path(REF_DIR, "kegg_term2gene.rds"))
kegg_names        <- readRDS(file.path(REF_DIR, "kegg_names.rds"))

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
full_files <- list.files(file.path(de_dir, "full"), pattern = "_full\\.tsv$",
                         full.names = TRUE)
if (length(full_files) == 0)
  stop("No ranked tables found in ", de_dir, "/full - run 08 first.")

cat("Ranked tables found:", length(full_files), "\n")
cat("Ranking metric:", RANK_METRIC, "(negated: NES > 0 = up in infected)\n\n")

results <- list(); logs <- list(); i_r <- 0L; i_l <- 0L

for (f in full_files) {

  res   <- read_tsv(f, show_col_types = FALSE)
  cl    <- unique(res$cluster)[1]
  stage <- unique(res$stage)[1]

  ranked <- make_ranked(res)

  cat(sprintf("%-26s %-24s ranked genes: %5d\n", cl, stage, length(ranked)))

  if (length(ranked) < MIN_RANKED) {
    cat("    too few ranked genes - skipping\n")
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
      error = function(e) {
        cat("    ", onto, "failed:", conditionMessage(e), "\n"); NULL }
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
}

tidy_gsea <- if (length(results)) bind_rows(results) else data.frame()
write_tsv(tidy_gsea, file.path(outdir, "GSEA_tidy.tsv"))
write_tsv(bind_rows(logs), file.path(outdir, "gsea_run_log.tsv"))

## Every comparison 08 wrote should appear in the run log, even with zero
## significant terms. A stage missing here was dropped at MIN_RANKED or
## failed to read, and that is not the same as "no enrichment".
cat("\n--- comparisons attempted ---\n")
print(as.data.frame(
  bind_rows(logs) %>% count(stage, name = "ontologies_run")))

cat("\nGSEA tidy table:", file.path(outdir, "GSEA_tidy.tsv"), "\n")

