# =========================================================
# Section 6: Dotplot matrix figures — term x (cluster x stage)
# Environment: go-enrich (local machine)
# =========================================================
# Input:  ORA_tidy.tsv  (from 05b)  and/or  GSEA_tidy.tsv (from 05c)
# Output: one PDF+PNG per method per ontology per direction, plus
#         candidate term lists for hand-curation and the supplementary
#         enrichment table.
#
# WHY A MATRIX AND NOT 24 SEPARATE DOTPLOTS
# 3 clusters x 4 stages x 2 directions = 24 panels. Reviewers who asked
# you to CUT text will not accept 24 new figures. One matrix — terms down
# the y-axis, stage across the x-axis, faceted by cluster — shows the
# shared core as a solid block and cluster-specific terms as sparse rows.
# That single figure is the argument for merging three results sections
# into two: the reader sees at a glance that most of the response is
# common across cell types.
#
# TWO-PASS WORKFLOW
#   Pass 1: no curated file exists -> top terms picked automatically and
#           term_selection_<method>_<ontology>.tsv written.
#   Pass 2: copy that to curated_terms_<method>_<ontology>.tsv, cut to the
#           ~15-20 terms that carry the story, re-run.
# Top-N by p-value alone returns five flavours of the same GO branch —
# a direct consequence of the ancestor propagation in 01b. Curate.

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# =========================================================
# CONFIG
# =========================================================
ora_file  <- "/home/pdewari/eggnog/results/enrichment/ora_across_stages/ORA_tidy.tsv"
gsea_file <- "/home/pdewari/eggnog/results/enrichment/gsea_across_stages/GSEA_tidy.tsv"
outdir    <- "/home/pdewari/eggnog/results/enrichment/figures"

# Clusters for the manuscript figure. NULL = all.
FOCUS_CLUSTERS <- c("Cluster_1", "Hepatopancreas_cells", "Gill_ciliary_cells")

STAGE_ORDER <- c("6hpi", "24hpi", "72hpi", "96hpi")

ONTOLOGIES <- c("GO_BP", "GO_CC", "KEGG")    # add "GO_MF" if wanted
DIRECTIONS <- c("up_in_target", "up_in_control")

DIRECTION_LABEL <- c(up_in_target  = "Upregulated in infected",
                     up_in_control = "Downregulated in infected")

TOP_N_PER_CELL <- 5     # pass-1 automatic selection

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

curated_file <- function(method, onto)
  file.path(outdir, paste0("curated_terms_", method, "_", onto, ".tsv"))

# =========================================================
# LOAD — whichever methods exist
# =========================================================
datasets <- list()

if (file.exists(ora_file)) {
  d <- read_tsv(ora_file, show_col_types = FALSE)
  if (nrow(d)) {
    d <- d %>% filter(background == "per_cluster") %>%
      mutate(size_stat = fold_enrichment, size_lab = "Fold\nenrichment")
    datasets[["ORA"]] <- d
  }
}

if (file.exists(gsea_file)) {
  d <- read_tsv(gsea_file, show_col_types = FALSE)
  if (nrow(d)) {
    d <- d %>% mutate(size_stat = abs(NES), size_lab = "|NES|")
    datasets[["GSEA"]] <- d
  }
}

if (length(datasets) == 0) stop("Neither ORA_tidy.tsv nor GSEA_tidy.tsv found.")
cat("Methods loaded:", paste(names(datasets), collapse = ", "), "\n")

prep <- function(d) {
  if (!is.null(FOCUS_CLUSTERS)) {
    missing <- setdiff(FOCUS_CLUSTERS, unique(d$cluster))
    if (length(missing))
      warning("FOCUS_CLUSTERS absent from results: ", paste(missing, collapse = ", "),
              "\nAvailable: ", paste(sort(unique(d$cluster)), collapse = ", "))
    d <- d %>% filter(cluster %in% FOCUS_CLUSTERS)
  }
  d %>%
    mutate(stage = factor(stage, levels = STAGE_ORDER),
           cluster = factor(cluster,
                            levels = if (is.null(FOCUS_CLUSTERS)) sort(unique(cluster))
                                     else intersect(FOCUS_CLUSTERS, unique(cluster))))
}

# =========================================================
# PLOT BUILDER
# =========================================================
build_plot <- function(d, method, onto, direction) {

  d <- d %>% filter(ontology == onto, direction == !!direction)
  if (nrow(d) == 0) return(NULL)

  cf <- curated_file(method, onto)
  if (file.exists(cf)) {
    keep_ids <- read_tsv(cf, show_col_types = FALSE)$ID
    d_sel <- d %>% filter(ID %in% keep_ids)
    sel_mode <- "curated"
  } else {
    keep_ids <- d %>%
      group_by(cluster, stage) %>%
      slice_min(p.adjust, n = TOP_N_PER_CELL, with_ties = FALSE) %>%
      ungroup() %>% pull(ID) %>% unique()
    d_sel <- d %>% filter(ID %in% keep_ids)
    sel_mode <- "auto"
  }
  if (nrow(d_sel) == 0) return(NULL)

  d_sel <- d_sel %>%
    mutate(term_label    = str_wrap(str_trunc(Description, 58), width = 42),
           neglog10_padj = -log10(p.adjust))

  # order by breadth first, then strength — makes the shared core read as
  # a block at the top and cluster-specific terms as sparse rows below
  term_order <- d_sel %>%
    group_by(term_label) %>%
    summarise(breadth = n_distinct(paste(cluster, stage)),
              best    = max(neglog10_padj), .groups = "drop") %>%
    arrange(breadth, best) %>% pull(term_label)

  d_sel <- d_sel %>% mutate(term_label = factor(term_label, levels = term_order))

  sub_txt <- if (method == "ORA")
    paste0("ORA, cluster-specific background; BH p.adj < 0.05; ", sel_mode, " selection")
  else
    paste0("GSEA on full ranked lists; BH p.adj < 0.05; ", sel_mode, " selection")

  p <- ggplot(d_sel, aes(x = stage, y = term_label)) +
    geom_point(aes(size = size_stat, colour = neglog10_padj)) +
    facet_grid(. ~ cluster, scales = "free_x", space = "free_x") +
    scale_colour_viridis_c(option = "plasma", end = 0.9,
                           name = expression(-log[10]~italic(p)[adj])) +
    scale_size_continuous(range = c(1.5, 6), name = unique(d_sel$size_lab)[1]) +
    labs(x = NULL, y = NULL,
         title = paste0(DIRECTION_LABEL[[direction]], " — ", sub("_", " ", onto)),
         subtitle = sub_txt) +
    theme_bw(base_size = 10) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1),
          axis.text.y      = element_text(size = 8, lineheight = 0.85),
          panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = "grey70"),
          strip.text       = element_text(face = "bold", size = 9),
          legend.position  = "right",
          plot.title       = element_text(face = "bold"))

  n_terms <- length(unique(d_sel$term_label))
  list(plot   = p,
       height = max(4, 0.32 * n_terms + 2.2),
       width  = 3.0 + 1.7 * length(unique(d_sel$cluster)),
       n_terms = n_terms,
       candidates = d %>% distinct(ID, Description) %>%
         left_join(d %>% group_by(ID) %>%
                     summarise(breadth   = n_distinct(paste(cluster, stage)),
                               best_padj = min(p.adjust),
                               max_size  = max(size_stat), .groups = "drop"),
                   by = "ID") %>%
         arrange(desc(breadth), best_padj))
}

# =========================================================
# RUN
# =========================================================
for (method in names(datasets)) {

  dat <- prep(datasets[[method]])
  cat("\n===", method, "===\n")

  for (onto in ONTOLOGIES) {
    wrote_candidates <- FALSE

    for (direction in DIRECTIONS) {
      built <- build_plot(dat, method, onto, direction)
      if (is.null(built)) { cat("  no results:", onto, direction, "\n"); next }

      f <- file.path(outdir, paste0("dotplot_", method, "_", onto, "_", direction, ".pdf"))
      ggsave(f, built$plot, width = built$width, height = built$height, limitsize = FALSE)
      ggsave(sub("\\.pdf$", ".png", f), built$plot,
             width = built$width, height = built$height, dpi = 300, limitsize = FALSE)
      cat("  wrote", basename(f), "(", built$n_terms, "terms )\n")

      if (!wrote_candidates) {
        write_tsv(built$candidates,
                  file.path(outdir, paste0("term_selection_", method, "_", onto, ".tsv")))
        wrote_candidates <- TRUE
      }
    }
  }

  # ---- supplementary table ----
  supp <- prep(datasets[[method]]) %>% arrange(cluster, direction, ontology, stage, p.adjust)
  write_tsv(supp, file.path(outdir, paste0("AdditionalFile_", method, "_full.tsv")))
}

# =========================================================
# ORA vs GSEA AGREEMENT — which method to lead with
# =========================================================
if (all(c("ORA", "GSEA") %in% names(datasets))) {

  cmp <- bind_rows(
    datasets$ORA  %>% transmute(method = "ORA",  cluster, stage, direction, ontology, ID, Description),
    datasets$GSEA %>% transmute(method = "GSEA", cluster, stage, direction, ontology, ID, Description)
  ) %>%
    distinct() %>%
    group_by(cluster, stage, direction, ontology, ID, Description) %>%
    summarise(methods = paste(sort(unique(method)), collapse = "+"), .groups = "drop")

  write_tsv(cmp, file.path(outdir, "ORA_vs_GSEA_agreement.tsv"))

  cat("\n--- ORA vs GSEA agreement ---\n")
  print(as.data.frame(count(cmp, methods)))
  cat("\nTerms found by BOTH are the safest to build the narrative on.\n")
  cat("GSEA-only terms are usually real signal that ORA missed because of\n")
  cat("the 2-fold cut. ORA-only terms are usually driven by a handful of\n")
  cat("very large fold changes — check the gene counts before claiming them.\n")
}

cat("\nFigures and tables in:", outdir, "\n")
cat("\nPASS 2: copy term_selection_<method>_<ontology>.tsv to\n")
cat("curated_terms_<method>_<ontology>.tsv, cut to ~15-20 terms, re-run.\n")

##########
supp_dir <- file.path(outdir, "per_cluster")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

d_all <- read_tsv(gsea_file, show_col_types = FALSE) %>%
  mutate(stage = factor(stage, levels = STAGE_ORDER))

for (cl in unique(d_all$cluster)) {
  d <- d_all %>% filter(cluster == cl)
  if (nrow(d) == 0) next
  
  d <- d %>%
    mutate(term_label = str_wrap(str_trunc(Description, 60), 45),
           neglog10_padj = -log10(p.adjust),
           panel = DIRECTION_LABEL[direction])
  
  p <- ggplot(d, aes(x = stage, y = reorder(term_label, neglog10_padj))) +
    geom_point(aes(size = abs(NES), colour = neglog10_padj)) +
    facet_wrap(~ panel, scales = "free_y", ncol = 1) +
    scale_colour_viridis_c(option = "plasma", end = 0.9,
                           name = expression(-log[10]~italic(p)[adj])) +
    scale_size_continuous(range = c(1.5, 5), name = "|NES|") +
    labs(title = gsub("_", " ", cl), x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7, lineheight = 0.85),
          strip.text = element_text(face = "bold"))
  
  n <- n_distinct(d$term_label)
  ggsave(file.path(supp_dir, paste0("GSEA_", cl, ".pdf")), p,
         width = 8, height = max(4, 0.22 * n + 2), limitsize = FALSE)
}
