# =========================================================
# 08_ora_panels.R — ORA dot-plot panels (supplementary)
# Environment: go-enrich. Recomputes nothing; reads ORA_tidy.tsv.
# =========================================================
# Same layout and same curation rules as 07_figure7_final.R, applied to
# over-representation analysis. Four panels: GO BP and KEGG, each direction.
#
# Two differences from the GSEA script, both forced by the data:
#   - redundancy is assessed on `geneID` (ORA's gene column) rather than
#     `core_enrichment`
#   - there is no NES, so colour encodes fold enrichment and direction is
#     carried by the panel title, as in the GSEA panels
#
# Expect sparse panels. Of 144 comparisons per collection, 34 were testable
# for GO BP and 17 returned terms; the thresholded lists carry a median of
# 5 annotated genes. See ORA_coverage_summary.tsv.

library(readr); library(dplyr); library(tidyr); library(ggplot2)

## ---- paths ----
ORA_IN  <- "/home/pdewari/eggnog/results/enrichment/ora_across_stages_20260923/ORA_tidy.tsv"
GSEA_IN <- "/home/pdewari/eggnog/results/enrichment_20260923/gsea_across_stages_20260923/GSEA_tidy.tsv"
EMAP    <- "/home/pdewari/eggnog/results/full_proteome_20260820_124651/full_proteome.emapper.annotations"
FIG_DIR <- "/home/pdewari/eggnog/results/enrichment/figures_20260923/ora_panels"

## ---- parameters ----
PADJ         <- 0.05
OVERLAP_CUT  <- 0.8      # overlap coefficient: intersection / smaller set
MIN_CLUSTERS <- 1        # set to 1 to see everything, then decide
STAGES       <- c("6hpi", "24hpi", "72hpi", "96hpi")
BACKGROUND   <- "per_cluster"

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# LOAD
# =========================================================
ora <- read_tsv(ORA_IN, show_col_types = FALSE)

cat("\n--- background values present ---\n")
print(as.data.frame(count(ora, background, direction)))

# ORA_tidy.tsv holds both the cluster-specific and the global runs. Plotting
# them together would mix two universes, so filter explicitly.
stopifnot(BACKGROUND %in% unique(ora$background))
ora <- ora %>% filter(background == BACKGROUND)
stopifnot(all(unique(ora$stage) %in% STAGES))

cat("\nrows after filtering to", BACKGROUND, ":", nrow(ora), "\n")

# =========================================================
# PROVENANCE
# =========================================================
prov <- data.frame(
  input       = ORA_IN,
  input_mtime = format(file.mtime(ORA_IN)),
  input_md5   = unname(tools::md5sum(ORA_IN)),
  annotation  = EMAP,
  built       = format(Sys.time()),
  background  = BACKGROUND,
  padj = PADJ, overlap_cut = OVERLAP_CUT, min_clusters = MIN_CLUSTERS)
write_tsv(prov, file.path(FIG_DIR, "ora_panels_provenance.tsv"))
print(t(prov))

# =========================================================
# GENE NAMES
# =========================================================
hdr  <- readLines(EMAP, n = 200)
skip <- grep("^#query", hdr)[1] - 1
egg  <- read_tsv(EMAP, skip = skip, comment = "##", show_col_types = FALSE) %>%
  rename(query = 1) %>%
  mutate(gene_id = sub("\\.\\d+$", "", sub("^transcript:", "", query)))

name_of <- function(ids) {
  i  <- match(ids, egg$gene_id)
  nm <- egg$Preferred_name[i]; pf <- egg$PFAMs[i]
  bad <- function(x) is.na(x) | x == "-" | x == ""
  nm[bad(nm)] <- pf[bad(nm)]
  nm[bad(nm)] <- ids[bad(nm)]
  paste(unique(nm), collapse = ", ")
}

# =========================================================
# CURATE
# =========================================================
collapse_redundant <- function(df, cut) {
  if (nrow(df) == 0) return(df)
  df   <- df[order(df$p.adjust), ]
  sets <- strsplit(df$geneID, "/")
  keep <- logical(nrow(df))
  for (i in seq_len(nrow(df))) {
    if (!any(keep)) { keep[i] <- TRUE; next }
    ov <- vapply(which(keep), function(k)
      length(intersect(sets[[i]], sets[[k]])) / min(length(sets[[i]]), length(sets[[k]])),
      numeric(1))
    if (max(ov) < cut) keep[i] <- TRUE
  }
  df[keep, ]
}

curate <- function(onto, dir_keep) {

  sig <- ora %>% filter(ontology == onto, direction == dir_keep, p.adjust < PADJ)
  if (!nrow(sig)) {
    cat(sprintf("\n%-6s %-14s  no terms\n", onto, dir_keep))
    return(list(display = sig, all = sig))
  }

  rep_terms <- sig %>%
    group_by(ID) %>%
    summarise(p.adjust = min(p.adjust),
              geneID = paste(unique(unlist(strsplit(geneID, "/"))), collapse = "/"),
              .groups = "drop")
  kept <- collapse_redundant(as.data.frame(rep_terms), OVERLAP_CUT)$ID

  all_terms <- sig %>%
    group_by(ID) %>% mutate(N_clusters = n_distinct(cluster)) %>% ungroup() %>%
    mutate(Kept_after_collapse = ID %in% kept,
           Displayed           = Kept_after_collapse & N_clusters >= MIN_CLUSTERS)

  cat(sprintf("\n%-6s %-14s  %4d rows | %3d terms | %3d collapsed | %3d shown (>=%d clusters)\n",
              onto, dir_keep, nrow(sig), n_distinct(sig$ID),
              length(kept), n_distinct(all_terms$ID[all_terms$Displayed]),
              MIN_CLUSTERS))

  # what a breadth of 1 would give, so the threshold choice is visible
  cat(sprintf("        breadth >=1: %3d terms | >=2: %3d | >=3: %3d\n",
              length(kept),
              n_distinct(all_terms$ID[all_terms$Kept_after_collapse & all_terms$N_clusters >= 2]),
              n_distinct(all_terms$ID[all_terms$Kept_after_collapse & all_terms$N_clusters >= 3])))

  list(display = filter(all_terms, Displayed), all = all_terms)
}

panels <- list(
  S1 = list(onto = "GO_BP", dir = "up_in_control",
            title = "ORA - GO biological process, genes lower in infected nuclei"),
  S2 = list(onto = "GO_BP", dir = "up_in_target",
            title = "ORA - GO biological process, genes higher in infected nuclei"),
  S3 = list(onto = "KEGG",  dir = "up_in_control",
            title = "ORA - KEGG, genes lower in infected nuclei"),
  S4 = list(onto = "KEGG",  dir = "up_in_target",
            title = "ORA - KEGG, genes higher in infected nuclei")
)

# =========================================================
# PLOT
# =========================================================
make_plot <- function(d, title, dir_keep) {

  ord <- d %>% group_by(Description) %>%
    summarise(stage_rank = mean(match(stage, STAGES)),
              nc = first(N_clusters), p = min(p.adjust), .groups = "drop") %>%
    arrange(desc(stage_rank), nc, desc(p)) %>% pull(Description)

  d <- d %>% mutate(
    stage       = factor(stage, levels = STAGES),
    cluster     = factor(gsub("_", " ", cluster)),
    Description = factor(Description, levels = ord))

  sc <- if (dir_keep == "up_in_control") {
    scale_colour_gradient(low = "#C6DBEF", high = "#08306B", name = "Fold\nenrichment")
  } else {
    scale_colour_gradient(low = "#FCBBA1", high = "#A50F15", name = "Fold\nenrichment")
  }

  ggplot(d, aes(cluster, Description,
                colour = fold_enrichment, size = -log10(p.adjust))) +
    geom_point() +
    facet_wrap(~ stage, nrow = 2, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    sc +
    scale_size_continuous(range = c(1.5, 5), name = "-log10 FDR") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey95", colour = NA),
          panel.spacing.x  = unit(0.6, "lines"),
          panel.grid.minor = element_blank(),
          legend.position  = "right")
}

as_table <- function(d, p) {
  d %>% arrange(desc(N_clusters), p.adjust, cluster, stage) %>%
    rowwise() %>%
    mutate(Genes = name_of(strsplit(geneID, "/")[[1]])) %>%
    ungroup() %>%
    transmute(Panel = p, Cluster = gsub("_", " ", cluster), Stage = stage,
              Ontology = ontology, Background = background,
              Term_ID = ID, Term = Description,
              N_clusters, Kept_after_collapse, Displayed,
              N_genes = Count, Fold_enrichment = round(fold_enrichment, 2),
              FDR = signif(p.adjust, 3),
              Genes_annotated = Genes, Gene_ids = geneID)
}

supp_all <- list()

for (p in names(panels)) {

  cfg <- panels[[p]]
  cur <- curate(cfg$onto, cfg$dir)
  if (!nrow(cur$all)) next
  supp_all[[p]] <- as_table(cur$all, p)

  d <- cur$display
  if (!nrow(d)) { cat("  -> panel", p, "empty at this breadth\n"); next }

  write_tsv(as_table(d, p), file.path(FIG_DIR, sprintf("ORA_%s_plotted_terms.tsv", p)))

  n_cl <- n_distinct(d$cluster); n_tm <- n_distinct(d$ID)

  ggsave(file.path(FIG_DIR, sprintf("ORA_%s.pdf", p)),
         make_plot(d, cfg$title, cfg$dir),
         width  = 5 + 2 * (0.30 * max(n_cl, 3)),
         height = max(5.5, 2.2 + 0.20 * n_tm),
         limitsize = FALSE)

  cat(sprintf("  -> ORA_%s.pdf : %d terms x %d clusters\n", p, n_tm, n_cl))
}

if (length(supp_all)) {
  write_tsv(bind_rows(supp_all),
            file.path(FIG_DIR, "AdditionalFile_ORA_all_terms.tsv"))
}

# =========================================================
# VERIFY — every plotted row must trace to an unmodified source row
# =========================================================
chk <- lapply(names(panels), function(p) {
  f <- file.path(FIG_DIR, sprintf("ORA_%s_plotted_terms.tsv", p))
  if (!file.exists(f)) return(NULL)
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(cl = gsub(" ", "_", Cluster)) %>%
    left_join(ora %>% select(cluster, stage, ID,
                             FE_s = fold_enrichment, p_s = p.adjust),
              by = c("cl" = "cluster", "Stage" = "stage", "Term_ID" = "ID")) %>%
    summarise(panel = p, rows = n(),
              FE_mismatch  = sum(abs(Fold_enrichment - round(FE_s, 2)) > 1e-6, na.rm = TRUE),
              FDR_mismatch = sum(abs(FDR - signif(p_s, 3)) > 1e-12, na.rm = TRUE),
              unmatched    = sum(is.na(FE_s)))
}) %>% bind_rows()

cat("\n--- verification (last three columns must be 0) ---\n")
print(as.data.frame(chk))

# =========================================================
# DOES ORA ADD GENES, OR THE SAME GENES UNDER MORE LABELS?
# =========================================================
# This is the question that decides whether these panels belong in the paper.

if (file.exists(GSEA_IN)) {

  gsea <- read_tsv(GSEA_IN, show_col_types = FALSE)

  for (dir_keep in c("up_in_control", "up_in_target")) {

    g <- gsea %>% filter(ontology == "GO_BP", direction == dir_keep)
    o <- ora  %>% filter(ontology == "GO_BP", direction == dir_keep, p.adjust < PADJ)

    if (!nrow(o)) { cat("\n--- GO_BP", dir_keep, ": no ORA terms ---\n"); next }

    g_ids   <- unique(g$ID)
    g_genes <- unique(unlist(strsplit(g$core_enrichment, "/")))

    extra       <- o %>% filter(!ID %in% g_ids)
    extra_genes <- unique(unlist(strsplit(extra$geneID, "/")))
    novel       <- setdiff(extra_genes, g_genes)

    cat(sprintf("\n--- GO_BP %s ---\n", dir_keep))
    cat("GSEA terms:", length(g_ids), "| ORA terms:", n_distinct(o$ID),
        "| ORA-only terms:", n_distinct(extra$ID), "\n")
    cat("genes behind ORA-only terms:", length(extra_genes),
        "| not in any GSEA leading edge:", length(novel),
        sprintf("(%.0f%%)\n", 100 * length(novel) / max(1, length(extra_genes))))

    if (length(novel)) {
      all_ids <- unlist(strsplit(extra$geneID, "/"))
      tb <- sort(table(all_ids[all_ids %in% novel]), decreasing = TRUE)
      cat("most frequent of those genes:\n")
      print(head(data.frame(gene    = names(tb),
                            name    = vapply(names(tb), name_of, character(1)),
                            n_terms = as.integer(tb)), 15), row.names = FALSE)
    }
  }
}

writeLines(capture.output(sessionInfo()), file.path(FIG_DIR, "sessionInfo.txt"))
cat("\nOutputs:", FIG_DIR, "\n")
