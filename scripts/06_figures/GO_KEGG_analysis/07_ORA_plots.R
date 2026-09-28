# =========================================================
# 08_ora_panels.R - ORA dot-plot panels (supplementary)
# Environment: go-enrich. Recomputes nothing; reads ORA_tidy.tsv.
# =========================================================
# ASCII only throughout. A non-ASCII character in a comment has previously
# lost its encoding on paste and been parsed as a unary minus.
#
# Same layout, curation rules and export standard as 07_figure7_final.R,
# applied to over-representation analysis. Four panels are computed: GO BP
# and KEGG, each direction. Only the panels given a `fig` reference below
# appear in the manuscript; the others are computed so the supplementary
# table is complete, but carry no figure reference.
#
# Two differences from the GSEA script, both forced by the data:
#   - redundancy is assessed on `geneID` (ORA's gene column) rather than
#     `core_enrichment`
#   - there is no NES, so colour encodes fold enrichment and direction is
#     carried by the panel title, as in the GSEA panels
#
# Expect sparse panels. Of 142 lists per collection, 34 were testable for
# GO BP and 17 returned terms; the thresholded lists carry a median of 5
# annotated genes. See ORA_coverage_summary.tsv.
# =========================================================
library(readr); library(dplyr); library(tidyr); library(ggplot2)
library(svglite); library(systemfonts); library(xml2)

## ---- paths ----
ENR     <- "/home/pdewari/eggnog/results/enrichment"
ORA_IN  <- file.path(ENR, "ora_across_stages_20260923",  "ORA_tidy.tsv")
GSEA_IN <- file.path(ENR, "gsea_across_stages_20260923", "GSEA_tidy.tsv")
EMAP    <- "/home/pdewari/eggnog/results/full_proteome_20260820_124651/full_proteome.emapper.annotations"
DE_DIR  <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_full_ranked_minpct01_20260923"
FIG_DIR <- file.path(ENR, "figures_20260923_27Sept", "ora_panels_27Sept")

## ---- parameters ----
PADJ          <- 0.05
OVERLAP_CUT   <- 0.8      # overlap coefficient: intersection / smaller set
MIN_CLUSTERS  <- 1        # 1 shows everything; raise once the panels are seen
STAGES        <- c("6hpi", "24hpi", "72hpi", "96hpi")
BACKGROUND    <- "per_cluster"
MARK_UNTESTED <- TRUE

# Display labels only, keyed by term ID so a relabel cannot drift onto
# another term. The Additional File keeps the official GO name.
LABEL_OVERRIDE <- c(
  "GO:0035872" = "NOD-like receptor signalling pathway"
)

FONT      <- "Arial"      # Liberation Sans substitutes on Linux
BASE_SIZE <- 9

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# FONT
# =========================================================
fm <- systemfonts::match_fonts(FONT)
cat("font requested:", FONT, "-> resolved:", basename(fm$path), "\n")

# =========================================================
# SVG POST-PROCESSING FOR INKSCAPE
# =========================================================
# Two svglite behaviours have to be undone before the file is editable.
#
# 1. Shared stroke and fill declarations live in a CSS block inside <defs>.
#    Inkscape drops that stylesheet when a group is ungrouped, and axis
#    lines, ticks and panel borders disappear. Copying the properties onto
#    every drawing element makes the file independent of the stylesheet.
#
# 2. Every <text> carries textLength and lengthAdjust="spacingAndGlyphs",
#    so that the rendered width matches R's font metrics whatever font the
#    viewer has. Inkscape honours it, so deleting a character stretches the
#    remaining ones to fill the original width. Removing it hands width
#    back to Inkscape's own Arial metrics, which is what makes a label
#    editable.
fix_svg_for_inkscape <- function(file) {
  
  x  <- xml2::read_xml(file)
  ns <- xml2::xml_ns(x)
  
  # 1. inline the stylesheet
  els <- xml2::xml_find_all(
    x, "//d1:line|//d1:polyline|//d1:polygon|//d1:path|//d1:rect|//d1:circle", ns)
  n <- 0
  for (e in els) {
    s <- xml2::xml_attr(e, "style"); if (is.na(s)) s <- ""
    changed <- FALSE
    if (!grepl("stroke\\s*:", s)) {
      s <- paste0(trimws(s), if (nzchar(trimws(s))) " " else "", "stroke: #000000;")
      changed <- TRUE
    }
    if (!grepl("fill\\s*:", s)) {
      s <- paste0(trimws(s), " fill: none;"); changed <- TRUE
    }
    if (changed) { xml2::xml_set_attr(e, "style", s); n <- n + 1 }
  }
  
  # 2. release the fixed text widths
  tx <- xml2::xml_find_all(x, "//d1:text|//d1:tspan", ns)
  for (e in tx) {
    xml2::xml_set_attr(e, "textLength",   NULL)
    xml2::xml_set_attr(e, "lengthAdjust", NULL)
  }
  
  xml2::write_xml(x, file)
  cat("     inlined styles on", n, "of", length(els), "elements;",
      "released", length(tx), "text widths\n")
}

# cairo_pdf keeps text as text and carries physical units unambiguously.
# The default pdf() device writes per-glyph positions and splits letters
# apart in Inkscape.
save_both <- function(plot, stem, width_mm, height_mm) {
  ggsave(paste0(stem, ".svg"), plot, width = width_mm, height = height_mm,
         units = "mm", device = svglite::svglite, limitsize = FALSE)
  fix_svg_for_inkscape(paste0(stem, ".svg"))
  ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, limitsize = FALSE)
}
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
  
  # what other breadth thresholds would give, so the choice stays visible
  cat(sprintf("        breadth >=1: %3d terms | >=2: %3d | >=3: %3d\n",
              length(kept),
              n_distinct(all_terms$ID[all_terms$Kept_after_collapse & all_terms$N_clusters >= 2]),
              n_distinct(all_terms$ID[all_terms$Kept_after_collapse & all_terms$N_clusters >= 3])))
  
  list(display = filter(all_terms, Displayed), all = all_terms)
}

# `fig` is the manuscript figure reference, or NA for a panel that is
# computed but not shown. Only S1 and S2 are in the supplementary figures,
# so KEGG terms reach the table with no figure reference rather than
# pointing at a figure that does not exist.
#
# w_mm / h_mm are the final printed size, set here rather than derived from
# the number of categories, so each panel imports into the A4 layout at the
# size it will occupy on the page. Tune after the first run: the console
# prints a height that would fit the terms of each panel.
panels <- list(
  S1 = list(onto = "GO_BP", dir = "up_in_control", fig = "Fig. S1",
            title = "ORA - GO biological process, genes lower in infected nuclei",
            w_mm = 180, h_mm = 200),
  S2 = list(onto = "GO_BP", dir = "up_in_target",  fig = "Fig. S2",
            title = "ORA - GO biological process, genes higher in infected nuclei",
            w_mm = 180, h_mm = 200),
  S3 = list(onto = "KEGG",  dir = "up_in_control", fig = NA_character_,
            title = "ORA - KEGG, genes lower in infected nuclei",
            w_mm = 160, h_mm = 110),
  S4 = list(onto = "KEGG",  dir = "up_in_target",  fig = NA_character_,
            title = "ORA - KEGG, genes higher in infected nuclei",
            w_mm = 160, h_mm = 110)
)

# =========================================================
# COMPARISONS THAT WERE NEVER RUN
# =========================================================
# A cluster-stage with no DE output was never tested. Mantle cell type 2
# has 19 nuclei at 96 hpi, too few for differential expression testing.
# Without shading, that cell is indistinguishable from a comparison that
# was tested and returned nothing.
NOT_TESTED <- data.frame(cluster = character(), stage = character(),
                         stringsAsFactors = FALSE)

if (MARK_UNTESTED) {
  if (dir.exists(file.path(DE_DIR, "full"))) {
    tested <- tibble(f = sub("_full\\.tsv$", "",
                             list.files(file.path(DE_DIR, "full"),
                                        pattern = "_full\\.tsv$"))) %>%
      separate(f, into = c("cluster", "stage"), sep = "_control_vs_")
    # grid built from every cluster with DE output, not from clusters that
    # happen to carry a significant term, or an untested comparison in a
    # cluster with no significant results would go unreported
    NOT_TESTED <- expand_grid(cluster = unique(tested$cluster), stage = STAGES) %>%
      anti_join(tested, by = c("cluster", "stage")) %>%
      mutate(cluster = gsub("_", " ", cluster)) %>%
      as.data.frame()
    cat("\nuntested cluster-stage comparisons:", nrow(NOT_TESTED), "\n")
    if (nrow(NOT_TESTED)) print(NOT_TESTED)
  } else {
    warning("DE_DIR/full not found; untested comparisons will not be shaded")
  }
}


# =========================================================
# PLOT
# =========================================================
make_plot <- function(d, title, dir_keep) {
  
  # display relabelling, before the term order is computed from Description
  if (length(LABEL_OVERRIDE)) {
    d <- d %>% mutate(
      Description = ifelse(ID %in% names(LABEL_OVERRIDE),
                           LABEL_OVERRIDE[ID], Description))
  }
  
  ord <- d %>% group_by(Description) %>%
    summarise(stage_rank = mean(match(stage, STAGES)),
              nc = first(N_clusters), p = min(p.adjust), .groups = "drop") %>%
    arrange(desc(stage_rank), nc, desc(p)) %>% pull(Description)
  
  d <- d %>% mutate(
    stage       = factor(stage, levels = STAGES),
    cluster     = factor(gsub("_", " ", cluster)),
    Description = factor(Description, levels = ord))
  
  nt <- NOT_TESTED %>%
    filter(cluster %in% levels(d$cluster), stage %in% STAGES) %>%
    tidyr::crossing(Description = levels(d$Description)) %>%
    mutate(stage       = factor(stage,       levels = STAGES),
           cluster     = factor(cluster,     levels = levels(d$cluster)),
           Description = factor(Description, levels = levels(d$Description)))
  
  # fold enrichment is above 1 for every over-represented set, so a
  # sequential scale is correct; direction is carried by the title
  sc <- if (dir_keep == "up_in_control") {
    scale_colour_gradient(low = "#C6DBEF", high = "#08306B", name = "Fold\nenrichment")
  } else {
    scale_colour_gradient(low = "#FCBBA1", high = "#A50F15", name = "Fold\nenrichment")
  }
  
  p <- ggplot(d, aes(cluster, Description,
                     colour = fold_enrichment, size = -log10(p.adjust)))
  
  if (nrow(nt)) {
    p <- p + geom_tile(data = nt, aes(x = cluster, y = Description),
                       fill = "grey88", colour = NA, inherit.aes = FALSE)
  }
  
  p +
    geom_point() +
    facet_wrap(~ stage, nrow = 2, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    sc +
    scale_size_continuous(range = c(1.5, 5), name = "-log10 FDR") +
    guides(colour = guide_colourbar(frame.colour = NA, ticks = FALSE)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      # every text element set explicitly; ggplot defaults several to grey30
      text              = element_text(size = BASE_SIZE, colour = "black", family = FONT),
      plot.title        = element_text(size = BASE_SIZE, colour = "black",
                                       family = FONT, face = "plain"),
      axis.text.x       = element_text(angle = 45, hjust = 1,
                                       size = BASE_SIZE, colour = "black"),
      axis.text.y       = element_text(size = BASE_SIZE, colour = "black"),
      strip.text        = element_text(size = BASE_SIZE, colour = "black"),
      legend.text       = element_text(size = BASE_SIZE, colour = "black"),
      legend.title      = element_text(size = BASE_SIZE, colour = "black"),
      legend.key.height = unit(12, "pt"),
      legend.key.width  = unit(7,  "pt"),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      panel.border      = element_rect(colour = "black", linewidth = 0.3, fill = NA),
      panel.grid.major  = element_line(colour = "grey94", linewidth = 0.25),
      panel.grid.minor  = element_blank(),
      axis.ticks        = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.5, "pt"),
      panel.spacing.x   = unit(0.6, "lines"),
      legend.position   = "right")
}

# =========================================================
# SUPPLEMENTARY TABLE
# =========================================================
# Method, Analysis and Direction describe the slice, so the file stands
# alone and cannot be confused with the GSEA additional file, which shares
# several column names. Figure_panel is filled only for displayed terms in
# panels that are actually in the manuscript.
as_table <- function(d, p, cfg) {
  d %>% arrange(desc(Displayed), desc(N_clusters), ID, p.adjust) %>%
    rowwise() %>%
    mutate(Genes = name_of(strsplit(geneID, "/")[[1]])) %>%
    ungroup() %>%
    transmute(
      Method       = "ORA",
      Analysis     = cfg$title,
      Direction    = cfg$dir,
      Figure_panel = if (is.na(cfg$fig)) NA_character_
      else ifelse(Displayed, cfg$fig, NA_character_),
      Cluster      = gsub("_", " ", cluster),
      Stage        = stage,
      Ontology     = ontology, Background = background,
      Term_ID      = ID, Term = Description,
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
  supp_all[[p]] <- as_table(cur$all, p, cfg)
  
  d <- cur$display
  if (!nrow(d)) { cat("  -> panel", p, "empty at this breadth\n"); next }
  
  write_tsv(as_table(d, p, cfg),
            file.path(FIG_DIR, sprintf("ORA_%s_plotted_terms.tsv", p)), na = "")
  
  n_cl <- n_distinct(d$cluster); n_tm <- n_distinct(d$ID)
  
  cat(sprintf("  -> ORA_%s : %d terms x %d clusters (set to %g x %g mm; about %.0f mm high would fit the terms)\n",
              p, n_tm, n_cl, cfg$w_mm, cfg$h_mm, 45 + 4.2 * n_tm))
  
  save_both(make_plot(d, cfg$title, cfg$dir),
            file.path(FIG_DIR, sprintf("ORA_%s", p)),
            width_mm = cfg$w_mm, height_mm = cfg$h_mm)
}

if (length(supp_all)) {
  write_tsv(bind_rows(supp_all),
            file.path(FIG_DIR, "AdditionalFile_ORA_all_terms.tsv"), na = "")
}

# =========================================================
# VERIFY - every plotted row must trace to an unmodified source row
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

# internal consistency of the deposited table
if (length(supp_all)) {
  a <- bind_rows(supp_all)
  stopifnot(
    sum(a$Displayed & !a$Kept_after_collapse) == 0,
    all(a$N_clusters[a$Displayed] >= MIN_CLUSTERS),
    all(a$FDR < PADJ),
    all(a$Fold_enrichment > 1),
    all(lengths(strsplit(a$Gene_ids, "/")) == a$N_genes),
    !any(duplicated(a[, c("Analysis", "Cluster", "Stage", "Term_ID")]))
  )
  cat("supplementary table checks passed:", nrow(a), "rows\n")
  print(as.data.frame(
    a %>% group_by(Analysis, Direction) %>%
      summarise(rows = n(), terms = n_distinct(Term_ID),
                kept = n_distinct(Term_ID[Kept_after_collapse]),
                shown = n_distinct(Term_ID[Displayed]), .groups = "drop")))
}

# =========================================================
# DOES ORA ADD GENES, OR THE SAME GENES UNDER MORE LABELS?
# =========================================================
# This is the question that decides whether these panels belong in the
# paper: a term ORA finds and GSEA does not is only worth reporting if the
# genes behind it are genes GSEA never surfaced.
if (file.exists(GSEA_IN)) {
  
  gsea <- read_tsv(GSEA_IN, show_col_types = FALSE)
  cat("\nGSEA rows:", nrow(gsea), "\n")
  
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
