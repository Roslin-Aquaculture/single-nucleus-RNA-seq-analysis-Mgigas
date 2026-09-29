###############################################################################
#  SUPPLEMENTARY FIGURE AND TABLE: SUB-CLUSTERING OF CLUSTER 1
#
#  Dewari et al., single-nucleus atlas of the Pacific oyster during OsHV-1
#  infection.
#
#  Reads the two Seurat objects written by reclustering_cluster1.R and produces
#  the files intended for the manuscript. The clustering is not repeated, so
#  figure styling can be iterated cheaply.
#
#  ASCII only throughout.
#
#  COLOUR
#  Categorical hues are assigned in one fixed order, shared between the
#  animal-of-origin panels and the composition bars, and validated for
#  colour-vision deficiency (worst adjacent pair: CVD dE 9.1, normal vision
#  dE 19.6, OKLab x100). A UMAP is an all-pairs case, in which no seven-colour
#  palette separates every pair, so the subcluster panels carry direct labels at
#  group centroids: identity never depends on colour alone.
###############################################################################

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(svglite); library(systemfonts); library(xml2)

IN_DIR   <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/cluster1_subclustering_final/full"
DIR_SUPP <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/cluster1_subclustering_final/supplementary"
dir.create(DIR_SUPP, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------

SAMPLE_LEVELS <- c("Uninfected", "Homogenate", "6-hpiA", "6-hpiD",
                   "24-hpiJ", "72-hpiJ", "96-hpiE")
COND_LEVELS   <- c("control", "6hpi", "24hpi", "72hpi", "96hpi")
RUN_LEVELS    <- c("Uncorrected", "Sample-corrected (Harmony)")

# The clustering script writes `condition`. Older objects carry `condition_new`;
# both are accepted so a previously saved object still loads.
COND_COL_CANDIDATES <- c("condition", "condition_new")

## Six panels: both runs coloured by subcluster, animal and infection stage.
## The uncorrected-by-stage panel is what shows that animals treated alike
## separate as completely as animals from different timepoints - the sentence
## that answers the reviewer. Set FALSE for the four-panel version.
SIX_PANEL <- TRUE

FONT       <- "Arial"
BASE_SIZE  <- 9
FIG_W_MM   <- 190
FIG_H_MM   <- if (SIX_PANEL) 230 else 200

## Fixed hue order, validated as a set. Never cycled, never reordered per panel.
PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
         "#e87ba4", "#008300", "#4a3aa7", "#e34948")

SAMPLE_COLS <- setNames(PAL[seq_along(SAMPLE_LEVELS)], SAMPLE_LEVELS)

## Infection stage is ordered and must not share hues with the animal key,
## or blue means Uninfected in one panel and control in the next.
COND_COLS <- setNames(
  c("#4d4d4d", "#8c6bb1", "#41ab5d", "#f16913", "#ce1256"), COND_LEVELS)

## Subclusters are labelled directly on the embedding, so colour is secondary
## here. The two runs are different partitions of the same nuclei, so they use
## different ramps: U_0 and H_0 are unrelated groups and should not look alike.
RAMP_UNCORR <- "Batlow"
RAMP_CORR   <- "Roma"

# ---------------------------------------------------------------------------
# EXPORT (same standard as the other manuscript figures)
# ---------------------------------------------------------------------------
fm <- systemfonts::match_fonts(FONT)
cat("font requested:", FONT, "-> resolved:", basename(fm$path), "\n")

# svglite puts shared stroke and fill in a CSS block inside <defs>; Inkscape
# drops it on ungroup and axis lines vanish. It also writes textLength on every
# <text>, so deleting a character stretches the rest. Both are undone here.
fix_svg_for_inkscape <- function(file) {
  x  <- xml2::read_xml(file); ns <- xml2::xml_ns(x)
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
    if (!grepl("fill\\s*:", s)) { s <- paste0(trimws(s), " fill: none;"); changed <- TRUE }
    if (changed) { xml2::xml_set_attr(e, "style", s); n <- n + 1 }
  }
  tx <- xml2::xml_find_all(x, "//d1:text|//d1:tspan", ns)
  for (e in tx) {
    xml2::xml_set_attr(e, "textLength",   NULL)
    xml2::xml_set_attr(e, "lengthAdjust", NULL)
  }
  xml2::write_xml(x, file)
  cat("     inlined styles on", n, "of", length(els), "elements;",
      "released", length(tx), "text widths\n")
}

# cairo_pdf keeps text as text and carries physical units unambiguously. The
# default pdf() device writes per-glyph positions and splits letters apart in
# Inkscape. The PNG is a quick-look copy, not the submission file.
save_all <- function(plot, stem, width_mm, height_mm) {
  ggsave(paste0(stem, ".svg"), plot, width = width_mm, height = height_mm,
         units = "mm", device = svglite::svglite, limitsize = FALSE)
  fix_svg_for_inkscape(paste0(stem, ".svg"))
  ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, limitsize = FALSE)
  ggsave(paste0(stem, ".png"), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 400, limitsize = FALSE)
}

theme_umap <- theme_minimal(base_size = BASE_SIZE, base_family = FONT) +
  theme(text            = element_text(size = BASE_SIZE, colour = "black",
                                       family = FONT),
        panel.grid      = element_blank(),
        axis.line       = element_line(colour = "black", linewidth = 0.3),
        axis.ticks      = element_line(colour = "black", linewidth = 0.3),
        axis.ticks.length = unit(1.5, "pt"),
        axis.text       = element_text(size = BASE_SIZE, colour = "black"),
        axis.title      = element_text(size = BASE_SIZE, colour = "black"),
        plot.title      = element_text(size = BASE_SIZE, colour = "black",
                                       family = FONT, face = "bold"),
        legend.text     = element_text(size = BASE_SIZE, colour = "black"),
        legend.title    = element_blank(),
        legend.key.size = unit(3, "mm"))

###############################################################################
# Load and relabel
###############################################################################

uncorr <- readRDS(file.path(IN_DIR, "cluster1_uncorrected.rds"))
corr   <- readRDS(file.path(IN_DIR, "cluster1_harmony.rds"))

cond_col <- intersect(COND_COL_CANDIDATES, colnames(uncorr@meta.data))[1]
stopifnot(!is.na(cond_col))
cat("infection-stage column:", cond_col, "\n")

## Underscore-separated labels, ordered numerically rather than alphabetically
## (which would place U_10 before U_2 at a higher resolution).
relabel <- function(obj, prefix) {
  k <- sub("^[UH]_?", "", as.character(obj$subcluster))
  obj$subcluster <- factor(
    paste0(prefix, "_", k),
    levels = paste0(prefix, "_", sort(as.numeric(unique(k)))))
  obj$sample <- factor(as.character(obj$sample), levels = SAMPLE_LEVELS)
  obj$cond   <- factor(as.character(obj@meta.data[[cond_col]]), levels = COND_LEVELS)
  obj
}

uncorr <- relabel(uncorr, "U")
corr   <- relabel(corr,   "H")

stopifnot(!anyNA(uncorr$sample), !anyNA(uncorr$cond))

###############################################################################
# Composition table
###############################################################################

## Recomputed here rather than read from the clustering script's CSV, so this
## script stands alone; the assertion below catches any divergence.

composition <- function(obj, label) {
  md <- obj@meta.data
  qc <- md %>% group_by(subcluster) %>%
    summarise(n_nuclei        = n(),
              median_nFeature = round(median(nFeature_RNA)), .groups = "drop")
  pct <- md %>%
    count(subcluster, sample, name = "n") %>%
    complete(subcluster, sample, fill = list(n = 0)) %>%
    group_by(subcluster) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    select(subcluster, sample, pct) %>%
    pivot_wider(names_from = sample, values_from = pct)
  qc %>% left_join(pct, by = "subcluster") %>%
    mutate(run = label, .before = 1) %>%
    mutate(across(where(is.numeric), ~round(.x, 1)))
}

comp <- bind_rows(composition(uncorr, RUN_LEVELS[1]),
                  composition(corr,   RUN_LEVELS[2])) %>%
  mutate(run = factor(run, levels = RUN_LEVELS))

print(as.data.frame(comp))
write_csv(comp, file.path(DIR_SUPP, "TableSX_cluster1_subcluster_composition.csv"))

## Expected contribution of each animal: its share of all cluster 1 nuclei.
## This is what "no subcluster exceeds 29%" has to be read against.
exp_file <- file.path(IN_DIR, "expected_share_per_animal.csv")
if (file.exists(exp_file)) {
  expected_share <- read_csv(exp_file, show_col_types = FALSE)
} else {
  expected_share <- uncorr@meta.data %>%
    count(sample, name = "n_nuclei") %>%
    mutate(pct_of_cluster1 = round(100 * n_nuclei / sum(n_nuclei), 1))
}
exp_lo <- min(expected_share$pct_of_cluster1)
exp_hi <- max(expected_share$pct_of_cluster1)
cat(sprintf("expected single-animal contribution: %.0f-%.0f%%\n", exp_lo, exp_hi))

###############################################################################
# Panel builders
###############################################################################

emb_df <- function(obj) {
  as.data.frame(Embeddings(obj, "umap")[, 1:2]) %>%
    setNames(c("UMAP_1", "UMAP_2")) %>%
    bind_cols(obj@meta.data[, c("subcluster", "sample", "cond")])
}

centroids <- function(df, col)
  df %>% group_by(grp = .data[[col]]) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

## label_groups = FALSE where groups are fully intermixed: centroids then
## coincide and a pile of overlapping labels reads as a rendering fault.
umap_panel <- function(df, col, cols, title,
                       label_groups = TRUE, show_legend = TRUE) {
  p <- ggplot(df, aes(UMAP_1, UMAP_2, colour = .data[[col]])) +
    geom_point(size = 0.35, alpha = 0.85) +
    scale_colour_manual(values = cols, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") + theme_umap

  if (label_groups) {
    ctr <- centroids(df, col)
    p <- p + geom_text_repel(
      data = ctr, aes(x = UMAP_1, y = UMAP_2, label = grp),
      inherit.aes = FALSE,                 # x and y must be given explicitly
      colour = "grey10", size = 2.6, fontface = "bold", family = FONT,
      max.overlaps = Inf,                  # never silently drop a label
      bg.color = "white", bg.r = 0.15,     # halo, so labels read over points
      box.padding = 0.2, min.segment.length = 0.3,
      segment.colour = "grey50", segment.size = 0.2)
  }
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

df_u <- emb_df(uncorr)
df_h <- emb_df(corr)

pA <- umap_panel(df_u, "subcluster",
                 hcl.colors(nlevels(df_u$subcluster), RAMP_UNCORR),
                 "A  Subclusters", show_legend = FALSE)
pB <- umap_panel(df_u, "sample", SAMPLE_COLS, "B  Animal of origin")
pC <- umap_panel(df_u, "cond",   COND_COLS,   "C  Infection stage")
pD <- umap_panel(df_h, "subcluster",
                 hcl.colors(nlevels(df_h$subcluster), RAMP_CORR),
                 "D  Subclusters", show_legend = FALSE)
pE <- umap_panel(df_h, "sample", SAMPLE_COLS, "E  Animal of origin",
                 label_groups = FALSE)
pF <- umap_panel(df_h, "cond",   COND_COLS,   "F  Infection stage",
                 label_groups = FALSE)

###############################################################################
# Composition panel
###############################################################################

## reverse = TRUE stacks in legend order, so the bars read top-to-bottom the way
## the animal legend reads down. The white separator keeps thin slivers visible
## in the uncorrected facet, where one animal occupies almost the whole bar.
## The dashed lines are the expected single-animal contribution.

bar_df <- comp %>%
  pivot_longer(-c(run, subcluster, n_nuclei, median_nFeature),
               names_to = "sample", values_to = "pct") %>%
  mutate(sample = factor(sample, levels = SAMPLE_LEVELS))

bar_title <- sprintf("%s  Contribution of each animal to each subcluster (dashed: expected %.0f-%.0f%%)",
                     if (SIX_PANEL) "G" else "E", exp_lo, exp_hi)

pBar <- ggplot(bar_df, aes(subcluster, pct, fill = sample)) +
  geom_col(position = position_stack(reverse = TRUE),
           colour = "white", linewidth = 0.25) +
  geom_hline(yintercept = c(exp_lo, exp_hi), linetype = "dashed",
             linewidth = 0.25, colour = "grey30") +
  scale_fill_manual(values = SAMPLE_COLS, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~ run, scales = "free_x") +
  labs(x = NULL, y = "% of subcluster", title = bar_title) +
  theme_minimal(base_size = BASE_SIZE, base_family = FONT) +
  theme(text               = element_text(size = BASE_SIZE, colour = "black",
                                          family = FONT),
        axis.text.x        = element_text(angle = 45, hjust = 1,
                                          size = BASE_SIZE, colour = "black"),
        axis.text.y        = element_text(size = BASE_SIZE, colour = "black"),
        axis.title         = element_text(size = BASE_SIZE, colour = "black"),
        plot.title         = element_text(size = BASE_SIZE, colour = "black",
                                          family = FONT, face = "bold"),
        strip.text         = element_text(size = BASE_SIZE, colour = "black",
                                          face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.25),
        panel.grid.minor   = element_blank(),
        legend.position    = "none")      # shares the scale shown in panel B

###############################################################################
# Assemble and save
###############################################################################

fig <- (pA | pB | pC) / (pD | pE | pF) / pBar +
  plot_layout(heights = c(1, 1, 0.85), guides = "collect")

save_all(fig, file.path(DIR_SUPP, "FigureSX_cluster1_subclustering"),
         FIG_W_MM, FIG_H_MM)

cat("\n  Written to ", normalizePath(DIR_SUPP), ":\n", sep = "")
cat("    FigureSX_cluster1_subclustering.svg / .pdf / .png\n")
cat("    TableSX_cluster1_subcluster_composition.csv\n\n")
cat("  Caption should note that the composition bars use the colours of the\n")
cat("  animal-of-origin panels, and that the dashed lines are the expected\n")
cat("  single-animal contribution.\n")
cat("  Open the PDF and check for label collisions before submitting.\n\n")

writeLines(capture.output(sessionInfo()), file.path(DIR_SUPP, "sessionInfo_figure.txt"))

