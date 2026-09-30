###############################################################################
#  UMAP SPLIT BY CONDITION, AND CLUSTER COMPOSITION ACROSS THE TIME COURSE
#
#  Dewari et al., single-nucleus atlas of the Pacific oyster during OsHV-1
#  infection. Reviewer request: show the UMAP per condition rather than pooled,
#  and report whether cluster proportions differ across infection time.
#
#  ASCII only throughout.
#
#  ALL EIGHT SAMPLES ARE INCLUDED HERE, matching the combined atlas figure.
#  Sample 24-hpiA is excluded from the differential expression and enrichment
#  analyses [state the reason in Methods] but is retained in the atlas, so the
#  sample set differs between this figure and those analyses. Say so once in
#  Methods rather than leaving a reader to notice it.
#
#  WHAT THIS CAN AND CANNOT SHOW
#  The later stages are represented by a single animal each, so a difference in
#  cluster proportion between stages cannot be separated from a difference
#  between individuals. The per-sample panel is therefore plotted alongside the
#  per-condition one: the spread within control, 6 hpi and 24 hpi, each of which
#  has two animals, is the scale of variation expected without any infection
#  effect, and is the yardstick for reading the rest.
#
#  OUTPUTS
#    FigureSX_umap_by_condition.svg / .pdf   UMAP facets + composition bars
#    cluster_proportions_by_condition.csv
#    cluster_proportions_by_sample.csv
###############################################################################

library(Seurat)
library(tidyverse)
library(patchwork)
library(svglite); library(systemfonts); library(xml2)

OUT <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/umap_by_condition"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

banner <- function(x) cat("\n\n", strrep("=", 74), "\n  ", x, "\n",
                          strrep("=", 74), "\n", sep = "")

# ---------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------

SEU_RDS     <- "/home/pdewari/Documents/parse_2025/seurat_2025/seu_obj_umap_18d_6r_3kRes.rds"
DROP_SAMPLE <- character(0)   # all eight samples, as in the combined atlas figure

COND_LEVELS   <- c("control", "6hpi", "24hpi", "72hpi", "96hpi")
SAMPLE_LEVELS <- c("Uninfected", "Homogenate", "6-hpiA", "6-hpiD",
                   "24-hpiA", "24-hpiJ", "72-hpiJ", "96-hpiE")

## Stages with two animals. Differences within these pairs carry no infection
## component and set the scale against which between-stage differences are read.
PAIRS <- list(
  "uninfected" = c("Uninfected", "Homogenate"),
  "6 hpi"      = c("6-hpiA", "6-hpiD"),
  "24 hpi"     = c("24-hpiA", "24-hpiJ")
)

CLUSTER_IDS <- c(
  "Cluster 0", "Cluster 1",
  "Gill ciliary cells", "Hepatopancreas cells",
  "Gill neuroepithelial cells", "Gill cell type 1",
  "Hyalinocytes", "Haemocyte cell type 1",
  "Mantle cell type 1", "Cluster 9",
  "Vesicular haemocytes", "Immature haemocytes",
  "Macrophage-like cells", "Adductor muscle cells",
  "Mantle cell type 2", "Mantle epithelial cells",
  "Gill cell type 2", "Small granule cells"
)

FONT      <- "Arial"
BASE_SIZE <- 9
FIG_W_MM  <- 270      # A4 landscape
FIG_H_MM  <- 200
PT_SIZE   <- 0.15     # UMAP points; small, because facets are small

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

save_both <- function(plot, stem, width_mm, height_mm) {
  ggsave(paste0(stem, ".svg"), plot, width = width_mm, height = height_mm,
         units = "mm", device = svglite::svglite, limitsize = FALSE)
  fix_svg_for_inkscape(paste0(stem, ".svg"))
  ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, limitsize = FALSE)
}

theme_ms <- function() {
  theme_bw(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      text              = element_text(size = BASE_SIZE, colour = "black", family = FONT),
      plot.title        = element_text(size = BASE_SIZE, colour = "black",
                                       family = FONT, face = "bold"),
      axis.text         = element_text(size = BASE_SIZE, colour = "black"),
      axis.title        = element_text(size = BASE_SIZE, colour = "black"),
      strip.text        = element_text(size = BASE_SIZE, colour = "black"),
      legend.text       = element_text(size = BASE_SIZE, colour = "black"),
      legend.title      = element_blank(),
      legend.key.size   = unit(3.5, "mm"),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      panel.border      = element_rect(colour = "black", linewidth = 0.3, fill = NA),
      panel.grid        = element_blank(),
      axis.ticks        = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.5, "pt"))
}

###############################################################################
banner("Load and label")
###############################################################################

stopifnot(file.exists(SEU_RDS))
seu <- readRDS(SEU_RDS)

names(CLUSTER_IDS) <- levels(seu)
seu <- RenameIdents(seu, CLUSTER_IDS)
seu$cell_type <- factor(as.character(Idents(seu)), levels = CLUSTER_IDS)

if (length(DROP_SAMPLE)) {
  seu <- subset(seu, subset = sample %in% setdiff(SAMPLE_LEVELS, DROP_SAMPLE))
  cat("\n  dropped:", paste(DROP_SAMPLE, collapse = ", "), "\n")
}
seu$sample <- factor(as.character(seu$sample),
                     levels = setdiff(SAMPLE_LEVELS, DROP_SAMPLE))

seu$condition <- factor(case_when(
  seu$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu$sample %in% c("6-hpiA", "6-hpiD")         ~ "6hpi",
  seu$sample %in% c("24-hpiA", "24-hpiJ")       ~ "24hpi",
  seu$sample == "72-hpiJ"                       ~ "72hpi",
  seu$sample == "96-hpiE"                       ~ "96hpi"
), levels = COND_LEVELS)
stopifnot(!anyNA(seu$condition), !anyNA(seu$sample), !anyNA(seu$cell_type))

cat("\n  Nuclei per condition:\n"); print(table(seu$condition))
cat("\n  Nuclei per sample:\n");    print(table(seu$sample))

## Palette: 18 clusters need a qualitative set that survives being shrunk into
## small facets. Seurat's polychrome is built for this; hcl.colors is a fallback.
CL_COLS <- tryCatch(
  setNames(Seurat::DiscretePalette(length(CLUSTER_IDS), palette = "polychrome"),
           CLUSTER_IDS),
  error = function(e)
    setNames(hcl.colors(length(CLUSTER_IDS), "Dark 3"), CLUSTER_IDS))

###############################################################################
banner("Cluster proportions")
###############################################################################

md <- seu@meta.data %>% select(cell_type, condition, sample)

prop_cond <- md %>%
  count(condition, cell_type, name = "n_nuclei") %>%
  group_by(condition) %>%
  mutate(pct = round(100 * n_nuclei / sum(n_nuclei), 2)) %>%
  ungroup()

prop_samp <- md %>%
  count(sample, cell_type, name = "n_nuclei") %>%
  group_by(sample) %>%
  mutate(pct = round(100 * n_nuclei / sum(n_nuclei), 2)) %>%
  ungroup()

write_csv(prop_cond, file.path(OUT, "cluster_proportions_by_condition.csv"))
write_csv(prop_samp, file.path(OUT, "cluster_proportions_by_sample.csv"))

cat("\n  Percentage of nuclei per cluster, by condition:\n\n")
print(as.data.frame(prop_cond %>% select(-n_nuclei) %>%
                      pivot_wider(names_from = condition, values_from = pct)))

## The yardstick: how much does a cluster's share differ between two animals
## that received the same treatment? Anything smaller than this between stages
## is not interpretable.
within_pairs <- map_dfr(names(PAIRS), function(nm) {
  p <- PAIRS[[nm]]
  if (!all(p %in% levels(prop_samp$sample))) return(NULL)
  prop_samp %>% filter(sample %in% p) %>%
    group_by(cell_type) %>%
    summarise(spread_pct_points = round(diff(range(pct)), 2), .groups = "drop") %>%
    mutate(pair = nm, .before = 1)
}) %>% arrange(desc(spread_pct_points))

cat("\n  Within-pair spread in cluster share (percentage points):\n\n")
print(as.data.frame(within_pairs), row.names = FALSE)
write_csv(within_pairs, file.path(OUT, "within_pair_proportion_spread.csv"))

cat(sprintf("\n  Median within-pair spread: %.2f points | maximum: %.2f points\n",
            median(within_pairs$spread_pct_points),
            max(within_pairs$spread_pct_points)))

###############################################################################
banner("Figure")
###############################################################################

emb <- as.data.frame(Embeddings(seu, "umap")[, 1:2]) %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  bind_cols(md)

## Every panel carries all nuclei in pale grey behind its own subset, so a
## reader compares the coloured cells against the same silhouette each time
## rather than against a differently shaped cloud.
bg <- emb %>% select(UMAP_1, UMAP_2)

pA <- ggplot(emb, aes(UMAP_1, UMAP_2)) +
  geom_point(data = bg, colour = "grey88", size = PT_SIZE) +
  geom_point(aes(colour = cell_type), size = PT_SIZE) +
  facet_wrap(~ condition, nrow = 1) +
  scale_colour_manual(values = CL_COLS, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = 1)) +
  labs(title = "A  Nuclei from each infection stage, on the shared embedding",
       x = "UMAP 1", y = "UMAP 2") +
  theme_ms()

pB <- ggplot(prop_cond, aes(condition, pct, fill = cell_type)) +
  geom_col(position = position_stack(reverse = FALSE),
           colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = CL_COLS, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "B  Cluster composition by stage", x = NULL, y = "% of nuclei") +
  theme_ms() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

pC <- ggplot(prop_samp, aes(sample, pct, fill = cell_type)) +
  geom_col(position = position_stack(reverse = FALSE),
           colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = CL_COLS, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "C  Cluster composition by animal", x = NULL, y = "% of nuclei") +
  theme_ms() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

fig <- pA / (pB | pC) +
  plot_layout(heights = c(1, 0.9), guides = "collect")

save_both(fig, file.path(OUT, "FigureSX_umap_by_condition"), FIG_W_MM, FIG_H_MM)

###############################################################################
banner("SUMMARY")
###############################################################################

cat("\n  Nuclei per condition:\n"); print(table(seu$condition))
cat(sprintf("\n  Within-pair spread in cluster share: median %.2f, max %.2f points\n",
            median(within_pairs$spread_pct_points),
            max(within_pairs$spread_pct_points)))
cat("  Any between-stage difference smaller than this is not interpretable,\n")
cat("  because each stage is one animal.\n")

cat("\n  Written to ", OUT, ":\n", sep = "")
cat("    FigureSX_umap_by_condition.svg / .pdf\n")
cat("    cluster_proportions_by_condition.csv\n")
cat("    cluster_proportions_by_sample.csv\n")
cat("    within_pair_proportion_spread.csv\n\n")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

