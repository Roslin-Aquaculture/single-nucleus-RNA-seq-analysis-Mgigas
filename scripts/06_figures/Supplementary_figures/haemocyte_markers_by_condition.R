############################
# Supplement to Figure 4:
# haemocyte marker expression by cell type x infection stage
#   1. subset haemocyte populations
#   2. cluster composition tables
#   3. dot plot (cell type x condition)
#   4. pheatmap of scaled average expression
# All outputs go to OUT_DIR (svg for Inkscape, cairo pdf, png, csv).
############################

library(Seurat)
library(tidyverse)
library(Matrix)
library(pheatmap)
library(svglite)
library(systemfonts)
library(xml2)

############################
# 0. Paths, style, export helpers
############################

setwd("/home/pdewari/Documents/parse_2025/seurat_2025/")
OUT_DIR <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/haemocyte_markers_by_cond"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(OUT_DIR, f)

# Folder holding the two gene-annotation tables (01_ORSON..., 02_Cg_gene_names).
# Defaults to the working directory; change if they live elsewhere.
SUPPORT_DIR <- getwd()

FONT      <- "Arial"
BASE_SIZE <- 9
FIG_W_MM  <- 190

# Infection stage: ordered, fixed colours
COND_LEVELS <- c("control", "6hpi", "24hpi", "72hpi", "96hpi")
COND_COLS <- setNames(
  c("#4d4d4d", "#8c6bb1", "#41ab5d", "#f16913", "#ce1256"), COND_LEVELS)

# Cell-type colours (same as Figure 4)
celltype_cols <- c(
  "Hyalinocytes"          = "#EDC948",
  "Haemocyte cell type 1" = "#59A14F",
  "Vesicular haemocytes"  = "#56B4E9",
  "Immature haemocytes"   = "#E15759",
  "Macrophage-like cells" = "#4E79A7",
  "Small granule cells"   = "#F28E2B"
)

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

# ggplot objects: svg (Inkscape), cairo pdf (text stays text), png quick-look
save_all <- function(plot, stem, width_mm, height_mm) {
  ggsave(paste0(stem, ".svg"), plot, width = width_mm, height = height_mm,
         units = "mm", device = svglite::svglite, limitsize = FALSE)
  fix_svg_for_inkscape(paste0(stem, ".svg"))
  ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, limitsize = FALSE)
  ggsave(paste0(stem, ".png"), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 400, limitsize = FALSE)
}

# grid objects (pheatmap gtable): same three formats
draw_grob <- function(g) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = FONT)))
  grid::grid.draw(g)
  grid::popViewport()
}
save_grob <- function(g, stem, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  draw_grob(g); dev.off()
  fix_svg_for_inkscape(paste0(stem, ".svg"))
  cairo_pdf(paste0(stem, ".pdf"), width = w, height = h)
  draw_grob(g); dev.off()
  png(paste0(stem, ".png"), width = width_mm, height = height_mm,
      units = "mm", res = 400)
  draw_grob(g); dev.off()
}

theme_fig <- theme_minimal(base_size = BASE_SIZE, base_family = FONT) +
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
        legend.title    = element_text(size = BASE_SIZE, colour = "black"),
        legend.key.size = unit(3, "mm"))

############################
# 1. Load object, final annotations, condition labels
############################

seu_obj <- read_rds("seu_obj_umap_18d_6r_3kRes.rds")

# NOTE: assumes cluster order in this object matches the order below.
new.cluster.ids <- c(
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
names(new.cluster.ids) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new.cluster.ids)
seu_obj$celltype <- Idents(seu_obj)

seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample == "24-hpiA" ~ "mid?",
  seu_obj$sample == "24-hpiJ" ~ "24hpi",
  seu_obj$sample == "72-hpiJ" ~ "72hpi",
  seu_obj$sample == "96-hpiE" ~ "96hpi"
)
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove 24-hpiA ("mid?")
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")
rm(seu_obj)
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

############################
# 2. Marker panel and gene labels
############################

# Marker groups (Divonne et al. 2025), same order as the Figure 4 panel.
# Names must match the final cluster annotations exactly.
marker_groups <- list(
  `Small granule cells`   = c("G12639", "G12733", "G22387", "G5864", "G32289"),
  `Macrophage-like cells` = c("G2310", "G6983", "G29373", "G28068", "G29966"),
  `Immature haemocytes`   = c("G4972", "G7023", "G1172", "G16773"),
  `Vesicular haemocytes`  = c("G28618", "G1202"),
  `Haemocyte cell type 1` = c("G21444", "G22986", "G384"),
  `Hyalinocytes`          = c("G2459", "G29330", "G32588", "G2457")
)
hae_types  <- names(marker_groups)
marker_ids <- unlist(marker_groups, use.names = FALSE)

missing_ids <- setdiff(marker_ids, rownames(seu_obj_clean))
if (length(missing_ids) > 0) {
  stop("Markers not found in object: ", paste(missing_ids, collapse = ", "))
}

# Readable labels for display only (object keeps original gene IDs as rownames)
ORSON <- read_csv(file.path(SUPPORT_DIR, "01_ORSON_french_group_2024_biorxiv_suppl2.csv")) %>%
  select(Gene.ID, Description = Sequence.Description)

cg_science <- read_tsv(file.path(SUPPORT_DIR, "02_Cg_gene_names.tsv"), col_names = FALSE) %>%
  rename(Gene.ID = X1, Gene.symbol.science = X2)

gene_map <- full_join(ORSON, cg_science, by = "Gene.ID") %>%
  mutate(
    gene_symbol = paste(
      coalesce(Gene.symbol.science, ""),
      coalesce(Description, ""),
      sep = " : "
    ) %>%
      trimws() %>%
      gsub("[^[:alnum:] :]", " ", x = .)
  ) %>%
  select(gene_id = Gene.ID, gene_symbol) %>%
  distinct(gene_id, .keep_all = TRUE)

label_df <- gene_map %>%
  filter(gene_id %in% marker_ids) %>%
  mutate(label = paste(gene_id, str_trunc(gene_symbol, 50)))

label_vec <- setNames(label_df$label, label_df$gene_id)
label_vec[setdiff(marker_ids, names(label_vec))] <- setdiff(marker_ids, names(label_vec))  # fallback: bare ID
label_vec <- label_vec[marker_ids]

############################
# 3. Subset haemocyte populations
############################

hae <- subset(seu_obj_clean, subset = celltype %in% hae_types)
hae$celltype      <- factor(as.character(hae$celltype), levels = hae_types)
hae$condition_new <- droplevels(hae$condition_new)

############################
# 4. Cluster composition tables (within haemocyte subset)
############################

comp_counts <- table(Cluster = hae$celltype, Condition = hae$condition_new)

comp_long <- as.data.frame(comp_counts, responseName = "n_cells") %>%
  group_by(Condition) %>%
  mutate(
    n_condition_total = sum(n_cells),
    pct_of_condition  = 100 * n_cells / n_condition_total   # composition of each condition
  ) %>%
  group_by(Cluster) %>%
  mutate(
    pct_of_cluster = 100 * n_cells / sum(n_cells),          # raw condition share of each cluster
    # condition share of each cluster after equalising condition sizes
    # (raw pct_of_cluster is biased toward conditions with more cells)
    pct_of_cluster_size_norm = 100 * (n_cells / n_condition_total) /
      sum(n_cells / n_condition_total)
  ) %>%
  ungroup() %>%
  mutate(across(starts_with("pct_"), ~ round(.x, 2)))

comp_wide_counts <- as.data.frame.matrix(comp_counts) %>%
  rownames_to_column("Cluster") %>%
  mutate(Total = rowSums(across(-Cluster)))

comp_wide_pct <- as.data.frame.matrix(round(100 * prop.table(comp_counts, margin = 2), 1)) %>%
  rownames_to_column("Cluster")

comp_sample <- as.data.frame.matrix(table(Cluster = hae$celltype, Sample = hae$sample)) %>%
  rownames_to_column("Cluster")

print(comp_wide_counts)
print(comp_wide_pct)

write_csv(comp_long,        out("haemocyte_composition_long.csv"))
write_csv(comp_wide_counts, out("haemocyte_composition_counts.csv"))
write_csv(comp_wide_pct,    out("haemocyte_composition_pct_of_condition.csv"))
write_csv(comp_sample,      out("haemocyte_composition_by_sample.csv"))

############################
# 5. Cell type x condition groups (drop small groups)
############################

min_cells <- 20   # groups with fewer cells give noisy averages and dots

hae$ct_cond <- paste(hae$celltype, hae$condition_new, sep = " | ")

grp_levels <- expand_grid(ct = hae_types, cond = levels(hae$condition_new)) %>%
  mutate(g = paste(ct, cond, sep = " | ")) %>%
  pull(g)

grp_n <- table(factor(hae$ct_cond, levels = grp_levels))
write_csv(
  tibble(group = names(grp_n), n_cells = as.integer(grp_n), kept = as.integer(grp_n) >= min_cells),
  out("haemocyte_group_cell_counts.csv")
)

grp_keep <- grp_levels[grp_n >= min_cells]
hae_f <- subset(hae, cells = colnames(hae)[hae$ct_cond %in% grp_keep])
hae_f$ct_cond <- factor(hae_f$ct_cond, levels = grp_keep)

############################
# 6. Dot plot
############################
# DotPlot(scale = TRUE) z-scores average expression across the displayed groups,
# so colour is comparable to the heatmap below.

dp <- DotPlot(
  hae_f,
  features = marker_ids,
  group.by = "ct_cond",
  cols = c("grey90", "#B2182B"),
  dot.scale = 5,
  cluster.idents = FALSE
) +
  scale_x_discrete(labels = label_vec) +
  scale_y_discrete(limits = rev(levels(hae_f$ct_cond))) +
  labs(x = NULL, y = NULL) +
  theme_fig +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.2)) +
  RotatedAxis()

print(dp)
save_all(dp, out("Figure_4_supp_dotplot_celltype_by_condition"),
         width_mm = FIG_W_MM, height_mm = 190)

############################
# 7. Pheatmap: scaled average expression
############################

# Average expression per ct_cond.
# Averaging is done on the non-log scale (expm1 of the log-normalised data),
# then log1p-transformed, then z-scored per gene across groups.
# Seurat v5 uses `layer`; for Seurat v4 use `slot = "data"`.
# If the assay still has split layers (v5), run JoinLayers() first.
expr <- GetAssayData(hae_f, assay = DefaultAssay(hae_f), layer = "data")[marker_ids, , drop = FALSE]

grp    <- hae_f$ct_cond
design <- Matrix::sparse.model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
n_per <- Matrix::colSums(design)

sums    <- as.matrix(expm1(expr) %*% design)
avg_log <- log1p(sweep(sums, 2, n_per, "/"))

# Per-gene z-score across groups; zero-variance genes -> 0; clip at +/-2
mat <- t(scale(t(avg_log)))
mat[is.nan(mat)] <- 0
mat <- pmax(pmin(mat, 2), -2)
rownames(mat) <- label_vec[rownames(mat)]

write.csv(avg_log, out("haemocyte_markers_avg_log_expression.csv"))
write.csv(mat,     out("haemocyte_markers_scaled_avg_expression.csv"))

# Column annotation: cell type and infection stage
ann_col <- hae_f@meta.data %>%
  as_tibble() %>%   # drops cell barcodes used as rownames
  transmute(ct_cond = as.character(ct_cond),
            Cell_type = celltype,
            Stage = condition_new) %>%
  distinct() %>%
  column_to_rownames("ct_cond")
ann_col <- ann_col[colnames(mat), , drop = FALSE]

# Row annotation: marker set
ann_row <- data.frame(
  Marker_set = factor(rep(hae_types, lengths(marker_groups)), levels = hae_types),
  row.names  = label_vec[marker_ids]
)

ann_colours <- list(
  Cell_type  = celltype_cols,
  Stage      = COND_COLS[levels(ann_col$Stage)],
  Marker_set = celltype_cols
)

# Gaps between cell types (columns) and marker sets (rows)
gaps_col <- head(cumsum(rle(as.character(ann_col$Cell_type))$lengths), -1)
gaps_row <- head(cumsum(rle(as.character(ann_row$Marker_set))$lengths), -1)

ph <- pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_row = ann_row,
  annotation_colors = ann_colours,
  gaps_col = gaps_col,
  gaps_row = gaps_row,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks = seq(-2, 2, length.out = 101),
  border_color = NA,
  cellwidth = 8,     # points
  cellheight = 9,    # points
  fontsize = BASE_SIZE,
  main = "Haemocyte markers: scaled average expression (cell type x stage)",
  silent = TRUE      # build only; drawn and saved below
)

draw_grob(ph$gtable)
save_grob(ph$gtable, out("Figure_4_supp_pheatmap_celltype_by_condition"),
          width_mm = FIG_W_MM + 50, height_mm = 160)


############################
# Supplementary: haemocyte composition by infection stage (stacked bar)
############################

comp_plot <- comp_long %>%
  mutate(Cluster   = factor(Cluster,   levels = hae_types),
         Condition = factor(Condition, levels = COND_LEVELS))

n_lab <- comp_plot %>% distinct(Condition, n_condition_total)

p_comp <- ggplot(comp_plot, aes(x = Condition, y = pct_of_condition, fill = Cluster)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  # % label inside segments, skipped when the segment is too small to read
  geom_text(aes(label = ifelse(pct_of_condition >= 5, sprintf("%.0f", pct_of_condition), "")),
            position = position_stack(vjust = 0.5),
            size = 7 / ggplot2::.pt, family = FONT, colour = "black") +
  # total haemocytes per stage above each bar
  geom_text(data = n_lab,
            aes(x = Condition, y = 100, label = paste0("n = ", n_condition_total)),
            inherit.aes = FALSE, vjust = -0.6,
            size = 7 / ggplot2::.pt, family = FONT, colour = "black") +
  scale_fill_manual(values = celltype_cols, breaks = hae_types) +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "% of haemocytes") +
  theme_fig +
  theme(legend.title = element_blank())

print(p_comp)
save_all(p_comp, out("Figure_S_haemocyte_composition_stacked_bar"),
         width_mm = 120, height_mm = 90)
# ------------------------------------------------------------------
# End of script
# ------------------------------------------------------------------
