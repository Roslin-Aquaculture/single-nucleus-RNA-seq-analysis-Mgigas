# =========================================================
# Supplementary figure: IAP paralogue expression across clusters
# =========================================================
# Control nuclei only. The panel shows baseline expression across cell
# types, not the infection response: pooling infected nuclei would let
# infection-responsive genes into the identity signature, and most
# nuclei in the dataset come from infected animals.
#
# Saves SVG (for Inkscape; text stays as text) and PDF via Cairo.

library(Seurat)
library(tidyverse)
library(Seurat.utils)   # RenameGenesSeurat

# =========================================================
# CONFIG
# =========================================================
OBJ  <- "/home/pdewari/Documents/parse_2025/seurat_2025/seu_obj_umap_18d_6r_3kRes.rds"
ANNO <- "/home/pdewari/Documents/parse_2025/seurat_2025"
OUT  <- "/home/pdewari/eggnog/results/figures"
FONT <- "Arial"         # "Liberation Sans" if Arial is unavailable

cluster_order <- c(
  "Cluster 0", "Cluster 1",
  "Gill ciliary cells", "Hepatopancreas cells",
  "Gill neuroepithelial cells", "Gill cell type 1",
  "Hyalinocytes", "Haemocyte cell type 1",
  "Mantle cell type 1", "Cluster 9",
  "Vesicular haemocytes", "Immature haemocytes",
  "Macrophage like cells", "Adductor muscle cells",
  "Mantle cell type 2", "Mantle epithelial cells",
  "Gill cell type 2", "Small granule cells"
)

gene_list_iap <- c(
  "G19384", "G19385", "G19421", "G19423", "G19424", "G19035",
  "G19414", "G20059", "G17608", "G15123", "G23947", "G25799"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# OBJECT: cluster names, then control nuclei only
# =========================================================
seu <- readRDS(OBJ)

names(cluster_order) <- levels(seu)
seu <- RenameIdents(seu, cluster_order)
cluster_order <- unname(cluster_order)

stopifnot(length(setdiff(cluster_order, as.character(unique(Idents(seu))))) == 0)

# 24-hpiA is excluded by not being a control
seu <- subset(seu, subset = sample %in% c("Homogenate", "Uninfected"))

Idents(seu) <- factor(as.character(Idents(seu)), levels = cluster_order)
stopifnot(!anyNA(Idents(seu)))

cat("control nuclei:", ncol(seu), "\n")

# =========================================================
# GENE LABELS: gene ID plus symbol and description
# =========================================================
# Skipped if the object has already been renamed (renamed rownames
# contain spaces; raw ones are bare IDs).
if (!any(grepl("\\s", head(rownames(seu), 100)))) {
  
  ORSON <- read_csv(file.path(ANNO, "01_ORSON_french_group_2024_biorxiv_suppl2.csv"),
                    show_col_types = FALSE) %>%
    select(Gene.ID, Description = Sequence.Description)
  
  cg_science <- read_tsv(file.path(ANNO, "02_Cg_gene_names.tsv"),
                         col_names = FALSE, show_col_types = FALSE) %>%
    rename(Gene.ID = X1, Gene.symbol.science = X2)
  
  gene_map <- full_join(ORSON, cg_science, by = "Gene.ID") %>%
    mutate(gene_symbol = paste(coalesce(Gene.symbol.science, ""),
                               coalesce(Description, ""), sep = " : ") %>%
             trimws() %>% gsub("[^[:alnum:] :]", " ", x = .)) %>%
    select(gene_id = Gene.ID, gene_symbol) %>%
    group_by(gene_symbol) %>%
    mutate(gene_symbol = if (n() > 1)
      paste0(gene_symbol, LETTERS[seq_along(gene_symbol)]) else gene_symbol) %>%
    ungroup() %>%
    mutate(gene_symbol = paste(gene_id, gene_symbol, sep = " ")) %>%
    filter(gene_id %in% rownames(seu))
  
  rename_vector <- setNames(
    gene_map$gene_symbol[match(rownames(seu), gene_map$gene_id)],
    rownames(seu))
  rename_vector[is.na(rename_vector)] <- names(rename_vector)[is.na(rename_vector)]
  rename_vector <- make.unique(rename_vector)
  
  seu <- RenameGenesSeurat(obj = seu, newnames = rename_vector)
  cat("genes renamed\n")
}

# =========================================================
# FEATURES: one rowname per gene ID
# =========================================================
genes_iap <- vapply(gene_list_iap, function(g) {
  m <- grep(paste0("^", g, "\\b"), rownames(seu), value = TRUE)
  if (length(m)) m[1] else NA_character_
}, character(1))

if (anyNA(genes_iap))
  cat("NOT FOUND:", paste(gene_list_iap[is.na(genes_iap)], collapse = ", "), "\n")

genes_iap     <- unname(genes_iap[!is.na(genes_iap)])
trimmed_names <- sub("^(\\S+)\\s+(.*)$", "\\1 (\\2)", sub(" *:.*", "", genes_iap))

stopifnot(!any(duplicated(genes_iap)))
print(data.frame(feature = genes_iap, label = trimmed_names))


# =========================================================
# PANEL SIZE
# =========================================================
# Exported at the size the panel should occupy on the A4 sheet, so it can
# be imported into Inkscape and positioned without rescaling. Rescaling
# after import changes the text size relative to the other panels, which
# is what breaks consistency across a composite figure.
PANEL_W_MM <- 180     # full page width, portrait A4 with margins
PANEL_H_MM <- 110     # adjust until the dot spacing looks right

# =========================================================
# SVG POST-PROCESSING
# =========================================================
# svglite writes shared stroke and fill properties into a single CSS block
# in <defs> and leaves them off the individual elements, which carry only
# stroke-width. Inkscape honours that stylesheet when rendering, but on
# ungroup it writes explicit style attributes that override the CSS rule,
# so any element relying on it loses its colour. Axis lines and ticks are
# drawn as polylines and are the ones affected: they survive the import
# and vanish on ungroup. This copies the inherited stroke and fill onto
# every drawing element so nothing depends on the stylesheet.
inline_svg_styles <- function(file) {
  x  <- xml2::read_xml(file)
  ns <- xml2::xml_ns(x)
  els <- xml2::xml_find_all(
    x, "//d1:line|//d1:polyline|//d1:polygon|//d1:path|//d1:rect|//d1:circle", ns)
  n <- 0
  for (e in els) {
    s <- xml2::xml_attr(e, "style")
    if (is.na(s)) s <- ""
    changed <- FALSE
    if (!grepl("stroke\\s*:", s)) {
      s <- paste0(trimws(s), if (nzchar(trimws(s))) " " else "", "stroke: #000000;")
      changed <- TRUE
    }
    if (!grepl("fill\\s*:", s)) {
      s <- paste0(trimws(s), " fill: none;")
      changed <- TRUE
    }
    if (changed) { xml2::xml_set_attr(e, "style", s); n <- n + 1 }
  }
  xml2::write_xml(x, file)
  cat("patched", n, "of", length(els), "elements\n")
}

# =========================================================
# PLOT
# =========================================================
# Text sizes are what they will be on the page: 8 pt axis labels, 7 pt
# legend. Most journals set a 5-7 pt minimum.
figS_iap <- DotPlot(seu, features = genes_iap) +
  coord_flip() +
  scale_x_discrete(labels = trimmed_names) +
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"),
    name = "Average\nexpression") +
  guides(colour = guide_colourbar(frame.colour = NA, ticks = FALSE)) +
  theme_minimal(base_size = 9, base_family = FONT) +
  theme(
    text              = element_text(size = 9, colour = "black", family = FONT),
    axis.line         = element_line(colour = "black", linewidth = 0.3),
    axis.ticks        = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(1.5, "pt"),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 9, colour = "black"),
    axis.text.y       = element_text(size = 9, colour = "black"),
    axis.title        = element_blank(),
    legend.text       = element_text(size = 9, colour = "black"),
    legend.title      = element_text(size = 9, colour = "black"),
    legend.key.height = unit(12, "pt"),
    legend.key.width  = unit(7, "pt"),
    panel.grid.major  = element_line(colour = "grey94", linewidth = 0.25),
    panel.grid.minor  = element_blank(),
    legend.position   = "right")

# =========================================================
# SAVE
# =========================================================
svg_path <- file.path(OUT, "FigS_IAP_paralogues.svg")
pdf_path <- file.path(OUT, "FigS_IAP_paralogues.pdf")

ggsave(svg_path, figS_iap, width = PANEL_W_MM, height = PANEL_H_MM,
       units = "mm", device = svglite::svglite, limitsize = FALSE)
inline_svg_styles(svg_path)

# PDF as a fallback import: svglite writes dimensions assuming 72 dpi while
# Inkscape assumes 96, so an SVG can arrive at 75% of the stated size. PDF
# carries physical units unambiguously. Check the imported width once and
# use whichever format lands at the right size.
ggsave(pdf_path, figS_iap, width = PANEL_W_MM, height = PANEL_H_MM,
       units = "mm", device = cairo_pdf, limitsize = FALSE)

cat("written to:", OUT, sprintf("(%.0f x %.0f mm)\n", PANEL_W_MM, PANEL_H_MM))

figS_iap

