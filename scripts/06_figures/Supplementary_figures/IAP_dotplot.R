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
# PLOT
# =========================================================
igS_iap <- DotPlot(seu, features = genes_iap) +
  coord_flip() +
  scale_x_discrete(labels = trimmed_names) +
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"),
    name = "Average\nexpression") +
  guides(colour = guide_colourbar(frame.colour = NA, ticks = FALSE)) +
  theme_minimal(base_size = 12, base_family = FONT) +
  theme(
    text              = element_text(colour = "black", family = FONT),
    axis.line         = element_line(colour = "black", linewidth = 0.3),
    axis.ticks        = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(2, "pt"),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 11, colour = "black"),
    axis.text.y       = element_text(size = 11, colour = "black"),
    axis.title        = element_blank(),
    legend.text       = element_text(size = 10, colour = "black"),
    legend.title      = element_text(size = 10, colour = "black"),
    legend.key.height = unit(12, "pt"),
    legend.key.width  = unit(8, "pt"),
    panel.grid.major  = element_line(colour = "grey94", linewidth = 0.25),
    panel.grid.minor  = element_blank(),
    legend.position   = "right")

W <- 3.5 + 0.34 * length(cluster_order)   
H <- 2 + 0.34 * length(genes_iap)       

ggsave(file.path(OUT, "FigS_IAP_paralogues.svg"), figS_iap,
       width = W, height = H, device = svglite::svglite, limitsize = FALSE)
ggsave(file.path(OUT, "FigS_IAP_paralogues.pdf"), figS_iap,
       width = W, height = H, device = cairo_pdf, limitsize = FALSE)

cat("written to:", OUT, sprintf("(%.1f x %.1f in)\n", W, H))

figS_iap
