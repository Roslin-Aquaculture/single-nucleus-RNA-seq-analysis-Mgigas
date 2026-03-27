# Figure 3: Marker gene overview across annotated oyster cell clusters

# ------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(Seurat.utils)

# ------------------------------------------------------------------
# Load processed Seurat object and assign final cluster identities
# ------------------------------------------------------------------

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object; this was created previously under scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds") 


# Rename clusters with final cell-type annotations
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

# ------------------------------------------------------------------
# Create simplified infection-stage groups and remove 24-hpiA sample
# ------------------------------------------------------------------

# Add simplified condition column
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample == "24-hpiA" ~ "mid?",
  seu_obj$sample == "24-hpiJ" ~ "24hpi",
  seu_obj$sample == "72-hpiJ" ~ "72hpi",
  seu_obj$sample == "96-hpiE" ~ "96hpi"
)

# Set factor order
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove 24-hpiA ("mid?")
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")
rm(seu_obj)

# Drop unused factor levels
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

head(seu_obj_clean)


# ------------------------------------------------------------------
# Add readable gene symbols to the Seurat object
# ------------------------------------------------------------------

# Replace gene IDs with combined gene symbol / description labels to improve readability in downstream plots

# Load annotation tables, these can be found in 06_figures/supporting_files
ORSON <- read_csv("01_ORSON_french_group_2024_biorxiv_suppl2.csv") %>%
  select(Gene.ID, Description = Sequence.Description)

cg_science <- read_tsv("02_Cg_gene_names.tsv", col_names = FALSE) %>%
  rename(Gene.ID = X1, Gene.symbol.science = X2)

# Merge annotations and build combined display labels
gene_map2 <- full_join(ORSON, cg_science, by = "Gene.ID") %>%
  mutate(
    gene_symbol = paste(
      coalesce(Gene.symbol.science, ""),
      coalesce(Description, ""),
      sep = " : "
    ) %>%
      trimws() %>%
      gsub("[^[:alnum:] :]", " ", x = .)
  ) %>%
  select(gene_id = Gene.ID, gene_symbol)

# Make duplicated gene labels unique
gene_map2 <- gene_map2 %>%
  group_by(gene_symbol) %>%
  mutate(
    gene_symbol = if (n() > 1) paste0(gene_symbol, LETTERS[seq_along(gene_symbol)]) else gene_symbol
  ) %>%
  ungroup()

# Prepend gene IDs to the display labels
gene_map2 <- gene_map2 %>%
  mutate(gene_symbol = paste(gene_id, gene_symbol, sep = " "))

# Keep only genes present in the Seurat object
gene_map_filtered <- gene_map2 %>%
  filter(gene_id %in% rownames(seu_obj_clean))

# Build rename vector for Seurat
gene_symbols <- gene_map_filtered$gene_symbol[match(rownames(seu_obj_clean), gene_map_filtered$gene_id)]
rename_vector <- setNames(gene_symbols, rownames(seu_obj_clean))
rename_vector[is.na(rename_vector)] <- names(rename_vector)[is.na(rename_vector)]
rename_vector <- make.unique(rename_vector)

# Rename genes in the Seurat object
seu_obj_clean <- RenameGenesSeurat(
  obj = seu_obj_clean,
  newnames = rename_vector
)

rownames(seu_obj_clean)

# ------------------------------------------------------------------
# Figure 3A: DotPlot of selected marker genes across all clusters
# ------------------------------------------------------------------

# Two key marker genes selected for each cluster
gene_list_dotplot <- c(
  "G21268", "G1208", "G11503", "G7179", "G21910", "G12025", "G3157", "G13652",
  "G6819", "G31525", "G2457", "G32588", "G384", "G22071", "G10190", "G25152",
  "G10388", "G4381", "G18439", "G25991", "G22661", "G26978", "G6983", "G5032",
  "G25637", "G3128", "G4180", "G7827", "G15965", "G15964", "G22673", "G27503",
  "G12733", "G5864"
)


# Match selected gene IDs to the renamed Seurat rownames
genes_dotplot <- unlist(
  lapply(gene_list_dotplot, function(g) {
    matches <- grep(paste0("^", g, "\\b"), rownames(seu_obj_clean), value = TRUE)
    if (length(matches) > 0) return(matches) else return(NA)
  })
)

print(genes_dotplot)

# Trim long gene labels for plotting
trimmed_names <- sub(" *:.*", "", genes_dotplot)

# Create DotPlot
fig3a <- DotPlot(seu_obj_clean, features = genes_dotplot) +
  scale_x_discrete(labels = trimmed_names) +
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  ggtitle("Figure 3A. Top marker genes for the 17 transcriptomic clusters")

ggsave("Figure_3A_DotPlot.pdf", plot = fig3a, width = 12, height = 6, units = "in", dpi = 300)

# ------------------------------------------------------------------
# Figure 3B: FeaturePlot of selected marker genes for representative cell types
# ------------------------------------------------------------------

# Selected markers for representative cell populations
gene_list_featurePlot <- c(
  "G21910", "G32588", "G384", 
  "G6983", "G25637", "G7827", 
  "G15965", "G22673", "G12733"
)

print(gene_list_featurePlot)

# Match selected gene IDs to the renamed Seurat rownames
full_gene_list_featurePlot <- unlist(
  lapply(gene_list_featurePlot, function(g) {
    matches <- grep(paste0("^", g, "\\b"), rownames(seu_obj_clean), value = TRUE)
    if (length(matches) > 0) return(matches) else return(NA)
  })
)

cols = c("#EDF4F8", "#B2182B")

fig3b <- FeaturePlot(
  seu_obj_clean,
  features = full_gene_list_featurePlot,
  ncol = 3
) &
  scale_color_gradientn(colours = cols) &
  theme(
    plot.title = element_text(size = 11)
  )

ggsave("Figure_3B_FeaturePlot.pdf", plot = fig3b, width = 12, height = 10, units = "in", dpi = 300)

# ------------------------------------------------------------------
# End of Figure 3 script
# ------------------------------------------------------------------
