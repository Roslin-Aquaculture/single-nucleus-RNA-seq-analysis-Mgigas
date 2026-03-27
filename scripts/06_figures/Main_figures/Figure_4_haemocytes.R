# Figure 4: stacked violin plots of haemocyte marker expression

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(Seurat.utils)

############################
# Load processed Seurat object and assign final cluster identities
############################

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object generated previously in scripts/05_seurat_processing/01_process_seurat_objects.R
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

############################
# Create simplified infection-stage groups and remove 24-hpiA sample
############################

# Add simplified condition labels for downstream comparisons
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

# Drop unused factor levels after subsetting
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

head(seu_obj_clean)

############################
# Add readable gene labels to the Seurat object
############################

# Load annotation tables from 06_figures/supporting_files
ORSON <- read_csv("01_ORSON_french_group_2024_biorxiv_suppl2.csv") %>%
  select(Gene.ID, Description = Sequence.Description)

cg_science <- read_tsv("02_Cg_gene_names.tsv", col_names = FALSE) %>%
  rename(Gene.ID = X1, Gene.symbol.science = X2)

# Merge annotation sources and build combined display labels
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

############################
# Figure 4: haemocyte marker panel
############################

# Marker set selected from Divonne et al 2025 used for the stacked violin plot
# Updated with macrophage-like markers from other sources
divonne_markers <- c(
  "G12639", "G12733", "G22387", "G5864", "G32289",
  "G2310", "G6983", "G29373", "G28068", "G29966",
  "G4972", "G7023", "G1172", "G16773",
  "G28618", "G1202",
  "G21444", "G22986", "G384",
  "G2459", "G29330", "G32588", "G2457"
)

# Marker grouping used during panel design
divonne_ann <- list(
  `Hyalinocytes` = c("G2459", "G29330", "G32588", "G2457"),
  `Haemocyte Cell Type 1` = c("G21444", "G22986", "G384"),
  `Vesicular haemocytes` = c("G28618", "G1202"),
  `Immature Haemocytes` = c("G4972", "G7023", "G1172", "G16773"),
  `Macrophage-like Cells` = c("G2310", "G6983", "G29373", "G28068", "G29966"),
  `Small Granule Cells` = c("G12639", "G12733", "G22387", "G5864", "G32289")
)

# Retrieve gene symbols for the marker panel and preserve the desired order
divonne_symbols <- gene_map_filtered %>%
  filter(gene_id %in% divonne_markers) %>%
  mutate(
    symbol_trimmed = sub("^\\S+\\s+", "", gene_symbol),
    gene_id = factor(gene_id, levels = divonne_markers)
  ) %>%
  arrange(gene_id)

# Prepare feature vector for plotting
features <- divonne_symbols$gene_symbol
names(features) <- divonne_symbols$symbol_trimmed

# Colour palette used for stacked violin plot groups
colour_palette <- c(
  "#EDC948", "#EDC948", "#EDC948", "#EDC948",
  "#59A14F", "#59A14F", "#59A14F",
  "#56B4E9", "#56B4E9",
  "#E15759", "#E15759", "#E15759", "#E15759",
  "#4E79A7", "#4E79A7", "#4E79A7", "#4E79A7", "#4E79A7",
  "#F28E2B", "#F28E2B", "#F28E2B", "#F28E2B", "#F28E2B"
)

# Create stacked violin plot
divonne_stackedvln <- VlnPlot(
  seu_obj_clean,
  features = rev(features),
  stack = TRUE,
  flip = TRUE,
  pt.size = 0,
  cols = colour_palette
) +
  NoLegend()

# Save figure
ggsave(
  "Figure_4_stacked_violin.pdf",
  plot = divonne_stackedvln,
  width = 12,
  height = 12,
  dpi = 600
)

# ------------------------------------------------------------------
# End of Figure 4 script
# ------------------------------------------------------------------
