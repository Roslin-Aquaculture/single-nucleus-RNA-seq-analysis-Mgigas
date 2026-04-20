## Suppl Figure S4 : Find conserved markers for haemocyte group


# This script identifies conserved marker genes for the haemocyte lineage ...
#..(six haemocyte clusters) compared to all other cell types.

# Haemocyte group included these six clusters:
# "Immature haemocytes", "Macrophage like cells","Vesicular haemocytes", 
# "Haemocyte cell type 1", "Hyalinocytes", "Small granule cells"




################## Seurat object preparation ########################################

# Load packages
library(Seurat)
library(tidyverse)
library(scCustomize)
library(Seurat.utils)
library(pheatmap)


##### Load Seurat object and prepare metadata

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object; this was created previously under scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds") 


# Assign cluster identities
new.cluster.ids <- c("Cluster 0", "Cluster 1", 
                     "Gill ciliary cells", "Hepatopancreas cells",
                     "Gill neuroepithelial cells", "Gill cell type 1",
                     "Hyalinocytes", "Haemocyte cell type 1",
                     "Mantle cell type 1", "Cluster 9", 
                     "Vesicular haemocytes","Immature haemocytes",
                     "Macrophage like cells","Adductor muscle cells",
                     "Mantle cell type 2", "Mantle epithelial cells",
                     "Gill cell type 2", "Small granule cells")

names(new.cluster.ids) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new.cluster.ids)

# Create simplified infection-stage metadata for grouping
# Note: 24-hpiA showed inconsistent qPCR signal and is flagged as "mid?"
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample %in% "24-hpiA" ~ "mid?",
  seu_obj$sample %in% "24-hpiJ" ~ "24hpi",
  seu_obj$sample %in% "72-hpiJ" ~ "72hpi",
  seu_obj$sample %in% "96-hpiE" ~ "96hpi",
  TRUE ~ "unknown"
)

# Set factor levels for consistent ordering
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove ambiguous 24-hpiA cells ("mid?")
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")

rm(seu_obj)

# Drop unused factor levels after subsetting
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

################## Add Gene annotation  #######################

# Load gene annotation sources
ORSON <- read_csv("01_ORSON_french_group_2024_biorxiv_suppl2.csv") %>%
  select(Gene.ID, Description = Sequence.Description)

cg_science <- read_tsv("02_Cg_gene_names.tsv", col_names = FALSE) %>%
  rename(Gene.ID = X1, Gene.symbol.science = X2)

# Merge annotations and construct combined gene labels
merged_df <- full_join(ORSON, cg_science, by = "Gene.ID") %>%
  mutate(
    gene_symbol = paste(
      coalesce(Gene.symbol.science, ""), 
      coalesce(Description, ""), 
      sep = " : "
    ) %>% trimws()
  ) %>%
  mutate(
    gene_symbol = gsub("[^[:alnum:] :]", " ", gene_symbol)
  )

gene_map2 <- merged_df %>%
  select(gene_id = Gene.ID, gene_symbol)

# Ensure unique gene labels to avoid downstream indexing errors
gene_map2 <- gene_map2 %>%
  group_by(gene_symbol) %>%
  mutate(
    gene_symbol = if(n() > 1) paste0(gene_symbol, LETTERS[seq_along(gene_symbol)]) else gene_symbol
  ) %>%
  ungroup()

# Prepend gene IDs to labels for traceability
gene_map2 <- gene_map2 %>%
  mutate(gene_symbol = paste(gene_id, gene_symbol, sep = " "))

# Keep only genes present in Seurat object
gene_map_filtered <- gene_map2 %>%
  filter(gene_id %in% rownames(seu_obj_clean))

# Build renaming vector
gene_symbols <- gene_map_filtered$gene_symbol[
  match(rownames(seu_obj_clean), gene_map_filtered$gene_id)
]

rename_vector <- setNames(gene_symbols, rownames(seu_obj_clean))
rename_vector[is.na(rename_vector)] <- names(rename_vector)[is.na(rename_vector)]
rename_vector <- make.unique(rename_vector)

# Apply gene renaming
seu_obj_clean <- RenameGenesSeurat(seu_obj_clean, newnames = rename_vector)



################## Cluster annotation for plotting ###################################

# Map numeric clusters to descriptive labels
cluster_names <- c("Cluster 0", "Cluster 1", 
                   "Gill ciliary cells", "Hepatopancreas cells",
                   "Gill neuroepithelial cells", "Gill cell type 1",
                   "Hyalinocytes", "Haemocyte cell type 1",
                   "Mantle cell type 1", "Cluster 9", 
                   "Vesicular haemocytes","Immature haemocytes",
                   "Macrophage like cells","Adductor muscle cells",
                   "Mantle cell type 2", "Mantle epithelial cells",
                   "Gill cell type 2", "Small granule cells")

cluster_index <- as.numeric(as.character(seu_obj_clean$seurat_clusters))
seu_obj_clean$cluster_names <- cluster_names[cluster_index + 1]

# Reorder clusters so haemocyte group appear last in plots
seu_obj_clean$cluster_names <- as.character(seu_obj_clean$cluster_names)

haemocyte_clusters <- c(
  "Immature haemocytes",
  "Macrophage like cells",
  "Vesicular haemocytes",
  "Haemocyte cell type 1",
  "Hyalinocytes",
  "Small granule cells"
)

non_haemocyte <- setdiff(unique(seu_obj_clean$cluster_names), haemocyte_clusters)
new_order <- c(non_haemocyte, haemocyte_clusters)

seu_obj_clean$cluster_names <- factor(seu_obj_clean$cluster_names, levels = new_order)

################## Differential expression ###########################################

# Collapse haemocyte clusters into a single identity class
haem_clusters <- haemocyte_clusters

seu_obj_clean$haem_group <- ifelse(
  Idents(seu_obj_clean) %in% haem_clusters,
  "haemocytes",
  "others"
)

Idents(seu_obj_clean) <- "haem_group"

# Identify conserved markers across samples
conserved_markers <- FindConservedMarkers(
  seu_obj_clean,
  ident.1 = "haemocytes",
  ident.2 = "others",
  grouping.var = "sample",
  only.pos = TRUE
)

# Format results and rank genes
conserved_markers_df <- conserved_markers %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  arrange(max_pval)

# Select top 50 conserved markers based on minimum p-value
top50_conserved <- conserved_markers_df %>%
  arrange(minimump_p_val) %>%
  slice_head(n = 50)

# Compute mean log2 fold change across conditions
logfc_cols <- grep("avg_log2FC", colnames(conserved_markers_df), value = TRUE)

conserved_markers_df <- conserved_markers_df %>%
  mutate(mean_log2FC = rowMeans(across(all_of(logfc_cols))))

################## Heatmap visualization #############################################

# Define diverging color palette
pd_palette <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#ABD9E9", "#E0F3F8", "#F7F7F7",
  "#FEE090", "#FDAE61", "#F46D43",
  "#D73027", "#A50026"
))

# Compute average expression per cluster
avg_exp <- AverageExpression(
  seu_obj_clean,
  group.by = "cluster_names",
  assays = "RNA",
  slot = "data"
)$RNA

# Subset matrix to selected genes
genes_present <- intersect(top50_conserved$gene, rownames(avg_exp))
mat <- avg_exp[genes_present, , drop = FALSE]

# Row-wise scaling (z-score per gene across clusters)
mat_scaled <- t(scale(t(mat)))

# Clean gene labels for display (remove gene IDs)
rownames(mat_scaled) <- sub("^G[0-9]+\\s*", "", rownames(mat_scaled))

# Save heatmap
pdf("conserved_markers_heatmap_top50.pdf",
    width = 14,
    height = 0.25 * nrow(mat_scaled) + 3)

pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = pd_palette(256),
  fontsize_row = 6,
  fontsize_col = 10,
  angle_col = 45,
  border_color = NA,
  cellwidth = 16,
  cellheight = 8,
  treeheight_row = 15,
  treeheight_col = 0,
  show_colnames = TRUE,
  show_rownames = TRUE,
  legend = TRUE,
  na_col = "grey95"
)

dev.off()

# ------------------------------------------------------------------
# End of Figure S4 script
# ------------------------------------------------------------------
