# Figure S6: Create Heatmaps For Differentially Expressed Genes in Each Cluster

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(Matrix)
library(dittoSeq)
library(Seurat.utils)


############################
# Load Seurat object and prepare metadata for group comparisons
############################

#...............................................................
# Original local working path used during analysis
#setwd("/home/pdewari/Documents/parse_2025/seurat_2025/")

# Original object load used during analysis
#seu_obj <- read_rds("seu_obj_umap_18d_6r_3kRes.rds")
#...............................................................

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object created previously in:
# scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds")

# Rename numeric clusters with final cell-type labels
new.cluster.ids <- c(
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

names(new.cluster.ids) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new.cluster.ids)

############################
# Create simplified infection-stage groups and remove 24-hpiA sample
############################

# Collapse samples into broader infection-stage groups for comparison
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample == "24-hpiA" ~ "mid?",
  seu_obj$sample == "24-hpiJ" ~ "24hpi",
  seu_obj$sample == "72-hpiJ" ~ "72hpi",
  seu_obj$sample == "96-hpiE" ~ "96hpi"
)

# Set the desired order of infection stages
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove the 24-hpiA sample ("mid?") from downstream analysis
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")
rm(seu_obj)

# Drop unused factor levels after subsetting
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

# Create a combined identity column: cluster + infection stage
# This is used to define matched cluster/stage groups
seu_obj_clean$sample_subtypes <- paste(Idents(seu_obj_clean), seu_obj_clean$condition_new, sep = "_")
Idents(seu_obj_clean) <- "sample_subtypes"

unique(seu_obj_clean$sample_subtypes)
unique(seu_obj_clean$condition_new)

############################
# Add readable gene symbols to the Seurat object
############################

# Replace gene IDs with combined gene symbol / description labels
# to improve readability in heatmaps and output tables

# Load annotation tables; these can be found in 06_figures/supporting_files
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

# Make duplicate display labels unique
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


################
# Generate aggregated expression data once
cluster_aggregate_exp <- AggregateExpression(seu_obj_clean, return.seurat = TRUE, group.by = "condition_new")
# Make condition_new a factor with the desired order
cluster_aggregate_exp$condition_new <- factor(
  cluster_aggregate_exp$condition_new,
  levels = c("control", "6hpi", "24hpi", "72hpi", "96hpi")
)

expr_mat <- GetAssayData(cluster_aggregate_exp, slot = "data")


# create directories to save heatmaps; and batch plot

base_dir <- "de_plots_full_ann_30032026/" # same as base_dir in scripts/06_figures/Main_figures/Figure_7_DEGs_heatmap.R

# Ensure necessary directories exist
unified_dir <- file.path(base_dir, "unified_gene_lists_v1")
dir.create(unified_dir, showWarnings = T)

heatmap_dir <- file.path(unified_dir, "heatmap_plots_v1")
dir.create(heatmap_dir, showWarnings = FALSE)


# List all DE gene files
de_files <- list.files(base_dir, pattern = "_de_genes_all.txt$", recursive = TRUE, full.names = TRUE)

# Initialize a list to store genes per cluster
cluster_gene_list <- list()

# Populate the list
for (file in de_files) {
  # Extract cluster name from the file path
  cluster_name <- basename(dirname(file))
  
  # Read genes from the file
  genes <- readLines(file)
  
  # Append genes to the corresponding cluster
  if (cluster_name %in% names(cluster_gene_list)) {
    cluster_gene_list[[cluster_name]] <- unique(c(cluster_gene_list[[cluster_name]], genes))
  } else {
    cluster_gene_list[[cluster_name]] <- unique(genes)
  }
}

# Loop through each cluster
processed_clusters <- unique(sub("_vs_.*", "", names(cluster_gene_list)))


pd_palette <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#ABD9E9", "#E0F3F8", "#F7F7F7",
  "#FEE090", "#FDAE61", "#F46D43",
  "#D73027", "#A50026"
))

# batch plotting
for (base_name in processed_clusters) {
  # Remove '_control' suffix if present
  base_name_clean <- sub("_control$", "", base_name)
  
  combined_genes <- character()
  
  for (condition in c("6hpi", "24hpi", "72hpi", "96hpi")) {
    full_cluster_name <- paste(base_name, condition, sep = "_vs_")
    if (full_cluster_name %in% names(cluster_gene_list)) {
      combined_genes <- unique(c(combined_genes, cluster_gene_list[[full_cluster_name]]))
    }
  }
  
  # Save unified gene list
  write.table(combined_genes,
              file = file.path(unified_dir, paste0(base_name_clean, "_unified_genes.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Filter for genes that exist in expression data
  valid_genes <- combined_genes[combined_genes %in% rownames(expr_mat)]
  if (length(valid_genes) == 0) {
    warning(paste("No valid genes found in expression matrix for:", base_name_clean))
    next
  }
  
  
  
  # Generate heatmap
  heatmap_plot <- dittoHeatmap(
    object = cluster_aggregate_exp,
    genes = valid_genes,
    annot.by = "condition_new",
    order.by = "condition_new",
    heatmap.colors = pd_palette(256),
    main = paste0("Unified DE genes heatmap: ", base_name_clean, " (", length(valid_genes), " genes)"),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    scale = "row",
    show_colnames = TRUE
  )
  
  
  # Save heatmap to PDF
  heatmap_file <- paste0(base_name_clean, "_unified_genes_heatmap.pdf")
  pdf(file.path(heatmap_dir, heatmap_file), width = 12, height = 10)
  print(heatmap_plot)
  dev.off()
}

# ------------------------------------------------------------------
# End of Figure S6 script
# ------------------------------------------------------------------
