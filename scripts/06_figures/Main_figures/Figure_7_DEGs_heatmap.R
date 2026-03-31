# Figure 7: top 15 DE gene heatmaps across infection stages


# This script generates heatmaps for each cluster (top 15 genes each for four conditions, 6, 24, 72, 96 hpi)
# Directory structure at the end of run should be:
# ├── heatmaps_avgexp_top15
#     ├── up_in_control
#     ├── up_in_target



############################### STEP 1: DATA PREPARATION ##############################

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(Matrix)
library(Seurat.utils)
library(pheatmap)

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



############################### STEP 2: PARSE DIRS AND FILES ##############################

############################
# Directories and setup
############################

# Base directory containing DE outputs generated in 06_figures/Main_figures/Figure_6_UpSetPlots.R
#
# Expected structure:
# de_plots_full_ann_30032026/
# ├── 6hpi
# ├── 24hpi
# ├── 72hpi
# └── 96hpi

base_dir <- "de_plots_full_ann_30032026/"

# Original absolute path used during analysis
#base_dir <- "/home/pdewari/Documents/parse_2025/seurat_2025/de_plots_full_ann_30032026/"

time_points <- c("6hpi", "24hpi", "72hpi", "96hpi")

# Detect cluster names from the first timepoint directory
first_tp_dir <- file.path(base_dir, time_points[1])
cluster_dirs <- list.dirs(first_tp_dir, full.names = FALSE, recursive = FALSE)
cluster_names <- unique(gsub("_control_vs_.*", "", cluster_dirs))

############################
# Prepare Seurat object for averaging
############################

# Set the final order of infection stages for average-expression heatmaps
desired_order <- c("control", "6hpi", "24hpi", "72hpi", "96hpi")
seu_obj_clean$condition_new <- factor(seu_obj_clean$condition_new, levels = desired_order)

############################
# Output directories
############################

# Separate output folders for genes up in control versus up in target
heatmap_dir_control <- file.path(base_dir, "heatmaps_avgexp_top15", "up_in_control")
heatmap_dir_target  <- file.path(base_dir, "heatmaps_avgexp_top15", "up_in_target")

dir.create(heatmap_dir_control, showWarnings = FALSE, recursive = TRUE)
dir.create(heatmap_dir_target, showWarnings = FALSE, recursive = TRUE)

############################
# Helper function
############################

# Convert directory-style cluster names back to space-separated labels
convert_cluster_name <- function(x) gsub("_", " ", x)

############################
# Colour palette
############################

# Diverging palette used for scaled average-expression heatmaps
pd_palette <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#ABD9E9", "#E0F3F8", "#F7F7F7",
  "#FEE090", "#FDAE61", "#F46D43",
  "#D73027", "#A50026"
))



############################### STEP 3: BATCH CREATE HEATMAPS ##############################


##############################################################
### MAIN LOOP: CLUSTERS × COMPARISON TYPES ###################
##############################################################

for (cluster_name in cluster_names) {
  cat("\n==============================\n")
  cat("Processing cluster:", cluster_name, "\n")
  cat("==============================\n")
  
  cluster_prefix <- paste0(cluster_name, "_control_vs_")
  cluster_actual <- convert_cluster_name(cluster_name)
  
  for (comparison in c("up_in_control", "up_in_target")) {
    
    cat("→ Comparison:", comparison, "\n")
    
    df_cluster <- data.frame(gene = character(), stage = character())
    
    # ---------------------------
    # Collect top 15 DE genes from each timepoint
    # ---------------------------
    for (tp in time_points) {
      file_path <- file.path(
        base_dir,
        tp,
        paste0(cluster_prefix, tp),
        paste0(cluster_prefix, tp, "_de_genes_", comparison, ".txt")
      )
      
      if (file.exists(file_path) && file.info(file_path)$size > 0) {
        genes <- read.table(file_path, header = FALSE, nrows = 15, sep = "\t")
        
        if (nrow(genes) > 0) {
          df_cluster <- rbind(df_cluster, data.frame(gene = genes[, 1], stage = tp))
        }
      }
    }
    
    if (nrow(df_cluster) == 0) {
      warning(paste("No DE genes for", cluster_name, comparison))
      next
    }
    
    # Keep unique genes across timepoints
    genes_to_plot <- unique(df_cluster$gene)
    
    # -----------------------------------------
    # Subset Seurat object to this cluster across all desired stages
    # -----------------------------------------
    idents_expected <- paste0(cluster_actual, "_", desired_order)
    idents_present <- idents_expected[idents_expected %in% Idents(seu_obj_clean)]
    
    if (length(idents_present) == 0) {
      cat("⚠ No matching identities for", cluster_name, "- skipping.\n")
      next
    }
    
    seu_sub <- subset(seu_obj_clean, idents = idents_present)
    
    # -----------------------------------------
    # Compute average expression by condition
    # -----------------------------------------
    avg_exp <- AverageExpression(
      seu_sub,
      group.by = "condition_new",
      assays = "RNA",
      slot = "data"
    )$RNA
    
    genes_present <- intersect(genes_to_plot, rownames(avg_exp))
    
    if (length(genes_present) == 0) {
      warning(paste("No DE genes found in object for", cluster_name))
      next
    }
    
    mat <- avg_exp[genes_present, , drop = FALSE]
    
    # Row-wise Z-score scaling for heatmap display
    mat_scaled <- t(scale(t(mat)))
    
    # -----------------------------------------
    # Save heatmap PDF
    # -----------------------------------------
    output_dir <- ifelse(comparison == "up_in_control", heatmap_dir_control, heatmap_dir_target)
    
    output_file <- file.path(
      output_dir,
      paste0(cluster_name, "_top15_", comparison, "_avgexp_heatmap.pdf")
    )
    
    pdf(output_file, width = 6, height = 0.25 * nrow(mat_scaled) + 2)
    
    p <- pheatmap(
      mat_scaled,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = pd_palette(256),
      main = paste(cluster_actual, "-", gsub("_", " ", comparison)),
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
    
    cat("✓ Saved heatmap:", output_file, "\n")
    
    # -----------------------------------------
    # Save clustered row-order gene list
    # -----------------------------------------
    # Genes are written out in the order returned by hierarchical clustering
    ordered_genes <- rownames(mat_scaled)[p$tree_row$order]
    
    ordered_gene_file <- file.path(
      output_dir,
      paste0(cluster_name, "_top15_", comparison, "_genes_clustered_order.txt")
    )
    
    write.table(
      ordered_genes,
      file = ordered_gene_file,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    cat("✓ Saved ordered gene list:", ordered_gene_file, "\n")
  }
}

# ------------------------------------------------------------------
# End of Figure 7 script
# ------------------------------------------------------------------
