# Supplementary Table S1.
# Distribution of nuclei expressing viral transcripts across cell subtypes and infection conditions

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(stringr)

############################
# Load Seurat object and prepare metadata for group comparisons
############################

# ------------------------------------------------------------------
# NOTE: Original analysis paths retained below for reproducibility
# setwd("/home/pdewari/Documents/parse_2025/seurat_2025/")
# seu_obj <- read_rds("seu_obj_umap_18d_6r_3kRes.rds")
# ------------------------------------------------------------------

# Set working directory to location of processed Seurat object
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object generated in:
# scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds")

# ------------------------------------------------------------------
# Assign biologically interpretable labels to clusters
# (replaces numeric cluster IDs with final annotated cell types)
# ------------------------------------------------------------------
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

# Map new labels to existing cluster levels
names(new.cluster.ids) <- levels(seu_obj)

# Apply new identities
seu_obj <- RenameIdents(seu_obj, new.cluster.ids)

############################
# Define infection-stage groups and filter samples
############################

# ------------------------------------------------------------------
# Collapse individual samples into broader infection-stage categories
# This simplifies downstream comparisons across timepoints
# ------------------------------------------------------------------
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample == "24-hpiA" ~ "mid?",
  seu_obj$sample == "24-hpiJ" ~ "24hpi",
  seu_obj$sample == "72-hpiJ" ~ "72hpi",
  seu_obj$sample == "96-hpiE" ~ "96hpi"
)

# Explicitly define ordering of infection stages
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# ------------------------------------------------------------------
# Exclude the 24-hpiA sample ("mid?") from analysis
# (treated as intermediate / ambiguous condition)
# ------------------------------------------------------------------
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")
rm(seu_obj)  # free memory

# Drop unused factor levels after filtering
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

############################
# Create combined identity: cell subtype × infection stage
############################

# ------------------------------------------------------------------
# Concatenate cluster identity and condition to generate a unique label
# for each subtype-condition combination (used for summarisation)
# ------------------------------------------------------------------
seu_obj_clean$sample_subtypes <- paste(
  Idents(seu_obj_clean),
  seu_obj_clean$condition_new,
  sep = "_"
)

# Set this combined label as the active identity class
Idents(seu_obj_clean) <- "sample_subtypes"

##########################################################

############################
# Viral transcript detection
############################

# ------------------------------------------------------------------
# Extract raw RNA count matrix
# ------------------------------------------------------------------
counts <- GetAssayData(seu_obj_clean, assay = "RNA", layer = "counts")

# ------------------------------------------------------------------
# Identify viral genes
# Criteria:
#   - Gene names starting with "ORF"
#   - Exclude entries containing "." (likely isoforms/artefacts)
# ------------------------------------------------------------------
viral_genes <- rownames(counts) %>%
  str_subset("^ORF") %>%
  str_subset("\\.", negate = TRUE)

# ------------------------------------------------------------------
# Quantify viral load per cell (sum of viral UMIs)
# ------------------------------------------------------------------
viral_umi_per_cell <- colSums(counts[viral_genes, , drop = FALSE])

# Identify cells with at least one viral UMI
cells_with_viral <- names(viral_umi_per_cell[viral_umi_per_cell > 0])

# Create dataframe of infected cells
df <- tibble(
  cell = cells_with_viral,
  viral_umi = viral_umi_per_cell[cells_with_viral],
  sample_subtype = seu_obj_clean$sample_subtypes[cells_with_viral]
)

# ------------------------------------------------------------------
# Summarise viral-positive cells per subtype-condition
# ------------------------------------------------------------------

# Total number of cells per subtype-condition
total_cells <- tibble(
  sample_subtype = names(table(seu_obj_clean$sample_subtypes)),
  total_cells = as.integer(table(seu_obj_clean$sample_subtypes))
)

# Count viral-positive cells and calculate proportions
df_summary <- df %>%
  count(sample_subtype, name = "cells_with_viral") %>%
  left_join(total_cells, by = "sample_subtype") %>%
  mutate(
    percent_viral = round((cells_with_viral / total_cells) * 100, 1),
    viral_cells = paste0(cells_with_viral, " (", percent_viral, "%)")
  ) %>%
  arrange(desc(percent_viral)) %>%
  select(sample_subtype, total_cells, viral_cells)

# View summary table
df_summary

# ------------------------------------------------------------------
# Export summary table for Supplementary Table S1
# ------------------------------------------------------------------
write_tsv(df_summary, "viral_cells_by_sample_subtype.tsv")

# ------------------------------------------------------------------
# End of script: Supplementary Table S1
# ------------------------------------------------------------------
