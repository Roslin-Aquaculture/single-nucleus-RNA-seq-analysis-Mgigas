# ------------------------------------------------------------------
# Script overview:
#
# This script evaluates whether the number of differentially expressed
# genes (DEGs) identified in each cell cluster is influenced by cluster
# size. It performs three main steps:
#
# 1. Summarises DEG results across infection stages (6, 24, 72 and 96 hpi)
#    for each annotated cell cluster, including the number of stage-specific
#    DEGs and the total number of unique DEGs across all stages.
#
# 2. Calculates the number of cells present in each cluster across infection
#    conditions using the processed Seurat object.
#
# 3. Combines DEG and cell abundance information to calculate DEG density
#    (unique DEGs per 100 cells) and ranks clusters based on this value.
#
# This normalised metric provides an additional assessment of transcriptional
# response magnitude that accounts for differences in cluster size.
# ------------------------------------------------------------------




################### PART 1: DEG SUMMARY ###################

# ------------------------------------------------------------------
# Summarise differential expression results across infection stages
#
# Output table:
# cluster | 6hpi | 24hpi | 72hpi | 96hpi | total_unique
#
# Each infection-stage column represents the number of DE genes
# identified in that cluster for the corresponding control vs
# infection comparison.
#
# total_unique represents the number of unique DE genes detected
# across all infection stages for each cluster.
# ------------------------------------------------------------------

library(tidyverse)
library(readr)

# Set working directory containing DGE output folders
setwd("/home/pdewari/Documents/parse_2025/seurat_2025/de_plots_full_ann_30032026")


# Infection stages included in the DEG summary
timepoints <- c("6hpi", "24hpi", "72hpi", "96hpi")


# Cluster order used for final output
# This matches the annotation order used in figures and manuscript
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


# Identify clusters automatically from DGE folder names
#
# Example folder:
# Gill_ciliary_cells_control_vs_6hpi
#
# The cluster name is extracted from the part before
# "control_vs"
clusters <- unique(unlist(lapply(timepoints, function(tp) {
  subdirs <- list.dirs(tp, recursive = FALSE, full.names = FALSE)
  sub("_control_vs.*$", "", subdirs)
})))


# Function to retrieve DE gene lists for a cluster and timepoint
#
# Returns:
# - unique gene names if file exists
# - empty vector if file is missing
get_genes <- function(tp, cl) {
  
  file <- file.path(
    tp,
    paste0(cl, "_control_vs_", tp),
    paste0(cl, "_control_vs_", tp, "_de_genes_all.txt")
  )
  
  if (!file.exists(file)) return(character(0))
  
  unique(read_lines(file))
}


# Generate DEG summary table
#
# For each cluster:
# - count DE genes at each infection stage
# - calculate total number of unique DE genes across stages
deg_summary <- map_dfr(clusters, function(cl) {
  
  g6  <- get_genes("6hpi", cl)
  g24 <- get_genes("24hpi", cl)
  g72 <- get_genes("72hpi", cl)
  g96 <- get_genes("96hpi", cl)
  
  tibble(
    cluster = gsub("_", " ", cl),
    `6hpi` = length(g6),
    `24hpi` = length(g24),
    `72hpi` = length(g72),
    `96hpi` = length(g96),
    total_unique = length(unique(c(g6, g24, g72, g96)))
  )
})


# Apply predefined cluster order
deg_summary <- deg_summary %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  arrange(cluster) %>%
  mutate(cluster = as.character(cluster))


deg_summary


################### PART 1 END ###################



################### PART 2: CELL COUNT SUMMARY ###################

library(Seurat)
library(tidyverse)
library(Seurat.utils)


# Load processed Seurat object
setwd("/home/pdewari/Documents/parse_2025/seurat_2025/")

seu_obj <- read_rds("seu_obj_umap_18d_6r_3kRes.rds")


# Rename clusters with final cell-type annotations
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



# Create simplified infection-stage groups
#
# 24-hpiA sample is excluded because it was not included in
# downstream comparisons
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample == "24-hpiA" ~ "mid?",
  seu_obj$sample == "24-hpiJ" ~ "24hpi",
  seu_obj$sample == "72-hpiJ" ~ "72hpi",
  seu_obj$sample == "96-hpiE" ~ "96hpi"
)


# Set condition order
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)


# Remove 24-hpiA sample
seu_obj_clean <- subset(
  seu_obj,
  subset = condition_new != "mid?"
)

rm(seu_obj)


# Remove unused factor levels
seu_obj_clean$condition_new <- droplevels(
  seu_obj_clean$condition_new
)


# Count cells per cluster and infection condition
cell_counts <- as.data.frame(
  table(
    Cluster = Idents(seu_obj_clean),
    Condition = seu_obj_clean$condition_new
  )
)


# Convert to wide format and calculate total cells per cluster
cell_counts_wide <- cell_counts %>%
  pivot_wider(
    names_from = Condition,
    values_from = Freq,
    values_fill = 0
  ) %>%
  mutate(
    Total = rowSums(across(where(is.numeric)))
  )


cell_counts_wide


################### PART 2 END ###################



################### PART 3: DEG NORMALISATION BY CELL NUMBER ###################

# Combine DEG summary with cell numbers
#
# Calculate DEG density as number of unique DE genes per 100 cells
# to account for differences in cluster size

deg_per_cell <- cell_counts_wide %>%
  rename(cluster = Cluster) %>%
  left_join(deg_summary, by = "cluster") %>%
  mutate(
    DEGs_per_100_cells = round((total_unique / Total) * 100, 2),
    Rank = rank(-DEGs_per_100_cells, ties.method = "min")
  ) %>%
  select(
    cluster,
    Total,
    total_unique,
    DEGs_per_100_cells,
    Rank
  )


deg_per_cell


################### PART 3 END ###################
