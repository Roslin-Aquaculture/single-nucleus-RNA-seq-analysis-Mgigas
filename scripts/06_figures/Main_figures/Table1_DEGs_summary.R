# ------------------------------------------------------------------
# Summarise differential expression results across infection stages
# Creates a table with one row per cluster and the following columns:
#   cluster | 6hpi | 24hpi | 72hpi | 96hpi | total_unique
#
# Stage-specific columns count the number of DE genes in each
# control-vs-infection comparison.
# total_unique counts the number of unique DE genes across all
# four infection stages for each cluster.
# ------------------------------------------------------------------

library(tidyverse)
library(readr)

# ------------------------------------------------------------------
# Set working directory to the parent folder containing the DGE outputs
# ------------------------------------------------------------------
setwd("/home/pdewari/Documents/parse_2025/seurat_2025/de_plots_full_ann_30032026")

# ------------------------------------------------------------------
# Define infection stages included in the DGE analysis
# ------------------------------------------------------------------
timepoints <- c("6hpi", "24hpi", "72hpi", "96hpi")

# ------------------------------------------------------------------
# Define the desired cluster order for the final summary table
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# Detect all clusters automatically from the comparison folder names
# Example folder name:
#   Gill_ciliary_cells_control_vs_6hpi
# This extracts just the cluster part before "_control_vs_"
# ------------------------------------------------------------------
clusters <- unique(unlist(lapply(timepoints, function(tp) {
  subdirs <- list.dirs(tp, recursive = FALSE, full.names = FALSE)
  sub("_control_vs_.*$", "", subdirs)
})))

# ------------------------------------------------------------------
# Helper function to read the DE gene list for one cluster/timepoint
# Returns a unique character vector of gene names
# If the expected file is missing, returns an empty character vector
# ------------------------------------------------------------------
get_genes <- function(tp, cl) {
  file <- file.path(
    tp,
    paste0(cl, "_control_vs_", tp),
    paste0(cl, "_control_vs_", tp, "_de_genes_all.txt")
  )
  
  if (!file.exists(file)) return(character(0))
  
  unique(read_lines(file))
}

# ------------------------------------------------------------------
# Build the summary table
# For each cluster:
#   - count DE genes at each infection stage
#   - count total unique DE genes across all stages
# Cluster names are converted from underscore format to readable labels
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# Apply the predefined cluster order so the output matches the figure
# and manuscript ordering rather than default alphabetical sorting
# ------------------------------------------------------------------
deg_summary <- deg_summary %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  arrange(cluster) %>%
  mutate(cluster = as.character(cluster))

deg_summary

# ------------------------------------------------------------------
# Save the summary table as a TSV file and print it to the console
# ------------------------------------------------------------------
write_tsv(deg_summary, "Table1_summary_DEGs_all_clusters.tsv")
