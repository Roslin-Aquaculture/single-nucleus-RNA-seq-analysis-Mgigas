# Figure 5: heatmap of viral transcript expression across infection stages

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(Seurat.utils)
library(dittoSeq)

############################
# Load processed Seurat object
############################

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object generated previously in scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds")


############################
# Rename clusters with final cell-type annotations
############################

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

# Add simplified infection-stage labels for downstream comparisons
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample == "24-hpiA" ~ "mid?",
  seu_obj$sample == "24-hpiJ" ~ "24hpi",
  seu_obj$sample == "72-hpiJ" ~ "72hpi",
  seu_obj$sample == "96-hpiE" ~ "96hpi"
)

# Set factor order for plotting and aggregation
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove 24-hpiA ("mid?") from downstream analysis
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")
rm(seu_obj)

# Drop unused factor levels after subsetting
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

head(seu_obj_clean)

############################
# Generate aggregated viral transcript heatmap
############################

# diverging colour palette for heatmap display
pd_palette <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#ABD9E9", "#E0F3F8", "#F7F7F7",
  "#FEE090", "#FDAE61", "#F46D43",
  "#D73027", "#A50026"
))

# Define infection-stage order for aggregation and plotting
desired_order <- c("control", "6hpi", "24hpi", "72hpi", "96hpi")

# Reset condition factor with the desired order
seu_obj_clean$condition_new <- factor(seu_obj_clean$condition_new, levels = desired_order)

# Aggregate expression by infection stage
cluster_aggregate_exp <- AggregateExpression(
  seu_obj_clean,
  return.seurat = TRUE,
  group.by = "condition_new"
)

# Remove the "g" prefix that Seurat may add to aggregated column names
colnames(cluster_aggregate_exp) <- gsub("^g", "", colnames(cluster_aggregate_exp))

# Re-apply condition factor ordering to the aggregated object
cluster_aggregate_exp$condition_new <- factor(
  cluster_aggregate_exp$condition_new,
  levels = c("control", "6hpi", "24hpi", "72hpi", "96hpi")
)

save_dir <- "./"

# Consistent per-gene height estimate for optional figure resizing
gene_height <- 0.25  # in inches per gene (adjust to taste)

# Optional: compute total figure height based on number of viral genes
# plot_height_viral <- length(viral$gene) * gene_height

# Extract viral genes from the Seurat object
viral <- rownames(seu_obj_clean) %>%
  as.data.frame() %>%
  filter(str_detect(., "^ORF")) %>%
  filter(!str_detect(., "\\.")) %>%
  rename(gene = ".")

############################
# Save heatmap as Cairo PDF with Arial font
############################

cairo_pdf(
  filename = paste0(save_dir, "Figure5_OsHV1_transcripts_dittoHeatmap.pdf"),
  width = 8,
  #height = plot_height_viral,
  family = "Arial"
)

viral_heatmap <- dittoHeatmap(
  object = cluster_aggregate_exp,
  genes = viral$gene,
  order.by = "condition_new",
  main = "Figure 5. OsHV-1 transcripts are expressed in the 72- and 96-hpi samples",
  heatmap.colors = pd_palette(256),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  scale = "row",
  show_colnames = TRUE,
  fontsize_row = 6,
  cellwidth = 16,
  cellheight = 8
) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(family = "Arial")
  )

dev.off()

# ------------------------------------------------------------------
# End of Figure 5 script
# ------------------------------------------------------------------
