# Figure 2: UMAP of the integrated Pacific oyster single-nucleus atlas

#.........................................................................................



library(Seurat)
library(tidyverse)

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

# Plot UMAP; Figure 2
p <- DimPlot(
  seu_obj,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5,
  repel = TRUE,
  label.box = FALSE,
  sizes.highlight = TRUE,
  label.size = 6
) +
  NoLegend() +
  ggtitle("Figure 2: A single-nucleus transcriptomic atlas of the\nPacific oyster following OsHV-1 infection") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

# Save as PDF
ggsave(
  filename = "Figure_2_UMAP.pdf",
  plot = p,
  device = cairo_pdf,
  width = 10,
  height = 8
)
