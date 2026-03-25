# Seurat processing workflow
# Updated: 24 March 2025
#
# Purpose:
#   1. Load and merge the eight per-sublibrary Seurat objects
#   2. Add sample metadata from Parse well identities
#   3. Filter nuclei by mitochondrial content and gene count
#   4. Remove ribosomal RNA genes from the expression matrix
#   5. Run normalisation, PCA, Harmony integration, clustering, and UMAP
#   6. Save the final processed Seurat object



################## STEP 1: LOAD AND MERGE SEURAT OBJECTS, ADD SAMPLE INFO ##################

# Load packages
library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(scCustomize)

# Set working directory
setwd("/home/pdewari/Documents/parse_2025/seurat_2025/")

# Load per-sublibrary Seurat objects
seu1 <- read_rds("seu_1.rds")
seu2 <- read_rds("seu_2.rds")
seu3 <- read_rds("seu_3.rds")
seu4 <- read_rds("seu_4.rds")
seu5 <- read_rds("seu_5.rds")
seu6 <- read_rds("seu_6.rds")
seu7 <- read_rds("seu_7.rds")
seu8 <- read_rds("seu_8.rds")

# Merge all eight objects
merged_8 <- merge(
  x = seu1,
  y = list(seu2, seu3, seu4, seu5, seu6, seu7, seu8),
  add.cell.ids = c("seu1", "seu2", "seu3", "seu4", "seu5", "seu6", "seu7", "seu8")
)

# Join layers into a single counts layer for downstream processing
seu_obj <- JoinLayers(merged_8)


# Map Parse well identities to sample labels
orig_to_sample <- c(
  "C9" = "Uninfected", "C10" = "Uninfected",
  "C11" = "Homogenate",  "C12" = "Homogenate",
  "D1" = "6-hpiA",       "D2" = "6-hpiA",
  "D3" = "6-hpiD",       "D4" = "6-hpiD",
  "D5" = "24-hpiA",      "D6" = "24-hpiA",
  "D7" = "24-hpiJ",      "D8" = "24-hpiJ",
  "D9" = "72-hpiJ",      "D10" = "72-hpiJ",
  "D11" = "96-hpiE",     "D12" = "96-hpiE"
)

# Add sample labels to metadata
seu_obj@meta.data$sample <- orig_to_sample[seu_obj@meta.data$orig.ident]

# Order sample labels for consistent downstream plotting
seu_obj$sample <- factor(
  seu_obj$sample,
  levels = c(
    "Uninfected", "Homogenate",
    "6-hpiA", "6-hpiD",
    "24-hpiA", "24-hpiJ",
    "72-hpiJ", "96-hpiE"
  )
)

head(seu_obj)

################## STEP 2: DATA FILTERING ##################

dim(seu_obj)
# [1] 33534 1179648

# ---------------------------
# Step 2.1: filter high-mitochondrial nuclei
# ---------------------------

# Identify additional MZ-prefixed mitochondrial genes present in the dataset
x <- as.data.frame(rownames(seu_obj))

x <- x %>%
  filter(str_detect(`rownames(seu_obj)`, "MZ"))

MZ_genes <- as.character(x$`rownames(seu_obj)`)

# Combined mitochondrial gene list
mt.genes <- c(
  "ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L",
  "ND5", "ND6", "COX1", "COX2", "COX3", "CYTB",
  MZ_genes
)

# Keep only mitochondrial genes present in the count matrix
mt.genes.present <- mt.genes[mt.genes %in% rownames(seu_obj)]

# Calculate mitochondrial percentage per nucleus
seu_obj <- PercentageFeatureSet(
  seu_obj,
  features = mt.genes.present,
  col.name = "percent.mt"
)

head(seu_obj)

# Filter out nuclei with >=5% mitochondrial reads
seu_obj_filt <- subset(seu_obj, subset = percent.mt < 5)

dim(seu_obj_filt)
# [1] 33534 967401  ~18% nuclei removed

# ---------------------------
# Step 2.2: filter by number of detected genes
# ---------------------------

seu_obj_filt <- subset(
  seu_obj_filt,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000
)

dim(seu_obj_filt)
# [1] 33534 23740

# ---------------------------
# Step 2.3: remove ribosomal RNA genes
# ---------------------------

# Ribosomal genes were extracted from the oyster GFF using:
# grep 'ribosomal RNA;' Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3 | awk '{print $9}' | awk -F";|:" '{print $2}' | grep -v "MZ*" > ribo_rRNA_genes
# grep '=28S' Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3 | awk '{print $9}' | awk -F";|:" '{print $2}' | grep -v "MZ*" >> ribo_rRNA_genes

ribo_genes <- read_table("ribo_rRNA_genes", col_names = FALSE)
ribo_genes_to_remove <- as.character(ribo_genes$X1)

# Optional additional highly abundant gene removed from the matrix
ribo_genes_to_remove <- c(ribo_genes_to_remove, "G32889")

dim(seu_obj_filt)
# [1] 33534 23740

seu_obj_filt <- seu_obj_filt[!(rownames(seu_obj_filt) %in% ribo_genes_to_remove), ]

dim(seu_obj_filt)
# [1] 33466 23740

# Check that the expected number of genes was removed
dim(seu_obj)[1] - dim(seu_obj_filt)[1] == length(ribo_genes_to_remove)

# Remove unfiltered object to avoid accidental reuse
rm(seu_obj)

# From this point onward, use seu_obj_filt

################## STEP 3: NORMALISATION, PCA, HARMONY, CLUSTERING, UMAP ##################

# Log-normalise counts
seu_obj_filt <- NormalizeData(
  seu_obj_filt,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Identify variable features
seu_obj_filt <- FindVariableFeatures(
  seu_obj_filt,
  selection.method = "vst",
  nfeatures = 3000
)

# Scale all genes
all.genes <- rownames(seu_obj_filt)
seu_obj_filt <- ScaleData(seu_obj_filt, features = all.genes)

# Run PCA on variable features
seu_obj_filt <- RunPCA(
  seu_obj_filt,
  features = VariableFeatures(object = seu_obj_filt)
)

DimPlot(seu_obj_filt, reduction = "pca") + NoLegend()

# Integrate across samples with Harmony
seu_obj_filt <- RunHarmony(
  seu_obj_filt,
  group.by.vars = "sample",
  plot_convergence = TRUE
)

# Construct graph, cluster nuclei, and run UMAP on Harmony dimensions
seu_obj_filt <- FindNeighbors(seu_obj_filt, dims = 1:18, reduction = "harmony")
seu_obj_filt <- FindClusters(seu_obj_filt, resolution = 0.6)
seu_obj_filt <- RunUMAP(seu_obj_filt, dims = 1:18, reduction = "harmony")

# Plot UMAP
DimPlot(
  seu_obj_filt,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5,
  repel = TRUE,
  label.box = TRUE,
  sizes.highlight = TRUE,
  label.size = 6
) +
  NoLegend() +
  ggtitle("Starsolo harmony umap: 18 dim x 6 res x 3k variable features") +
  theme(plot.title = element_text(hjust = 0.5))

# Save final processed Seurat object
saveRDS(seu_obj_filt, file = "seu_obj_filt_umap_18d_6r_3kRes_2026.rds")
