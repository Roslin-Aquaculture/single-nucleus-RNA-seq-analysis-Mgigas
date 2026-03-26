# Figure S2: QC violin plots and summary table for raw and filtered nuclei

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(scCustomize)
library(patchwork)
library(scales)
library(tibble)

############################
# Load Seurat objects and build the raw merged object
############################

# Set working directory; this would be output of 04_seurat_object/02_merge_parse_hexamer_polya_barcodes.R
setwd("/path/to/seurat_objects/")

# Load all sublibrary Seurat objects
seu1 <- read_rds("seu_1.rds")
seu2 <- read_rds("seu_2.rds")
seu3 <- read_rds("seu_3.rds")
seu4 <- read_rds("seu_4.rds")
seu5 <- read_rds("seu_5.rds")
seu6 <- read_rds("seu_6.rds")
seu7 <- read_rds("seu_7.rds")
seu8 <- read_rds("seu_8.rds")

# Merge all eight sublibraries
merged_8 <- merge(
  x = seu1,
  y = list(seu2, seu3, seu4, seu5, seu6, seu7, seu8),
  add.cell.ids = c("seu1", "seu2", "seu3", "seu4", "seu5", "seu6", "seu7", "seu8")
)

# Join layers into a single counts layer for downstream processing
seu_obj <- JoinLayers(merged_8)

############################
# Add sample metadata
############################

# Map Parse well identities to biological sample labels
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

# Set sample order for consistent plotting and summaries
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

############################
# Calculate mitochondrial percentage
############################

# Identify all MZ-prefixed mitochondrial genes present in the matrix
x <- as.data.frame(rownames(seu_obj))

x <- x %>%
  filter(str_detect(`rownames(seu_obj)`, "MZ"))

MZ_genes <- as.character(x$`rownames(seu_obj)`)

# Combine standard mitochondrial genes with MZ-prefixed genes
mt.genes = c(
  "ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L",
  "ND5", "ND6", "COX1", "COX2", "COX3", "CYTB", MZ_genes
) # got the MZ ones from a cluster that looks very weird!

# Keep only mitochondrial genes present in the count matrix
mt.genes.present <- mt.genes[mt.genes %in% rownames(seu_obj)]

# Calculate mitochondrial percentage per nucleus
seu_obj <- PercentageFeatureSet(seu_obj, features = mt.genes.present, col.name = "percent.mt")

head(seu_obj)
tail(seu_obj)

dim(seu_obj)

############################
# Figure S2A: QC violin plots for the raw dataset
############################

# Subsample nuclei for percent.mt plotting because the full dataset is too large
set.seed(123)
cells_use <- sample(colnames(seu_obj), 50000)
seu_obj_p1 <- subset(seu_obj, cells = cells_use)

# Violin plot: mitochondrial percentage
p1 <- VlnPlot(
  seu_obj_p1,
  features = "percent.mt",
  group.by = "sample",
  pt.size = 0.05,
  y.max = 20,
  raster = TRUE
) +
  NoLegend() +
  xlab("") +
  ylab("percent.mt") +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )

# Violin plot: UMI counts
p2 <- VlnPlot(
  seu_obj,
  features = "nCount_RNA",
  group.by = "sample",
  pt.size = 0.02,
  y.max = 40000
) +
  NoLegend() +
  xlab("") +
  ylab("nUMI") +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )

# Violin plot: number of detected genes
p3 <- VlnPlot(
  seu_obj,
  features = "nFeature_RNA",
  group.by = "sample",
  pt.size = 0.02
) +
  NoLegend() +
  xlab("") +
  ylab("nFeatures") +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )

# Combine and save QC violin plots
combined_plot <- p1 / p2 / p3 + plot_layout(ncol = 1, heights = c(1, 1, 1))

ggsave("Figure_S2A_QC_violin.pdf", combined_plot, width = 8, height = 10, device = cairo_pdf)

############################
# Apply filtering to generate the filtered object
############################

# Filter out nuclei with >=5% mitochondrial reads
seu_obj_filt <- subset(seu_obj, subset = percent.mt < 5)

dim(seu_obj_filt)
# [1] 33534 967401  ~18% nuclei removed

# Filter nuclei by number of detected genes
seu_obj_filt <- subset(
  seu_obj_filt,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000
)

dim(seu_obj_filt)
# [1] 33534 23740

# Remove ribosomal RNA genes
# Ribosomal genes were extracted from the oyster GFF using:
# grep 'ribosomal RNA;' Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3 | awk '{print $9}' | awk -F";|:" '{print $2}' | grep -v "MZ*" > ribo_rRNA_genes
# grep '=28S' Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3 | awk '{print $9}' | awk -F";|:" '{print $2}' | grep -v "MZ*" >> ribo_rRNA_genes

ribo_genes <- read_table("ribo_rRNA_genes", col_names = FALSE)
ribo_genes_to_remove <- as.character(ribo_genes$X1)

# Optional additional highly abundant gene removed from the matrix
ribo_genes_to_remove <- c(ribo_genes_to_remove, "G32889")

# remove rRNA genes
seu_obj_filt <- seu_obj_filt[!(rownames(seu_obj_filt) %in% ribo_genes_to_remove), ]

dim(seu_obj_filt)
# [1] 33466 23740

# Check that the expected number of genes was removed
dim(seu_obj)[1] - dim(seu_obj_filt)[1] == length(ribo_genes_to_remove)

############################
# Create QC summary table for raw and filtered data
############################

# Function to summarise QC metrics from a Seurat object
make_qc_summary <- function(seu) {
  sample_levels <- c(
    "Uninfected", "Homogenate",
    "6-hpiA", "6-hpiD",
    "24-hpiA", "24-hpiJ",
    "72-hpiJ", "96-hpiE"
  )
  
  qc_summary <- seu@meta.data %>%
    as_tibble(rownames = "cell_id") %>%
    mutate(sample = factor(sample, levels = sample_levels)) %>%
    group_by(sample) %>%
    summarise(
      Number_of_nuclei = n(),
      Total_UMI_count = sum(nCount_RNA, na.rm = TRUE),
      Mean_UMI_per_nucleus = round(mean(nCount_RNA, na.rm = TRUE), 0),
      Mean_genes_per_nucleus = round(mean(nFeature_RNA, na.rm = TRUE), 0),
      Mean_mitochondrial_percent = round(mean(percent.mt, na.rm = TRUE), 1),
      .groups = "drop"
    )
  
  total_row <- seu@meta.data %>%
    as_tibble() %>%
    summarise(
      sample = "Total",
      Number_of_nuclei = n(),
      Total_UMI_count = sum(nCount_RNA, na.rm = TRUE),
      Mean_UMI_per_nucleus = round(mean(nCount_RNA, na.rm = TRUE), 0),
      Mean_genes_per_nucleus = round(mean(nFeature_RNA, na.rm = TRUE), 0),
      Mean_mitochondrial_percent = round(mean(percent.mt, na.rm = TRUE), 1)
    )
  
  bind_rows(qc_summary, total_row) %>%
    mutate(
      Number_of_nuclei = format(Number_of_nuclei, big.mark = ",", scientific = FALSE),
      Total_UMI_count = format(Total_UMI_count, big.mark = ",", scientific = FALSE),
      Mean_UMI_per_nucleus = format(Mean_UMI_per_nucleus, big.mark = ",", scientific = FALSE),
      Mean_genes_per_nucleus = format(Mean_genes_per_nucleus, big.mark = ",", scientific = FALSE)
    ) %>%
    rename(Sample = sample)
}

# Generate summary tables for raw and filtered objects
qc_raw <- make_qc_summary(seu_obj)
qc_filt <- make_qc_summary(seu_obj_filt)

# Combine raw and filtered values into a single table as raw (filtered)
qc_combined <- qc_raw %>%
  rename(
    Number_of_nuclei_raw = Number_of_nuclei,
    Total_UMI_count_raw = Total_UMI_count,
    Mean_UMI_per_nucleus_raw = Mean_UMI_per_nucleus,
    Mean_genes_per_nucleus_raw = Mean_genes_per_nucleus,
    Mean_mitochondrial_percent_raw = Mean_mitochondrial_percent
  ) %>%
  left_join(
    qc_filt %>%
      rename(
        Number_of_nuclei_filt = Number_of_nuclei,
        Total_UMI_count_filt = Total_UMI_count,
        Mean_UMI_per_nucleus_filt = Mean_UMI_per_nucleus,
        Mean_genes_per_nucleus_filt = Mean_genes_per_nucleus,
        Mean_mitochondrial_percent_filt = Mean_mitochondrial_percent
      ),
    by = "Sample"
  ) %>%
  mutate(
    Number_of_nuclei = paste0(Number_of_nuclei_raw, " (", Number_of_nuclei_filt, ")"),
    Total_UMI_count = paste0(Total_UMI_count_raw, " (", Total_UMI_count_filt, ")"),
    Mean_UMI_per_nucleus = paste0(Mean_UMI_per_nucleus_raw, " (", Mean_UMI_per_nucleus_filt, ")"),
    Mean_genes_per_nucleus = paste0(Mean_genes_per_nucleus_raw, " (", Mean_genes_per_nucleus_filt, ")"),
    Mean_mitochondrial_percent = paste0(
      Mean_mitochondrial_percent_raw, " (", Mean_mitochondrial_percent_filt, ")"
    )
  ) %>%
  select(
    Sample,
    Number_of_nuclei,
    Total_UMI_count,
    Mean_UMI_per_nucleus,
    Mean_genes_per_nucleus,
    Mean_mitochondrial_percent
  )


# Export combined QC summary table
write.table(
  qc_combined,
  file = "Figure_S2B_qc_summary_26032026.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############### ENDS ###################
