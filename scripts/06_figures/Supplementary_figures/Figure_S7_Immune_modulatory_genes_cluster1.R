# Figure S7: Immune modulatory genes upregulated in Cluster 1 in response to OsHV-1 infection

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(Seurat.utils)
library(pheatmap)

################## Seurat object preparation ################################

#...............................................................
# Original local working path used during analysis
#setwd("/home/pdewari/Documents/parse_2025/seurat_2025/")

# Original object load used during analysis
#seu_obj <- read_rds("seu_obj_umap_18d_6r_3kRes.rds")
#...............................................................

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object created previously in scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds")


# Assign cluster annotations
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

# ------------------------------------------------------------------
# Define infection conditions
# Note: 24-hpiA excluded due to poor viral qPCR signal
# ------------------------------------------------------------------
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample %in% "24-hpiA" ~ "mid?",   # low confidence sample
  seu_obj$sample %in% "24-hpiJ" ~ "24hpi",
  seu_obj$sample %in% "72-hpiJ" ~ "72hpi",
  seu_obj$sample %in% "96-hpiE" ~ "96hpi",
  TRUE ~ "unknown"
)

# Ensure correct temporal ordering
seu_obj$condition_new <- factor(
  seu_obj$condition_new, 
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove uncertain condition
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")
rm(seu_obj)

# Drop unused factor levels
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

# Create combined identity (cluster + condition) for downstream grouping
seu_obj_clean$sample_subtypes <- paste(
  Idents(seu_obj_clean), 
  seu_obj_clean$condition_new, 
  sep = "_"
)
Idents(seu_obj_clean) <- "sample_subtypes"

################## Gene annotation harmonisation #############################

# Combine ORSON and CG annotation to improve gene interpretability
ORSON <- read_csv("01_ORSON_french_group_2024_biorxiv_suppl2.csv") %>% 
  select(Gene.ID, Description = Sequence.Description)

cg_science <- read_tsv("02_Cg_gene_names.tsv", col_names = FALSE) %>% 
  rename(Gene.ID = X1, Gene.symbol.science = X2)

# Merge annotations and construct readable gene labels
merged_df <- full_join(ORSON, cg_science, by = "Gene.ID") %>%
  mutate(
    gene_symbol = paste(
      coalesce(Gene.symbol.science, ""), 
      coalesce(Description, ""), 
      sep = " : "
    ) %>% trimws()
  ) %>%
  mutate(
    gene_symbol = gsub("[^[:alnum:] :]", " ", gene_symbol)  # clean special characters
  )

gene_map2 <- merged_df %>%
  select(gene_id = Gene.ID, gene_symbol)

# Ensure unique gene names to avoid downstream indexing errors
gene_map2 <- gene_map2 %>%
  group_by(gene_symbol) %>%
  mutate(
    gene_symbol = if (n() > 1) {
      paste0(gene_symbol, LETTERS[seq_along(gene_symbol)])
    } else {
      gene_symbol
    }
  ) %>%
  ungroup()

# Retain gene ID in label for traceability
gene_map2 <- gene_map2 %>%
  mutate(gene_symbol = paste(gene_id, gene_symbol, sep = " "))

# Restrict to genes present in Seurat object
gene_map_filtered <- gene_map2 %>% 
  filter(gene_id %in% rownames(seu_obj_clean))

# Create renaming vector
gene_symbols <- gene_map_filtered$gene_symbol[
  match(rownames(seu_obj_clean), gene_map_filtered$gene_id)
]

rename_vector <- setNames(gene_symbols, rownames(seu_obj_clean))
rename_vector[is.na(rename_vector)] <- names(rename_vector)[is.na(rename_vector)]
rename_vector <- make.unique(rename_vector)

# Apply gene renaming
seu_obj_clean <- RenameGenesSeurat(
  obj = seu_obj_clean,
  newnames = rename_vector
)

################## Expression summarisation ##################################

# Custom diverging colour palette for heatmap
pd_palette <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#ABD9E9", "#E0F3F8", "#F7F7F7",
  "#FEE090", "#FDAE61", "#F46D43",
  "#D73027", "#A50026"
))

# Compute average expression per condition
avg_exp <- AverageExpression(
  seu_obj_clean,
  assays = "RNA",
  group.by = "condition_new",
  slot = "data"
)$RNA

################## Immune gene selection #####################################

# Load pre-computed cluster-specific DEG list, generated in scripts/06_figures/Supplementary_figures/Figure_S6_Heatmaps_DEGs.R
cluster1_deg <- read_csv(
  "path_to_file/de_plots_full_ann_30032026/unified_gene_lists_v1/Cluster_1_unified_genes.txt",
  col_names = FALSE
) %>% 
  rename(gene = X1)

# Keyword sets capturing conserved immune functions across timepoints
# (RNA sensing, ubiquitination, apoptosis, signalling, ECM/adhesion)
hpi_6 <- c(
  "Malt1", "Jun", "Crebl2", "Fer", "Abl",
  "Epha2", "Fgfr3 2", "Ced3", "VWFA",
  "fibronectin", "Ig like", "Znfx1",
  "Ddx58", "Dsrad", "Helz2", "Traf3",
  "trim3", "Ubp", "Dzip3", "E3",
  "ubiquitin"
)

hpi24 <- c(
  "Znfx1", "Ddx58", "Dsrad",
  "Helz2", "Traf3", "trim", "Ubp",
  "Dzip3", "E3", "ubiquitin"
)

hpi_72_96 <- c(
  "Birc", "Diap", "Ced3", "Traf3",
  "trim3", "RING", "E3", "Mxd1",
  "Crebl2", "Hmcn1", "Fibrillin",  
  "Cadherin", "Protocadherin"
)

# Combine and deduplicate keywords
cluster1_unique <- unique(c(hpi_6, hpi24, hpi_72_96))

# Pattern matching against gene annotations (case-insensitive)
pattern <- paste0("\\b(", paste(cluster1_unique, collapse = "|"), ")\\b")

cluster1_deg_filtered <- cluster1_deg[
  grepl(pattern, cluster1_deg$gene, ignore.case = TRUE),
]

# Intersect with expressed genes and scale expression (row-wise)
cluster1_deg_genes <- intersect(cluster1_deg_filtered$gene, rownames(avg_exp))

mat_cluster1 <- avg_exp[cluster1_deg_genes, , drop = FALSE]
mat_cluster1_scaled <- t(scale(t(mat_cluster1)))

# Remove gene IDs (G####) for cleaner visualisation
rownames(mat_cluster1_scaled) <- sub("^G\\d+\\s+", "", rownames(mat_cluster1_scaled))

################## Heatmap visualisation #####################################

# Set font for publication consistency
par(family = "Arial")

pdf(
  "Figure S7. Immune modulatory genes upregulated in cluster 1 in response to OsHV-1 infection.pdf",
  width = 10,
  height = 0.25 * nrow(mat_cluster1_scaled) + 2
)

pheatmap(
  mat_cluster1_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = pd_palette(256),
  main = "Figure S7: Expression of immune modulatory genes in cluster 1",
  fontsize_row = 8,
  fontsize_col = 9,
  angle_col = 45,
  border_color = NA,
  cellwidth = 20,
  cellheight = 9,
  treeheight_row = 15,
  treeheight_col = 0,
  show_colnames = TRUE,
  show_rownames = TRUE,
  legend = TRUE,
  na_col = "grey95"
)

dev.off()

# ------------------------------------------------------------------
# End of Figure S7 script
# ------------------------------------------------------------------
