# =========================================================
# Section 9: Cluster 1 supplementary figure and table
# Environment: SEURAT session
# =========================================================
# Builds the two Cluster 1 identity outputs from results already
# generated upstream:
#
#   Supplementary Figure S_  dot plot of the 12 IAP paralogues identified
#                            as cluster markers, across all 18 clusters,
#                            control nuclei only
#   Supplementary Table S_   the 59 GO/KEGG terms enriched among Cluster 1
#                            markers and enriched in no other cluster
#
# Inputs:
#   07a output — cluster1_identity/markers_control_only.tsv
#   07b output — cluster1_identity/cluster1_discriminating_terms.tsv
#                cluster1_identity/cluster1_identity_tidy.tsv
#   eggNOG     — full_proteome.emapper.annotations (for gene symbols)
#
# The 12 IAP accessions below are hardcoded for reproducibility. They were
# derived in 07b as: genes carrying a BIR domain in the eggNOG annotation
# (n = 47), intersected with the control-only cluster marker sets (n = 12).
# =========================================================

library(Seurat)
library(Matrix)
library(dplyr)
library(readr)

# =========================================================
# CONFIG
# =========================================================
seurat_obj <- readRDS("/home/pdewari/Documents/parse_2025/seurat_2025/seu_obj_umap_18d_6r_3kRes.rds")

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

names(new.cluster.ids) <- levels(seurat_obj)

seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)


###############
# Add simplified condition column
seurat_obj$condition_new <- case_when(
  seurat_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seurat_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seurat_obj$sample == "24-hpiA" ~ "mid?",
  seurat_obj$sample == "24-hpiJ" ~ "24hpi",
  seurat_obj$sample == "72-hpiJ" ~ "72hpi",
  seurat_obj$sample == "96-hpiE" ~ "96hpi"
)

# Set factor order
seurat_obj$condition_new <- factor(
  seurat_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove 24-hpiA ("mid?")
seurat_obj_clean <- subset(seurat_obj, subset = condition_new != "mid?")
rm(seurat_obj)

# Drop unused factor levels
seurat_obj_clean$condition_new <- droplevels(seurat_obj_clean$condition_new)

seurat_obj_clean$cluster_annotation <- as.character(Idents(seurat_obj_clean))

head(seurat_obj_clean)



#
library(Seurat); library(ggplot2)

# the 12 cell-type-restricted IAPs, grouped by the cluster they mark
iap12 <- c("G19384","G19385","G19421","G19423","G19424",   # Cluster 1 tandem arrays
           "G19035","G19414","G20059","G17608",            # gill ciliary
           "G15123",                                        # haemocyte type 1
           "G23947",                                        # vesicular haemocytes
           "G25799")                                        # gill type 2 / mantle epithelial

rn  <- rownames(seurat_obj_clean)
map <- setNames(rn, sub("^(\\S+)\\s.*$", "\\1", rn))

# ---- 8A: identity, control nuclei only ----
ctrl <- subset(seurat_obj_clean, subset = condition_new == "control")

ord <- c("Cluster 0", "Cluster 1",
         "Gill ciliary cells", "Hepatopancreas cells",
         "Gill neuroepithelial cells", "Gill cell type 1",
         "Hyalinocytes", "Haemocyte cell type 1",
         "Mantle cell type 1", "Cluster 9",
         "Vesicular haemocytes", "Immature haemocytes",
         "Macrophage like cells", "Adductor muscle cells",
         "Mantle cell type 2", "Mantle epithelial cells",
         "Gill cell type 2", "Small granule cells")

Idents(ctrl) <- factor(ctrl$cluster_annotation, levels = ord)

#Idents(ctrl) <- ctrl$cluster_annotation

p8a <- DotPlot(ctrl, features = unname(map[iap12])) +
  coord_flip() +
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(x = NULL, y = NULL)



ggsave("/home/pdewari/eggnog/results/plots/cluster1_identity/Figure8A_IAP_across_clusters.pdf", p8a, width = 9, height = 5)

############################################################################

# table suppl for clsuter 1 identity

id_dir <- "/home/pdewari/eggnog/results/plots/cluster1_identity"
disc <- read_tsv(file.path(id_dir, "cluster1_discriminating_terms.tsv"), show_col_types = FALSE)
tidy_id <- read_tsv(file.path(id_dir, "cluster1_identity_tidy.tsv"), show_col_types = FALSE)

eggnog_full <- read_tsv(
  "/home/pdewari/eggnog/results/full_proteome_20260820_124651/full_proteome.emapper.annotations",
  comment = "##", show_col_types = FALSE) %>%
  mutate(gene_id = sub("\\..*$", "", sub("^transcript:", "", `#query`)))

ann <- eggnog_full %>%
  dplyr::select(gene_id, Preferred_name, Description, PFAMs, COG_category) %>%
  distinct(gene_id, .keep_all = TRUE)



gene_names <- function(ids) {
  v <- strsplit(ids, "/")[[1]]
  nm <- ann$Preferred_name[match(v, ann$gene_id)]
  paste(ifelse(is.na(nm) | nm == "-", v, paste0(v, " (", nm, ")")), collapse = ", ")
}

tab_full <- disc %>%
  filter(n_other_clusters_sharing == 0) %>%
  left_join(tidy_id %>%
              filter(analysis == "one_vs_all", cluster == "Cluster_1") %>%
              dplyr::select(ID, geneID) %>% distinct(),
            by = "ID") %>%
  rowwise() %>%
  mutate(Genes_annotated = gene_names(geneID)) %>%
  ungroup() %>%
  arrange(ontology, p.adjust) %>%
  transmute(Ontology        = sub("GO_", "GO ", ontology),
            Term_ID         = ID,
            Term            = Description,
            N_genes         = Count,
            Fold_enrichment = round(fold_enrichment, 1),
            P_adjusted      = signif(p.adjust, 2),
            Genes           = Genes_annotated)

write_tsv(tab_full, file.path(id_dir, "TableS_cluster1_terms.tsv"))
nrow(tab_full)

nrow(tab_full)                      # 59
table(tab_full$Ontology)            # GO BP 51ish, GO MF 9ish, KEGG 1
sum(grepl("\\(", tab_full$Genes))   # how many rows got at least one gene symbol
