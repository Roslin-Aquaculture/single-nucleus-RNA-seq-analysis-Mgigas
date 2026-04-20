## Suppl Figure S3 : Heatmap plotting of key markers for each of the 17 transcriptomic clusters

# load packages
library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(scCustomize)
library(dittoSeq)
library(Seurat.utils)
library(pheatmap)



#......................................................................................

##### step 1: load Seurat object and add annotation, i.e. gene names/symbol

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object; this was created previously under scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds") 


# spaces and hyphens replaced with underscores
new.cluster.ids <- c("Cluster_0", "Cluster_1", 
                     "Gill_ciliary_cells", "Hepatopancreas_cells",
                     "Gill_neuroepithelial_cells", "Gill_cell_type_1",
                     "Hyalinocytes", "Haemocyte_cell_type_1",
                     "Mantle_cell_type_1", "Cluster_9", 
                     "Vesicular_haemocytes", "Immature_haemocytes",
                     "Macrophage_like_cells", "Adductor_muscle_cells",
                     "Mantle_cell_type_2", "Mantle_epithelial_cells",
                     "Gill_cell_type_2", "Small_granule_cells")

names(new.cluster.ids) <- levels(seu_obj)

seu_obj <- RenameIdents(seu_obj, new.cluster.ids) # add cluster info to Idents


# Load gene annotation sources
ORSON <- read_csv("01_ORSON_french_group_2024_biorxiv_suppl2.csv") %>%
  select(Gene.ID, Description = Sequence.Description)

cg_science <- read_tsv("02_Cg_gene_names.tsv", col_names = FALSE) %>%
  rename(Gene.ID = X1, Gene.symbol.science = X2)

# Perform a full join to keep all rows
merged_df <- full_join(ORSON, cg_science, by = "Gene.ID")


# Perform a full join and concatenate Gene.symbol.science first
merged_df <- full_join(ORSON, cg_science, by = "Gene.ID") %>%
  mutate(
    gene_symbol = paste(
      coalesce(Gene.symbol.science, ""), 
      coalesce(Description, ""), 
      sep = " : "
    ) %>% trimws()  # Remove extra spaces
  ) %>% 
  mutate(
    gene_symbol = gsub("[^[:alnum:] :]", " ", gene_symbol)  # Keep "-" but replace other special characters
  )

gene_map2 <- merged_df %>%
  select(gene_id = Gene.ID, gene_symbol)



# duplicates are a PAIN and will generate error below
#Error in .which_data(assay, slot, object)[genes, cells.use] : 
#subscript out of bounds
#Called from: as.matrix(.which_data(assay, slot, object)[genes, cells.use])

#Custom suffix for duplicates
gene_map2 <- gene_map2 %>%
  group_by(gene_symbol) %>%
  mutate(
    gene_symbol = if(n() > 1) {
      paste0(gene_symbol, LETTERS[seq_along(gene_symbol)])
    } else {
      gene_symbol
    }
  ) %>%
  ungroup()


# add gene ID to gene symbol so that we can plot later using gene id
gene_map2 <- gene_map2 %>%
  mutate(gene_symbol = paste(gene_id, gene_symbol, sep = " "))

# Filter the mapping to include only genes present in Seurat object
gene_map_filtered <- gene_map2 %>% 
  filter(gene_id %in% rownames(seu_obj))

# Create a named vector for renaming
gene_symbols <- gene_map_filtered$gene_symbol[match(rownames(seu_obj), gene_map_filtered$gene_id)]
rename_vector <- setNames(gene_symbols, rownames(seu_obj))

# Handle missing gene symbols by keeping their original IDs
rename_vector[is.na(rename_vector)] <- names(rename_vector)[is.na(rename_vector)]

# Ensure uniqueness of gene symbols
#rename_vector <- make.unique(rename_vector)

# Ensure uniqueness of gene symbols, use _1, _2 suffixes instead of .1, .2
rename_vector <- make.unique(rename_vector, sep = "_")

# Rename genes in Seurat object
seu_obj <- RenameGenesSeurat(
  obj = seu_obj,
  newnames = rename_vector
)

rownames(seu_obj)

#......................................................................................




##### step 2: Find marker genes for each trasncriptomic cluster

all_markers <- FindAllMarkers(object = seu_obj, min.pct = 0.20, log2fc.threshold = 0.25)

seu_obj

# convert all markers to df
df_all_markers <- as.data.frame(all_markers) %>% 
  filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = "Gene.ID")

#reorder
df_all_markers <- df_all_markers %>%
  relocate(cluster, .after = Gene.ID)


# save as tsv
write_tsv(df_all_markers, "df_all_markers_20Jan26.tsv" )

#......................................................................................




#### step 3: make a list of key marker genes for each transcriptomic cluster
#### these genes are described in Results section of the manuscript

########### marker genes list: gills and mantle *****************************************

# cluster = gill_ciliary_genes
gill_ciliary_cells <- c("G1821 Kif19 : Kinesin like protein",
                        "G26785 Trpv6 1 : ANK REP REGION domain containing protein",
                        "G7179 Trpm8 4 : LSDAT euk domain containing protein",
                        "G31230 Unchar 7516 : Monocarboxylate transporter 14")


# cluster = Gill_neuroepithelial_cells
Gill_neuroepithelial_cells <- c(
  "G13652 Kcnae : Potassium voltage gated channel protein eag",
  "G2482 Kcnab : BTB domain containing protein",
  "G29679 Kcnal 2 : BTB domain containing protein",
  "G3874 Kcnaw 10 : BTB domain containing protein",
  "G3157 Ccg5 2 : Voltage dependent calcium channel gamma 7 subunit",
  "G12270 Pk1l2 3 : Organic cation transporter protein",
  "G2196 Pk1l2 7 : REJ domain containing protein",
  "G2278 Pkdre 4 : Polycystic kidney disease and receptor for egg jelly related protein",
  "G2772 Pkd1 5 : Polycystic kidney disease protein 1 like 2",
  "G3220 Pkd2 4 : Polycystin 2",
  "G30689 Pclo : Tripartite motif containing protein 45",
  "G25486 Dscam 4 : Down syndrome cell adhesion molecule",
  "G25489 Dscam 3 : Down syndrome cell adhesion molecule",
  "G473 Fstl5 : Agrin")


# cluster = Gill_cell_type_1
Gill_cell_type_1 <- c(
  "G32686 Ctp5b : Mucin 5B",
  "G32688 Nrp2 4 : mucin 5AC like isoform X2",
  "G13460 MEIS2 MEIS1 TF : Homeobox domain containing protein",
  "G15847 Ece 2 : Endothelin converting enzyme 1",
  "G6819 Egfl8 2 : Protein tyrosine phosphatase",
  "G24730 Sc6a5 7 : Transporter",
  "G11316 Notc2 4 : Neurogenic locus notch like protein 2")

# cluster = Gill_cell_type_2
Gill_cell_type_2 <- c(
  "G24730 Sc6a5 7 : Transporter",
  "G24590 Sc6a5 3 : Transporter",
  "G9706 Catl 4 : Cathepsin L",
  "G9707 Catl 5 : Cathepsin L",
  "G9704 Catl1 1 : CTSL",
  "G22673 Bgbp 1 : GH16 domain containing protein",
  "G27166 Chs2 8 : Chitin synthase",
  "G26475 Cd109 2 : PUA domain containing protein",
  "G1009 Cr1 : Complement component receptor 1 like protein",
  "G27503 Ca3a1 : VWFA domain containing protein",
  "G24966 Xpp2 : Xaa Pro aminopeptidase 2")

# cluster = Mantle_cell_type_1
Mantle_cell_type_1 <- c(
  "G3104 Chi10 2 : Putative chitinase 3",
  "G3105 Chia 4 : Chitotriosidase 1",
  "G28856 Unchar 8225 : Peritrophin",
  "G25153 Co6a4 4 : Collagen alpha 4 VI  chain",
  "G26943 Matn3 : Collagen alpha 1 XIV  chain",
  "G505 Cad23 3 : Cadherin 23",
  "G3229 Fat4 14 : protocadherin Fat 4 isoform X3",
  "G3231 Ds 2 : Hedgehog protein",
  "G5010 GATA6 1 TF : GATA 4")

# cluster = Mantle_cell_type_2
Mantle_cell_type_2 <- c(
  "G4180 Pif 7 : von Willebrand factor",
  "G9685 Mp : Mantle protein 10",
  "G15658 Notc1 16 : Neurogenic locus Notch protein",
  "G15659 Notc2 9 : Neurogenic locus notch like protein 2",
  "G15660 Unchar 14453 : Neurogenic locus notch like protein 2",
  "G15661 Notc1 15 : Neurogenic locus notch homolog protein 1",
  "G7827 Unchar 2395 : Chitin binding type 4 domain containing protein",
  "G16369 Unchar 13418 : Chitin binding protein",
  "G10716 Chs2 14 : Chitin synthase",
  "G5250 Hex 7 : Beta N acetylhexosaminidase")

# cluster = Mantle_epithelial_cells
Mantle_epithelial_cells <- c(
  "G25288 Cml6 : Putative insoluble matrix shell protein 5 like protein",
  "G12306 Epdr1 8 : Mammalian ependymin related protein 1",
  "G3105 Chia 4 : Chitotriosidase 1",
  "G25681 Unchar 4648 : Teneurin 2",
  "G21247 Sas 2 : Collagen alpha 4 VI  chain",
  "G21248 Matn1 7 : Collagen alpha 4 VI  chain",
  "G21316 Sas 1 : Collagen alpha 4 VI  chain",
  "G21249 Co6a5 8 : Collagen alpha 5 VI  chain",
  "G13467 Co6a4 7 : Collagen alpha 6 VI  chain",
  "G23525 Unchar 6399 : Collagen alpha 6 VI  chain",
  "G8881 Co6a6 1 : Collagen alpha 6 VI  chain",
  "G22023 Unchar 5987 : WISP1",
  "G22691 Pgca 15 : C type lectin domain containing protein",
  "G20091 Plcl 6 : C type lectin domain containing protein")


########### marker genes list: heamocytes ******************************************

#Hyalinocytes
Hyalinocytes <- c(
  "G2457 Bgh3 20 : Transforming growth factor beta induced protein ig h3",
  "G32588 Pzp : Murinoglobulin 2"
)


# haemocyte_type1
haemocyte_type1 <- c(
  "G22986 Tll2 2 : Metalloendopeptidase",
  "G384 Ptpra 14 : Protein tyrosine phosphatase",
  "G21444 Unchar 5724 : Uncharacterized protein")


#Immature_haemocytes
Immature_haemocytes <- c(
  "G26978 Plcl 29 : C type lectin domain containing protein",
  "G26288 Yat7 1 : H lectin domain containing protein",
  "G31903 Isp2 2 : C type lectin domain containing protein",
  "G17173 Pgca 2 : C type lectin domain containing protein",
  "G26289 Yat7 3 : H lectin domain containing protein",
  "G26963 Plcl 34 : C type lectin domain containing protein",
  "G31904 Plcl 75 : C type lectin domain containing protein",
  "G32748 Unchar 15019 : Bactericidal permeability increasing protein lipopolysaccharide binding protein",
  "G24125 Lamp1 1 : Lysosome associated membrane glycoprotein 1",
  "G2853 Rs12 : 40S ribosomal protein S12",
  "G23358 Rla2 : 60S acidic ribosomal protein P2",
  "G5649 Rla1 : 60S acidic ribosomal protein P1",
  "G7023 Rl19 : Ribosomal protein L19",
  "G31071 Rs23 : 40S ribosomal protein S23",
  "G26894 Rl31 : 60S ribosomal protein L31",
  "G1172 Rs2 : 40S ribosomal protein S2",
  "G4972 Rl12 : 60S ribosomal protein L12"
)

#small_granule_cells
small_granule_cells <- c(
  "G12733 Leg4 4 : Galectin",
  "G5864 S39aa 1 : Zinc transporter ZIP10",
  "G4920 Znt2 4 : Zinc transporter 8",
  "G32289 Yl446 8 : PNPLA domain containing protein",
  "G22387 Tgm1 1 : TGc domain containing protein",
  "G12639 Pygm : Alpha 1 4 glucan phosphorylase")

#Vesicular_haemocytes
Vesicular_haemocytes <- c(
  "G28618 Mylk 4 : Sporozoite and liver stage asparagine rich protein",
  "G1202 Unchar 11539 : DBB domain containing protein",
  "G18439 Wasf3 2 : Wiskott Aldrich syndrome protein family member",
  "G28165 Plsp 1 : Myeloperoxidase",
  "G32837 Pat3 1 : Integrin beta",
  "G31976 Ita9 1 : Integrin alpha2 domain containing protein",
  "G32780 Pat3 3 : Integrin beta",
  "G32483 Ita4 3 : Integrin alpha2 domain containing protein",
  "G32776 Itb6 1 : Integrin beta",
  "G25996 Unchar 5029 : Alpha 2 macroglobulin family protein",
  "G9588 Srec 29 : Ig like domain containing protein",
  "G9563 Unchar 3437 : Ig like domain containing protein",
  "G6683 Ceam5 2 : Ig like domain containing protein",
  "G6695 Nect4 1 : Ig like domain containing protein")

#Macrophage_like_cells
Macrophage_like_cells <- c(
  "G6983 Mrc1 7 : Macrophage mannose receptor 1",
  "G28068 Mrc1 2 : Macrophage mannose receptor 1",
  "G2310 Sap 1 : Proactivator polypeptide",
  "G29373 Gst8 1 : Glutathione S transferase sigma class protein",
  "G29966 Klf15 : Krueppel like factor 15",
  "G27980 Dmbt1 10 : Elastin b")

########### marker genes list: other clusters *************************************

#cluster = Adductor_muscle_cells
Adductor_muscle_cells <- c(
  "G25637 Obscn 2 : Muscle M line assembly protein unc 89",
  "G3128 Mys 2 : Myosin heavy chain  striated muscle",
  "G24796 Mysp : Myosin tail 1 domain containing protein",
  "G27652 Cnn3 : Transgelin",
  "G5908 Tnnt : Troponin T  skeletal muscle",
  "G5252 Co1a2 2 : Collagen alpha 2 I  chain",
  "G5254 Co3a1 : Fibrillar collagen NC1 domain containing protein",
  "G25218 Co6a3 1 : Collagen alpha 3 VI  chain",
  "G24458 Co1a2 1 : Fibrillar collagen NC1 domain containing protein",
  "G25215 Coca1 3 : Collagen alpha 3 VI  chain",
  "G4997 Unc22 3 : Titin",
  "G4999 Titin 21 : Titin"
)

#cluster= Hepatopancreas_cells
Hepatopancreas_cells <- c(
  "G3629 Unchar 10692 : Salivary glue protein Sgs 4",
  "G5009 Unchar 12157 : Death domain containing protein",
  "G24428 Co6a5 5 : Programmed cell death protein 5",
  "G20476 Unchar 5930 : Suppressor of tumorigenicity 14 protein homolog",
  "G15473 Slit 3 : Protein slit",
  "G28009 Pkhl1 1 : Fibrocystin L",
  "G26105 Pkhl1 6 : Fibrocystin L",
  "G26108 Pkhl1 3 : Fibrocystin L",
  "G21910 Jag2 : SAP domain containing protein",
  "G12025 Dner 1 : Apple domain containing protein")

#cluster= cluster_1
cluster_1 <- c(
  "G5430 Mlrp1 15 : Metalloendopeptidase",
  "G25111 Uvs2 : Metalloendopeptidase",
  "G6764 Mlrp2 23 : Metalloendopeptidase",
  "G16354 Unchar 14539 : Enterin neuropeptide",
  "G1208 Bcap 5 : DBB domain containing protein",
  "G1206 Bcap 4 : Phosphoinositide 3 kinase adapter protein 1",
  "G25110 Mif 1 : Macrophage migration inhibitory factor",
  "G19423 Diap2 3 : RING type domain containing protein",
  "G19421 Birc3 5 : RING type domain containing protein",
  "G19385 Diap2 6 : RING type domain containing protein",
  "G17327 Birc7 5 : RING type domain containing protein",
  "G19424 Birc3 3 : RING type domain containing protein",
  "G21268 Cahd1 2 : VWFA domain containing protein",
  "G3801 Cahd1 8 : VWFA domain containing protein",
  "G4113 Chs2 16 : Chitin synthase",
  "G21268 Cahd1 2 : VWFA domain containing protein",
  "G3801 Cahd1 8 : VWFA domain containing protein",
  "G7431 Vwde 6 : von Willebrand factor D and EGF domain containing protein",
  "G6266 Unchar 11767 : Hemicentin 1")

#cluster = cluster_9
cluster_9 <- c(
  "G4381 Zan 5 : Zonadhesin",
  "G4383 Unchar 12364 : Zonadhesin",
  "G5291 Rims2 2 : Regulating synaptic membrane exocytosis protein 2",
  "G15649 Kcip4 : Kv channel interacting protein 4",
  "G29367 Cac1m : Voltage dependent L type calcium channel subunit alpha",
  "G415 Kcma1 : BK channel",
  "G3492 Cac1h : Voltage dependent T type calcium channel subunit alpha")


########### combine into one list ************************************************
# we need to split them into three groups for manageable figures

combined_gill_mantle <- unique(
  c(
    gill_ciliary_cells,
    Gill_neuroepithelial_cells,
    Gill_cell_type_1,
    Gill_cell_type_2,
    Mantle_cell_type_1,
    Mantle_cell_type_2,
    Mantle_epithelial_cells
  ))

combined_haemocyte_group <- unique(
  c(
    Hyalinocytes,
    haemocyte_type1,
    Immature_haemocytes,
    small_granule_cells,
    Vesicular_haemocytes,
    Macrophage_like_cells
  ))

combined_other <- unique(
  c(
    Adductor_muscle_cells,
    Hepatopancreas_cells,
    cluster_1,
    cluster_9
  ))

#......................................................................................


#### step 4: compute AverageExpression and prepare matrix

# manual colour palette 
pd_palette <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#ABD9E9", "#E0F3F8", "#F7F7F7",
  "#FEE090", "#FDAE61", "#F46D43",
  "#D73027", "#A50026"))

# compute average expression per cluster
avg_exp <- AverageExpression(
  seu_obj,
  assays = "RNA",
  slot = "data"
)$RNA

rownames(avg_exp)


############## make matrix for each group
# keep gill mantle together
genes_gill_mantle <- intersect(combined_gill_mantle, rownames(avg_exp))
mat_gill_mantle <- avg_exp[genes_gill_mantle, , drop = FALSE]
mat_gill_mantle_scaled <- t(scale(t(mat_gill_mantle)))
# Remove the initial G#### and space
rownames(mat_gill_mantle_scaled) <- sub("^G\\d+\\s+", "", rownames(mat_gill_mantle_scaled))

# haemocyte group
genes_haemocyte <- intersect(combined_haemocyte_group, rownames(avg_exp))
mat_haemocyte <- avg_exp[genes_haemocyte, , drop = FALSE]
mat_haemocyte_scaled <- t(scale(t(mat_haemocyte)))
rownames(mat_haemocyte_scaled) <- sub("^G\\d+\\s+", "", rownames(mat_haemocyte_scaled))

# all remaining clusters 
genes_other <- intersect(combined_other, rownames(avg_exp))
mat_other <- avg_exp[genes_other, , drop = FALSE]
mat_other_scaled <- t(scale(t(mat_other)))
rownames(mat_other_scaled) <- sub("^G\\d+\\s+", "", rownames(mat_other_scaled))


#......................................................................................


#### step 5: now plot each of the three gene lists here

# Set global font family to Arial
par(family = "Arial")

# gill mantle
pdf("key_markers_gill_mantle_heatmap.pdf", width = 12, height = 0.25 * nrow(mat_gill_mantle_scaled) + 2)
p <- pheatmap(
  mat_gill_mantle_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = pd_palette(256),
  main = "average expression of key marker genes for gill and mantle groups",
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  border_color = NA,
  cellwidth = 15,
  cellheight = 10,
  treeheight_row = 15,
  treeheight_col = 0,
  show_colnames = TRUE,
  show_rownames = TRUE,
  legend = TRUE,
  na_col = "grey95")
dev.off()

# heamocyte
pdf("key_markers_haemocytes_heatmap.pdf", width = 12, height = 0.25 * nrow(mat_haemocyte_scaled) + 2)
p <- pheatmap(
  mat_haemocyte_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = pd_palette(256),
  main = "average expression of key markers for haemocyte group",
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  border_color = NA,
  cellwidth = 15,
  cellheight = 10,
  treeheight_row = 15,
  treeheight_col = 0,
  show_colnames = TRUE,
  show_rownames = TRUE,
  legend = TRUE,
  na_col = "grey95")
dev.off()

### other
pdf("key_markers_all_other_heatmap.pdf", width = 12, height = 0.25 * nrow(mat_other_scaled) + 2)
p <- pheatmap(
  mat_other_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = pd_palette(256),
  main = "average expression of key marker genes: others",
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  border_color = NA,
  cellwidth = 15,
  cellheight = 10,
  treeheight_row = 15,
  treeheight_col = 0,
  show_colnames = TRUE,
  show_rownames = TRUE,
  legend = TRUE,
  na_col = "grey95")
dev.off()



# ------------------------------------------------------------------
# End of Figure S3 script
# ------------------------------------------------------------------
