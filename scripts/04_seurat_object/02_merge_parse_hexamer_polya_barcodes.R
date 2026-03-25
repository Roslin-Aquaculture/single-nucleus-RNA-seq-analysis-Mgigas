
setwd("/path/to/starsolo_outputs/")

#load Libraries
library(dplyr)
library(Seurat)
library (ggplot2)

#scripts for merging barcodes and reading in the matrix/feature/barcode outputs from Star
source("merge_parse_hexamers_polyA_captures/merge_parse_hexamers_polyA_captures.R")

mat_1 <-  Read10X(data.dir ="star_outputs_Uni_Mult_EM_collated/a1_subdir/")
merge_1 <- merge_hexamer_polyA_columns(mtx = mat_1, kit = "WT_mega", version = "v2", bc_directory = "merge_parse_hexamers_polyA_captures")
seu_1 <- CreateSeuratObject(counts = merge_1, project = "seu_1")
saveRDS(seu_1, "seu_1.rds")

mat_2 <-  Read10X(data.dir ="star_outputs_Uni_Mult_EM_collated/a2_subdir/")
merge_2 <- merge_hexamer_polyA_columns(mtx = mat_2, kit = "WT_mega", version = "v2", bc_directory = "merge_parse_hexamers_polyA_captures")
seu_2 <- CreateSeuratObject(counts = merge_2, project = "seu_2")
saveRDS(seu_2, "seu_2.rds")

mat_3 <-  Read10X(data.dir ="star_outputs_Uni_Mult_EM_collated/a3_subdir/")
merge_3 <- merge_hexamer_polyA_columns(mtx = mat_3, kit = "WT_mega", version = "v2", bc_directory = "merge_parse_hexamers_polyA_captures")
seu_3 <- CreateSeuratObject(counts = merge_3, project = "seu_3")
saveRDS(seu_3, "seu_3.rds")

mat_4 <-  Read10X(data.dir ="star_outputs_Uni_Mult_EM_collated/a4_subdir/")
merge_4 <- merge_hexamer_polyA_columns(mtx = mat_4, kit = "WT_mega", version = "v2", bc_directory = "merge_parse_hexamers_polyA_captures")
seu_4 <- CreateSeuratObject(counts = merge_4, project = "seu_4")
saveRDS(seu_4, "seu_4.rds")

mat_5 <-  Read10X(data.dir ="star_outputs_Uni_Mult_EM_collated/a5_subdir/")
merge_5 <- merge_hexamer_polyA_columns(mtx = mat_5, kit = "WT_mega", version = "v2", bc_directory = "merge_parse_hexamers_polyA_captures")
seu_5 <- CreateSeuratObject(counts = merge_5, project = "seu_5")
saveRDS(seu_5, "seu_5.rds")

mat_6 <-  Read10X(data.dir ="star_outputs_Uni_Mult_EM_collated/a6_subdir/")
merge_6 <- merge_hexamer_polyA_columns(mtx = mat_6, kit = "WT_mega", version = "v2", bc_directory = "merge_parse_hexamers_polyA_captures")
seu_6 <- CreateSeuratObject(counts = merge_6, project = "seu_6")
saveRDS(seu_6, "seu_6.rds")

mat_7 <-  Read10X(data.dir ="star_outputs_Uni_Mult_EM_collated/a7_subdir/")
merge_7 <- merge_hexamer_polyA_columns(mtx = mat_7, kit = "WT_mega", version = "v2", bc_directory = "merge_parse_hexamers_polyA_captures")
seu_7 <- CreateSeuratObject(counts = merge_7, project = "seu_7")
saveRDS(seu_7, "seu_7.rds")

mat_8 <-  Read10X(data.dir ="star_outputs_Uni_Mult_EM_collated/a8_subdir/")
merge_8 <- merge_hexamer_polyA_columns(mtx = mat_8, kit = "WT_mega", version = "v2", bc_directory = "merge_parse_hexamers_polyA_captures")
seu_8 <- CreateSeuratObject(counts = merge_8, project = "seu_8")
saveRDS(seu_8, "seu_8.rds")
