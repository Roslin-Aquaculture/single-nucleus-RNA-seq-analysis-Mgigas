## Single-nucleus RNA-seq analysis in *Magallana/Crassostrea gigas*

This repository contains the analysis pipeline for single-nucleus RNA sequencing (snRNA-seq) data from the Pacific oyster (*Magallana/Crassostrea gigas*), developed to characterise cell-type diversity and transcriptional responses to OsHV-1 infection.

---

### Overview

This workflow enables:

- Construction of a whole-organism cell atlas  
- Identification and annotation of cell types  
- Cluster-specific differential gene expression analysis  
- Investigation of temporal responses to viral infection  
- Targeted analysis of immune-related genes  

The pipeline is designed for non-model organisms with limited gene annotation.

---

### ⚙️ Workflow

**01_preprocessing/**  
- Quality control and filtering of nuclei  
- Preparation of input data  

**02_reference_preparation/**  
- Reference genome/transcriptome setup  
- Annotation formatting  

**03_mapping_quantification/**  
- Read alignment
- Generation of count matrices  

**04_seurat_object/**  
- Creation of Seurat objects  
- Initial data structuring  

**05_seurat_processing/**  
- Normalisation and scaling  
- PCA, UMAP, and Harmony integration  
- Clustering and cell type annotation  

**06_figures/**  
- **Main_figures/**: scripts for main manuscript figures  
- **Supplementary_figures/**: scripts for supplementary analyses  
- **supporting_files/**: curated gene lists and annotation inputs  

---
