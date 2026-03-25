### Seurat processing

This directory contains the main Seurat workflow used to generate the final integrated object for downstream analysis and figure generation.

#### Overview

Starting from the eight sublibrary Seurat objects generated in the previous step, this stage:

1. **Loads and merges all sublibrary Seurat objects**
2. **Adds sample metadata based on well identity**
3. **Applies nucleus-level filtering**
4. **Removes ribosomal RNA genes from the expression matrix**
5. **Performs normalisation, variable feature selection, scaling, PCA, Harmony integration, clustering, and UMAP**
6. **Saves the final processed Seurat object**

The main output from this step is the integrated Seurat object used for downstream analyses.

#### Main script

Files in this directory:

- `01_process_seurat_objects.R`  
  Loads the eight per-sublibrary Seurat objects, merges them, adds metadata, performs filtering and dimensionality reduction, and saves the final processed object.
- `02_ribo_rRNA_genes` list of rRNA genes in M gigas

---

#### Summary of workflow

##### 1. Load and merge Seurat objects

The eight sublibrary Seurat objects (`seu_1.rds` to `seu_8.rds`) are loaded and merged into a single object. The layers are then joined so that all nuclei are stored in a single counts layer for downstream processing.

##### 2. Add sample metadata

Sample labels are assigned from the original Parse well identities (`orig.ident`) and stored in the metadata as a new `sample` column. Sample levels are then ordered for consistent downstream plotting and analysis.

##### 3. Filter nuclei

Quality control filtering is applied in three stages:

- remove nuclei with high mitochondrial content (`percent.mt >= 5`)
- retain nuclei with between 200 and 3000 detected genes
- remove ribosomal RNA genes from the matrix

A custom mitochondrial gene list is used, including both standard mitochondrial genes and additional `MZ`-prefixed genes detected in the dataset.

##### 4. Normalisation and dimensionality reduction

The filtered object is processed using the following standard Seurat workflow:

- `NormalizeData()` with LogNormalize and scale factor 10,000
- `FindVariableFeatures()` using 3000 variable genes
- `ScaleData()`
- `RunPCA()`
- `RunHarmony()` using `sample` as the grouping variable
- `FindNeighbors()` using Harmony dimensions 1:18
- `FindClusters()` at resolution 0.6
- `RunUMAP()` using Harmony dimensions 1:18

##### 5. Save final object

The final processed object is saved as an `.rds` file for downstream annotation, marker analysis, differential expression, viral transcript analysis, and figure generation.

---

#### Input files

Typical inputs for this step are:

- `seu_1.rds` to `seu_8.rds`
- `02_ribo_rRNA_genes`

#### Output files

Typical outputs from this step are:

- final processed Seurat object  
  e.g. `seu_obj_filt_umap_18d_6r_3kRes_2026.rds`

#### Notes

- This script performs the main merge and processing workflow for all sublibraries.
- Harmony integration is performed using sample identity as the batch variable.
- The filtered and integrated Seurat object generated here is the main object used in all subsequent analyses.

#### Downstream use

The final processed Seurat object is used for cluster annotation, marker gene analysis, differential expression testing, viral transcript summaries, and manuscript figure generation.
