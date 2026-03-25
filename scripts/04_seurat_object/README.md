### Seurat object generation

This directory contains scripts and notes used to generate Seurat objects from the STARsolo count matrices.

#### Overview

In this project, Seurat object generation follows the STARsolo mapping step and uses the **virus-inclusive oyster + OsHV-1 reference-based count matrices** generated earlier in the pipeline.

The recommended approach for this project is:

1. **Use the collated STARsolo EM matrices as input**
2. **Merge oligo-dT and random-hexamer capture barcodes before nucleus calling**
3. **Create one Seurat object per sublibrary**
4. **Save those per-sublibrary Seurat objects for downstream merging and analysis**

#### Important note

For this project, the correct approach is to **combine oligo-dT and hexamer capture barcodes before creating Seurat objects**.

Parse libraries use a mixture of oligo-dT and random-hexamer priming. As a result, reads derived from the same nucleus can appear under separate barcode representations if these capture types are not merged first. To avoid splitting the same biological nucleus across multiple barcodes, the oligo-dT and hexamer-derived barcodes should be merged prior to Seurat object creation.


#### Background

This workflow uses STARsolo outputs for all eight sublibraries. From each STARsolo output directory, the files used for Seurat input were:

- `barcodes.tsv`
- `features.tsv`
- `UniqueAndMult-EM.mtx`

To make these compatible with Seurat, the EM matrix was copied and renamed as `matrix.mtx`, and all three files were gzipped. These files were then collated into one directory per sublibrary.

The next step was to merge Parse hexamer and poly(A) capture barcodes using the barcode-merging functions stored in `merge_parse_hexamers_polyA_captures/`. The merged matrix was then used to create one Seurat object for each sublibrary.

#### Directory contents

Suggested files for this directory may include:

- `README.md`  
  Instructions for preparing Seurat input files and creating Seurat objects.

- `01_collate_starsolo_matrices.sh`  
  Copy `barcodes.tsv`, `features.tsv`, and `UniqueAndMult-EM.mtx` from each STARsolo output directory, rename the EM matrix to `matrix.mtx`, and gzip the files.

- `02_merge_parse_hexamer_polya_barcodes.R`  
  Merge oligo-dT and random-hexamer capture barcodes before Seurat object creation. This requires `merge_parse_hexamers_polyA_captures/`, see below.   

- `merge_parse_hexamers_polyA_captures/`  
  Helper files and barcode-definition tables required for barcode merging.

#### Recommended workflow

##### 1. Prepare Seurat-compatible input files from STARsolo output

The first step is to collate the required STARsolo files from each sublibrary output directory. For this project, we have eight sub-libraries a1-a8.    
This collation can be done either manually or using a bash script `01_collate_starsolo_matrices.sh`.  
##### 1a: Manually collate.
From each STARsolo `Solo.out/GeneFull/raw/` directory, use:

- `barcodes.tsv`
- `features.tsv`
- `UniqueAndMult-EM.mtx`

The original workflow used the **EM matrix** for downstream analysis.

Because Seurat expects the count matrix file to be named `matrix.mtx`, rename:

- `UniqueAndMult-EM.mtx` → `matrix.mtx`

Then gzip all three files after copying them into a new sublibrary-specific directory.  

##### 1b: Use bash script to collate.
Run the `01_collate_starsolo_matrices.sh` bash script to collate all relevant output files from the 8 sub-libraries into star_outputs_collated dir.


###### Example expected collated structure

```text
star_outputs_Uni_Mult_EM_collated/
├── a1_subdir/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── a2_subdir/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── a3_subdir/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── a4_subdir/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── a5_subdir/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── a6_subdir/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── a7_subdir/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── a8_subdir/
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

##### 2. Merge barcodes and create seurat objects

This is the key step in this part of the pipeline.

Before creating Seurat objects, merge oligo-dT and random-hexamer capture barcodes so that reads from the same nucleus are assigned together rather than being treated as separate nuclei.


The barcode merge step used the helper files in:

```text
merge_parse_hexamers_polyA_captures/
```

These include the barcode-definition CSV files and the barcode-merging R script required for Parse barcode reconciliation.
To merge the barcodes and create one seurat object for each sub-library, download the helper files `merge_parse_hexamers_polyA_captures/` and run the R script `02_merge_parse_hexamer_polya_barcodes.R`.  
An example run is shown below.  

```
#!/bin/bash

#$ -V -cwd
#$ -l h_rt=20:00:00
#$ -l h_vmem=60G
#$ -pe sharedmem 6

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate seurat5
Rscript 02_merge_parse_hexamer_polya_barcodes.R

```

This will create one Seurat object for each sublibrary and save it as an `.rds` file.


Typical outputs at this stage are:

- `seu_1.rds`
- `seu_2.rds`
- `seu_3.rds`
- `seu_4.rds`
- `seu_5.rds`
- `seu_6.rds`
- `seu_7.rds`
- `seu_8.rds`


#### Input files

Typical inputs for this step are:

- collated STARsolo EM matrices
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`
- barcode merge helper files from `merge_parse_hexamers_polyA_captures/`

#### Output files

This step should produce:

- one merged expression matrix per sublibrary
- one Seurat object per sublibrary
- `.rds` files ready for downstream merging, QC, clustering, and annotation

#### Notes

- We used **STARsolo EM matrix** as the count matrix input for Seurat object generation.
- Rename `UniqueAndMult-EM.mtx` to `matrix.mtx` before loading with Seurat.
- Merge oligo-dT and hexamer capture barcodes **before** nucleus calling and Seurat object creation.
- Keep the same barcode-merging settings across all sublibraries.
- Save one Seurat object per sublibrary rather than merging all sublibraries immediately at this stage.

#### Suggested checks after Seurat object generation

Before moving to the next stage, confirm that:

- all eight sublibraries were collated correctly
- `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` exist for each sublibrary
- barcode merging completed without errors
- one Seurat object was created for each sublibrary
- the number of nuclei is consistent with expectations after barcode merging

#### Downstream use

The per-sublibrary Seurat objects generated here are used in the next stage of the workflow for merging sublibraries, adding metadata, quality control, clustering, annotation, and downstream host-virus expression analysis.
