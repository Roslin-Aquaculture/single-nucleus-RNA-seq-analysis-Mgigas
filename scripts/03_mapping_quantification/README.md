### Mapping and quantification

This directory contains scripts and notes used to generate barcode whitelist files and perform read mapping and quantification for the Pacific oyster OsHV-1 single-nucleus RNA-seq analysis.

#### Overview

In this project, mapping and quantification followed a **two-step workflow**:

1. **Run Parse split-pipe on each sublibrary to generate barcode whitelist information**
2. **Run STARsolo using the Parse-derived whitelist files and the virus-inclusive oyster + OsHV-1 reference**

Both steps are required in the final pipeline.

Parse split-pipe was used to obtain barcode whitelist information compatible with the Parse library structure. STARsolo was then used as the main aligner and quantification tool, while always mapping against the **composite Pacific oyster + OsHV-1 genome and annotation**.

#### Background

The libraries in this project were generated using Parse Biosciences chemistry, so barcode structure and whitelist information first needed to be recovered in a Parse-compatible way. In the original workflow, each sublibrary was therefore processed separately with Parse split-pipe.

The resulting Parse output, particularly the barcode information in `process/barcode_data.csv`, was then used to prepare the whitelist files required for STARsolo. STARsolo was subsequently used for alignment and quantification against the final virus-inclusive oyster + OsHV-1 reference, allowing more flexible control over mapping while retaining Parse-compatible barcode handling.

#### Suggested directory contents

Suggested files for this directory may include:

- `README.md`  
  Instructions for whitelist generation, read mapping, and count matrix generation.

- `01_run_parse_for_whitelists.sh`  
  Run Parse split-pipe on each sublibrary to generate barcode information.

- `02_extract_whitelists.sh`  
  Extract barcode whitelist files from Parse output, typically from `process/barcode_data.csv`.

- `03_run_starsolo.sh`  
  Run STARsolo on each sublibrary using the virus-inclusive oyster + OsHV-1 reference and Parse-derived whitelist files.

- `04_collate_starsolo_outputs.sh`  
  Collect Seurat-ready output files from each sublibrary into a single structured output directory.

- `05_combine_parse_sublibraries.sh`  
  Optional Parse command to combine sublibrary outputs for record keeping or QC.

#### Workflow

##### 1. Confirm inputs before starting

Before running this step, ensure that you have:

- preprocessed FASTQ files from `01_preprocessing`
- the final **combined oyster + OsHV-1 reference** from `02_reference_preparation`
- the correct Parse chemistry and kit settings
- sufficient compute resources for both Parse and STARsolo runs, see below
- For processing a single sublibrary:
  - 100M reads or less =  64GB memory, 8 CPU
  - 100M-500M reads = 128GB memory, 116 CPU
  - 500M-1B  reads = 256GB memory, 24-32 CPU

##### 2. Run Parse split-pipe on each sublibrary

The first stage is to process each sublibrary with Parse split-pipe so that barcode whitelist information can be recovered from the Parse output.

###### Example Parse command for one sublibrary

```bash
#$ -V -cwd
#$ -l h_rt=24:10:00
#$ -l h_vmem=20G
#$ -pe sharedmem 18

module load anaconda/2024.02
conda activate spipe

split-pipe \
  --mode all \
  --chemistry v2 \
  --kit WT_mega \
  --kit_score_skip \
  --genome_dir /path/to/virus_inclusive_parse_reference \
  --fq1 /path/to/A1_group_oyster_R1.fastq.gz \
  --fq2 /path/to/A1_group_oyster_R2.fastq.gz \
  --output_dir /path/to/a1_results \
  --sample Ambre1 C9-C10 \
  --sample Ambre2 C11-C12 \
  --sample Ambre3 D1-D2 \
  --sample Ambre4 D3-D4 \
  --sample Ambre5 D5-D6 \
  --sample Ambre6 D7-D8 \
  --sample Ambre7 D9-D10 \
  --sample Ambre8 D11-D12
```

Run this separately for each sublibrary.

Combine all the sub-library results into one.

```bash
#$ -V -cwd
#$ -l h_rt=02:10:00 ###HH:MM:SS
#$ -l h_vmem=20G
#$ -pe sharedmem 8

module load anaconda/2024.02
conda activate spipe

find "$(pwd)" -maxdepth 1 -type d -name 'a[12345678]_*' | sort > sublibs.lis

split-pipe \
    --mode comb \
    --sublib_list sublibs.lis \
    --output_dir ./
```

##### 3. Extract whitelist files from Parse output

After Parse sub-library results are combined, obtain the whitelist files from the Parse output directory, specifically from:

```text
combine output > process > barcode_data.csv
```

The whitelist preparation follows this logic:

- **barcode 1** corresponds to sample identity or well-level sample assignment, copy barcode 1 int owhitelist1.txt
- **barcode 2 and barcode 3** and the same files and correspond to the cellular barcode structure; copy barcode 2 into whitelist2.txt
- only wells that actually contain samples should be included
- do not include all 96 wells unless the entire plate was used


These whitelist files are then used as input for STARsolo.

##### 4. Prepare the STAR index using the virus-inclusive oyster + OsHV-1 reference

STARsolo mapping was always run against the final composite oyster + OsHV-1 reference.

Before mapping, generate a STAR index from the virus-inclusive FASTA and GTF files.

```bash
STAR --runMode genomeGenerate --runThreadN 4 \
  --genomeDir /path/to/cgigas_star_oshv \
  --genomeFastaFiles Crassostrea_gigas_OsHV1_combined.fa \
  --sjdbGTFfile Crassostrea_gigas_OsHV1_combined.gtf
```

This step only needs to be run once per reference build.

##### 5. Run STARsolo on each sublibrary using Parse-derived whitelist files

Once the whitelist files are ready (step 3 above), run STARsolo separately for each sublibrary.

The original STARsolo workflow used:

- gzipped paired FASTQ files
- the virus-inclusive oyster + OsHV-1 genome index
- whitelist files derived from Parse output

###### Example STARsolo command
Note: Do not create temp dir beforehand!

```bash
#!/bin/bash

#$ -V -cwd
#$ -l rl9=false
#$ -l h_rt=16:00:00
#$ -l h_vmem=12G
#$ -pe sharedmem 12

. /etc/profile.d/modules.sh
module load roslin/star/2.7.10b

read1="/path/to/A3_group_oyster_R1.fastq.gz"
read2="/path/to/A3_group_oyster_R2.fastq.gz"
sample="3"
temp="temp"$sample
out="a"$sample"_results"

STAR \
  --runThreadN 32 \
  --genomeDir /path/to/cgigas_star_oshv \
  --readFilesCommand zcat \
  --readFilesIn $read1 $read2 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --runDirPerm All_RWX \
  --soloType CB_UMI_Complex \
  --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
  --soloUMIposition 0_0_0_9 \
  --soloCBwhitelist /path/to/whitelist2.txt /path/to/whitelist2.txt /path/to/whitelist1.txt \
  --soloCBmatchWLtype EditDist_2 \
  --soloBarcodeReadLength 0 \
  --soloFeatures GeneFull Gene \
  --soloCellReadStats Standard \
  --outTmpDir /path/to/$temp \
  --limitBAMsortRAM 230854492160 \
  --outFileNamePrefix /path/to/$out \
  --soloMultiMappers Rescue EM Uniform \
  --outFilterScoreMinOverLread 0.50 \
  --outFilterMatchNminOverLread 0.50
```

###### Important notes

- Use the **virus-inclusive composite oyster + OsHV-1 STAR index** for all STARsolo mapping.
- Keep the same barcode layout and whitelist logic across all sublibraries.
- Run each sublibrary separately and keep output directories clearly named.
- Avoid creating inconsistent temporary output directories or changing barcode settings between runs.

##### 6. Collect STARsolo outputs for downstream analysis

For downstream Seurat analysis, the main files needed from each STARsolo run are:

- `matrix.mtx`
- `features.tsv` or `genes.tsv`
- `barcodes.tsv`

The STARsolo workflow generates multiple matrix outputs, including:

- `matrix.mtx`
- `UniqueAndMult-EM.mtx`
- `UniqueAndMult-Rescue.mtx`
- `UniqueAndMult-Uniform.mtx`


A collation script should gather the relevant output files from all sublibraries into a single results directory, with one subfolder per sublibrary or one clearly organised combined directory.
Copy and rename UniqueAndMult-Rescue.mtx to matrix.mtx; rename features.tsv to genes.tsv.


#### Input files

Typical inputs for this step are:

- split FASTQ files from preprocessing
- Parse-compatible virus-inclusive reference
- STAR-compatible virus-inclusive reference
- Parse-derived barcode whitelist files

#### Output files

This step produces:

- Parse output directories for each sublibrary
- barcode whitelist files derived from Parse output
- STARsolo-aligned BAM files
- per-sublibrary count matrices
- barcode and feature tables
- STARsolo count matrix variants
- collated matrices ready for Seurat import

#### Notes

- Both Parse and STARsolo are part of the final workflow and should be documented together.
- Parse is used here to recover barcode whitelist information, not as the final quantification source.
- STARsolo is the main aligner and quantification tool for downstream analysis.
- Always use the **virus-inclusive oyster + OsHV-1 genome and annotation** for both Parse reference preparation and STARsolo mapping.
- Keep sublibrary names, output directories, and sample labels consistent across all jobs.
- Retain Parse and STARsolo logs, as these are useful for troubleshooting barcode, memory, and mapping issues.


#### Downstream use

The STARsolo matrices generated here are used to construct Seurat objects, perform quality control, quantify host and viral transcript abundance, and carry out all downstream clustering and differential expression analyses.
