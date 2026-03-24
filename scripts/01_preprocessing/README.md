## Preprocessing

This directory contains scripts used to prepare raw FASTQ files before alignment and downstream single-nucleus RNA-seq analysis.

### Overview

In this project, some sequencing data required preprocessing before mapping. This step has two main purposes:

1. **Concatenate FASTQ files from multiple sequencing lanes** for the same sublibrary.
2. **Split FASTQ files into predefined biological groups/species** before alignment.

This approach is useful when reads from different biological sources are present in the same sequencing run and need to be separated before downstream analysis.

### Directory contents

- `split_fastq_by_group.py`  
  Python [script ](https://www.dropbox.com/scl/fi/z8u9cj4rngoqd9087b5xp/fastq_sep_groups.py?rlkey=ah2v8p1qqz4ji21uazecee58o&e=1&st=clpvhb6q&dl=0) used to split FASTQ files into predefined groups based on well positions or sample layout.

### Recommended workflow

#### 1. Concatenate lane-level FASTQ files

If a sublibrary was sequenced across more than one lane, first concatenate the FASTQ files for that same sublibrary.

**Important:**  
Concatenate lanes belonging to the same sublibrary only. Do **not** merge different sublibraries together.

##### Example

```bash
# Example for sublibrary A1
# for sublibrary A1
cat A1_EKDL240002473-1A_223M7CLT4_L7_1.fq.gz A1_EKDL240002473-1A_223M7CLT4_L8_1.fq.gz > A1_EKDL240002473-1A_223M7CLT4_cat_1.fq.gz
cat A1_EKDL240002473-1A_223M7CLT4_L7_2.fq.gz A1_EKDL240002473-1A_223M7CLT4_L8_2.fq.gz > A1_EKDL240002473-1A_223M7CLT4_cat_2.fq.gz
```

After this step, each sublibrary should have one merged `R1` file and one merged `R2` file.

#### 2. Split FASTQ files by group/species

After concatenation, run the Python script to split the FASTQ files into the relevant biological groups.

This is the main preprocessing step when a sequencing pool contains more than one species or experimental group.

##### Before running the script

Check and update the following:

- input FASTQ file names
- output directory
- chemistry or kit settings
- group names
- well ranges corresponding to each group

##### Example usage

###### Example 1: splitting by multiple groups

```bash
python $SCRIPTPATH \
--chemistry v2 \
--fq1 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_cat_1.fq.gz \
--fq2 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_cat_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12
```

### Input files

Typical inputs for this step are:

- gzip-compressed FASTQ files (`*.fastq.gz` or `*.fq.gz`)
- one `R1` file and one matching `R2` file
- predefined group-to-well assignments

### Output files

This step produces:

- group-specific FASTQ files
- one output set per defined group
- files ready for alignment in the next stage of the pipeline

### Notes

- Input FASTQ files should be gzip-compressed.
- If lane-level FASTQs exist, concatenate them before splitting.
- Keep sublibraries separate throughout preprocessing.
- Double-check group definitions carefully before running the script.
- A small test run or dry run is recommended before processing the full dataset.

### Suggested checks after preprocessing

After splitting, confirm that:

- output files were created for all expected groups
- file sizes are non-zero
- group names are correct
- the total number of reads appears reasonable

### Downstream use

The FASTQ files generated here are used as input for the alignment and quantification steps in the main analysis pipeline.
