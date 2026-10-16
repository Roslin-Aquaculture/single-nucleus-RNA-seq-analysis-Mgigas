# Old *C. gigas* gene IDs → current *M. gigas* annotation (Liftoff)

Maps gene IDs from the original Pacific oyster annotation used for our snRNA-seq
alignment (Ensembl Metazoa release 58 / Roslin gene models, `G#####`, on
**cgigas_uk_roslin_v1**) to genes in the current NCBI annotation of the new
reference genome **xbMagGiga1.1** (GCF_963853765.1, annotation RS_2024_06).

The output is a single table (TSV + Excel) with one row per gene: the old gene ID,
its old symbol/description, and the matching current gene with its NCBI GeneID,
symbol, description, biotype and coordinates.

## How it works

```
old FASTA + old GFF3 + new FASTA ──Liftoff──▶ old gene models placed on the new genome
                                                         │
                       new NCBI GFF3 ──── R script ──────┘──▶ old_gene_id → current gene (TSV / XLSX)
```

1. **Liftoff** (`run_liftoff_oyster.sh`) extracts each old gene's sequence, aligns
   it to the new genome with minimap2, and rebuilds the gene model (exons) at its
   new position. Old gene IDs are kept; each lifted gene gets coverage and
   sequence-identity scores. Genes that cannot be placed (< 50% coverage or
   identity) are listed in `unmapped_genes.txt`.
2. **Mapping** (`map_seurat_genes_liftoff.R`) merges the exons of each lifted gene
   and of each current NCBI gene, counts shared exonic bases on the same strand,
   and assigns each old gene to the current gene with the best reciprocal overlap:
   - **high**: reciprocal exonic overlap ≥ 0.5
   - **medium**: reciprocal best match, with ≥ 50% of the lifted exons inside the gene
   - weaker overlaps are left unassigned (the nearest gene is still reported in `best_GeneID`)

   It also flags old genes that share one current gene (`many_to_one`) or span two
   (`spans_2+_genes`), and adds old gene symbols/descriptions from our annotation
   tables.

## Prerequisites

### Software

- **Liftoff** (includes minimap2), installed with mamba:

  ```bash
  mamba create -n liftoff -c conda-forge -c bioconda liftoff
  mamba activate liftoff        # or: conda activate liftoff
  liftoff --version
  ```

- **R ≥ 4.1** with `data.table`, `openxlsx` (Excel output) and `Seurat` (only if
  reporting genes from a Seurat object):

  ```r
  install.packages(c("data.table", "openxlsx", "Seurat"))
  ```

- **Hardware**: ~16 GB RAM is enough (peak ~12 GB). Runs in under an hour on a
  modern laptop; also runs on an SGE cluster (e.g. Eddie) via `qsub`.

### Input files

| File | Source |
|---|---|
| Old genome FASTA (uncompressed) | `Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa` (Ensembl Metazoa) |
| Old annotation | `Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3` (Ensembl Metazoa release 58) |
| New genome FASTA (uncompressed) | `GCF_963853765.1_xbMagGiga1.1_genomic.fna` (NCBI) |
| New annotation | `GCF_963853765.1_xbMagGiga1.1_genomic.gff.gz` (NCBI) |
| Old gene symbols *(optional)* | `02_Cg_gene_names.tsv` (no header: gene ID, symbol) |
| Old gene descriptions *(optional)* | `01_ORSON_french_group_2024_biorxiv_suppl2.csv` (`Gene.ID`, `Sequence.Description`) |
| Seurat object *(optional)* | `.rds` whose rownames are old `G#####` IDs |

Download the new genome and annotation from NCBI:

```bash
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1
wget $BASE/GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz
wget $BASE/GCF_963853765.1_xbMagGiga1.1_genomic.gff.gz
gunzip GCF_963853765.1_xbMagGiga1.1_genomic.fna.gz      # Liftoff needs an uncompressed FASTA
```

The old FASTA's sequence names must match the old GFF3 (both use `LR761634.1`, …);
the script prints both so you can check.

## Usage

### 1. Run Liftoff

Edit the input paths at the top of `run_liftoff_oyster.sh`, then:

```bash
conda activate liftoff
bash run_liftoff_oyster.sh        # laptop / workstation
# qsub run_liftoff_oyster.sh      # SGE cluster (Eddie)
```

Outputs: `Gmodels_on_xbMagGiga1.gff3` (lifted gene models) and `unmapped_genes.txt`.

### 2. Build the mapping table

Edit the settings at the top of `map_seurat_genes_liftoff.R`:

```r
SEURAT_RDS   <- NA     # NA = all old genes; or path to a Seurat .rds to report its features
LIFTED_GFF   <- "Gmodels_on_xbMagGiga1.gff3"
UNMAPPED_TXT <- "unmapped_genes.txt"
NEW_GFF      <- "GCF_963853765.1_xbMagGiga1.1_genomic.gff.gz"
SCIENCE_TSV  <- "02_Cg_gene_names.tsv"                            # NA to skip
ORSON_CSV    <- "01_ORSON_french_group_2024_biorxiv_suppl2.csv"   # NA to skip
```

Then run it **outside** the Liftoff environment (so your usual R with Seurat is used):

```bash
conda deactivate
Rscript map_seurat_genes_liftoff.R
```

## Output

`TableS_gene_ID_mapping_liftoff_Cgigas_to_Mgigas.tsv` and `.xlsx` (sheets: mapping,
summary, README with every column described, methods).

| Column | Description |
|---|---|
| `seurat_rowname` | Seurat feature name (only when a Seurat object is given) |
| `old_gene_id` | Old gene ID (`G#####`) |
| `old_symbol_science`, `old_description_ORSON` | Old gene symbol and description |
| `liftoff_status` | `lifted` / `unmapped_by_liftoff` / `not_in_liftoff_output` |
| `lift_*` | Lifted position, Liftoff coverage and sequence identity |
| `mapping_status`, `match_confidence`, `one_to_one` | Assignment result and confidence |
| `recip_overlap`, `frac_lifted_exons_in_gene`, `reciprocal_best` | Overlap scores behind the assignment |
| `new_GeneID`, `new_symbol`, `new_description`, `new_gene_biotype` | Current NCBI gene |
| `new_seqid`, `new_start`, `new_end`, `new_strand` | Current gene coordinates |

For downstream use, filter on `one_to_one == TRUE`; `match_confidence` distinguishes
high from medium assignments.

### Results for our data

| | Genes |
|---|---|
| Old genes in Seurat object | 33,466 |
| Lifted to xbMagGiga1.1 | 29,907 |
| Assigned to a current gene | 19,258 |
| — of which one-to-one | 18,589 |

## Citation

Shumate A, Salzberg SL (2021). Liftoff: accurate mapping of gene annotation.
*Bioinformatics* 37(12):1639–1643. https://doi.org/10.1093/bioinformatics/btaa1016
