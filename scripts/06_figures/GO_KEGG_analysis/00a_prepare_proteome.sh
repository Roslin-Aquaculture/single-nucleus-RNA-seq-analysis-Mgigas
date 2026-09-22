#!/bin/bash
set -euo pipefail

# =========================================================
# Section 2.3: Generate protein FASTA (genome-wide) and pick
# longest isoform per gene, then fix invalid characters and
# stage the cleaned FASTA for the eggNOG-mapper job (00b).
# Environment: eggnog (HPC / Eddie)
# =========================================================

# =========================================================
# PROJECT PATHS
# =========================================================
PROJECT="/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/eggnog"
GENOME_DIR="$PROJECT/genome"
GENOME="$GENOME_DIR/Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa"
GFF="$GENOME_DIR/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3"
FASTA="$PROJECT/input/cgigas_longest_clean.fa"
# generic (uncomment and edit if replicating elsewhere):
# PROJECT="/path/to/your/hpc/project/eggnog"
# GENOME_DIR="$PROJECT/genome"
# GENOME="$GENOME_DIR/your_genome.dna_sm.primary_assembly.fa"
# GFF="$GENOME_DIR/your_annotation.gff3"
# FASTA="$PROJECT/input/cgigas_longest_clean.fa"

cd "$GENOME_DIR"

# =========================================================
# GENERATE PROTEIN FASTA (ALL TRANSCRIPTS, GENOME-WIDE)
# =========================================================
gffread "$GFF" -g "$GENOME" -y cgigas_proteins.fa

# =========================================================
# GENE -> TRANSCRIPT -> LENGTH TABLE
# =========================================================
seqkit fx2tab -n -l cgigas_proteins.fa \
| awk '{
    header=$1;
    len=$2;
    gene=header;
    sub(/transcript:/,"",gene);
    sub(/\..*/,"",gene);
    print gene"\t"len"\t"header
}' > id_len_map.txt

# =========================================================
# RANK ISOFORMS PER GENE (1 = LONGEST, 2 = SECOND LONGEST, ...)
# Makes the "longest transcript per gene" selection auditable,
# rather than just trusting a sort+awk chain blindly.
# =========================================================
sort -k1,1 -k2,2nr id_len_map.txt > id_len_map_sorted.txt

awk '{
    if ($1 != prev_gene) { rank = 1 }
    else { rank++ }
    prev_gene = $1
    print $0"\t"rank
}' id_len_map_sorted.txt > id_len_map_ranked.txt

awk '$4 == 1 {print $3}' id_len_map_ranked.txt > longest_ids.txt

echo "========================================"
echo "Total transcripts in table:  $(wc -l < id_len_map.txt)"
echo "Total unique genes (rank 1): $(awk '$4==1' id_len_map_ranked.txt | wc -l)"
echo "Total longest_ids.txt lines: $(wc -l < longest_ids.txt)"
echo "========================================"

# spot-check a gene's ranked isoforms, e.g.:
grep -w "G1001" id_len_map_ranked.txt

# =========================================================
# BUILD ONE-PROTEIN-PER-GENE PROTEOME
# =========================================================
seqkit grep -f longest_ids.txt cgigas_proteins.fa > cgigas_longest.fa

# sanity check — these two numbers should match
grep -c ">" cgigas_longest.fa
wc -l longest_ids.txt

# =========================================================
# FIX INVALID CHARACTERS BEFORE RUNNING DIAMOND
# gffread -y can leave "." placeholder characters in
# incompletely-translated sequences, commonly in mitochondrial
# genes (different genetic code than nuclear genes). Diamond
# errors on these ("Invalid character (.) in sequence").
# Recommended fix: exclude mitochondrial genes entirely
# (typically a tiny number, not usually detected/relevant in
# snRNA-seq nuclear transcriptome data).
# =========================================================
echo "Sequences with invalid characters: $(grep -v '^>' cgigas_longest.fa | grep -c '\.')"
echo "Mitochondrial gene sequences:      $(grep -c '_df_mr' cgigas_longest.fa)"

awk '/^>/{p = ($0 !~ /_df_mr/)} p' cgigas_longest.fa > cgigas_longest_clean.fa

echo "Invalid characters remaining: $(grep -v '^>' cgigas_longest_clean.fa | grep -c '\.')"

# =========================================================
# STAGE CLEANED FASTA FOR THE EGGNOG-MAPPER JOB (00b)
# =========================================================
mkdir -p "$(dirname "$FASTA")"
cp cgigas_longest_clean.fa "$FASTA"

echo "Staged for eggNOG-mapper:"
ls -lh "$FASTA"
