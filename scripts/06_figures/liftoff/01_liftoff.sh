OLD_FA=/home/pdewari/Documents/parse_2025/seurat_2025/genome/Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa
OLD_GFF=/home/pdewari/Documents/parse_2025/seurat_2025/genome/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3
NEW_FA=/home/pdewari/genomes/xbMagGiga1.1/GCF_963853765.1_xbMagGiga1.1_genomic.fna

# FASTAs must be uncompressed
ls -lh $OLD_FA $NEW_FA

# old FASTA names must match old GFF3 seqids (expect LR761634.1 in both)
grep -m3 '>' $OLD_FA
grep -v '^#' $OLD_GFF | cut -f1 | sort -u | head -3

# new FASTA should be NC_0888xx
grep -m2 '>' $NEW_FA

printf "gene\nncRNA_gene\npseudogene\n" > feature_types.txt

nohup liftoff \
  -g $OLD_GFF \
  -o Gmodels_on_xbMagGiga1.gff3 \
  -u unmapped_genes.txt \
  -f feature_types.txt \
  -dir liftoff_intermediate \
  -p $(( $(nproc) - 1 )) \
  $NEW_FA $OLD_FA > liftoff.log 2>&1 &
