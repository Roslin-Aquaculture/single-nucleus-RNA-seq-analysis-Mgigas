### Reference preparation

This directory contains scripts and notes used to prepare the reference files required for alignment and quantification in the Pacific oyster OsHV-1 single-nucleus RNA-seq analysis.

---

#### Overview

Reference preparation in this project involved four main tasks:

1. **Prepare the Pacific oyster genome FASTA**
2. **Obtain and convert the oyster gene annotation to GTF format**
3. **Check and refine the GTF so that it is compatible with the Parse pipeline**
4. **Generate a composite oyster + OsHV-1 reference for the final viral-inclusive analysis**

The final goal of this step is to create reference FASTA and GTF files that can be used for genome indexing and read mapping.

#### Background

The original oyster reference preparation in the source repository used the **Crassostrea gigas cgigas_uk_roslin_v1** assembly from Ensembl Metazoa release 58, downloaded as individual soft-masked chromosome FASTA files and combined into a single genome FASTA. The corresponding annotation was available as **GFF3**, which was converted to **GTF** using AGAT because the Parse pipeline requires GTF input.

In our workflow, the oyster genome and annotation were combined with an OsHV-1 genome FASTA and converted viral annotation to generate a composite reference for Parse indexing.

#### Directory contents

- `oyster_genome_files.txt`
  This file has fasta filenames for all the individual chromosomes for Pacific oyster. Make a copy of this .txt file in your current working directory to download oyster genome.
- `viv46-2-m_assembly_NR_2017.fasta` Genome file for OsHV-1
- `viv46-2-m_assembly_NR_final_ok.gff` GFF annotation file for OsHV-1
- 
---

#### Recommended workflow

##### 1. Download the oyster genome FASTA

The original workflow used the soft-masked Pacific oyster genome from [Ensmebl Metazoa release 58](http://ftp.ensemblgenomes.org/pub/metazoa/release-58/fasta/crassostrea_gigas/dna/). Individual chromosome FASTA files were downloaded and combined into one reference FASTA.

Soft-masked genome sequence was used rather than hard-masked sequence so that the full genome sequence remained available for alignment, while repetitive reads could still be filtered downstream if needed.

###### Example approach

```bash
# Run the code below to download chromosome-level FASTA files listed in oyster_genome_files.txt
while read LINE;
do
  LINK=$(echo "$LINE" | awk '{print $3}')
  echo "downloading $LINK"
  wget -c http://ftp.ensemblgenomes.org/pub/metazoa/release-58/fasta/crassostrea_gigas/dna/$LINK
  echo "unzipping $LINK"
  gunzip $LINK
done < genome_files.txt

# Now combine into one FASTA
cat *.fa > Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa
```

After this step, you should have a single oyster genome FASTA file ready for annotation matching and indexing.

##### 2. Obtain and convert the oyster annotation to GTF

The Parse pipeline requires a **GTF** annotation file. 
Download Pacific oyster [Ensembl GFF3 annotation](https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-58/gff3/crassostrea_gigas/) and then convert to GTF using [AGAT tool](https://agat.readthedocs.io/en/latest/gff_to_gtf.html). For AGAT documentation and download instructions, see [github page](https://github.com/NBISweden/AGAT).


###### Download and convert annotation

```bash
# download
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-58/gff3/crassostrea_gigas/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3.gz

# unzip
gunzip Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3.gz

# convert to gtf
agat_convert_sp_gff2gtf.pl \
  --gff Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3 \
  -o Crassostrea_gigas_uk_roslin_v1.gtf
```

A conversion log and summary statistics are useful to keep for troubleshooting and record keeping.

##### 3. Check and fix the GTF for Parse compatibility

Before indexing, the GTF should be checked for compliance with Parse formatting requirements.

Change the attribute name `biotype` to `gene_biotype`, because Parse expects the latter.

###### Example fix

```bash
sed -i 's/\bbiotype\b/gene_biotype/g' Crassostrea_gigas_uk_roslin_v1.gtf
```

This step is important because non-compliant GTF formatting can cause `split-pipe --mode mkref` to fail.


##### 4. Prepare the OsHV-1 annotation and genome files

We used OsHV-1 genome FASTA file and GFF annotation provided by Aurélie Dotto-Maurel.

###### Convert viral GFF to GTF

```bash
# convert viv46-2-m_assembly_NR_final_ok.gff -> GTF format
agat_convert_sp_gff2gtf.pl --gff viv46-2-m_assembly_NR_final_ok.gff -o viv46-2-m_assembly_NR_final_ok.gtf

# manually remove the first two lines with ## in the viral GFT file
```


##### 5. Build the composite oyster + OsHV-1 reference

In the final workflow, the oyster and viral references were combined into a single FASTA and a single GTF.

###### Example approach

```bash

# Combine GTF files
cat Crassostrea_gigas_uk_roslin_v1.gtf viv46-2-m_assembly_NR_final_ok.gtf >  Crassostrea_gigas_OsHV1_combined.gtf

# Manually ensure the viral sequence-region header is present if needed
##sequence-region   viv46-2-m_assembly_NR_final 1 186279

# Combine FASTA files
cat Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa viv46-2-m_assembly_NR_2017.fasta > Crassostrea_gigas_OsHV1_combined.fa
```

Before combining files, check that sequence names in the viral FASTA and viral annotation match exactly.

##### 6. Index the final reference for Parse

Once the final FASTA and GTF are ready, build the Parse reference index.

###### Example command

```bash
split-pipe --mode mkref \
  --genome_name cgigas_parse \
  --fasta /path/to/Crassostrea_gigas_OsHV1_combined.fa \
  --genes /path/to/Crassostrea_gigas_OsHV1_combined.gtf \
  --output_dir /path/to/cgigas_parse
```

This indexed reference can then be used for read mapping and quantification in the downstream pipeline.

#### Input files

Typical inputs for this step are:

- oyster genome FASTA files
- oyster annotation in GFF3 format
- OsHV-1 genome FASTA
- OsHV-1 annotation in GFF format

#### Output files

This step should produce:

- combined oyster genome FASTA
- oyster annotation in Parse-compatible GTF format
- composite oyster + OsHV-1 FASTA file
- composite oyster + OsHV-1 GTF file
- indexed Parse reference directory

#### Notes

- Keep the FASTA and annotation files matched to the same reference assembly wherever possible.
- The Parse pipeline requires a valid GTF and may fail if expected attributes are missing.
- Check that chromosome or contig names are identical between FASTA and GTF files.
- Keep logs from AGAT conversion and Parse indexing, as these are useful for troubleshooting.

#### Suggested checks after reference preparation

Before moving to mapping, confirm that:

- the composite oyster + OsHV-1 reference files were generated successfully
- Parse indexing completed without errors

#### Downstream use

The reference files generated here are used in the alignment and quantification step, including construction of the final host-virus composite reference used for viral transcript detection.
