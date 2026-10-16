#!/usr/bin/env Rscript
# =============================================================================
# Map Seurat rownames (old C. gigas genes, Ensembl/Roslin G##### models on
# cgigas_uk_roslin_v1) to current M. gigas genes (xbMagGiga1.1, GCF_963853765.1)
# using Liftoff.
#
#   1. Liftoff (run_liftoff_oyster.sh) places each old gene model on the new genome
#   2. Lifted exons are overlapped with the current NCBI annotation (same strand)
#   3. Each old gene gets its best current gene, scored by reciprocal exonic overlap
#
# Output: one row per Seurat feature -> supplementary table (CSV + XLSX)
# =============================================================================
suppressPackageStartupMessages(library(data.table))

setwd("/home/pdewari/genomes/xbMagGiga1.1/")

## ---- settings ----------------------------------------------------------------
SEURAT_RDS   <- "/home/pdewari/Documents/parse_2025/seurat_2025/seu_obj_umap_18d_6r_3kRes.rds"
LIFTED_GFF   <- "Gmodels_on_xbMagGiga1.gff3"                     # Liftoff output
UNMAPPED_TXT <- "unmapped_genes.txt"                             # Liftoff -u output
NEW_GFF      <- "GCF_963853765.1_xbMagGiga1.1_genomic.gff.gz"    # current NCBI annotation
OUT_PREFIX   <- "TableS_gene_ID_mapping_liftoff_Cgigas_to_Mgigas"
# old-gene annotation added next to old_gene_id (set either to NA to skip)
SCIENCE_TSV  <- "/home/pdewari/Documents/parse_2025/seurat_2025/02_Cg_gene_names.tsv"                                # no header: gene ID, symbol
ORSON_CSV    <- "/home/pdewari/Documents/parse_2025/seurat_2025/01_ORSON_french_group_2024_biorxiv_suppl2.csv"       # columns Gene.ID, Sequence.Description

MIN_RECIP    <- 0.5   # HIGH: reciprocal exonic overlap >= this
MIN_FRAC_OLD <- 0.5   # MEDIUM: reciprocal best pair AND >= this fraction of lifted exons in the gene
LOW_SEQ_ID   <- 0.9   # flag lifts whose Liftoff sequence identity is below this
NEW_ASSEMBLY <- "xbMagGiga1.1 (GCF_963853765.1), NCBI annotation RS_2024_06"

## ---- helpers -----------------------------------------------------------------
cat_cmd  <- function(p) sprintf(if (grepl("\\.gz$", p)) "gzip -dc '%s'" else "cat '%s'", p)
gff_cols <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr")
attr_val <- function(attr, key) {
  out <- rep(NA_character_, length(attr))
  hit <- grepl(paste0("(^|;)", key, "="), attr)
  out[hit] <- sub(paste0(".*(^|;)", key, "=([^;]+).*"), "\\2", attr[hit])
  out
}
url_decode <- function(x) {
  i <- !is.na(x) & grepl("%", x, fixed = TRUE)
  x[i] <- vapply(x[i], utils::URLdecode, character(1), USE.NAMES = FALSE)
  x
}
read_gff <- function(p) {
  g <- fread(cmd = sprintf("%s | awk -F'\\t' '!/^#/ && NF==9'", cat_cmd(p)),
             sep = "\t", header = FALSE, quote = "", col.names = gff_cols)
  g[, `:=`(start = as.integer(start), end = as.integer(end))]
  g
}
# union of exons per gene (isoforms merged so no base is counted twice)
reduce_exons <- function(e) {
  e <- unique(e)
  setorder(e, gene_id, seqid, strand, start, end)
  e[, grp := cumsum(c(TRUE, start[-1] > head(cummax(end), -1) + 1L)),
    by = .(gene_id, seqid, strand)]
  e[, .(start = min(start), end = max(end)), by = .(gene_id, seqid, strand, grp)][, grp := NULL][]
}

## ---- 1. which genes to report -------------------------------------------------
use_seurat <- !is.na(SEURAT_RDS)
if (use_seurat) {
  suppressPackageStartupMessages(library(Seurat))
  obj  <- readRDS(SEURAT_RDS)
  feat <- data.table(seurat_rowname = rownames(obj))
  rm(obj); invisible(gc())
  message(nrow(feat), " features in Seurat object")
}

## ---- 2. lifted old gene models -----------------------------------------------
lg <- read_gff(LIFTED_GFF)
lg[, `:=`(ID = attr_val(attr, "ID"), Parent = attr_val(attr, "Parent"))]

genes_l <- lg[grepl("gene$", type) & is.na(Parent)]
genes_l[, gene_id := attr_val(attr, "gene_id")]
genes_l[is.na(gene_id), gene_id := sub("^gene[:-]", "", ID)]
genes_l[, `:=`(lift_coverage   = as.numeric(attr_val(attr, "coverage")),
               lift_seq_id     = as.numeric(attr_val(attr, "sequence_ID")),
               lift_extra_copies = as.integer(attr_val(attr, "extra_copy_number")))]
if (anyDuplicated(genes_l$gene_id)) {
  message(sum(duplicated(genes_l$gene_id)), " duplicated lifted gene IDs (extra copies) - keeping first")
}
genes_l <- unique(genes_l, by = "gene_id")

tx <- lg[!is.na(ID) & Parent %in% genes_l$ID, .(tx = ID, gene_raw = Parent)]
ex <- lg[type == "exon" & !is.na(Parent), .(seqid, start, end, strand, par = Parent)]
multi <- grepl(",", ex$par, fixed = TRUE)
if (any(multi)) ex <- rbind(ex[!multi],
                            ex[multi, .(par = unlist(strsplit(par, ",", fixed = TRUE))), by = .(seqid, start, end, strand)])
ex[, gene_raw := tx$gene_raw[match(par, tx$tx)]]
ex[is.na(gene_raw) & par %in% genes_l$ID, gene_raw := par]
ex_old <- ex[!is.na(gene_raw), .(gene_id = genes_l$gene_id[match(gene_raw, genes_l$ID)],
                                 seqid, start, end, strand)]
ex_old <- rbind(ex_old, genes_l[!gene_id %in% ex_old$gene_id, .(gene_id, seqid, start, end, strand)])
lifted <- genes_l[, .(gene_id, lift_seqid = seqid, lift_start = start, lift_end = end,
                      lift_strand = strand, lift_coverage, lift_seq_id, lift_extra_copies)]
rm(lg, genes_l, tx, ex); invisible(gc())

unmapped <- if (file.exists(UNMAPPED_TXT)) {
  u <- trimws(readLines(UNMAPPED_TXT)); sub("^gene[:-]", "", u[u != ""])
} else character(0)
message(nrow(lifted), " lifted genes; ", length(unmapped), " unmapped by Liftoff")
message("example lifted IDs: ", paste(head(lifted$gene_id, 5), collapse = ", "),
        " | seqids: ", paste(head(unique(lifted$lift_seqid), 3), collapse = ", "))

## ---- 3. current NCBI annotation ----------------------------------------------
ng <- read_gff(NEW_GFF)
new_genes <- ng[type %in% c("gene", "pseudogene") & grepl("GeneID:", attr)]
new_genes[, GeneID := sub(".*GeneID:([0-9]+).*", "\\1", attr)]
new_genes <- unique(new_genes, by = "GeneID")
new_info <- new_genes[, .(new_GeneID = GeneID,
                          new_symbol = attr_val(attr, "Name"),
                          new_description = url_decode(attr_val(attr, "description")),
                          new_gene_biotype = attr_val(attr, "gene_biotype"),
                          new_seqid = seqid, new_start = start, new_end = end, new_strand = strand)]
ex_new <- ng[type == "exon" & grepl("GeneID:", attr),
             .(gene_id = sub(".*GeneID:([0-9]+).*", "\\1", attr), seqid, start, end, strand)]
ex_new <- rbind(ex_new, new_genes[!GeneID %in% ex_new$gene_id,
                                  .(gene_id = GeneID, seqid, start, end, strand)])
rm(ng, new_genes); invisible(gc())
message(nrow(new_info), " genes in current annotation")

shared_seq <- intersect(unique(ex_old$seqid), unique(ex_new$seqid))
if (!length(shared_seq)) stop("Lifted GFF and NEW_GFF share no sequence names - ",
                              "was Liftoff run against the NCBI (GCF) FASTA?")

## ---- 4. exonic overlap: lifted old model vs current gene ---------------------
xo <- reduce_exons(ex_old); xn <- reduce_exons(ex_new)
rm(ex_old, ex_new); invisible(gc())
len_old <- xo[, .(len = sum(end - start + 1)), by = gene_id]
len_new <- xn[, .(len = sum(end - start + 1)), by = gene_id]
setnames(xn, "gene_id", "GeneID")
setkey(xn, seqid, strand, start, end)
ov <- foverlaps(xo, xn, by.x = c("seqid", "strand", "start", "end"), type = "any", nomatch = NULL)
ov[, bp := pmin(end, i.end) - pmax(start, i.start) + 1L]
sh <- ov[, .(shared_bp = sum(bp)), by = .(gene_id, GeneID)]
sh[, `:=`(frac_old = shared_bp / len_old$len[match(gene_id, len_old$gene_id)],
          frac_new = shared_bp / len_new$len[match(GeneID, len_new$gene_id)])]
sh[, recip := pmin(frac_old, frac_new)]
setorder(sh, gene_id, -recip, -shared_bp)
best <- sh[, .(best_GeneID = GeneID[1],
               recip_overlap = round(recip[1], 3),
               frac_lifted_exons_in_gene = round(frac_old[1], 3),
               n_overlapping_genes = .N,
               second_recip_overlap = if (.N > 1) round(recip[2], 3) else NA_real_),
           by = gene_id]
new_best <- sh[order(GeneID, -recip, -shared_bp)][, .(best_old = gene_id[1]), by = GeneID]
best[, reciprocal_best := new_best$best_old[match(best_GeneID, new_best$GeneID)] == gene_id]
best[, split_across_genes := !is.na(second_recip_overlap) & second_recip_overlap >= MIN_RECIP]
best[, match_confidence := fifelse(recip_overlap >= MIN_RECIP, "high",
                                   fifelse(reciprocal_best & frac_lifted_exons_in_gene >= MIN_FRAC_OLD,
                                           "medium", NA_character_))]
rm(ov, sh, xo, xn, new_best); invisible(gc())

## ---- 5. assemble one row per gene / Seurat feature ---------------------------
known <- unique(c(lifted$gene_id, unmapped))
if (!use_seurat) feat <- data.table(seurat_rowname = sort(known))
feat[, old_gene_id := fifelse(seurat_rowname %in% known, seurat_rowname,
                              fifelse(gsub("-", "_", seurat_rowname) %in% known,
                                      gsub("-", "_", seurat_rowname), NA_character_))]
res <- merge(feat, lifted, by.x = "old_gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
res <- merge(res, best, by.x = "old_gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
res[, liftoff_status := fifelse(!is.na(lift_seqid), "lifted",
                                fifelse(old_gene_id %in% unmapped, "unmapped_by_liftoff",
                                        "not_in_liftoff_output"))]
res[, new_GeneID := fifelse(!is.na(match_confidence), best_GeneID, NA_character_)]
res[, mapping_status := fcase(
  !is.na(new_GeneID),                           "mapped",
  liftoff_status == "lifted" & !is.na(best_GeneID), "weak_overlap_only",
  liftoff_status == "lifted",                   "lifted_no_overlapping_gene",
  rep(TRUE, .N),                                liftoff_status)]

res[!is.na(new_GeneID), n_features_per_new_gene := .N, by = new_GeneID]
res[!is.na(new_GeneID) & n_features_per_new_gene > 1L,
    mapping_status := paste0(mapping_status, ";many_to_one")]
res[!is.na(new_GeneID) & split_across_genes %in% TRUE,
    mapping_status := paste0(mapping_status, ";spans_2+_genes")]
res[, low_identity_lift := !is.na(lift_seq_id) & lift_seq_id < LOW_SEQ_ID]
res[, one_to_one := !is.na(new_GeneID) & n_features_per_new_gene == 1L &
      !(split_across_genes %in% TRUE)]

res <- merge(res, new_info, by = "new_GeneID", all.x = TRUE, sort = FALSE)
res[, new_assembly := fifelse(is.na(new_GeneID), NA_character_, NEW_ASSEMBLY)]
res <- res[match(feat$seurat_rowname, seurat_rowname)]

## ---- 5b. old-gene symbol / description ---------------------------------------
collapse_ann <- function(id, val) {
  d <- data.table(id = as.character(id), val = as.character(val))[!is.na(id) & !is.na(val) & val != ""]
  d[, .(val = paste(unique(val), collapse = "; ")), by = id]
}
res[, .key := fifelse(is.na(old_gene_id), seurat_rowname, old_gene_id)]
res[, `:=`(old_symbol_science = NA_character_, old_description_ORSON = NA_character_)]
if (!is.na(SCIENCE_TSV)) {
  sci <- fread(SCIENCE_TSV, header = FALSE, sep = "\t", quote = "")
  a <- collapse_ann(sci[[1]], sci[[2]])
  res[, old_symbol_science := a$val[match(.key, a$id)]]
  message("symbols (", basename(SCIENCE_TSV), "): ", sum(!is.na(res$old_symbol_science)), " of ", nrow(res))
}
if (!is.na(ORSON_CSV)) {
  ors <- fread(ORSON_CSV)
  stopifnot("ORSON file needs columns Gene.ID and Sequence.Description" =
              all(c("Gene.ID", "Sequence.Description") %in% names(ors)))
  a <- collapse_ann(ors$Gene.ID, ors$Sequence.Description)
  res[, old_description_ORSON := a$val[match(.key, a$id)]]
  message("descriptions (", basename(ORSON_CSV), "): ", sum(!is.na(res$old_description_ORSON)), " of ", nrow(res))
}
res[, .key := NULL]

cols <- c(
  seurat_rowname = "Feature name in the Seurat object",
  old_gene_id = "Gene ID in the old annotation (Ensembl release 58 / Roslin, cgigas_uk_roslin_v1)",
  old_symbol_science = "Gene symbol for the old gene ID (02_Cg_gene_names.tsv)",
  old_description_ORSON = "Sequence description for the old gene ID (ORSON, bioRxiv 2024, supplementary table 2)",
  liftoff_status = "lifted / unmapped_by_liftoff / not_in_liftoff_output",
  lift_seqid = "Sequence on the new assembly where Liftoff placed the old gene",
  lift_start = "Lifted gene start", lift_end = "Lifted gene end", lift_strand = "Lifted gene strand",
  lift_coverage = "Liftoff: fraction of the old gene's exonic sequence aligned to the new genome",
  lift_seq_id = "Liftoff: sequence identity of the aligned exons",
  lift_extra_copies = "Liftoff: number of additional copies found",
  low_identity_lift = paste0("TRUE if Liftoff sequence identity < ", LOW_SEQ_ID),
  mapping_status = "mapped / weak_overlap_only / lifted_no_overlapping_gene / unmapped_by_liftoff / not_in_liftoff_output; flags ';many_to_one' and ';spans_2+_genes'",
  match_confidence = paste0("high = reciprocal exonic overlap >= ", MIN_RECIP,
                            "; medium = reciprocal best match with >= ", MIN_FRAC_OLD * 100,
                            "% of lifted exons inside the gene"),
  one_to_one = "TRUE if the feature maps to exactly one current gene, not shared with another feature and not split across genes",
  n_features_per_new_gene = "Number of Seurat features mapping to the same current gene",
  recip_overlap = "Reciprocal exonic overlap with the best current gene: shared bp / larger of the two merged exon lengths",
  frac_lifted_exons_in_gene = "Fraction of the lifted gene's exonic bp inside the best current gene",
  n_overlapping_genes = "Number of current genes sharing exonic sequence with the lifted gene (same strand)",
  second_recip_overlap = "Reciprocal exonic overlap of the second-best current gene",
  reciprocal_best = "TRUE if the lifted gene and the current gene are each other's best match",
  split_across_genes = paste0("TRUE if a second current gene also has reciprocal overlap >= ", MIN_RECIP),
  best_GeneID = "Best-overlapping current gene, reported even when below the confidence thresholds",
  new_GeneID = "Assigned current NCBI GeneID (blank unless high or medium confidence)",
  new_symbol = "Current gene symbol", new_description = "Current gene description",
  new_gene_biotype = "Current gene biotype",
  new_seqid = "RefSeq sequence accession (new assembly)", new_start = "Gene start (new assembly)",
  new_end = "Gene end (new assembly)", new_strand = "Strand (new assembly)",
  new_assembly = "Assembly and annotation release")
res <- res[, names(cols), with = FALSE]
stopifnot(nrow(res) == nrow(feat), !anyDuplicated(res$seurat_rowname))

## ---- 6. summary + save -------------------------------------------------------
summ <- res[, .N, by = .(liftoff_status, match_confidence, mapping_status)][order(-N)]
summ <- rbind(summ, data.table(liftoff_status = "TOTAL", match_confidence = NA_character_,
                               mapping_status = "one_to_one", N = sum(res$one_to_one)))
print(summ)

if (!use_seurat) { res[, seurat_rowname := NULL]; cols <- cols[names(cols) != "seurat_rowname"] }
fwrite(res, paste0(OUT_PREFIX, ".tsv"), sep = "\t")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  readme <- data.table(column = names(cols), description = unname(cols))
  methods <- data.table(note = c(
    "Old gene models (Ensembl Metazoa release 58, cgigas_uk_roslin_v1) were lifted onto xbMagGiga1.1 (GCF_963853765.1) with Liftoff.",
    "Lifted exons were merged per gene and compared with exons of the current NCBI annotation (RS_2024_06) on the same strand.",
    paste0("High confidence: reciprocal exonic overlap >= ", MIN_RECIP,
           ". Medium: reciprocal best match with >= ", MIN_FRAC_OLD * 100, "% of lifted exons inside the gene.")))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "gene_ID_mapping"); openxlsx::writeData(wb, 1, res)
  openxlsx::addWorksheet(wb, "summary");         openxlsx::writeData(wb, 2, summ)
  openxlsx::addWorksheet(wb, "README");          openxlsx::writeData(wb, 3, readme)
  openxlsx::addWorksheet(wb, "methods");         openxlsx::writeData(wb, 4, methods)
  openxlsx::freezePane(wb, 1, firstRow = TRUE)
  openxlsx::saveWorkbook(wb, paste0(OUT_PREFIX, ".xlsx"), overwrite = TRUE)
  message("Saved ", OUT_PREFIX, ".tsv and .xlsx")
} else {
  message("Saved ", OUT_PREFIX, ".tsv (install 'openxlsx' for the .xlsx)")
}
