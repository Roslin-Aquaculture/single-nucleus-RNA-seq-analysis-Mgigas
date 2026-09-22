# =========================================================
# Section 3 (v2): Building the Enrichment Reference Tables
# Environment: go-enrich (local machine)
# =========================================================
# Drop-in replacement for 01_build_reference_tables.R.
#
# WHAT CHANGED FROM v1 AND WHY
#
#   1. GO ANNOTATIONS ARE NOW PROPAGATED UP THE GO DAG.
#      eggNOG's GOs column holds directly-assigned terms only. enricher()
#      with a custom TERM2GENE does NO ancestor propagation (unlike
#      enrichGO() with an OrgDb, which uses the full ancestor closure
#      internally). Without propagation a gene annotated to
#      "GO:0006955 immune response" does not count toward
#      "GO:0002376 immune system process", parent terms are starved, hits
#      fragment across sparse child terms, and broad themes under-enrich.
#
#   2. GO IS SPLIT BY ONTOLOGY (BP / MF / CC).
#      v1 mixed all three in one TERM2GENE, so BH correction pooled three
#      ontologies and output interleaved "cilium" (CC) with "cilium
#      movement" (BP). Three separate runs is standard practice. It also
#      isolates CC, which is where cell-identity signal concentrates.
#
#   3. KEGG map##### IDs ARE FOLDED INTO ko#####.
#      eggNOG lists BOTH namespaces for the same pathway
#      ("ko00010,map00010"), so v1 carried every KEGG pathway twice —
#      doubling the multiple-testing burden with exact duplicates and
#      putting paired rows in every dotplot.
#
#   4. KEGG GLOBAL/OVERVIEW MAPS ARE DROPPED.
#      ko01100 "Metabolic pathways" contains thousands of genes and is
#      significant in essentially any list. Left in, it tops every panel
#      and says nothing.
#
#   5. THE GLOBAL BACKGROUND IS KEPT, BUT ITS SCOPE IS NARROWED.
#      intersect(seurat_genes, eggnog_genes) is the CORRECT universe for
#      the cluster-marker enrichment in 02/03 (any gene in the object
#      could have been a marker). It is the WRONG universe for the
#      infection-stage DE enrichment in 04/05, where a gene was only
#      eligible if it was testable IN THAT CLUSTER. Per-cluster universes
#      come from 01c_export_cluster_universes.R.
#
# Run once per R session before any enrichment analysis.
#
# Before running: transfer full_proteome.emapper.annotations from Eddie to
# the local machine (e.g. scp or rsync).

library(readr)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)

# Update the setwd() path below to match your own eggNOG output directory.
setwd("/home/pdewari/eggnog/results/full_proteome_20260820_124651")
# generic (uncomment and edit if replicating elsewhere):
# setwd("/path/to/your/eggnog/results/full_proteome_<timestamp>")

REF_DIR <- "reference_tables_v2"
dir.create(REF_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 3.1 Load + clean eggNOG annotation  (unchanged from v1)
# ---------------------------------------------------------
eggnog_full <- read_tsv("full_proteome.emapper.annotations", comment = "##",
                        show_col_types = FALSE)

# eggNOG's #query column looks like "transcript:G3710.1" — strip down to
# gene-level ID ("G3710") to match Seurat rownames / marker table gene IDs
eggnog_full <- eggnog_full %>%
  mutate(gene_id = sub("^transcript:", "", `#query`)) %>%
  mutate(gene_id = sub("\\..*$", "", gene_id))

cat("eggNOG genes annotated:", nrow(eggnog_full), "\n")

# ---------------------------------------------------------
# 3.2 Global background — for the CLUSTER MARKER analysis only
# ---------------------------------------------------------
# Every gene eligible to be a marker = all genes in the Seurat matrix
# after QC, restricted to those with an eggNOG annotation.
#
# Do NOT use this for the infection-stage DE enrichment. See 01c.

# Run once in the Seurat session (wherever seurat_obj lives):
#   raw_names <- rownames(seurat_obj[["RNA"]])
#   seurat_genes_vec <- sub("^(\\S+)\\s.*$", "\\1", raw_names)
#   write.csv(data.frame(gene_id = seurat_genes_vec), "seurat_genes.csv", row.names = FALSE)

seurat_genes <- read_csv("/home/pdewari/Documents/parse_2025/seurat_2025/seurat_genes.csv",
                         show_col_types = FALSE)
seurat_genes_vec <- seurat_genes$gene_id

background_global <- intersect(seurat_genes_vec, eggnog_full$gene_id)
cat("Global background gene count (markers only):", length(background_global), "\n")

# kept under the old name so 02/03 run unchanged
background <- background_global

write.csv(data.frame(gene_id = background_global),
          file.path(REF_DIR, "background_global.csv"), row.names = FALSE)

# ---------------------------------------------------------
# 3.3 GO term2gene — raw (directly assigned), then PROPAGATED
# ---------------------------------------------------------
# IMPORTANT: enricher() reads TERM2GENE by column POSITION, not name.
# Column 1 = term ID, column 2 = gene ID. Reversing this causes a silent
# "No gene can be mapped" failure.

term2gene_raw <- eggnog_full %>%
  dplyr::select(gene_id, GOs) %>%
  filter(GOs != "-", !is.na(GOs)) %>%
  separate_rows(GOs, sep = ",") %>%
  dplyr::rename(GO_ID = GOs) %>%
  distinct() %>%
  dplyr::select(GO_ID, gene_id)

cat("GO term2gene rows (direct, as in v1):", nrow(term2gene_raw), "\n")

# --- ancestor closure -------------------------------------------------
# GOBPANCESTOR / GOMFANCESTOR / GOCCANCESTOR give every ancestor of a term
# within its own ontology. Obsolete or non-GO.db IDs coming out of eggNOG
# simply map to themselves.
anc_list <- c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR), as.list(GOCCANCESTOR))
anc_names <- names(anc_list)

uniq_go <- unique(term2gene_raw$GO_ID)
expand_map <- vector("list", length(uniq_go))
names(expand_map) <- uniq_go
for (g in uniq_go) {
  a <- if (g %in% anc_names) anc_list[[g]] else character(0)
  a <- a[!is.na(a) & a != "all"]
  expand_map[[g]] <- unique(c(g, a))   # always includes the term itself
}

term2gene_prop <- term2gene_raw %>%
  mutate(anc = expand_map[GO_ID]) %>%
  tidyr::unnest(anc) %>%
  transmute(GO_ID = anc, gene_id = gene_id) %>%
  distinct()

cat("GO term2gene rows (propagated):", nrow(term2gene_prop), "\n")
cat("  expansion factor:", round(nrow(term2gene_prop) / nrow(term2gene_raw), 1), "x\n")

# SANITY CHECK: if the expansion factor is ~1x, eggNOG already propagated
# for you and this step was a no-op (fine). If it is 5-40x, it did not,
# and v1 was losing every parent term.

# drop the three ontology roots — they contain everything and are useless
GO_ROOTS <- c("GO:0008150", "GO:0003674", "GO:0005575")
term2gene_prop <- term2gene_prop %>% filter(!GO_ID %in% GO_ROOTS)

# ---------------------------------------------------------
# 3.4 GO term names + ontology split
# ---------------------------------------------------------
go_info <- AnnotationDbi::select(
  GO.db,
  keys     = unique(term2gene_prop$GO_ID),
  columns  = c("TERM", "ONTOLOGY"),
  keytype  = "GOID"
) %>% filter(!is.na(TERM))

term2name_go <- go_info %>% transmute(GO_ID = GOID, name = TERM)

bp_ids <- go_info$GOID[go_info$ONTOLOGY == "BP"]
mf_ids <- go_info$GOID[go_info$ONTOLOGY == "MF"]
cc_ids <- go_info$GOID[go_info$ONTOLOGY == "CC"]

t2g_BP <- term2gene_prop %>% filter(GO_ID %in% bp_ids)
t2g_MF <- term2gene_prop %>% filter(GO_ID %in% mf_ids)
t2g_CC <- term2gene_prop %>% filter(GO_ID %in% cc_ids)

cat("GO terms after propagation — BP:", length(unique(t2g_BP$GO_ID)),
    "| MF:", length(unique(t2g_MF$GO_ID)),
    "| CC:", length(unique(t2g_CC$GO_ID)), "\n")

# ---------------------------------------------------------
# 3.5 KEGG term2gene — deduplicated, global maps removed
# ---------------------------------------------------------
# KEGG "Global and overview maps" (BRITE category 1.0). These are
# umbrella maps containing thousands of genes; they are significant in
# almost any gene list and carry no interpretive content.
KEGG_GLOBAL <- c("ko01100", "ko01110", "ko01120", "ko01200", "ko01210",
                 "ko01212", "ko01230", "ko01232", "ko01240", "ko01250")

kegg_term2gene <- eggnog_full %>%
  dplyr::select(gene_id, KEGG_Pathway) %>%
  filter(KEGG_Pathway != "-", !is.na(KEGG_Pathway)) %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  # fold map##### into ko##### so the same pathway is not counted twice
  mutate(pathway_id = sub("^map", "ko", KEGG_Pathway)) %>%
  filter(grepl("^ko[0-9]{5}$", pathway_id)) %>%
  filter(!pathway_id %in% KEGG_GLOBAL) %>%
  distinct(pathway_id, gene_id) %>%
  dplyr::select(pathway_id, gene_id)   # term first, gene second

cat("KEGG term2gene rows (deduplicated):", nrow(kegg_term2gene), "\n")
cat("KEGG pathways retained:", length(unique(kegg_term2gene$pathway_id)), "\n")

# --- KEGG pathway names, with defensive ID normalisation --------------
# download_KEGG("ko")$KEGGPATHID2NAME returns two columns, but whether the
# IDs arrive as "ko00010" or bare "00010" varies between clusterProfiler /
# KEGG REST versions. Normalise both to "ko#####".
kegg_names_raw <- download_KEGG("ko")$KEGGPATHID2NAME
stopifnot(ncol(kegg_names_raw) == 2)
names(kegg_names_raw) <- c("from", "to")

kegg_names <- kegg_names_raw %>%
  mutate(pathway_id = paste0("ko", sub("^[A-Za-z]+", "", from))) %>%
  transmute(pathway_id, name = to) %>%
  distinct()

# confirm the join will actually work before you run 800 enrichments
n_matched <- length(intersect(unique(kegg_term2gene$pathway_id), kegg_names$pathway_id))
cat("KEGG pathways with a matching name:", n_matched,
    "/", length(unique(kegg_term2gene$pathway_id)), "\n")
if (n_matched == 0) stop("KEGG ID formats do not match — inspect kegg_names_raw$from")

# ---------------------------------------------------------
# 3.6 Save all reference tables
# ---------------------------------------------------------
# Propagated GO tables are large; .rds is far faster to reload than .csv.
saveRDS(t2g_BP,        file.path(REF_DIR, "t2g_BP.rds"))
saveRDS(t2g_MF,        file.path(REF_DIR, "t2g_MF.rds"))
saveRDS(t2g_CC,        file.path(REF_DIR, "t2g_CC.rds"))
saveRDS(term2name_go,  file.path(REF_DIR, "term2name_go.rds"))
saveRDS(kegg_term2gene, file.path(REF_DIR, "kegg_term2gene.rds"))
saveRDS(kegg_names,    file.path(REF_DIR, "kegg_names.rds"))
saveRDS(background_global, file.path(REF_DIR, "background_global.rds"))

# also keep the v1-style unpropagated table, for a background/propagation
# sensitivity comparison if a reviewer asks
saveRDS(term2gene_raw, file.path(REF_DIR, "t2g_GO_unpropagated.rds"))

cat("\nAll reference tables built and saved to:", REF_DIR, "\n")

# To reload in a future session instead of rebuilding:
#   REF_DIR        <- "reference_tables_v2"
#   t2g_BP         <- readRDS(file.path(REF_DIR, "t2g_BP.rds"))
#   t2g_MF         <- readRDS(file.path(REF_DIR, "t2g_MF.rds"))
#   t2g_CC         <- readRDS(file.path(REF_DIR, "t2g_CC.rds"))
#   term2name_go   <- readRDS(file.path(REF_DIR, "term2name_go.rds"))
#   kegg_term2gene <- readRDS(file.path(REF_DIR, "kegg_term2gene.rds"))
#   kegg_names     <- readRDS(file.path(REF_DIR, "kegg_names.rds"))
#   background_global <- readRDS(file.path(REF_DIR, "background_global.rds"))

# ---------------------------------------------------------
# Sanity-check numbers (record these for the Methods section)
# ---------------------------------------------------------
cat("\n--- record these for Methods ---\n")
cat("Total Seurat genes:                 ", length(seurat_genes_vec), "\n")
cat("With eggNOG annotation (global bg): ", length(background_global), "\n")
cat("With >=1 GO term (propagated):      ", length(unique(term2gene_prop$gene_id)), "\n")
cat("  BP:", length(unique(t2g_BP$gene_id)),
    "| MF:", length(unique(t2g_MF$gene_id)),
    "| CC:", length(unique(t2g_CC$gene_id)), "\n")
cat("With >=1 KEGG pathway:              ", length(unique(kegg_term2gene$gene_id)), "\n")


###
for (nm in c("BP","MF","CC")) {
s <- table(get(paste0("t2g_", nm))$GO_ID)
cat(nm, "— terms total:", length(s),
    "| testable (size 10-500):", sum(s >= 10 & s <= 500), "\n")
}
s <- table(kegg_term2gene$pathway_id)
cat("KEGG — pathways total:", length(s),
    "| testable (size 10-500):", sum(s >= 10 & s <= 500), "\n")

