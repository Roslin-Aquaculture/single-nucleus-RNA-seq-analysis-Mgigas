###############################################################################
#  CLUSTER 1 SUB-CLUSTERING, WITH AND WITHOUT CORRECTION FOR SAMPLE OF ORIGIN
#
#  Dewari et al., single-nucleus atlas of the Pacific oyster during OsHV-1
#  infection. Reviewer request: re-analyse cluster 1 over infection time, and
#  cluster it on its own.
#
#  ASCII only throughout.
#
#  This script computes and writes tables, and saves the two Seurat objects.
#  The supplementary figure is drawn by reclustering_cluster1_figure.R, which
#  reads those objects and recomputes nothing - so styling can be iterated
#  without re-running the clustering.
#
#  WHAT THIS DESIGN CAN AND CANNOT ANSWER
#  Each infection stage is a single animal, so sample and timepoint are the
#  same variable. Correcting for sample therefore removes any genuine
#  infection-stage effect along with the individual effect: the corrected run
#  bounds residual structure, it does not test for infection. The uncorrected
#  run is where the question is actually answered.
#
#  Two comparisons bear on this. The 6 hpi animals are the strict case:
#  identical treatment, identical timepoint, two animals. The two controls are
#  the near case: one unchallenged, one mock-challenged with uninfected
#  homogenate, so they differ in handling but neither received virus, and their
#  separation therefore cannot be an infection effect. Reported as two separate
#  comparisons rather than pooled, because they support slightly different
#  claims.
#
#  The analysis is run twice on the same nuclei:
#    "U"  uncorrected  -- structure driven by sample is left in place
#    "H"  Harmony      -- sample of origin regressed out before clustering
#
#  EVERY NUMBER IN THE RESPONSE TO REVIEWERS IS COMPUTED AND WRITTEN HERE:
#    between-animal vs between-stage DE        -> between_animal_vs_between_stage_de
#    N subclusters, M of them >90% one animal  -> subcluster_purity
#    where each animal's nuclei went           -> animal_destination
#    expected single-animal contribution       -> expected_share_per_animal
#    purity at every resolution tested         -> resolution_sweep
#    detection limit (smallest subcluster)     -> printed in the summary
#    viral-positive nuclei of the late animal  -> viral_focus_late_animal
###############################################################################

library(Seurat)
library(tidyverse)

OUT <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/cluster1_subclustering_final"
DIR_FULL <- file.path(OUT, "full")
dir.create(DIR_FULL, showWarnings = FALSE, recursive = TRUE)

banner <- function(x) cat("\n\n", strrep("=", 74), "\n  ", x, "\n",
                          strrep("=", 74), "\n", sep = "")

# ---------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------

SEU_RDS     <- "/home/pdewari/Documents/parse_2025/seurat_2025/seu_obj_umap_18d_6r_3kRes.rds"
TARGET      <- "Cluster 1"
N_HVG       <- 2000
N_DIMS      <- 15
RES_FINAL   <- 0.4
RESOLUTIONS <- c(0.2, 0.3, 0.4, 0.6, 0.8)
DROP_SAMPLE <- "24-hpiA"        # excluded from the atlas
SEED        <- 42

## Pairs of animals whose separation cannot be an infection effect.
##   uninfected: neither received virus (unchallenged vs mock-challenge)
##   6 hpi     : identical treatment at the same timepoint - the strict case
PAIRS_MATCHED <- list(
  "uninfected" = c("Uninfected", "Homogenate"),
  "6 hpi"      = c("6-hpiA", "6-hpiD")
)
VIRAL_ANIMAL  <- "72-hpiJ"

## Thresholds matched to the DE settings used elsewhere in the paper.
DE_MIN_PCT <- 0.1
DE_LOGFC   <- 0.25
DE_PADJ    <- 0.05

CLUSTER_IDS <- c(
  "Cluster 0", "Cluster 1",
  "Gill ciliary cells", "Hepatopancreas cells",
  "Gill neuroepithelial cells", "Gill cell type 1",
  "Hyalinocytes", "Haemocyte cell type 1",
  "Mantle cell type 1", "Cluster 9",
  "Vesicular haemocytes", "Immature haemocytes",
  "Macrophage-like cells", "Adductor muscle cells",
  "Mantle cell type 2", "Mantle epithelial cells",
  "Gill cell type 2", "Small granule cells"
)

STAGES_INF <- c("6hpi", "24hpi", "72hpi", "96hpi")
RUN_U <- "Uncorrected"
RUN_H <- "Sample-corrected (Harmony)"

set.seed(SEED)

###############################################################################
banner("Subset cluster 1")
###############################################################################

stopifnot(file.exists(SEU_RDS))
seu <- readRDS(SEU_RDS)

names(CLUSTER_IDS) <- levels(seu)
seu <- RenameIdents(seu, CLUSTER_IDS)
seu$cell_type <- Idents(seu)

# 24-hpiA is excluded from the atlas, so it is dropped by name rather than by
# deriving a second condition variable that duplicates `sample`.
seu <- subset(seu, subset = sample != DROP_SAMPLE)
seu$sample <- droplevels(seu$sample)

seu$condition <- factor(case_when(
  seu$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu$sample %in% c("6-hpiA", "6-hpiD")         ~ "6hpi",
  seu$sample == "24-hpiJ"                       ~ "24hpi",
  seu$sample == "72-hpiJ"                       ~ "72hpi",
  seu$sample == "96-hpiE"                       ~ "96hpi"
), levels = c("control", STAGES_INF))
stopifnot(!anyNA(seu$condition))

sub0 <- subset(seu, idents = TARGET)
rm(seu); gc()
try(sub0 <- JoinLayers(sub0), silent = TRUE)

cat(sprintf("\n  Nuclei in %s: %s\n", TARGET, format(ncol(sub0), big.mark = ",")))

## Each animal's share of the cluster: the expectation any subcluster should
## meet if nuclei were distributed at random, and the power statement.
expected_share <- sub0@meta.data %>%
  count(sample, name = "n_nuclei") %>%
  mutate(pct_of_cluster1 = round(100 * n_nuclei / sum(n_nuclei), 1)) %>%
  arrange(desc(n_nuclei))
print(as.data.frame(expected_share))
write_csv(expected_share, file.path(DIR_FULL, "expected_share_per_animal.csv"))

exp_lo <- min(expected_share$pct_of_cluster1)
exp_hi <- max(expected_share$pct_of_cluster1)
cat(sprintf("\n  Expected contribution of any one animal: %.0f-%.0f%%\n", exp_lo, exp_hi))
cat(sprintf("  Nuclei per animal: %d-%d\n",
            min(expected_share$n_nuclei), max(expected_share$n_nuclei)))

## Is the animal effect biological or a sublibrary batch effect? If each animal
## spans more than one well, the effect is the animal.
if ("orig.ident" %in% colnames(sub0@meta.data)) {
  cat("\n  Animal by sequencing well:\n")
  print(table(sub0$sample, sub0$orig.ident))
  wells <- sub0@meta.data %>% group_by(sample) %>%
    summarise(n_wells = n_distinct(orig.ident), .groups = "drop")
  write_csv(wells, file.path(DIR_FULL, "wells_per_animal.csv"))
  cat(sprintf("  Wells per animal: %d-%d\n", min(wells$n_wells), max(wells$n_wells)))
}

###############################################################################
banner("Between-animal variation at fixed treatment")
###############################################################################

## The decisive comparison, and the one that needs no correction. Animals within
## each pair cannot differ by infection, so differences between them are
## individual variation. Set against control-versus-stage differences, this
## measures the nuisance rather than describing it.
##
## Note these are nucleus-level tests: with hundreds of nuclei per group, any
## comparison at this depth returns several hundred genes. The comparison is
## internally calibrated because all six were made the same way on the same
## nuclei.

de_count <- function(obj, group_col, a, b) {
  Idents(obj) <- group_col
  if (!all(c(a, b) %in% as.character(unique(Idents(obj))))) return(NA_integer_)
  nrow(FindMarkers(obj, ident.1 = a, ident.2 = b,
                   min.pct = DE_MIN_PCT, logfc.threshold = DE_LOGFC,
                   verbose = FALSE) %>% filter(p_val_adj < DE_PADJ))
}

de_compare <- bind_rows(
  map_dfr(names(PAIRS_MATCHED), function(nm) {
    p <- PAIRS_MATCHED[[nm]]
    tibble(contrast_type = "between animals, no infection difference",
           comparison = paste(p[1], "vs", p[2]),
           n_de = de_count(sub0, "sample", p[1], p[2]))
  }),
  map_dfr(STAGES_INF, function(s)
    tibble(contrast_type = "between stages",
           comparison = paste("control vs", s),
           n_de = de_count(sub0, "condition", "control", s)))
)
print(as.data.frame(de_compare), row.names = FALSE)
write_csv(de_compare, file.path(DIR_FULL, "between_animal_vs_between_stage_de.csv"))

###############################################################################
banner("Re-cluster, twice")
###############################################################################

## Variable features are re-selected inside the subset: genes that vary across
## the whole atlas carry information about differences between cell types, not
## about structure within one of them.

run_subclustering <- function(obj, use_harmony, prefix) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = N_HVG, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, seed.use = SEED, verbose = FALSE)

  red <- "pca"
  if (use_harmony) {
    if (!requireNamespace("harmony", quietly = TRUE))
      stop("harmony is not installed but USE_HARMONY was requested.")
    set.seed(SEED)
    obj <- harmony::RunHarmony(obj, group.by.vars = "sample",
                               reduction.use = "pca", verbose = FALSE)
    red <- "harmony"
  }

  obj <- FindNeighbors(obj, reduction = red, dims = 1:N_DIMS, verbose = FALSE)

  ## Resolution sweep reporting purity, not just cluster counts. A null at one
  ## resolution is weak; a null at every resolution tested is the claim.
  sweep <- map_dfr(RESOLUTIONS, function(r) {
    s  <- FindClusters(obj, resolution = r, random.seed = SEED, verbose = FALSE)
    md <- s@meta.data
    pur <- md %>% count(seurat_clusters, sample, name = "n") %>%
      group_by(seurat_clusters) %>%
      summarise(n_nuclei = sum(n), max_animal = 100 * max(n) / sum(n),
                .groups = "drop")
    stg <- md %>% count(seurat_clusters, condition, name = "n") %>%
      group_by(seurat_clusters) %>%
      summarise(max_stage = 100 * max(n) / sum(n), .groups = "drop")
    tibble(run = prefix, resolution = r,
           n_subclusters         = nrow(pur),
           smallest_subcluster   = min(pur$n_nuclei),
           max_single_animal_pct = round(max(pur$max_animal), 1),
           max_single_stage_pct  = round(max(stg$max_stage),  1))
  })

  obj <- FindClusters(obj, resolution = RES_FINAL, random.seed = SEED, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = red, dims = 1:N_DIMS,
                 seed.use = SEED, verbose = FALSE)
  obj$subcluster <- factor(paste0(prefix, obj$seurat_clusters))

  list(obj = obj, sweep = sweep, reduction = red)
}

uncorr <- run_subclustering(sub0, use_harmony = FALSE, prefix = "U")
corr   <- run_subclustering(sub0, use_harmony = TRUE,  prefix = "H")

sweep <- bind_rows(uncorr$sweep, corr$sweep)
print(as.data.frame(sweep))
write_csv(sweep, file.path(DIR_FULL, "resolution_sweep.csv"))

cat(sprintf("\n  Uncorrected: %d subclusters\n  Harmony    : %d subclusters\n",
            nlevels(uncorr$obj$subcluster), nlevels(corr$obj$subcluster)))

###############################################################################
banner("Composition and QC")
###############################################################################

composition <- function(obj, label) {
  md <- obj@meta.data
  qc <- md %>% group_by(subcluster) %>%
    summarise(n_nuclei = n(),
              median_nFeature = round(median(nFeature_RNA)), .groups = "drop")
  pct <- md %>%
    count(subcluster, sample, name = "n") %>%
    complete(subcluster, sample, fill = list(n = 0)) %>%
    group_by(subcluster) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    select(subcluster, sample, pct) %>%
    pivot_wider(names_from = sample, values_from = pct)
  qc %>% left_join(pct, by = "subcluster") %>%
    mutate(run = label, .before = 1) %>%
    mutate(across(where(is.numeric), ~round(.x, 1)))
}

comp <- bind_rows(composition(uncorr$obj, RUN_U),
                  composition(corr$obj,   RUN_H))

cat("\n  Per cent of each subcluster's nuclei contributed by each animal:\n\n")
print(as.data.frame(comp))
write_csv(comp, file.path(DIR_FULL, "subcluster_composition.csv"))

purity <- comp %>%
  rowwise() %>%
  mutate(max_pct_one_animal =
           max(c_across(-c(run, subcluster, n_nuclei, median_nFeature)))) %>%
  ungroup() %>%
  select(run, subcluster, n_nuclei, max_pct_one_animal)
print(as.data.frame(purity))
write_csv(purity, file.path(DIR_FULL, "subcluster_purity.csv"))

for (r in unique(purity$run)) {
  p  <- purity %>% filter(run == r)
  hi <- p %>% filter(max_pct_one_animal > 90)
  cat(sprintf("\n  %-28s %d of %d subclusters >90%% one animal | smallest %d nuclei\n",
              r, nrow(hi), nrow(p), min(p$n_nuclei)))
  if (nrow(hi)) cat(sprintf("     those %d span %.0f-%.0f%%\n", nrow(hi),
                            min(hi$max_pct_one_animal), max(hi$max_pct_one_animal)))
}

###############################################################################
banner("Where each animal's nuclei went")
###############################################################################

animal_destination <- function(obj, label) {
  obj@meta.data %>%
    count(sample, subcluster, name = "n") %>%
    group_by(sample) %>%
    mutate(pct_of_this_animal = 100 * n / sum(n)) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(purity %>% filter(run == label) %>%
                select(subcluster, subcluster_purity = max_pct_one_animal),
              by = "subcluster") %>%
    mutate(run = label, .before = 1) %>%
    mutate(across(where(is.numeric), ~round(.x, 1)))
}

dest <- bind_rows(animal_destination(uncorr$obj, RUN_U),
                  animal_destination(corr$obj,   RUN_H))
cat("\n  Dominant subcluster of each animal:\n\n")
print(as.data.frame(dest))
write_csv(dest, file.path(DIR_FULL, "animal_destination.csv"))

###############################################################################
banner("Viral transcripts per subcluster")
###############################################################################

viral_of <- function(obj) {
  vg <- grep("^ORF", rownames(obj), value = TRUE)
  if (!length(vg)) return(NULL)
  as.numeric(Matrix::colSums(
    GetAssayData(obj, assay = "RNA", layer = "counts")[vg, , drop = FALSE]))
}

viral_table <- function(obj, label) {
  v <- viral_of(obj); if (is.null(v)) return(NULL)
  obj@meta.data %>%
    mutate(viral = v) %>%
    group_by(subcluster, condition) %>%
    summarise(n_nuclei = n(),
              viral_umi_total = sum(viral),
              n_pos1 = sum(viral >= 1),
              n_pos2 = sum(viral >= 2), .groups = "drop") %>%
    mutate(run = label, .before = 1)
}

# Zero-viral subclusters are retained: dropping them removes the denominator.
vir <- bind_rows(viral_table(uncorr$obj, RUN_U), viral_table(corr$obj, RUN_H))
print(as.data.frame(vir))
write_csv(vir, file.path(DIR_FULL, "subcluster_viral.csv"))

## Counts, not percentages. The late-timepoint animal carries only a handful of
## viral-positive nuclei, and a percentage would imply precision they cannot
## support.
viral_focus <- function(obj, label, animal = VIRAL_ANIMAL) {
  v <- viral_of(obj); if (is.null(v)) return(NULL)
  md <- obj@meta.data %>% mutate(viral = v) %>% filter(sample == animal)
  if (!nrow(md) || sum(md$viral >= 1) == 0) return(NULL)
  md %>% group_by(subcluster) %>%
    summarise(n_nuclei = n(), n_viral_pos = sum(viral >= 1), .groups = "drop") %>%
    mutate(run = label, animal = animal,
           pct_of_animal_nuclei = round(100 * n_nuclei / sum(n_nuclei), 1),
           total_viral_pos      = sum(n_viral_pos), .before = 1) %>%
    arrange(desc(pct_of_animal_nuclei))
}

vf <- bind_rows(viral_focus(uncorr$obj, RUN_U), viral_focus(corr$obj, RUN_H))
if (!is.null(vf) && nrow(vf)) {
  cat(sprintf("\n  Where %s's nuclei and its viral-positive nuclei sit:\n\n", VIRAL_ANIMAL))
  print(as.data.frame(vf))
  write_csv(vf, file.path(DIR_FULL, "viral_focus_late_animal.csv"))
}

###############################################################################
banner("Markers")
###############################################################################

## Deposited with the code rather than placed in the supplementary material:
## no subcluster is interpreted in the manuscript, so these are reference
## material for a reader who wants to look, not evidence for a claim.

markers_for <- function(obj, label, min_cells = 20) {
  Idents(obj) <- "subcluster"
  keep <- names(which(table(obj$subcluster) >= min_cells))
  FindAllMarkers(subset(obj, idents = keep), only.pos = TRUE,
                 min.pct = 0.20, logfc.threshold = DE_LOGFC, verbose = FALSE) %>%
    mutate(run = label, .before = 1)
}

mk <- bind_rows(markers_for(uncorr$obj, RUN_U), markers_for(corr$obj, RUN_H))
write_csv(mk, file.path(DIR_FULL, "subcluster_markers_all.csv"))
cat(sprintf("\n  Markers written: %s rows\n", format(nrow(mk), big.mark = ",")))

###############################################################################
banner("Save objects for the figure script")
###############################################################################

saveRDS(uncorr$obj, file.path(DIR_FULL, "cluster1_uncorrected.rds"))
saveRDS(corr$obj,   file.path(DIR_FULL, "cluster1_harmony.rds"))
cat("\n  cluster1_uncorrected.rds and cluster1_harmony.rds written.\n")
cat("  Run reclustering_cluster1_figure.R next for the supplementary figure.\n")

###############################################################################
banner("NUMBERS FOR THE RESPONSE TO REVIEWERS")
###############################################################################

pu <- purity %>% filter(run == RUN_U)
ph <- purity %>% filter(run == RUN_H)
sw <- sweep  %>% filter(run == "H")

cat("\n  1. DE genes, between animals with no infection difference vs between stages:\n")
print(as.data.frame(de_compare), row.names = FALSE)

cat(sprintf("\n  2. Uncorrected: %d subclusters; %d are >90%% one animal\n",
            nrow(pu), sum(pu$max_pct_one_animal > 90)))
if (any(pu$max_pct_one_animal > 90))
  cat(sprintf("     those span %.0f-%.0f%%\n",
              min(pu$max_pct_one_animal[pu$max_pct_one_animal > 90]),
              max(pu$max_pct_one_animal[pu$max_pct_one_animal > 90])))

for (nm in names(PAIRS_MATCHED)) {
  d <- dest %>% filter(run == RUN_U, sample %in% PAIRS_MATCHED[[nm]])
  if (nrow(d) == 2)
    cat(sprintf("\n  3. %s pair: %s (%.0f%%) and %s (%.0f%%) - %s subclusters\n",
                nm, d$sample[1], d$subcluster_purity[1],
                d$sample[2], d$subcluster_purity[2],
                ifelse(d$subcluster[1] != d$subcluster[2], "different", "SAME")))
}

cat(sprintf("\n  4. Corrected: %d subclusters; highest single-animal share %.0f%%, expected %.0f-%.0f%%\n",
            nrow(ph), max(ph$max_pct_one_animal), exp_lo, exp_hi))

cat(sprintf("\n  5. Corrected run across resolutions %s: highest single-stage share %.0f-%.0f%%\n",
            paste(RESOLUTIONS, collapse = ", "),
            min(sw$max_single_stage_pct), max(sw$max_single_stage_pct)))
cat("     (expected share of each stage:)\n")
print(as.data.frame(sub0@meta.data %>% count(condition) %>%
                      mutate(pct = round(100 * n / sum(n), 1))), row.names = FALSE)

cat(sprintf("\n  6. Detection limit: smallest subcluster recovered %d nuclei (%.0f%% of cluster 1)\n",
            min(pu$n_nuclei), 100 * min(pu$n_nuclei) / ncol(sub0)))
cat(sprintf("     Nuclei per animal: %d-%d\n",
            min(expected_share$n_nuclei), max(expected_share$n_nuclei)))

if (!is.null(vf) && nrow(vf)) {
  top <- vf %>% filter(run == RUN_U) %>% slice_max(pct_of_animal_nuclei, n = 1)
  cat(sprintf("\n  7. %s: %d of %d viral-positive nuclei fell in %s, which holds %.0f%% of its nuclei\n",
              top$animal, top$n_viral_pos, top$total_viral_pos,
              top$subcluster, top$pct_of_animal_nuclei))
}

cat("\n  Written to ", DIR_FULL, ":\n", sep = "")
cat("    between_animal_vs_between_stage_de.csv, expected_share_per_animal.csv,\n")
cat("    wells_per_animal.csv, resolution_sweep.csv, subcluster_composition.csv,\n")
cat("    subcluster_purity.csv, animal_destination.csv, subcluster_viral.csv,\n")
cat("    viral_focus_late_animal.csv, subcluster_markers_all.csv,\n")
cat("    and both Seurat objects\n\n")

writeLines(capture.output(sessionInfo()), file.path(DIR_FULL, "sessionInfo.txt"))
print(sessionInfo())
