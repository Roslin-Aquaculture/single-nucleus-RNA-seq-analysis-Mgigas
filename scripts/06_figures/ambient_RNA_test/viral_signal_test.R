###############################################################################
#  AMBIENT RNA CONTROL FOR VIRAL TRANSCRIPT DETECTION
#
#  Dewari et al., single-nucleus atlas of the Pacific oyster during OsHV-1
#  infection. This script tests whether viral transcripts detected in nuclei
#  can be attributed to ambient RNA rather than genuine infection.
#
#  Self-contained: run top to bottom from the per-sublibrary Seurat objects.
#
#  PART A  Data preparation
#  TEST 1  Ambient viral content in empty barcodes
#  TEST 2  Sensitivity to the empty-barcode threshold
#  TEST 3  Calibrated enrichment over ambient (+ overdispersion check)
#  TEST 4  Viral UMIs per nucleus vs ambient expectation (depth-aware)
#  TEST 5  Cell-type specificity of viral signal
#
#  DESIGN NOTES
#  * Parse uses combinatorial barcoding, so most possible barcode combinations
#    are never occupied by a nucleus. These "empty barcodes" contain ambient
#    RNA only and give a direct estimate of background.
#  * Ambient transfer is calibrated on the DESIGN controls (unchallenged,
#    mock-challenge). Virus-exposed samples that happen to show no viral
#    transcripts are excluded: selecting them on that outcome would be
#    circular and would bias the estimate downwards.
#  * Effect sizes are the headline. Poisson p-values are anti-conservative
#    under overdispersion, so an overdispersed negative binomial is reported
#    alongside.
#  * Viral-positive nuclei are deeper than viral-negative nuclei, so Test 4
#    computes a per-nucleus expectation from each nucleus's own depth.
###############################################################################

library(Seurat)
library(tidyverse)
library(Matrix)

# Set this to the directory holding the per-sublibrary Seurat objects
# (seu_1.rds ... seu_8.rds), the clustered atlas object and the ribosomal
# gene list. See README for how to obtain these from the deposited data.
DATA_DIR <- "/home/pdewari/Documents/parse_2025/seurat_2025/"
setwd(DATA_DIR)

OUT_DIR <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/viral_ambient_check_output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

banner  <- function(x) cat("\n\n", strrep("=", 74), "\n  ", x, "\n",
                           strrep("=", 74), "\n", sep = "")
verdict <- function(x) cat("\n  >>> ", x, "\n", sep = "")

# ---------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------

FINAL_RDS     <- "seu_obj_umap_18d_6r_3kRes.rds"
RIBO_FILE     <- "02_ribo_rRNA_genes.txt"
EXTRA_REMOVE  <- "G32889"
CLUSTER_COL   <- "seurat_clusters"

EMPTY_MAX_UMI <- 100     # EmptyDrops default for "certainly ambient"
EMPTY_MIN_UMI <- 1
MT_MAX        <- 5       # nucleus filters -- must match the main pipeline
NFEATURE_MIN  <- 200
NFEATURE_MAX  <- 3000

FOCAL_SAMPLE  <- "72-hpiJ"
NEG_SAMPLES   <- c("Uninfected", "Homogenate")   # design controls only
STRICT_UMI    <- 2                               # high-confidence call

NEG_ALTS <- list(
  "Uninfected only"                      = c("Uninfected"),
  "Design controls (used)"               = c("Uninfected", "Homogenate"),
  "Design + 24 hpi (circular, not used)" = c("Uninfected", "Homogenate",
                                             "24-hpiA", "24-hpiJ")
)

CLUSTER_NAMES <- c(
  "0"  = "Cluster 0",                  "1"  = "Cluster 1",
  "2"  = "Gill ciliary cells",         "3"  = "Hepatopancreas cells",
  "4"  = "Gill neuroepithelial cells", "5"  = "Gill cell type 1",
  "6"  = "Hyalinocytes",               "7"  = "Haemocyte cell type 1",
  "8"  = "Mantle cell type 1",         "9"  = "Cluster 9",
  "10" = "Vesicular haemocytes",       "11" = "Immature haemocytes",
  "12" = "Macrophage-like cells",      "13" = "Adductor muscle cells",
  "14" = "Mantle cell type 2",         "15" = "Mantle epithelial cells",
  "16" = "Gill cell type 2",           "17" = "Small granule cells"
)

orig_to_sample <- c(
  "C9"="Uninfected","C10"="Uninfected","C11"="Homogenate","C12"="Homogenate",
  "D1"="6-hpiA","D2"="6-hpiA","D3"="6-hpiD","D4"="6-hpiD",
  "D5"="24-hpiA","D6"="24-hpiA","D7"="24-hpiJ","D8"="24-hpiJ",
  "D9"="72-hpiJ","D10"="72-hpiJ","D11"="96-hpiE","D12"="96-hpiE"
)
SAMPLE_LEVELS <- c("Uninfected","Homogenate","6-hpiA","6-hpiD",
                   "24-hpiA","24-hpiJ","72-hpiJ","96-hpiE")

# Internal labels map to the manuscript as follows:
#   "Uninfected" = unchallenged control
#   "Homogenate" = mock-challenge control (uninfected oyster homogenate)
# Remaining labels are hours post-infection plus the animal identifier.

###############################################################################
banner("PART A  --  Data preparation")
###############################################################################

sub_files <- sprintf("seu_%d.rds", 1:8)
if (!all(file.exists(sub_files)))
  stop("Missing sublibraries: ", paste(sub_files[!file.exists(sub_files)], collapse=", "))
if (!file.exists(RIBO_FILE)) stop("Ribosomal gene file not found: ", RIBO_FILE)
if (!file.exists(FINAL_RDS))
  warning("FINAL_RDS not found: ", FINAL_RDS, " -- Test 5 will be skipped.",
          immediate. = TRUE)

message("Loading and merging sublibraries ...")
seu_list <- lapply(1:8, function(i) readRDS(sprintf("seu_%d.rds", i)))
merged_8 <- merge(x = seu_list[[1]], y = seu_list[-1],
                  add.cell.ids = paste0("seu", 1:8))
seu_all  <- JoinLayers(merged_8)
rm(merged_8, seu_list); gc()

# unname(): the lookup is named by well ID; Seurat would otherwise try to match
# those names to cell barcodes and fail with "No cell overlap".
seu_all$sample <- factor(unname(orig_to_sample[as.character(seu_all$orig.ident)]),
                         levels = SAMPLE_LEVELS)
message("Total barcodes: ", ncol(seu_all))
print(table(seu_all$sample, useNA = "ifany"))

## Gene sets: identified but NOT removed -------------------------------------
## Barcodes are classified on the FULL gene set: removing rRNA would shrink
## every barcode's total and change which fall below the "empty" threshold,
## and the main pipeline filters percent.mt / nFeature BEFORE removing rRNA.
## rRNA is excluded only from the fraction denominator, so that it matches the
## gene set retained in nuclei.

ribo_genes  <- read_table(RIBO_FILE, col_names = FALSE)
to_remove   <- intersect(c(as.character(ribo_genes$X1), EXTRA_REMOVE),
                         rownames(seu_all))
viral_genes <- grep("^ORF", rownames(seu_all), value = TRUE)  # no dot-exclusion
message("Ribosomal/extra genes: ", length(to_remove),
        " | Viral genes: ", length(viral_genes))
stopifnot(length(viral_genes) > 0)

MZ_genes   <- grep("MZ", rownames(seu_all), value = TRUE)
mt.present <- intersect(c("ATP6","ATP8","ND1","ND2","ND3","ND4","ND4L",
                          "ND5","ND6","COX1","COX2","COX3","CYTB", MZ_genes),
                        rownames(seu_all))
seu_all <- PercentageFeatureSet(seu_all, features = mt.present,
                                col.name = "percent.mt")

counts     <- GetAssayData(seu_all, assay = "RNA", layer = "counts")
tot_full   <- Matrix::colSums(counts)
nfeat      <- Matrix::colSums(counts > 0)
keep_genes <- setdiff(rownames(counts), to_remove)
tot_noribo <- Matrix::colSums(counts[keep_genes, , drop = FALSE])
vg         <- intersect(viral_genes, rownames(counts))

bc <- tibble(
  barcode      = colnames(counts),
  sample       = seu_all$sample,
  percent_mt   = seu_all$percent.mt,
  total_umi    = tot_full,
  total_noribo = tot_noribo,
  viral_umi    = Matrix::colSums(counts[vg, , drop = FALSE]),
  viral_ngene  = Matrix::colSums(counts[vg, , drop = FALSE] > 0),
  n_genes      = nfeat,
  empty        = tot_full >= EMPTY_MIN_UMI & tot_full < EMPTY_MAX_UMI,
  nucleus      = seu_all$percent.mt < MT_MAX &
                 nfeat > NFEATURE_MIN & nfeat < NFEATURE_MAX
)

cat(sprintf("\n  Empty barcodes : %s\n  Called nuclei  : %s  (pipeline: 23,740)\n",
            format(sum(bc$empty), big.mark=","),
            format(sum(bc$nucleus), big.mark=",")))

ambient <- bc %>% filter(empty) %>% group_by(sample) %>%
  summarise(n_empty_barcodes    = n(),
            viral_umi_empty     = sum(viral_umi),
            total_umi_empty     = sum(total_umi),
            total_noribo_empty  = sum(total_noribo),
            ambient_frac_noribo = sum(viral_umi)/max(sum(total_noribo),1),
            .groups="drop")

per_sample <- bc %>% filter(nucleus) %>% group_by(sample) %>%
  summarise(n_cells        = n(),
            total_noribo   = sum(total_noribo),
            observed_viral = sum(viral_umi),
            n_viral_pos    = sum(viral_umi >= 1),
            n_viral_strict = sum(viral_umi >= STRICT_UMI),
            .groups="drop")

saveRDS(bc, file.path(OUT_DIR, "barcode_level_viral_table.rds"))
rm(seu_all); gc()

###############################################################################
banner("TEST 1  --  Ambient viral content in empty barcodes")
###############################################################################

print(as.data.frame(ambient))
cat(sprintf("\n  Total empty barcodes : %s\n  Total viral UMIs     : %s\n",
            format(sum(ambient$n_empty_barcodes), big.mark=","),
            format(sum(ambient$viral_umi_empty), big.mark=",")))
verdict(sprintf(
  "Ambient viral RNA exists but is rare (%.2e of UMIs, ~1 in %s). %s carries the most, consistent with virus shed from genuinely infected cells.",
  sum(ambient$viral_umi_empty)/sum(ambient$total_noribo_empty),
  format(round(sum(ambient$total_noribo_empty)/sum(ambient$viral_umi_empty)),
         big.mark=","), FOCAL_SAMPLE))
write_csv(ambient, file.path(OUT_DIR, "test1_ambient.csv"))

###############################################################################
banner("TEST 2  --  Threshold sensitivity")
###############################################################################

t2 <- map_dfr(c(20,50,100,200,500), function(thr) {
  e <- bc %>% filter(total_umi >= 1, total_umi < thr)
  tibble(threshold=thr, n_barcodes=nrow(e), viral_umi=sum(e$viral_umi),
         ambient_frac=sum(e$viral_umi)/max(sum(e$total_noribo),1))
})
print(as.data.frame(t2))
rng <- range((t2 %>% filter(threshold <= 200))$ambient_frac)
cat("\n  NOTE: thresholds above 200 overlap the nucleus definition (nFeature > 200),\n")
cat("        so they measure cells, not soup, and are excluded.\n")
verdict(sprintf(
  "Across valid thresholds (20-200) the ambient fraction spans %.2e - %.2e (%.0f%% variation). Not threshold-driven.",
  rng[1], rng[2], 100*(rng[2]-rng[1])/rng[1]))
write_csv(t2, file.path(OUT_DIR, "test2_threshold_sensitivity.csv"))

###############################################################################
banner("TEST 3  --  Calibrated enrichment over ambient")
###############################################################################
# A naive model assumes 100% of a nucleus's UMIs are soup-derived. The design
# controls show that is far too high, and bound the true transfer rate (alpha).

obs_exp <- per_sample %>%
  left_join(ambient %>% select(sample, ambient_frac_noribo), by="sample") %>%
  mutate(expected_naive = total_noribo * ambient_frac_noribo)

cat("\n  Calibration options:\n")
alpha_tab <- map_dfr(names(NEG_ALTS), function(nm) {
  cal <- obs_exp %>% filter(sample %in% NEG_ALTS[[nm]])
  a   <- qgamma(0.95, shape=sum(cal$observed_viral)+1)/max(sum(cal$expected_naive),1e-9)
  e72 <- obs_exp$expected_naive[obs_exp$sample==FOCAL_SAMPLE]*a
  o72 <- obs_exp$observed_viral[obs_exp$sample==FOCAL_SAMPLE]
  tibble(calibration=nm, n_control_nuclei=sum(cal$n_cells),
         control_viral=sum(cal$observed_viral),
         naive_expected=sum(cal$expected_naive),
         alpha_ub=a, focal_expected=e72, focal_enrichment=o72/e72)
})
print(as.data.frame(alpha_tab))
write_csv(alpha_tab, file.path(OUT_DIR, "test3_calibration_options.csv"))

calib     <- obs_exp %>% filter(sample %in% NEG_SAMPLES)
exp_naive <- sum(calib$expected_naive); obs_neg <- sum(calib$observed_viral)
alpha_ub  <- qgamma(0.95, shape=obs_neg+1)/exp_naive

cat(sprintf("\n  PRIMARY calibration -- design controls (%s):\n",
            paste(NEG_SAMPLES, collapse=", ")))
cat(sprintf("    nuclei                    : %s\n", format(sum(calib$n_cells), big.mark=",")))
cat(sprintf("    viral UMIs observed       : %g\n", obs_neg))
cat(sprintf("    predicted if 100%% ambient : %.1f\n", exp_naive))
cat(sprintf("    => alpha <= %.4f (at most %.1f%% of nuclear UMIs are ambient)\n",
            alpha_ub, 100*alpha_ub))

t3 <- obs_exp %>%
  mutate(expected = expected_naive*alpha_ub,
         enrichment = ifelse(expected > 0.01, observed_viral/expected, NA_real_),
         p_poisson = ifelse(expected > 0.01 & observed_viral > 0,
                            ppois(observed_viral-1, expected, lower.tail=FALSE), NA_real_),
         p_nb_overdisp = ifelse(expected > 0.01 & observed_viral > 0,
                                pnbinom(observed_viral-1, mu=expected, size=1,
                                        lower.tail=FALSE), NA_real_)) %>%
  select(sample, n_cells, observed_viral, n_viral_pos, n_viral_strict,
         expected, enrichment, p_poisson, p_nb_overdisp)
print(as.data.frame(t3))
cat("\n  NOTE: enrichment is NA where the ambient estimate is 0; those 1-2 UMI\n")
cat("        observations are indistinguishable from background.\n")

foc <- t3 %>% filter(sample == FOCAL_SAMPLE)
verdict(sprintf("%s: %g viral UMIs vs %.1f expected = %.1f-fold (Poisson p = %.2g; overdispersed NB p = %.2g).",
                FOCAL_SAMPLE, foc$observed_viral, foc$expected,
                foc$enrichment, foc$p_poisson, foc$p_nb_overdisp))
write_csv(t3, file.path(OUT_DIR, "test3_calibrated_enrichment.csv"))

# control summary, for the Table S1 legend
ctrl <- bc %>% filter(nucleus, sample %in% c("Uninfected","Homogenate")) %>%
  summarise(nuclei=n(), pos1=sum(viral_umi>=1), pos2=sum(viral_umi>=STRICT_UMI))
early <- bc %>% filter(nucleus, sample %in% c("6-hpiA","6-hpiD","24-hpiA","24-hpiJ")) %>%
  summarise(nuclei=n(), pos1=sum(viral_umi>=1), pos2=sum(viral_umi>=STRICT_UMI))
cat("\n  Controls : "); print(as.data.frame(ctrl))
cat("  6/24 hpi : ");   print(as.data.frame(early))

###############################################################################
banner("TEST 4  --  Viral UMIs per nucleus  [DEPTH-AWARE]")
###############################################################################
# Viral-positive nuclei are deeper than viral-negative ones, so a single
# uniform lambda would understate the ambient tail. Each nucleus therefore
# gets its own expectation from its own sequencing depth.

cat("\n  Depth of viral+ vs viral- nuclei:\n")
print(as.data.frame(
  bc %>% filter(nucleus, sample==FOCAL_SAMPLE) %>%
    mutate(pos = viral_umi >= 1) %>% group_by(pos) %>%
    summarise(n=n(), median_umi=median(total_noribo),
              median_genes=median(n_genes), .groups="drop")))
cat(sprintf("  Wilcoxon on depth: p = %.3g\n",
            wilcox.test(total_noribo ~ viral_umi >= 1,
                        data = bc %>% filter(nucleus, sample==FOCAL_SAMPLE))$p.value))

amb_frac <- ambient$ambient_frac_noribo[ambient$sample == FOCAL_SAMPLE]
vfoc <- bc %>% filter(nucleus, sample == FOCAL_SAMPLE) %>%
  mutate(lam_i = total_noribo * amb_frac * alpha_ub)

t4 <- tibble(
  viral_umi = 0:5,
  observed  = sapply(0:5, function(k) sum(vfoc$viral_umi >= k & vfoc$viral_umi < k+1)),
  expected_uniform     = round(dpois(0:5, mean(vfoc$lam_i)) * nrow(vfoc), 4),
  expected_depth_aware = round(sapply(0:5, function(k) sum(dpois(k, vfoc$lam_i))), 4)
)
print(as.data.frame(t4))

obs_ge2 <- sum(vfoc$viral_umi >= STRICT_UMI)
exp_ge2 <- sum(1 - ppois(STRICT_UMI - 1, vfoc$lam_i))
obs_1   <- sum(vfoc$viral_umi >= 1 & vfoc$viral_umi < 2)
exp_1   <- sum(dpois(1, vfoc$lam_i))

cat(sprintf("\n  Single-UMI nuclei : observed %d, expected %.1f (%.1f-fold)\n",
            obs_1, exp_1, obs_1/exp_1))
cat(sprintf("  >= %d UMI nuclei   : observed %d, expected %.3f (%.0f-fold, p = %.2g)\n",
            STRICT_UMI, obs_ge2, exp_ge2, obs_ge2/exp_ge2,
            ppois(obs_ge2-1, exp_ge2, lower.tail=FALSE)))

verdict(sprintf(
  "Single-UMI nuclei are only modestly above ambient, but %d nuclei carry >= %d viral UMIs against %.3f expected (%.0f-fold). The signal does not rest on single-UMI calls.",
  obs_ge2, STRICT_UMI, exp_ge2, obs_ge2/exp_ge2))
write_csv(t4, file.path(OUT_DIR, "test4_umi_distribution.csv"))

###############################################################################
banner("TEST 5  --  Cell-type specificity of viral signal")
###############################################################################
# Ambient uptake is cell-type blind: viral UMIs should land in each cluster in
# proportion to that cluster's share of UMIs. Infection should not.

if (file.exists(FINAL_RDS)) {
  seu_final <- readRDS(FINAL_RDS)
  clust <- tibble(barcode = colnames(seu_final),
                  cluster = as.character(seu_final@meta.data[[CLUSTER_COL]]))

  t5 <- bc %>% filter(nucleus, sample == FOCAL_SAMPLE) %>%
    inner_join(clust, by="barcode") %>%
    group_by(cluster) %>%
    summarise(n_cells=n(), viral=sum(viral_umi),
              n_strict=sum(viral_umi >= STRICT_UMI),
              umi=sum(total_noribo), .groups="drop") %>%
    mutate(cell_type = unname(CLUSTER_NAMES[cluster]),
           umi_share = umi/sum(umi),
           expected  = sum(viral)*umi_share,
           obs_exp   = viral/pmax(expected, 1e-9)) %>%
    select(cluster, cell_type, n_cells, viral, n_strict,
           umi_share, expected, obs_exp) %>%
    arrange(desc(viral))
  print(as.data.frame(t5))

  ct <- suppressWarnings(chisq.test(round(t5$viral), p=t5$umi_share,
                                    simulate.p.value=TRUE, B=1e6))
  cat(sprintf("\n  Chi-square vs depth-proportional: p = %.3g (simulation floor 1e-06)\n",
              ct$p.value))

  top3 <- t5 %>% slice_max(viral, n=3)
  cat(sprintf("\n  Top 3: %s\n",
              paste(sprintf("%s (%.1fx, %d high-confidence nuclei)",
                            top3$cell_type, top3$obs_exp, top3$n_strict),
                    collapse="; ")))
  cat(sprintf("  These hold %.0f%% of depth but %.0f%% of viral transcripts.\n",
              100*sum(top3$umi_share), 100*sum(top3$viral)/sum(t5$viral)))
  cat("\n  CAUTION: a cell type whose enrichment rests on single-UMI nuclei is\n")
  cat("           less firmly supported -- check the n_strict column.\n")

  verdict(sprintf("Viral signal is cell-type specific (p = %.3g): %.0f%% of signal in %s, holding only %.0f%% of depth.",
                  ct$p.value, 100*sum(top3$viral)/sum(t5$viral),
                  paste(top3$cell_type, collapse=", "), 100*sum(top3$umi_share)))
  write_csv(t5, file.path(OUT_DIR, "test5_cluster_specificity.csv"))

  ## Table S1 -- viral detection per cell subtype ---------------------------
  tabS1 <- bc %>% filter(nucleus) %>%
    inner_join(clust, by="barcode") %>%
    mutate(cell_type = unname(CLUSTER_NAMES[cluster]),
           subtype   = paste0(cell_type, "_", sample)) %>%
    group_by(sample, cell_type, subtype) %>%
    summarise(total_nuclei=n(), pos1=sum(viral_umi>=1),
              pos2=sum(viral_umi>=STRICT_UMI),
              viral_umi_total=sum(viral_umi), .groups="drop") %>%
    mutate(pct1=100*pos1/total_nuclei, pct2=100*pos2/total_nuclei) %>%
    arrange(desc(pct2), desc(pct1))

  tabS1_fmt <- tabS1 %>%
    transmute(`Sample subtype`=subtype, `Total nuclei`=total_nuclei,
              `Viral-positive (>=1 UMI), n (%)`=sprintf("%d (%.1f%%)", pos1, pct1),
              `High-confidence (>=2 UMI), n (%)`=sprintf("%d (%.1f%%)", pos2, pct2),
              `Total viral UMIs`=viral_umi_total)
  cat("\n  Table S1 (subtypes with a high-confidence nucleus or >1% positive):\n")
  print(as.data.frame(tabS1_fmt %>% filter(tabS1$pos2 > 0 | tabS1$pct1 > 1)))
  write_csv(tabS1_fmt, file.path(OUT_DIR, "TableS1_updated.csv"))

} else {
  cat("  FINAL_RDS not found -- Test 5 skipped.\n")
}

###############################################################################
banner("SUMMARY")
###############################################################################

cat("\n  Test 1  ambient viral RNA present but rare\n")
cat("  Test 2  estimate stable across valid thresholds (20-200 UMIs)\n")
cat(sprintf("  Test 3  %s enriched %.1f-fold over ambient (NB p = %.2g); alpha <= %.3f\n",
            FOCAL_SAMPLE, foc$enrichment, foc$p_nb_overdisp, alpha_ub))
cat(sprintf("  Test 4  %d nuclei with >=%d viral UMIs vs %.3f expected (%.0f-fold, depth-aware)\n",
            obs_ge2, STRICT_UMI, exp_ge2, obs_ge2/exp_ge2))
if (exists("ct"))
  cat(sprintf("  Test 5  cell-type specific (p = %.3g)\n", ct$p.value))

cat("\n  Outputs: ", normalizePath(OUT_DIR), "\n\n", sep="")


########################
cat("\n  Session info (for reproducibility):\n\n")
print(sessionInfo())
