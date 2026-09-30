###############################################################################
#  SUPPLEMENTARY TABLES FOR THE AMBIENT RNA CONTROL
#
#  Dewari et al., single-nucleus atlas of the Pacific oyster during OsHV-1
#  infection. Reviewer comment: a single viral UMI cannot separate infection
#  from ambient capture in a late-stage infection.
#
#  Reads the outputs of viral_signal_tests.R and produces only the tables cited
#  in the manuscript:
#
#    TableSX_a_ambient_and_viral_detection_per_sample.csv
#    TableSX_b_detection_threshold_vs_ambient.csv
#    TableSX_c_viral_transcripts_per_cell_type.csv
#    TableSX_d_viral_detection_per_cell_subtype.csv   (updates existing Table S1)
#
#  Everything else the analysis produces -- threshold sensitivity, the
#  alternative calibrations, the depth comparison, the raw UMI distribution --
#  stays in the analysis output directory and is deposited with the code. It
#  supports the Methods rather than any statement in the text.
###############################################################################

library(tidyverse)

setwd("/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/viral_ambient_check_output")
IN_DIR   <- "/home/pdewari/eggnog/results/enrichment/figures_20260923_27Sept/viral_ambient_check_output"
DIR_SUPP <- file.path(IN_DIR, "supplementary")
dir.create(DIR_SUPP, showWarnings = FALSE, recursive = TRUE)

FOCAL      <- "72-hpiJ"
STRICT_UMI <- 2
EARLY      <- c("Uninfected", "Homogenate", "6-hpiA", "6-hpiD",
                "24-hpiA", "24-hpiJ")

###############################################################################
# Inputs
###############################################################################

need <- c("barcode_level_viral_table.rds", "test1_ambient.csv",
          "test3_calibration_options.csv", "test3_calibrated_enrichment.csv",
          "test5_cluster_specificity.csv", "TableS1_updated.csv")
missing <- need[!file.exists(file.path(IN_DIR, need))]
if (length(missing))
  stop("Run viral_signal_tests.R first; missing: ", paste(missing, collapse = ", "))

bc      <- readRDS(file.path(IN_DIR, "barcode_level_viral_table.rds"))
ambient <- read_csv(file.path(IN_DIR, "test1_ambient.csv"), show_col_types = FALSE)
alpha   <- read_csv(file.path(IN_DIR, "test3_calibration_options.csv"), show_col_types = FALSE)
t3      <- read_csv(file.path(IN_DIR, "test3_calibrated_enrichment.csv"), show_col_types = FALSE)
t5      <- read_csv(file.path(IN_DIR, "test5_cluster_specificity.csv"), show_col_types = FALSE)
tabS1   <- read_csv(file.path(IN_DIR, "TableS1_updated.csv"), show_col_types = FALSE)

## The calibration actually used: design controls only. Selecting virus-exposed
## samples that happen to show no viral transcripts would be circular.
alpha_ub <- alpha$alpha_ub[alpha$calibration == "Design controls (used)"]
stopifnot(length(alpha_ub) == 1)

## Blank rather than NA or 0 where a quantity is undefined: "0-fold enrichment"
## and "NA" both invite a second reading, and NA would otherwise mean two
## different things in the same table.
blank <- function(x, keep) ifelse(keep, format(x, trim = TRUE), "—")

###############################################################################
# Table a: ambient and viral detection per sample
###############################################################################

## The ambient rate is given as "one viral UMI per N" rather than a fraction:
## 1 in 27,000 against 1 in 149,000 is easier to read than 3.7e-5 against
## 6.7e-6. Larger values therefore mean LESS ambient virus.

tab_a <- ambient %>%
  select(sample, n_empty_barcodes, viral_umi_empty,
         ambient_frac_noribo, total_noribo_empty) %>%
  left_join(t3 %>% select(sample, n_cells, observed_viral, n_viral_pos,
                          n_viral_strict, expected, enrichment),
            by = "sample") %>%
  transmute(
    Sample                             = as.character(sample),
    `Empty barcodes`                   = format(n_empty_barcodes, big.mark = ","),
    `Viral UMIs in empty barcodes`     = viral_umi_empty,
    `Ambient rate (1 viral UMI per N total UMIs)` =
      blank(round(1 / ambient_frac_noribo, -2), viral_umi_empty > 0),
    Nuclei                             = format(n_cells, big.mark = ","),
    `Viral UMIs in nuclei`             = observed_viral,
    `Viral-positive nuclei (>=1 UMI)`  = n_viral_pos,
    `High-confidence nuclei (>=2 UMI)` = n_viral_strict,
    `Viral UMIs expected from ambient` = blank(round(expected, 1), expected > 0.01),
    `Enrichment over ambient`          =
      blank(round(enrichment, 1), expected > 0.01 & observed_viral > 0))

totals <- tibble(
  Sample = "All samples",
  `Empty barcodes` = format(sum(ambient$n_empty_barcodes), big.mark = ","),
  `Viral UMIs in empty barcodes` = sum(ambient$viral_umi_empty),
  `Ambient rate (1 viral UMI per N total UMIs)` =
    format(round(sum(ambient$total_noribo_empty) /
                   sum(ambient$viral_umi_empty), -2), big.mark = ","),
  Nuclei = format(sum(t3$n_cells), big.mark = ","),
  `Viral UMIs in nuclei` = sum(t3$observed_viral),
  `Viral-positive nuclei (>=1 UMI)`  = sum(t3$n_viral_pos),
  `High-confidence nuclei (>=2 UMI)` = sum(t3$n_viral_strict),
  `Viral UMIs expected from ambient` = "—",
  `Enrichment over ambient`          = "—")

tab_a <- bind_rows(tab_a, totals)
cat("\n--- Table a: per sample ---\n")
print(as.data.frame(tab_a))
write_csv(tab_a, file.path(DIR_SUPP,
                           "TableSX_a_ambient_and_viral_detection_per_sample.csv"))

###############################################################################
# Table b: detection threshold against a depth-aware ambient expectation
###############################################################################

## This is a different quantity from the "expected" column in Table a, which is
## expected viral UMIs from pooled depth. Here the expectation is the number of
## NUCLEI that would carry that many viral UMIs by chance, computed for each
## nucleus from its own non-ribosomal depth. Viral-positive nuclei are deeper
## than viral-negative ones, so a pooled expectation would understate the
## background in precisely the nuclei that matter.

amb_frac <- ambient$ambient_frac_noribo[ambient$sample == FOCAL]

vfoc <- bc %>% filter(nucleus, sample == FOCAL) %>%
  mutate(lambda = total_noribo * amb_frac * alpha_ub)

ctrl <- bc %>% filter(nucleus, sample %in% EARLY)

obs_ge1 <- sum(vfoc$viral_umi >= 1)
exp_ge1 <- sum(1 - dpois(0, vfoc$lambda))
obs_eq1 <- sum(vfoc$viral_umi == 1)
exp_eq1 <- sum(dpois(1, vfoc$lambda))
obs_ge2 <- sum(vfoc$viral_umi >= STRICT_UMI)
exp_ge2 <- sum(1 - ppois(STRICT_UMI - 1, vfoc$lambda))

tab_b <- tibble(
  `Detection threshold` = c(
    sprintf("At least 1 viral UMI (%s)", FOCAL),
    sprintf("Exactly 1 viral UMI (%s)", FOCAL),
    sprintf("At least %d viral UMIs (%s)", STRICT_UMI, FOCAL),
    sprintf("At least %d viral UMIs (control and early samples)", STRICT_UMI)),
  `Nuclei examined` = format(c(nrow(vfoc), nrow(vfoc), nrow(vfoc), nrow(ctrl)),
                             big.mark = ","),
  `Nuclei observed` = c(obs_ge1, obs_eq1, obs_ge2,
                        sum(ctrl$viral_umi >= STRICT_UMI)),
  `Nuclei expected from ambient` = c(sprintf("%.1f", exp_ge1),
                                     sprintf("%.1f", exp_eq1),
                                     sprintf("%.2f", exp_ge2), "—"),
  `Fold over ambient` = c(sprintf("%.1f", obs_ge1 / exp_ge1),
                          sprintf("%.1f", obs_eq1 / exp_eq1),
                          sprintf("%.0f", obs_ge2 / exp_ge2), "—"),
  `Poisson P` = c(
    sprintf("%.2g", ppois(obs_ge1 - 1, exp_ge1, lower.tail = FALSE)),
    "—",
    sprintf("%.2g", ppois(obs_ge2 - 1, exp_ge2, lower.tail = FALSE)),
    "—"))

cat("\n--- Table b: detection threshold ---\n")
print(as.data.frame(tab_b))
write_csv(tab_b, file.path(DIR_SUPP,
                           "TableSX_b_detection_threshold_vs_ambient.csv"))

###############################################################################
# Table c: viral transcripts per cell type at the focal timepoint
###############################################################################

## Ambient uptake is cell-type agnostic, so under an ambient explanation viral
## transcripts should follow each cell type's share of sequencing depth. The
## high-confidence column is kept so that a cell type whose enrichment rests
## only on single-UMI calls is visible rather than implied.

tab_c <- t5 %>%
  transmute(
    `Cell type`                        = cell_type,
    Nuclei                             = n_cells,
    `Viral UMIs`                       = viral,
    `High-confidence nuclei (>=2 UMI)` = n_strict,
    `Share of sequencing depth (%)`    = round(100 * umi_share, 1),
    `Viral UMIs expected from depth`   = round(expected, 1),
    `Observed / expected`              = round(obs_exp, 2)) %>%
  arrange(desc(`Viral UMIs`))

cat("\n--- Table c: per cell type ---\n")
print(as.data.frame(tab_c))
write_csv(tab_c, file.path(DIR_SUPP,
                           "TableSX_c_viral_transcripts_per_cell_type.csv"))

###############################################################################
# Table d: viral detection per cell subtype (updates the existing Table S1)
###############################################################################

write_csv(tabS1, file.path(DIR_SUPP,
                           "TableSX_d_viral_detection_per_cell_subtype.csv"))

###############################################################################

cat("\n  Written to ", normalizePath(DIR_SUPP), ":\n", sep = "")
cat("    TableSX_a_ambient_and_viral_detection_per_sample.csv\n")
cat("    TableSX_b_detection_threshold_vs_ambient.csv\n")
cat("    TableSX_c_viral_transcripts_per_cell_type.csv\n")
cat("    TableSX_d_viral_detection_per_cell_subtype.csv\n\n")
cat(sprintf("  Headline: %d nuclei with >=%d viral UMIs at %s against %.2f expected (%.0f-fold)\n",
            obs_ge2, STRICT_UMI, FOCAL, exp_ge2, obs_ge2 / exp_ge2))
cat(sprintf("  and %d such nuclei among %s control and early-timepoint nuclei.\n\n",
            sum(ctrl$viral_umi >= STRICT_UMI), format(nrow(ctrl), big.mark = ",")))

