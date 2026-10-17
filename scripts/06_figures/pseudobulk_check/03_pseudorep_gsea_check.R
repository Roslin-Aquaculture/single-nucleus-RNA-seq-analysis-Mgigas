# =========================================================
# Section 5d: GSEA on the pseudoreplication diagnostic comparisons
# Environment: go-enrich (local machine)
# =========================================================
# Runs the same GSEA pipeline as 05c, but on the output of
# 08_pseudoreplication_checks.R rather than on the main control-vs-stage
# comparisons, then computes the summary statistics reported in Methods.
#
# RUN ORDER
#   1. Seurat session:     08_pseudoreplication_checks.R
#   2. go-enrich session:  source("00_reload_go_enrich_session.R")
#   3. go-enrich session:  source("05d_gsea_pseudorep_checks.R")  <- this file
#
# This file sources 05c_gsea_pseudorep_engine.R, a dedicated copy of the
# GSEA engine with the diagnostic paths fixed. Nothing in the main
# pipeline (04b, 05b, 05c, 06) is read or written.
#
# Requires in this session:
#   - reference tables from 01b (t2g_BP/MF/CC, term2name_go,
#     kegg_term2gene, kegg_names) - see 00_reload_go_enrich_session.R
# Requires on disk:
#   - 08_pseudoreplication_checks.R output
#
# BOTH REPLICATED TIMEPOINTS. 24-hpiA is back in the analysis, so 6 and
# 24 hpi each have two animals and each supports B, C and D. Stages are
# matched by PREFIX rather than by exact name, so adding a third pair to
# TIMEPOINTS in 08 needs no change here.
#
#   A_ctrl_vs_ctrl            two uninfected animals against one another
#   B_<tp>_between_animals    the two animals of one timepoint
#   C_ctrl_vs_<sample>        each animal against the pooled controls
#   D_pseudobulk_<tp>         animal-level profiles, 2 controls vs 2 infected
#
# READ B WITH THE qPCR TABLE OPEN. Host-normalised qPCR shows one animal
# of each early pair carried little or no virus (6-hpiD undetectable;
# 24-hpiA at 9.7e-7 relative to host against 0.313 in 24-hpiJ), so B mixes
# individual variation with a real infection difference and is an upper
# bound rather than a measurement. A, the two uninfected controls, is the
# only comparison where no infection difference is possible, and it is
# the one the pseudoreplication argument rests on.

library(readr); library(dplyr); library(clusterProfiler)

# =========================================================
# The GSEA engine is a dedicated copy reserved for these diagnostics, so
# the paths live there and cannot collide with the main analysis. The main
# 05c_gsea_across_stages.R is not touched by this script and needs no
# guards.
# =========================================================
source("05c_gsea_pseudorep_engine.R")   # sets de_dir, outdir; writes GSEA_tidy.tsv

# =========================================================
# Summary statistics reported in Methods
# =========================================================
tg <- read_tsv(file.path(outdir, "GSEA_tidy.tsv"), show_col_types = FALSE)

# terms matching the translational programme
TRANS <- "ribosom|translat|peptide biosynth|peptide metab|amide biosynth|tRNA|aminoacyl"

cat("\n=========================================================\n")
cat("PSEUDOREPLICATION DIAGNOSTICS - numbers reported in Methods\n")
cat("=========================================================\n")

cat("\n  stages present:", paste(sort(unique(tg$stage)), collapse = ", "), "\n")

cat("\n--- terms per comparison ---\n")
print(as.data.frame(count(tg, stage, direction)))

cat("\n--- translational terms per comparison ---\n")
print(as.data.frame(
  tg %>% filter(grepl(TRANS, Description, ignore.case = TRUE)) %>%
    count(stage, direction, name = "n_translational")))

# A - the negative control. Few or no ribosomal terms here means the
#     translational signature does not separate two uninfected animals.
cat("\n--- A: control vs control ---\n")
a_all  <- sum(tg$stage == "A_ctrl_vs_ctrl")
a_ribo <- sum(tg$stage == "A_ctrl_vs_ctrl" &
              grepl("ribosom", tg$Description, ignore.case = TRUE))
cat("total terms:", a_all, "| ribosomal terms:", a_ribo, "\n")

# B - upper bound on between-animal variation, at each timepoint
cat("\n--- B: animal vs animal within a timepoint ---\n")
b <- tg %>% filter(grepl("^B_", stage),
                   grepl(TRANS, Description, ignore.case = TRUE)) %>%
  count(stage, direction, name = "n_translational")
if (nrow(b)) print(as.data.frame(b)) else cat("  no translational terms\n")

# C - reproducibility across the two animals of each pair
cat("\n--- C: each animal vs pooled control ---\n")
cc <- tg %>% filter(grepl("^C_", stage),
                    grepl(TRANS, Description, ignore.case = TRUE)) %>%
  count(stage, cluster, direction)
if (nrow(cc)) print(as.data.frame(cc)) else cat("  no translational terms\n")

# D - does the signature survive animal-level aggregation, at each
#     replicated timepoint? Built from the stages actually present, so a
#     comparison that produced no GSEA results at all is reported as
#     absent rather than as zero - those mean different things.
cat("\n--- D: pseudobulk, 2 controls vs 2 infected ---\n")
d_stages <- sort(unique(tg$stage[grepl("^D_", tg$stage)]))
if (!length(d_stages)) {
  cat("  no D_ stages in GSEA_tidy.tsv - check that 08 wrote the\n")
  cat("  D_pseudobulk_* files and that 05c read them\n")
}
for (s in d_stages) {
  ds <- tg %>% filter(stage == s, grepl(TRANS, Description, ignore.case = TRUE))
  cat(sprintf("  %-22s %3d translational terms", s, nrow(ds)))
  if (nrow(ds))
    cat(sprintf(" | NES %.2f to %.2f | clusters: %s",
                min(ds$NES), max(ds$NES),
                paste(sort(unique(ds$cluster)), collapse = ", ")))
  cat("\n")
}

d <- tg %>% filter(grepl("^D_", stage), grepl(TRANS, Description, ignore.case = TRUE))
if (nrow(d)) write_tsv(d, file.path(outdir, "pseudobulk_translational_terms.tsv"))

# =========================================================
cat("\n--- HOW TO READ THIS ---\n")
cat("
  A with few or no translational terms, and C showing them in both
  animals of a pair independently, is the result that answers the
  pseudoreplication criticism: the signature is not what separates two
  uninfected animals, and it reproduces in individuals rather than
  depending on one.

  D surviving at either timepoint is the stronger claim, because the
  fold change there is animal-level. Note that at each timepoint one of
  the two infected animals carried little virus, so D is diluted - a
  signature that survives it is not an artefact of aggregation, and the
  dilution is worth stating rather than leaving a reader to find it.

  B is an upper bound inflated by that same infection difference. Do not
  report it as the magnitude of between-animal variation.
")

cat("\nOutputs in:", outdir, "\n")

