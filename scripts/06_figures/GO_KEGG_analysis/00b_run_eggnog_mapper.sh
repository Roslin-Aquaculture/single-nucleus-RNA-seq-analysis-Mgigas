#!/bin/bash
#$ -V -cwd
#$ -l h_rt=6:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 32

# generic (uncomment and edit if replicating elsewhere):
# #$ -P your_project_allocation
set -euo pipefail

# =========================================================
# Section 2.4: Run eggNOG-mapper on the full proteome
# Environment: eggnog (HPC / Eddie)
# Run once, genome-wide — not per cluster. Submit with:
#   qsub 00b_run_eggnog_mapper.sh
# Monitor with:
#   qstat -u <username>
#   tail -f "$OUT/logs/stdout.log"
#
# Output: full_proteome.emapper.annotations — a TSV with #query, GOs,
# KEGG_Pathway, KEGG_ko, and other annotation columns for every gene in
# the genome. This single file is reused for all downstream analyses in
# Sections 3-4.
# =========================================================

# =========================================================
# PROJECT PATHS
# =========================================================
PROJECT="/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/eggnog"
# generic:
# PROJECT="/path/to/your/hpc/project/eggnog"
DB="$PROJECT/eggnog_db"
RESULTS="$PROJECT/results"
FASTA="$PROJECT/input/cgigas_longest_clean.fa"

# =========================================================
# CONDA (HPC SAFE FOR EDDIE)
# =========================================================
set +u
source /exports/applications/apps/SL7/anaconda/5.3.1/etc/profile.d/conda.sh
# generic: source /path/to/your/anaconda/etc/profile.d/conda.sh
conda activate eggnog
set -u

# =========================================================
# TIMESTAMPED OUTPUT DIRECTORY
# =========================================================
TS=$(date +"%Y%m%d_%H%M%S")
OUT="$RESULTS/full_proteome_$TS"
LOGDIR="$OUT/logs"
mkdir -p "$OUT" "$LOGDIR"

# =========================================================
# ARCHIVE SCRIPT + ENVIRONMENT (REPRODUCIBILITY)
# =========================================================
cp "$0" "$OUT/run_script.sh"
env > "$OUT/env.txt"

# =========================================================
# REDIRECT ALL OUTPUT TO LOG FILES
# =========================================================
exec > >(tee -a "$LOGDIR/stdout.log") 2> >(tee -a "$LOGDIR/stderr.log" >&2)
echo "========================================"
echo "EggNOG RUN STARTED"
echo "Time   : $(date)"
echo "Input  : $FASTA"
echo "Output : $OUT"
echo "DB     : $DB"
echo "========================================"

# =========================================================
# RUN EGGNOG-MAPPER
# =========================================================
emapper.py \
  -i "$FASTA" \
  --itype proteins \
  -m diamond \
  --cpu 32 \
  --tax_scope Metazoa \
  --target_orthologs one2one \
  --data_dir "$DB" \
  -o "$OUT/full_proteome"

# =========================================================
# FINAL SUMMARY
# =========================================================
echo "========================================"
echo "EggNOG RUN COMPLETED"
echo "Time   : $(date)"
echo "Results:"
ls -lh "$OUT"
echo "Logs:"
ls -lh "$LOGDIR"
echo "========================================"
