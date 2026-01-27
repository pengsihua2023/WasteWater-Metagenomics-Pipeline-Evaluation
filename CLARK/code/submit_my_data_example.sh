#!/bin/bash
# Usage example: Submit your data for analysis
# 
# Usage:
# 1. Modify the parameters below (if needed)
# 2. Run: bash submit_my_data_example.sh
# Or directly use sbatch to submit the corresponding script

# ============================================
# Data Path Configuration
# ============================================

# Long-read data path
LONG_READ_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/long_reads"
LONG_READ_FILE="llnl_66d1047e.fastq.gz"
LONG_READ_PATH="$LONG_READ_DIR/$LONG_READ_FILE"

# Short-read data path (paired-end)
SHORT_READ_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/short_reads"
SHORT_READ_R1="llnl_66ce4dde_R1.fastq.gz"
SHORT_READ_R2="llnl_66ce4dde_R2.fastq.gz"
SHORT_READ_R1_PATH="$SHORT_READ_DIR/$SHORT_READ_R1"
SHORT_READ_R2_PATH="$SHORT_READ_DIR/$SHORT_READ_R2"

# Output directories
OUTPUT_DIR_SHORT="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/run_CLARK/results_short"
OUTPUT_DIR_LONG="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/run_CLARK/results_long"

# ============================================
# Option 1: Process short-read data (paired-end sequencing)
# ============================================
echo "=========================================="
echo "Option 1: Submit short-read data (paired-end) analysis job"
echo "=========================================="
echo ""
echo "Set environment variables and submit job:"
echo ""
echo "export INPUT_R1=$SHORT_READ_R1_PATH"
echo "export INPUT_R2=$SHORT_READ_R2_PATH"
echo "export OUTPUT_PREFIX=llnl_66ce4dde_virus"
echo "export OUTPUT_DIR=$OUTPUT_DIR_SHORT"
echo "export K_MER=20"
echo "export MODE=0"
echo "sbatch slurm_virus_classification_paired.sh"
echo ""
echo "Or use complete pipeline (including abundance estimation):"
echo ""
echo "export INPUT_FILE=$SHORT_READ_R1_PATH"
echo "export INPUT_R2=$SHORT_READ_R2_PATH"
echo "export OUTPUT_PREFIX=llnl_66ce4dde_virus"
echo "export OUTPUT_DIR=$OUTPUT_DIR_SHORT"
echo "export K_MER=20"
echo "export MODE=0"
echo "export ESTIMATE_ABUNDANCE=true"
echo "export GAMMA_THRESHOLD=0.03"
echo "sbatch slurm_virus_pipeline.sh"
echo ""

# ============================================
# Option 2: Process long-read data (single-end sequencing)
# ============================================
echo "=========================================="
echo "Option 2: Submit long-read data (single-end) analysis job"
echo "=========================================="
echo ""
echo "Set environment variables and submit job:"
echo ""
echo "export INPUT_FILE=$LONG_READ_PATH"
echo "export OUTPUT_PREFIX=llnl_66d1047e_virus"
echo "export OUTPUT_DIR=$OUTPUT_DIR_LONG"
echo "export K_MER=20"
echo "export MODE=0"
echo "sbatch slurm_virus_classification.sh"
echo ""
echo "Or use complete pipeline (including abundance estimation):"
echo ""
echo "export INPUT_FILE=$LONG_READ_PATH"
echo "export OUTPUT_PREFIX=llnl_66d1047e_virus"
echo "export OUTPUT_DIR=$OUTPUT_DIR_LONG"
echo "export K_MER=20"
echo "export MODE=0"
echo "export ESTIMATE_ABUNDANCE=true"
echo "export GAMMA_THRESHOLD=0.03"
echo "sbatch slurm_virus_pipeline.sh"
echo ""

# ============================================
# Quick submission commands (uncomment to use)
# ============================================

# Submit short-read data job (paired-end, complete pipeline)
# export INPUT_FILE=$SHORT_READ_R1_PATH
# export INPUT_R2=$SHORT_READ_R2_PATH
# export OUTPUT_PREFIX=llnl_66ce4dde_virus
# export OUTPUT_DIR=$OUTPUT_DIR_SHORT
# export K_MER=20
# export MODE=0
# export ESTIMATE_ABUNDANCE=true
# export GAMMA_THRESHOLD=0.03
# sbatch slurm_virus_pipeline.sh

# Submit long-read data job (single-end, complete pipeline)
# export INPUT_FILE=$LONG_READ_PATH
# export OUTPUT_PREFIX=llnl_66d1047e_virus
# export OUTPUT_DIR=$OUTPUT_DIR_LONG
# export K_MER=20
# export MODE=0
# export ESTIMATE_ABUNDANCE=true
# export GAMMA_THRESHOLD=0.03
# sbatch slurm_virus_pipeline.sh
