#!/bin/bash
#SBATCH --job-name=Viral_Classification_ShortRead
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Viral_Classification_ShortRead_%j.out
#SBATCH --error=Viral_Classification_ShortRead_%j.err

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "🦠  Metagenome Viral Classification Workflow (Short-Read Mode)"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0

# User's conda environment path
USER_CONDA_ENV="/home/sp96859/.conda/envs/gottcha2_env"

# Activate conda environment
source activate "$USER_CONDA_ENV"

# Run GOTTCHA2 classification analysis
echo "🔬 2. Running GOTTCHA2 taxonomic profiling..."
# Database path: directory containing .mmi index files
DB_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/gottcha2_db"
# Database prefix: determined by index filename gottcha_db.species.fna.mmi
DB_PREFIX="${DB_DIR}/gottcha_db.species.fna"

gottcha2 profile \
  -i llnl_66ce4dde_R1.fastq.gz llnl_66ce4dde_R2.fastq.gz \
  -d "$DB_PREFIX" \
  -t 32 \
  -o result_shortread
