#!/bin/bash
#SBATCH --job-name=Viral_Classification_LongRead
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Viral_Classification_LongRead_%j.out
#SBATCH --error=Viral_Classification_LongRead_%j.err

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "🦠  Metagenome Viral Classification Workflow (Long-Read Mode)"
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

# Run GOTTCHA2 classification analysis (Long-read mode)
echo "🔬 2. Running GOTTCHA2 taxonomic profiling (Long-Read Mode)..."
# Database path: directory containing .mmi index files
DB_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/gottcha2_db"
# Database prefix: determined by index filename gottcha_db.species.fna.mmi
DB_PREFIX="${DB_DIR}/gottcha_db.species.fna"

# Long-read data input (single file)
# Supported formats: .fastq.gz, .fq.gz, .fastq, .fq
# Select appropriate parameters based on your data type:
# - Oxford Nanopore: use --nanopore parameter
# - PacBio: use --pacbio parameter

gottcha2 profile \
  -i llnl_66d1047e.fastq.gz \
  -d "$DB_PREFIX" \
  -t 32 \
  --nanopore \
  -o result_longread

echo ""
echo "=========================================="
echo "✅ GOTTCHA2 analysis completed!"
echo "End time: $(date)"
echo "=========================================="
echo ""
echo "📊 Output files are in: result_longread/"
echo "   - result_longread.tsv: Classification results table"
echo "   - result_longread.summary.tsv: Classification summary"
echo "   - result_longread.full.tsv: Full classification information"
