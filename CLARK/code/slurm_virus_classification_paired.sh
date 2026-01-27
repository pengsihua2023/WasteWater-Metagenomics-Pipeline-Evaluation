#!/bin/bash
#SBATCH --job-name=CLARK_Virus_Classification_Paired
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=CLARK_Virus_Classification_Paired_%j.out
#SBATCH --error=CLARK_Virus_Classification_Paired_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# ============================================
# Configuration Parameters (modify as needed)
# ============================================
# CLARK_DB_DIR: CLARK database directory
export CLARK_DB_DIR=${CLARK_DB_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/CLARK_db"}
# CLARK_DIR: CLARK software installation directory, should contain classify_metagenome.sh, set_targets.sh, estimate_abundance.sh, etc.
export CLARK_DIR=${CLARK_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/CLARK"}
INPUT_R1=${INPUT_R1:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/short_reads/llnl_66ce4dde_R1.fastq.gz"}  # Paired-end R1 file
INPUT_R2=${INPUT_R2:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/short_reads/llnl_66ce4dde_R2.fastq.gz"}  # Paired-end R2 file
OUTPUT_PREFIX=${OUTPUT_PREFIX:-"virus_results"}  # Output file prefix
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/run_CLARK/results_short"}  # Output directory
K_MER=${K_MER:-20}  # k-mer size (recommended 20-21 for virus detection, higher sensitivity)
MODE=${MODE:-0}  # Execution mode: 0=full, 1=default, 2=express
THREADS=$SLURM_CPUS_PER_TASK  # Use CPU count allocated by SLURM

# ============================================
# Start Classification
# ============================================
echo "=========================================="
echo "CLARK Virus Classification Job (Paired-end)"
echo "Time: $(date)"
echo "Input file R1: $INPUT_R1"
echo "Input file R2: $INPUT_R2"
echo "Output directory: $OUTPUT_DIR"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Database directory: $CLARK_DB_DIR"
echo "k-mer size: $K_MER"
echo "Execution mode: $MODE"
echo "Threads: $THREADS"
echo "=========================================="

# Check input files
if [ ! -f "$INPUT_R1" ]; then
    echo "Error: R1 file does not exist: $INPUT_R1"
    exit 1
fi

if [ ! -f "$INPUT_R2" ]; then
    echo "Error: R2 file does not exist: $INPUT_R2"
    exit 1
fi

# Check database
if [ ! -d "$CLARK_DB_DIR/Viruses" ]; then
    echo "Error: Virus database does not exist: $CLARK_DB_DIR/Viruses"
    echo "Please run slurm_build_virus_db.sh first to build the database"
    exit 1
fi

# Check CLARK
if [ ! -f "$CLARK_DIR/classify_metagenome.sh" ]; then
    echo "Error: Cannot find classify_metagenome.sh script!"
    echo "Please set CLARK_DIR environment variable"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Enter CLARK directory
cd "$CLARK_DIR" || exit 1

# Build classification command (paired-end data)
CLASSIFY_CMD="./classify_metagenome.sh"
CLASSIFY_CMD="$CLASSIFY_CMD -O $INPUT_R1 $INPUT_R2"
CLASSIFY_CMD="$CLASSIFY_CMD -R $OUTPUT_DIR/$OUTPUT_PREFIX"
CLASSIFY_CMD="$CLASSIFY_CMD -m $MODE"
CLASSIFY_CMD="$CLASSIFY_CMD -k $K_MER"
CLASSIFY_CMD="$CLASSIFY_CMD -n $THREADS"

# If compressed file, add gzipped option
if [[ "$INPUT_R1" == *.gz ]] || [[ "$INPUT_R2" == *.gz ]]; then
    CLASSIFY_CMD="$CLASSIFY_CMD --gzipped"
    echo "Detected gzip compressed file, adding --gzipped option"
fi

# Run classification
echo "Executing command: $CLASSIFY_CMD"
echo "Start time: $(date)"
echo "----------------------------------------"

eval $CLASSIFY_CMD

CLASSIFY_EXIT_CODE=$?

echo "----------------------------------------"
echo "End time: $(date)"

if [ $CLASSIFY_EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "Virus classification completed successfully!"
    echo "Result file: $OUTPUT_DIR/${OUTPUT_PREFIX}.csv"
    echo "=========================================="
    
    # Display result file information
    if [ -f "$OUTPUT_DIR/${OUTPUT_PREFIX}.csv" ]; then
        echo "Result file statistics:"
        wc -l "$OUTPUT_DIR/${OUTPUT_PREFIX}.csv"
    fi
else
    echo "=========================================="
    echo "Error: Virus classification failed! Exit code: $CLASSIFY_EXIT_CODE"
    echo "=========================================="
    exit $CLASSIFY_EXIT_CODE
fi
