#!/bin/bash
#SBATCH --job-name=CLARK_Virus_Pipeline
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=CLARK_Virus_Pipeline_%j.out
#SBATCH --error=CLARK_Virus_Pipeline_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# ============================================
# Configuration Parameters (modify as needed)
# ============================================
# CLARK_DB_DIR: CLARK database directory
export CLARK_DB_DIR=${CLARK_DB_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/CLARK_db"}
# CLARK_DIR: CLARK software installation directory, should contain classify_metagenome.sh, estimate_abundance.sh, etc.
export CLARK_DIR=${CLARK_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/CLARK"}

# Input files (can be single-end or paired-end)
INPUT_FILE=${INPUT_FILE:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/long_reads/llnl_66d1047e.fastq.gz"}  # Single-end file, or R1 file
INPUT_R2=${INPUT_R2:-""}  # Paired-end R2 file (if empty use single-end mode, set R2 path to use paired-end mode)
OUTPUT_PREFIX=${OUTPUT_PREFIX:-"virus_results"}  # Output file prefix
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/run_CLARK/results_long"}  # Output directory (default long-read, can be overridden by environment variable)

# CLARK parameters
K_MER=${K_MER:-20}  # k-mer size (recommended 20-21 for virus detection, higher sensitivity)
MODE=${MODE:-0}  # Execution mode: 0=full, 1=default, 2=express
THREADS=$SLURM_CPUS_PER_TASK  # Use CPU count allocated by SLURM

# Abundance estimation parameters
ESTIMATE_ABUNDANCE=${ESTIMATE_ABUNDANCE:-"true"}  # Whether to perform abundance estimation
GAMMA_THRESHOLD=${GAMMA_THRESHOLD:-"0.03"}  # Gamma threshold
MPA_FORMAT=${MPA_FORMAT:-"false"}  # Whether to output MetaPhlAn format

# ============================================
# Start Processing
# ============================================
echo "=========================================="
echo "CLARK Virus Detection Complete Pipeline"
echo "Start time: $(date)"
echo "Input file: $INPUT_FILE"
if [ -n "$INPUT_R2" ]; then
    echo "Input file R2: $INPUT_R2"
    echo "Mode: Paired-end sequencing"
else
    echo "Mode: Single-end sequencing"
fi
echo "Output directory: $OUTPUT_DIR"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Database directory: $CLARK_DB_DIR"
echo "k-mer size: $K_MER"
echo "Execution mode: $MODE"
echo "Threads: $THREADS"
echo "=========================================="

# Check input files
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file does not exist: $INPUT_FILE"
    exit 1
fi

if [ -n "$INPUT_R2" ] && [ ! -f "$INPUT_R2" ]; then
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

# ============================================
# Step 1: Virus Classification
# ============================================
echo ""
echo "=========================================="
echo "Step 1: Start Virus Classification"
echo "=========================================="

cd "$CLARK_DIR" || exit 1

# Build classification command
CLASSIFY_CMD="./classify_metagenome.sh"
if [ -n "$INPUT_R2" ]; then
    # Paired-end mode
    CLASSIFY_CMD="$CLASSIFY_CMD -O $INPUT_FILE $INPUT_R2"
else
    # Single-end mode
    CLASSIFY_CMD="$CLASSIFY_CMD -O $INPUT_FILE"
fi

CLASSIFY_CMD="$CLASSIFY_CMD -R $OUTPUT_DIR/$OUTPUT_PREFIX"
CLASSIFY_CMD="$CLASSIFY_CMD -m $MODE"
CLASSIFY_CMD="$CLASSIFY_CMD -k $K_MER"
CLASSIFY_CMD="$CLASSIFY_CMD -n $THREADS"

# If compressed file, add gzipped option
if [[ "$INPUT_FILE" == *.gz ]] || [[ "$INPUT_R2" == *.gz ]]; then
    CLASSIFY_CMD="$CLASSIFY_CMD --gzipped"
    echo "Detected gzip compressed file, adding --gzipped option"
fi

# Run classification
echo "Executing command: $CLASSIFY_CMD"
echo "Start time: $(date)"
echo "----------------------------------------"

eval $CLASSIFY_CMD

CLASSIFY_EXIT_CODE=$?

if [ $CLASSIFY_EXIT_CODE -ne 0 ]; then
    echo "Error: Virus classification failed! Exit code: $CLASSIFY_EXIT_CODE"
    exit $CLASSIFY_EXIT_CODE
fi

echo "Classification completion time: $(date)"
CLASSIFY_RESULT="$OUTPUT_DIR/${OUTPUT_PREFIX}.csv"

if [ ! -f "$CLASSIFY_RESULT" ]; then
    echo "Error: Classification result file does not exist: $CLASSIFY_RESULT"
    exit 1
fi

echo "Classification result file: $CLASSIFY_RESULT"
echo "Result file statistics:"
wc -l "$CLASSIFY_RESULT"

# ============================================
# Step 2: Abundance Estimation
# ============================================
if [ "$ESTIMATE_ABUNDANCE" = "true" ]; then
    echo ""
    echo "=========================================="
    echo "Step 2: Start Abundance Estimation"
    echo "=========================================="
    
    if [ ! -f "$CLARK_DIR/estimate_abundance.sh" ]; then
        echo "Warning: Cannot find estimate_abundance.sh script, skipping abundance estimation"
    else
        # Build abundance estimation command
        ABUNDANCE_CMD="./estimate_abundance.sh"
        ABUNDANCE_CMD="$ABUNDANCE_CMD -F $CLASSIFY_RESULT"
        ABUNDANCE_CMD="$ABUNDANCE_CMD -D $CLARK_DB_DIR"
        
        if [ -n "$GAMMA_THRESHOLD" ]; then
            ABUNDANCE_CMD="$ABUNDANCE_CMD -g $GAMMA_THRESHOLD"
        fi
        
        if [ "$MPA_FORMAT" = "true" ]; then
            ABUNDANCE_CMD="$ABUNDANCE_CMD --mpa"
        fi
        
        echo "Executing command: $ABUNDANCE_CMD"
        echo "Start time: $(date)"
        echo "----------------------------------------"
        
        eval $ABUNDANCE_CMD
        
        ABUNDANCE_EXIT_CODE=$?
        
        if [ $ABUNDANCE_EXIT_CODE -eq 0 ]; then
            echo "Abundance estimation completion time: $(date)"
            echo "Abundance estimation completed successfully!"
        else
            echo "Warning: Abundance estimation failed, but classification completed. Exit code: $ABUNDANCE_EXIT_CODE"
        fi
    fi
else
    echo ""
    echo "Skip abundance estimation step"
fi

# ============================================
# Completion
# ============================================
echo ""
echo "=========================================="
echo "CLARK Virus Detection Pipeline Completed!"
echo "Completion time: $(date)"
echo "Result file: $CLASSIFY_RESULT"
echo "=========================================="
