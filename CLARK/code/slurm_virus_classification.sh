#!/bin/bash
#SBATCH --job-name=CLARK_Virus_Classification
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=CLARK_Virus_Classification_%j.out
#SBATCH --error=CLARK_Virus_Classification_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# ============================================
# Configuration Parameters (modify as needed)
# ============================================
# CLARK_DB_DIR: CLARK database directory
export CLARK_DB_DIR=${CLARK_DB_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/CLARK_db"}
# CLARK_DIR: CLARK software installation directory, should contain classify_metagenome.sh, set_targets.sh, estimate_abundance.sh, etc.
export CLARK_DIR=${CLARK_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/CLARK"}
INPUT_FILE=${INPUT_FILE:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/long_reads/llnl_66d1047e.fastq.gz"}  # Input file path
OUTPUT_PREFIX=${OUTPUT_PREFIX:-"virus_results"}  # Output file prefix
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/run_CLARK/results_long"}  # Output directory
K_MER=${K_MER:-20}  # k-mer size (recommended 20-21 for virus detection, higher sensitivity)
MODE=${MODE:-0}  # Execution mode: 0=full, 1=default, 2=express
THREADS=$SLURM_CPUS_PER_TASK  # Use CPU count allocated by SLURM

# ============================================
# Start Classification
# ============================================
echo "=========================================="
echo "CLARK Virus Classification Job"
echo "Time: $(date)"
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Database directory: $CLARK_DB_DIR"
echo "k-mer size: $K_MER"
echo "Execution mode: $MODE"
echo "Threads: $THREADS"
echo "=========================================="

# Check input file
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file does not exist: $INPUT_FILE"
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

# Build classification command
CLASSIFY_CMD="./classify_metagenome.sh"
CLASSIFY_CMD="$CLASSIFY_CMD -O $INPUT_FILE"
CLASSIFY_CMD="$CLASSIFY_CMD -R $OUTPUT_DIR/$OUTPUT_PREFIX"
CLASSIFY_CMD="$CLASSIFY_CMD -m $MODE"
CLASSIFY_CMD="$CLASSIFY_CMD -k $K_MER"
CLASSIFY_CMD="$CLASSIFY_CMD -n $THREADS"

# If compressed file, add gzipped option
if [[ "$INPUT_FILE" == *.gz ]]; then
    CLASSIFY_CMD="$CLASSIFY_CMD --gzipped"
    echo "Detected gzip compressed file, adding --gzipped option"
fi

# If FASTA file, may need special handling
if [[ "$INPUT_FILE" == *.fa ]] || [[ "$INPUT_FILE" == *.fasta ]]; then
    echo "Detected FASTA format file"
fi

# Run classification
echo "Executing command: $CLASSIFY_CMD"
echo "Start time: $(date)"
echo "----------------------------------------"

# Run classification command, capture output and exit code
eval $CLASSIFY_CMD 2>&1 | tee /tmp/clark_classify_output.log
CLASSIFY_EXIT_CODE=${PIPESTATUS[0]}

echo "----------------------------------------"
echo "End time: $(date)"

# Check for segmentation fault or other serious errors
if grep -q "Segmentation fault\|core dumped\|Killed\|Out of memory" /tmp/clark_classify_output.log 2>/dev/null || [ $CLASSIFY_EXIT_CODE -ne 0 ]; then
    echo "=========================================="
    echo "Warning: Error or segmentation fault detected!"
    echo "Exit code: $CLASSIFY_EXIT_CODE"
    echo "Please check error log for details"
    echo "=========================================="
    
    # Even with errors, check if result file is valid
    RESULT_FILE="$OUTPUT_DIR/${OUTPUT_PREFIX}.csv"
    if [ -f "$RESULT_FILE" ]; then
        FILE_LINES=$(wc -l < "$RESULT_FILE")
        if [ "$FILE_LINES" -gt 1 ]; then
            echo ""
            echo "Note: Result file generated ($FILE_LINES lines), but may be incomplete"
            echo "Recommend checking result file integrity"
            echo "Result file: $RESULT_FILE"
            echo ""
        else
            echo "Error: Result file exists but content is invalid"
            exit 1
        fi
    else
        echo "Error: Result file not generated"
        exit 1
    fi
elif [ $CLASSIFY_EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "Virus classification completed successfully!"
    echo "Result file: $OUTPUT_DIR/${OUTPUT_PREFIX}.csv"
    echo "=========================================="
    
    # Display result file information
    if [ -f "$OUTPUT_DIR/${OUTPUT_PREFIX}.csv" ]; then
        echo "Result file statistics:"
        wc -l "$OUTPUT_DIR/${OUTPUT_PREFIX}.csv"
        
        # Check if file has valid content (at least header and data rows)
        FILE_LINES=$(wc -l < "$OUTPUT_DIR/${OUTPUT_PREFIX}.csv")
        if [ "$FILE_LINES" -lt 2 ]; then
            echo "Warning: Result file has too few lines, processing may be incomplete"
        fi
    fi
else
    echo "=========================================="
    echo "Error: Virus classification failed! Exit code: $CLASSIFY_EXIT_CODE"
    echo "=========================================="
    exit $CLASSIFY_EXIT_CODE
fi
