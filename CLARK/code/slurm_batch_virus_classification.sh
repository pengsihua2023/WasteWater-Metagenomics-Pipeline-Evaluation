#!/bin/bash
#SBATCH --job-name=CLARK_Batch_Virus_Classification
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=CLARK_Batch_Virus_Classification_%j.out
#SBATCH --error=CLARK_Batch_Virus_Classification_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# ============================================
# Configuration Parameters (modify as needed)
# ============================================
# CLARK_DB_DIR: CLARK database directory
export CLARK_DB_DIR=${CLARK_DB_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/CLARK_db"}
# CLARK_DIR: CLARK software installation directory, should contain classify_metagenome.sh, estimate_abundance.sh, etc.
export CLARK_DIR=${CLARK_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/CLARK"}
INPUT_DIR=${INPUT_DIR:-"$(pwd)/input"}  # Input file directory
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/run_CLARK/results"}  # Output directory

# CLARK parameters
K_MER=${K_MER:-20}  # k-mer size (recommended 20-21 for virus detection)
MODE=${MODE:-0}  # Execution mode: 0=full, 1=default, 2=express
THREADS=$SLURM_CPUS_PER_TASK  # Use CPU count allocated by SLURM

# Whether to perform abundance estimation
ESTIMATE_ABUNDANCE=${ESTIMATE_ABUNDANCE:-"true"}

# ============================================
# Start Batch Processing
# ============================================
echo "=========================================="
echo "CLARK Batch Virus Classification Job"
echo "Start time: $(date)"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Database directory: $CLARK_DB_DIR"
echo "k-mer size: $K_MER"
echo "Execution mode: $MODE"
echo "Threads: $THREADS"
echo "=========================================="

# Check input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
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

# Find all input files
if [ -n "$(find "$INPUT_DIR" -maxdepth 1 -name "*.fastq" -o -name "*.fq" -o -name "*.fa" -o -name "*.fasta" 2>/dev/null)" ]; then
    # Single-end files
    FILES=($(find "$INPUT_DIR" -maxdepth 1 \( -name "*.fastq" -o -name "*.fq" -o -name "*.fa" -o -name "*.fasta" \) -type f | sort))
    PAIRED_MODE=false
elif [ -n "$(find "$INPUT_DIR" -maxdepth 1 -name "*_R1*" 2>/dev/null)" ]; then
    # Paired-end files
    R1_FILES=($(find "$INPUT_DIR" -maxdepth 1 -name "*_R1*" -type f | sort))
    PAIRED_MODE=true
else
    echo "Error: No input files found in $INPUT_DIR!"
    exit 1
fi

# Process files
if [ "$PAIRED_MODE" = "true" ]; then
    echo "Detected paired-end sequencing files"
    TOTAL_SAMPLES=${#R1_FILES[@]}
    echo "Found $TOTAL_SAMPLES samples"
    
    for R1_FILE in "${R1_FILES[@]}"; do
        # Find corresponding R2 file
        R2_FILE=$(echo "$R1_FILE" | sed 's/_R1/_R2/g' | sed 's/_1\./_2./g')
        
        if [ ! -f "$R2_FILE" ]; then
            echo "Warning: R2 file not found: $R2_FILE, skipping sample"
            continue
        fi
        
        SAMPLE_NAME=$(basename "$R1_FILE")
        SAMPLE_NAME=${SAMPLE_NAME%_R1*}
        SAMPLE_NAME=${SAMPLE_NAME%_1.*}
        SAMPLE_NAME=${SAMPLE_NAME%.*}
        
        echo ""
        echo "=========================================="
        echo "Processing sample: $SAMPLE_NAME"
        echo "R1 file: $R1_FILE"
        echo "R2 file: $R2_FILE"
        echo "=========================================="
        
        # Build classification command
        CLASSIFY_CMD="./classify_metagenome.sh"
        CLASSIFY_CMD="$CLASSIFY_CMD -O $R1_FILE $R2_FILE"
        CLASSIFY_CMD="$CLASSIFY_CMD -R $OUTPUT_DIR/${SAMPLE_NAME}_virus"
        CLASSIFY_CMD="$CLASSIFY_CMD -m $MODE"
        CLASSIFY_CMD="$CLASSIFY_CMD -k $K_MER"
        CLASSIFY_CMD="$CLASSIFY_CMD -n $THREADS"
        
        if [[ "$R1_FILE" == *.gz ]] || [[ "$R2_FILE" == *.gz ]]; then
            CLASSIFY_CMD="$CLASSIFY_CMD --gzipped"
        fi
        
        # Run classification
        eval $CLASSIFY_CMD
        
        if [ $? -eq 0 ]; then
            echo "Sample $SAMPLE_NAME classification completed"
            
            # Abundance estimation
            if [ "$ESTIMATE_ABUNDANCE" = "true" ] && [ -f "$CLARK_DIR/estimate_abundance.sh" ]; then
                RESULT_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_virus.csv"
                if [ -f "$RESULT_FILE" ]; then
                    echo "Calculating abundance for sample $SAMPLE_NAME..."
                    ./estimate_abundance.sh -F "$RESULT_FILE" -D "$CLARK_DB_DIR" -g 0.03
                fi
            fi
        else
            echo "Error: Sample $SAMPLE_NAME classification failed"
        fi
    done
else
    echo "Detected single-end sequencing files"
    TOTAL_SAMPLES=${#FILES[@]}
    echo "Found $TOTAL_SAMPLES samples"
    
    for INPUT_FILE in "${FILES[@]}"; do
        SAMPLE_NAME=$(basename "$INPUT_FILE")
        SAMPLE_NAME=${SAMPLE_NAME%.*}
        SAMPLE_NAME=${SAMPLE_NAME%.gz}
        
        echo ""
        echo "=========================================="
        echo "Processing sample: $SAMPLE_NAME"
        echo "Input file: $INPUT_FILE"
        echo "=========================================="
        
        # Build classification command
        CLASSIFY_CMD="./classify_metagenome.sh"
        CLASSIFY_CMD="$CLASSIFY_CMD -O $INPUT_FILE"
        CLASSIFY_CMD="$CLASSIFY_CMD -R $OUTPUT_DIR/${SAMPLE_NAME}_virus"
        CLASSIFY_CMD="$CLASSIFY_CMD -m $MODE"
        CLASSIFY_CMD="$CLASSIFY_CMD -k $K_MER"
        CLASSIFY_CMD="$CLASSIFY_CMD -n $THREADS"
        
        if [[ "$INPUT_FILE" == *.gz ]]; then
            CLASSIFY_CMD="$CLASSIFY_CMD --gzipped"
        fi
        
        # Run classification
        eval $CLASSIFY_CMD
        
        if [ $? -eq 0 ]; then
            echo "Sample $SAMPLE_NAME classification completed"
            
            # Abundance estimation
            if [ "$ESTIMATE_ABUNDANCE" = "true" ] && [ -f "$CLARK_DIR/estimate_abundance.sh" ]; then
                RESULT_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_virus.csv"
                if [ -f "$RESULT_FILE" ]; then
                    echo "Calculating abundance for sample $SAMPLE_NAME..."
                    ./estimate_abundance.sh -F "$RESULT_FILE" -D "$CLARK_DB_DIR" -g 0.03
                fi
            fi
        else
            echo "Error: Sample $SAMPLE_NAME classification failed"
        fi
    done
fi

echo ""
echo "=========================================="
echo "Batch processing completed!"
echo "Completion time: $(date)"
echo "Results saved in: $OUTPUT_DIR"
echo "=========================================="
