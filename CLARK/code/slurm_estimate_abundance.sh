#!/bin/bash
#SBATCH --job-name=CLARK_Estimate_Abundance
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=CLARK_Estimate_Abundance_%j.out
#SBATCH --error=CLARK_Estimate_Abundance_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# ============================================
# Configuration Parameters (modify as needed)
# ============================================
# CLARK_DB_DIR: CLARK database directory
export CLARK_DB_DIR=${CLARK_DB_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/CLARK_db"}
# CLARK_DIR: CLARK software installation directory, should contain estimate_abundance.sh, classify_metagenome.sh, etc.
export CLARK_DIR=${CLARK_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/CLARK"}
INPUT_CSV=${INPUT_CSV:-"results/virus_results.csv"}  # CLARK classification result CSV file
OUTPUT_PREFIX=${OUTPUT_PREFIX:-"abundance"}  # Output file prefix
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/run_CLARK/results_long"}  # Output directory (default: results_long)

# Filtering parameters (optional)
CONFIDENCE_THRESHOLD=${CONFIDENCE_THRESHOLD:-""}  # Confidence threshold, e.g.: 0.80
GAMMA_THRESHOLD=${GAMMA_THRESHOLD:-""}  # Gamma threshold, e.g.: 0.03
ABUNDANCE_THRESHOLD=${ABUNDANCE_THRESHOLD:-""}  # Abundance threshold (%), e.g.: 2
HIGH_CONFIDENCE=${HIGH_CONFIDENCE:-"false"}  # Whether to use only high confidence assignments
MPA_FORMAT=${MPA_FORMAT:-"false"}  # Whether to output MetaPhlAn format

# ============================================
# Start Abundance Estimation
# ============================================
echo "=========================================="
echo "CLARK Abundance Estimation Job"
echo "Time: $(date)"
echo "Input file: $INPUT_CSV"
echo "Output directory: $OUTPUT_DIR"
echo "Database directory: $CLARK_DB_DIR"
echo "=========================================="

# Convert input file path to absolute path
if [ ! -f "$INPUT_CSV" ]; then
    echo "Error: Input file does not exist: $INPUT_CSV"
    exit 1
fi
# Get absolute path of input file
INPUT_CSV_ABS=$(readlink -f "$INPUT_CSV" 2>/dev/null || realpath "$INPUT_CSV" 2>/dev/null || echo "$(cd "$(dirname "$INPUT_CSV")" && pwd)/$(basename "$INPUT_CSV")")
INPUT_CSV="$INPUT_CSV_ABS"

# Check database
if [ ! -d "$CLARK_DB_DIR" ]; then
    echo "Error: Database directory does not exist: $CLARK_DB_DIR"
    exit 1
fi

# Check CLARK
if [ ! -f "$CLARK_DIR/estimate_abundance.sh" ]; then
    echo "Error: Cannot find estimate_abundance.sh script!"
    echo "Please set CLARK_DIR environment variable"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Enter CLARK directory
cd "$CLARK_DIR" || exit 1

# Build abundance estimation command
ABUNDANCE_CMD="./estimate_abundance.sh"
ABUNDANCE_CMD="$ABUNDANCE_CMD -F $INPUT_CSV"
ABUNDANCE_CMD="$ABUNDANCE_CMD -D $CLARK_DB_DIR"

# Add filtering options
if [ "$HIGH_CONFIDENCE" = "true" ]; then
    ABUNDANCE_CMD="$ABUNDANCE_CMD --highconfidence"
    echo "Using high confidence filtering"
fi

if [ -n "$CONFIDENCE_THRESHOLD" ]; then
    ABUNDANCE_CMD="$ABUNDANCE_CMD -c $CONFIDENCE_THRESHOLD"
    echo "Confidence threshold: $CONFIDENCE_THRESHOLD"
fi

if [ -n "$GAMMA_THRESHOLD" ]; then
    ABUNDANCE_CMD="$ABUNDANCE_CMD -g $GAMMA_THRESHOLD"
    echo "Gamma threshold: $GAMMA_THRESHOLD"
fi

if [ -n "$ABUNDANCE_THRESHOLD" ]; then
    ABUNDANCE_CMD="$ABUNDANCE_CMD -a $ABUNDANCE_THRESHOLD"
    echo "Abundance threshold: $ABUNDANCE_THRESHOLD%"
fi

if [ "$MPA_FORMAT" = "true" ]; then
    ABUNDANCE_CMD="$ABUNDANCE_CMD --mpa"
    echo "Output MetaPhlAn format"
fi

# Run abundance estimation
echo "Executing command: $ABUNDANCE_CMD"
echo "Start time: $(date)"
echo "----------------------------------------"

# Create temporary log file to capture output
TEMP_LOG="/tmp/clark_abundance_$$.log"
eval $ABUNDANCE_CMD 2>&1 | tee "$TEMP_LOG"

ABUNDANCE_EXIT_CODE=${PIPESTATUS[0]}

echo "----------------------------------------"
echo "End time: $(date)"

# Prepare output file path
INPUT_DIR=$(dirname "$INPUT_CSV")
INPUT_BASE=$(basename "$INPUT_CSV" .csv)
OUTPUT_DIR_ABS=$(readlink -f "$OUTPUT_DIR" 2>/dev/null || realpath "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")
mkdir -p "$OUTPUT_DIR_ABS"
OUTPUT_CSV="$OUTPUT_DIR_ABS/${INPUT_BASE}_abundance.csv"

if [ $ABUNDANCE_EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "Abundance estimation completed successfully!"
    echo "=========================================="
    
    # Try to extract CSV format abundance data from log file
    echo "Extracting abundance data from output..."
    echo "Output file path: $OUTPUT_CSV"
    
    # Find CSV data start line (match various possible header formats)
    CSV_START_LINE=""
    
    # Method 1: Find exact header "Name,TaxID,Lineage,Count,Proportion_All(%),Proportion_Classified(%)"
    CSV_START_LINE=$(grep -n "^Name,TaxID,Lineage,Count,Proportion_All" "$TEMP_LOG" | head -1 | cut -d: -f1)
    
    # Method 2: Find line starting with "Name,TaxID"
    if [ -z "$CSV_START_LINE" ]; then
        CSV_START_LINE=$(grep -n "^Name,TaxID" "$TEMP_LOG" | head -1 | cut -d: -f1)
    fi
    
    # Method 3: Find line containing "Name,TaxID,Lineage"
    if [ -z "$CSV_START_LINE" ]; then
        CSV_START_LINE=$(grep -n "Name,TaxID,Lineage" "$TEMP_LOG" | head -1 | cut -d: -f1)
    fi
    
    # Method 4: Find line containing "Name.*TaxID" (loose match)
    if [ -z "$CSV_START_LINE" ]; then
        CSV_START_LINE=$(grep -n -E "^Name[^,]*,[^,]*TaxID" "$TEMP_LOG" | head -1 | cut -d: -f1)
    fi
    
    if [ -n "$CSV_START_LINE" ]; then
        echo "Found CSV data start position at line $CSV_START_LINE in log"
        # Extract all content from CSV start line to end of file
        tail -n +$CSV_START_LINE "$TEMP_LOG" > "$OUTPUT_CSV"
        
        # Clean file: remove possible empty lines and non-data lines
        # Keep header and data lines (lines containing at least one comma)
        awk '/^Name,.*TaxID|^[^,]*,[0-9]+/' "$OUTPUT_CSV" > "$OUTPUT_CSV.tmp" 2>/dev/null || cat "$OUTPUT_CSV" > "$OUTPUT_CSV.tmp"
        mv "$OUTPUT_CSV.tmp" "$OUTPUT_CSV"
        
        # Verify if extracted file is valid
        if [ -f "$OUTPUT_CSV" ] && [ -s "$OUTPUT_CSV" ]; then
            # Check if file contains valid CSV data (at least header)
            FIRST_LINE=$(head -1 "$OUTPUT_CSV")
            if echo "$FIRST_LINE" | grep -q -E "Name.*TaxID|TaxID.*Lineage"; then
                echo "✓ Successfully extracted abundance data and saved to: $OUTPUT_CSV"
                echo "File size: $(ls -lh "$OUTPUT_CSV" | awk '{print $5}')"
                echo "Total lines: $(wc -l < "$OUTPUT_CSV")"
                echo ""
                echo "Header: $FIRST_LINE"
                echo ""
                echo "Top 10 most abundant taxa:"
                head -11 "$OUTPUT_CSV"
                
                # Clean temporary log file
                [ -f "$TEMP_LOG" ] && rm -f "$TEMP_LOG"
                RESULT_FILE="$OUTPUT_CSV"
            else
                echo "Warning: Extracted file format is incorrect, header is: $FIRST_LINE"
                CSV_START_LINE=""
            fi
        else
            echo "Warning: Failed to create output file"
            CSV_START_LINE=""
        fi
    fi
    
    # If still not successful, try to find existing files
    if [ -z "$RESULT_FILE" ] || [ ! -f "$RESULT_FILE" ]; then
        echo ""
        echo "Trying to find CLARK auto-created result files..."
        CURRENT_DIR=$(pwd)
        CLARK_DIR_ABS=$(readlink -f "$CLARK_DIR" 2>/dev/null || realpath "$CLARK_DIR" 2>/dev/null || echo "$CLARK_DIR")
        
        # Possible output file locations (ordered by priority)
        POSSIBLE_FILES=(
            "$INPUT_DIR/${INPUT_BASE}_abundance.csv"
            "$INPUT_DIR/${INPUT_BASE}.csv_abundance.csv"
            "$INPUT_DIR/abundance.csv"
            "$OUTPUT_DIR_ABS/${INPUT_BASE}_abundance.csv"
            "$OUTPUT_DIR_ABS/${INPUT_BASE}.csv_abundance.csv"
            "$OUTPUT_DIR_ABS/abundance.csv"
            "$CURRENT_DIR/${INPUT_BASE}_abundance.csv"
            "$CURRENT_DIR/${INPUT_BASE}.csv_abundance.csv"
            "$CURRENT_DIR/abundance.csv"
            "$CLARK_DIR_ABS/${INPUT_BASE}_abundance.csv"
            "$CLARK_DIR_ABS/${INPUT_BASE}.csv_abundance.csv"
            "$SLURM_SUBMIT_DIR/${INPUT_BASE}_abundance.csv"
        )
        
        for file in "${POSSIBLE_FILES[@]}"; do
            if [ -f "$file" ]; then
                # If file found, copy to target output directory
                if [ "$file" != "$OUTPUT_CSV" ]; then
                    cp "$file" "$OUTPUT_CSV"
                    echo "✓ Found result file and copied to: $OUTPUT_CSV"
                else
                    echo "✓ Found result file: $OUTPUT_CSV"
                fi
                echo "File size: $(ls -lh "$OUTPUT_CSV" | awk '{print $5}')"
                echo ""
                echo "Top 10 most abundant taxa:"
                head -11 "$OUTPUT_CSV"
                RESULT_FILE="$OUTPUT_CSV"
                break
            fi
        done
        
        # If still not found, finally try to extract all CSV-like lines from log
        if [ -z "$RESULT_FILE" ] || [ ! -f "$RESULT_FILE" ]; then
            echo ""
            echo "Trying to extract all CSV format data from log..."
            
            # Find all lines containing multiple commas (CSV format)
            # Match pattern: contains at least 2 commas, and second field is number (TaxID)
            grep -E "^[^,]+,[0-9]+," "$TEMP_LOG" > "$OUTPUT_CSV" 2>/dev/null || true
            
            if [ -f "$OUTPUT_CSV" ] && [ -s "$OUTPUT_CSV" ]; then
                # Check if first line is header
                FIRST_LINE=$(head -1 "$OUTPUT_CSV")
                if ! echo "$FIRST_LINE" | grep -q -E "Name.*TaxID|TaxID.*Lineage"; then
                    # If no header, add standard header
                    echo "Name,TaxID,Lineage,Count,Proportion_All(%),Proportion_Classified(%)" > "$OUTPUT_CSV.tmp"
                    cat "$OUTPUT_CSV" >> "$OUTPUT_CSV.tmp"
                    mv "$OUTPUT_CSV.tmp" "$OUTPUT_CSV"
                    echo "Added standard CSV header"
                fi
                
                # Verify file format
                if head -1 "$OUTPUT_CSV" | grep -q -E "Name|TaxID"; then
                    echo "✓ Extracted data from log and saved to: $OUTPUT_CSV"
                    echo "File size: $(ls -lh "$OUTPUT_CSV" | awk '{print $5}')"
                    echo "Total lines: $(wc -l < "$OUTPUT_CSV")"
                    echo ""
                    echo "Top 10 most abundant taxa:"
                    head -11 "$OUTPUT_CSV"
                    RESULT_FILE="$OUTPUT_CSV"
                else
                    echo "Warning: Extracted data format is incorrect"
                    rm -f "$OUTPUT_CSV"
                fi
            else
                echo "Error: Unable to extract valid CSV data from log"
                echo ""
                echo "Log file content preview (searching for lines containing 'Name' or 'TaxID'):"
                grep -i -E "Name|TaxID" "$TEMP_LOG" | head -10
                echo ""
                echo "Log file content preview (last 50 lines):"
                tail -50 "$TEMP_LOG"
            fi
        fi
    fi
    
    # Final confirmation
    if [ -n "$RESULT_FILE" ] && [ -f "$RESULT_FILE" ]; then
        echo ""
        echo "=========================================="
        echo "Abundance estimation result file: $RESULT_FILE"
        echo "=========================================="
    else
        echo ""
        echo "Warning: Failed to extract or create result file"
        echo "Please check log file $TEMP_LOG for more information"
    fi
    
    # Clean temporary log file
    [ -f "$TEMP_LOG" ] && rm -f "$TEMP_LOG"
else
    echo "=========================================="
    echo "Error: Abundance estimation failed! Exit code: $ABUNDANCE_EXIT_CODE"
    echo "=========================================="
    # Display last few lines of error log
    if [ -f "$TEMP_LOG" ]; then
        echo ""
        echo "Error output summary:"
        tail -20 "$TEMP_LOG"
        rm -f "$TEMP_LOG"
    fi
    exit $ABUNDANCE_EXIT_CODE
fi
