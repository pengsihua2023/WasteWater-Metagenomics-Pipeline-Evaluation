#!/bin/bash
# Virus detection and classification complete pipeline script

# Configuration parameters
CLARK_DB_DIR=${CLARK_DB_DIR:-"/path/to/CLARK_DB"}
INPUT_DIR=${INPUT_DIR:-"./input"}
OUTPUT_DIR=${OUTPUT_DIR:-"./virus_results"}
THREADS=${THREADS:-8}
THRESHOLD=${THRESHOLD:-0.7}
K_MER=${K_MER:-31}

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if CLARK is installed
if [ ! -f "./CLARK" ]; then
    print_error "CLARK executable file does not exist! Please compile CLARK first."
    exit 1
fi

# Check virus database
if [ ! -d "$CLARK_DB_DIR/viruses" ]; then
    print_warning "Virus database does not exist, starting build..."
    if [ -f "./set_targets.sh" ]; then
        ./set_targets.sh viruses --db $CLARK_DB_DIR
        if [ $? -ne 0 ]; then
            print_error "Database build failed!"
            exit 1
        fi
    else
        print_error "set_targets.sh script does not exist!"
        exit 1
    fi
else
    print_info "Virus database already exists: $CLARK_DB_DIR/viruses"
fi

# Create output directory
mkdir -p $OUTPUT_DIR

# Process input files
if [ -d "$INPUT_DIR" ]; then
    # Process all FASTQ files in directory
    FASTQ_FILES=($INPUT_DIR/*.fastq $INPUT_DIR/*.fastq.gz)
elif [ -f "$INPUT_DIR" ]; then
    # Process single file
    FASTQ_FILES=("$INPUT_DIR")
else
    print_error "Input directory or file does not exist: $INPUT_DIR"
    exit 1
fi

# Check if there are input files
if [ ${#FASTQ_FILES[@]} -eq 0 ] || [ ! -f "${FASTQ_FILES[0]}" ]; then
    print_error "No input files found!"
    exit 1
fi

print_info "Found ${#FASTQ_FILES[@]} input files"

# Virus classification for each sample
for fastq_file in "${FASTQ_FILES[@]}"; do
    if [ ! -f "$fastq_file" ]; then
        continue
    fi
    
    sample_name=$(basename "$fastq_file")
    sample_name=${sample_name%.fastq}
    sample_name=${sample_name%.fq}
    sample_name=${sample_name%.gz}
    
    print_info "Processing sample: $sample_name"
    
    # Build CLARK command
    CLARK_CMD="./CLARK -T $fastq_file"
    CLARK_CMD="$CLARK_CMD -O $OUTPUT_DIR"
    CLARK_CMD="$CLARK_CMD -R ${sample_name}_virus"
    CLARK_CMD="$CLARK_CMD -D $CLARK_DB_DIR"
    CLARK_CMD="$CLARK_CMD -m 0"
    CLARK_CMD="$CLARK_CMD -k $K_MER"
    CLARK_CMD="$CLARK_CMD -t $THREADS"
    CLARK_CMD="$CLARK_CMD --threshold $THRESHOLD"
    
    # If compressed file, add gzipped option
    if [[ "$fastq_file" == *.gz ]]; then
        CLARK_CMD="$CLARK_CMD --gzipped"
    fi
    
    # Run CLARK
    print_info "Running command: $CLARK_CMD"
    eval $CLARK_CMD
    
    if [ $? -eq 0 ]; then
        print_info "Sample $sample_name classification completed"
        
        # Calculate abundance (if Python script exists)
        result_file="$OUTPUT_DIR/${sample_name}_virus.csv"
        if [ -f "$result_file" ] && [ -f "calculate_virus_abundance.py" ]; then
            print_info "Calculating virus abundance for sample $sample_name..."
            python calculate_virus_abundance.py \
                "$result_file" \
                "$OUTPUT_DIR/${sample_name}_abundance.csv"
        fi
    else
        print_error "Sample $sample_name classification failed!"
    fi
    
    echo ""
done

print_info "All samples processed! Results saved in: $OUTPUT_DIR"


