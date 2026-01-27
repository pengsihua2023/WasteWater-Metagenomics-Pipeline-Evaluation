#!/bin/bash
# Quick submission script for long-read data classification
# Usage: bash run_long_read_classification.sh

# Set environment variables
export INPUT_FILE="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/long_reads/llnl_66d1047e.fastq.gz"
export OUTPUT_PREFIX="llnl_66d1047e_virus"
export OUTPUT_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/run_CLARK/results_long"
export K_MER=20
export MODE=0

# Display configuration information
echo "=========================================="
echo "Preparing to submit CLARK virus classification job"
echo "=========================================="
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Output prefix: $OUTPUT_PREFIX"
echo "k-mer size: $K_MER"
echo "Execution mode: $MODE"
echo "=========================================="
echo ""

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file does not exist: $INPUT_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Submit job to SLURM
echo "Submitting job to SLURM..."
sbatch slurm_virus_classification.sh

if [ $? -eq 0 ]; then
    echo ""
    echo "Job submitted successfully!"
    echo "Use 'squeue -u $USER' to check job status"
    echo "Results will be saved in: $OUTPUT_DIR"
else
    echo ""
    echo "Error: Job submission failed"
    exit 1
fi
