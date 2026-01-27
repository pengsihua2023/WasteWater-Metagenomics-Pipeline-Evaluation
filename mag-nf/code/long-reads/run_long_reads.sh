#!/bin/bash
#SBATCH --job-name=long_reads_flye
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=long_reads_flye_%j.out
#SBATCH --error=long_reads_flye_%j.err

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "🦠 Long Reads Analysis - Direct Pipeline"
echo "=========================================="
echo "Pipeline: Flye Assembly → Kraken2 Classification"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Set paths
LONG_READS="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/long_reads/llnl_66d1047e.fastq.gz"
KRAKEN2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref"
OUTPUT_DIR="results_long"
SAMPLE="llnl_66d1047e"

mkdir -p "$OUTPUT_DIR"

# Set Apptainer cache to /scratch (avoid home directory quota)
export APPTAINER_CACHEDIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/cache"
export APPTAINER_TMPDIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/tmp"
mkdir -p "$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_TMPDIR"

echo "✅ Apptainer cache directory: $APPTAINER_CACHEDIR"

# Container paths (use verified versions)
FLYE_CONTAINER="docker://quay.io/biocontainers/flye:2.9.1--py39h6935b12_0"
KRAKEN2_CONTAINER="docker://quay.io/biocontainers/kraken2:2.1.2--pl5262h7d875b9_0"

echo "📁 Input: $LONG_READS"
echo "🗄️  Kraken2 DB: $KRAKEN2_DB"
echo "📊 Output: $OUTPUT_DIR"
echo ""

# ==============================
# Step 1: Flye metagenomic assembly
# ==============================
echo "=========================================="
echo "📊 Step 1: Flye Metagenomic Assembly"
echo "=========================================="

apptainer run \
    --no-mount /lscratch \
    --bind /scratch:/scratch \
    $FLYE_CONTAINER \
    flye \
        --nano-raw "$LONG_READS" \
        --out-dir "${OUTPUT_DIR}/flye_assembly" \
        --threads 32 \
        --meta \
        --genome-size 5m

FLYE_EXIT=$?

if [ $FLYE_EXIT -eq 0 ]; then
    echo "✅ Flye assembly completed"
    
    # Check assembly stats
    if [ -f "${OUTPUT_DIR}/flye_assembly/assembly.fasta" ]; then
        CONTIGS=$(grep -c "^>" "${OUTPUT_DIR}/flye_assembly/assembly.fasta")
        SIZE=$(stat -f%z "${OUTPUT_DIR}/flye_assembly/assembly.fasta" 2>/dev/null || stat -c%s "${OUTPUT_DIR}/flye_assembly/assembly.fasta" 2>/dev/null)
        echo "   - Contigs: $CONTIGS"
        echo "   - Assembly size: $SIZE bytes"
    fi
else
    echo "❌ Flye assembly failed with exit code: $FLYE_EXIT"
    exit 1
fi
echo ""

# ==============================
# Step 2: Kraken2 classification on contigs
# ==============================
echo "=========================================="
echo "📊 Step 2: Kraken2 Viral Classification (Contigs)"
echo "=========================================="

apptainer run \
    --no-mount /lscratch \
    --bind /scratch:/scratch \
    $KRAKEN2_CONTAINER \
    kraken2 \
        --db "$KRAKEN2_DB" \
        --threads 16 \
        --output "${OUTPUT_DIR}/${SAMPLE}_contigs_classification.txt" \
        --report "${OUTPUT_DIR}/${SAMPLE}_contigs_kraken2_report.txt" \
        "${OUTPUT_DIR}/flye_assembly/assembly.fasta"

echo "✅ Kraken2 classification on contigs completed"
echo ""

# ==============================
# Step 3: Kraken2 classification on reads
# ==============================
echo "=========================================="
echo "📊 Step 3: Kraken2 Viral Classification (Reads)"
echo "=========================================="

apptainer run \
    --no-mount /lscratch \
    --bind /scratch:/scratch \
    $KRAKEN2_CONTAINER \
    kraken2 \
        --db "$KRAKEN2_DB" \
        --threads 16 \
        --output "${OUTPUT_DIR}/${SAMPLE}_reads_classification.txt" \
        --report "${OUTPUT_DIR}/${SAMPLE}_reads_kraken2_report.txt" \
        "$LONG_READS"

echo "✅ Kraken2 classification on reads completed"
echo ""

# ==============================
# Results Summary
# ==============================
echo "=========================================="
echo "📊 Analysis Summary"
echo "=========================================="
echo ""
echo "✅ All steps completed successfully!"
echo ""
echo "📁 Results directory: $OUTPUT_DIR"
echo ""
echo "Generated files:"
echo "  1. Assembly:"
echo "     - ${OUTPUT_DIR}/flye_assembly/assembly.fasta"
echo "     - ${OUTPUT_DIR}/flye_assembly/assembly_info.txt"
echo ""
echo "  2. Viral Classification (Contigs):"
echo "     - ${OUTPUT_DIR}/${SAMPLE}_contigs_kraken2_report.txt"
echo ""
echo "  3. Viral Classification (Reads):"
echo "     - ${OUTPUT_DIR}/${SAMPLE}_reads_kraken2_report.txt"
echo ""

# Show top viral classifications
echo "📈 Top viral classifications (from contigs):"
awk '$4=="D" || $4=="P" || $4=="F" {printf "  %6.2f%% - %s\n", $1, $6}' \
    "${OUTPUT_DIR}/${SAMPLE}_contigs_kraken2_report.txt" | head -10

echo ""
echo "End time: $(date)"
echo "=========================================="
