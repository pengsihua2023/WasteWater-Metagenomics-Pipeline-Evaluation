#!/bin/bash
#SBATCH --job-name=nfcore_mag_short
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=9-00:00:00
#SBATCH --output=nfcore_mag_short_%j.out
#SBATCH --error=nfcore_mag_short_%j.err

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "🦠 nf-core/mag - Short Reads Analysis"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load environment
module load Miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nextflow_env

# Set environment variables
export NXF_OPTS="-Xms1g -Xmx4g"
export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases:/databases"
export APPTAINER_NO_MOUNT="/lscratch"
export SINGULARITY_NO_MOUNT="/lscratch"
export NXF_SINGULARITY_CACHEDIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/singularity"

KRAKEN2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref"

echo "=========================================="
echo "📊 Short Reads Analysis with nf-core/mag"
echo "=========================================="

# Create short reads samplesheet
cat > samplesheet_short.csv << 'EOF'
sample,group,short_reads_1,short_reads_2,long_reads
llnl_66ce4dde,viral_group1,/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/short_reads/llnl_66ce4dde_R1.fastq.gz,/scratch/sp96859/Meta-genome-data-analysis/Apptainer/data/short_reads/llnl_66ce4dde_R2.fastq.gz,
EOF

echo "✅ Short reads samplesheet created"
echo ""

# Run nf-core/mag v2.5.4
nextflow run nf-core/mag \
    -r 2.5.4 \
    -profile apptainer \
    -c nfcore_mag_v3_no_validation.config \
    --input samplesheet_short.csv \
    --outdir results_short \
    --max_cpus 32 \
    --max_memory '256.GB' \
    --max_time '72.h' \
    --coassemble_group \
    --megahit_options '--min-contig-len 1000' \
    --spades_options '--meta --only-assembler' \
    --fastp_save_trimmed_fail \
    --fastp_qualified_quality 20 \
    --save_clipped_reads \
    --kraken2_db "$KRAKEN2_DB" \
    --skip_binning \
    --skip_binqc \
    --skip_gtdbtk \
    --skip_prokka \
    --skip_quast \
    -work-dir work_short \
    -resume \
    -with-report results_short/execution_report.html \
    -with-timeline results_short/execution_timeline.html

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "📊 Results Summary"
echo "=========================================="

if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Short reads analysis completed!"
    echo ""
    echo "Results: results_short/"
    
    if [ -d "results_short/Taxonomy" ]; then
        echo ""
        echo "📊 Viral classification results:"
        REPORTS=$(find results_short/Taxonomy -name "*.kraken2.report.txt" 2>/dev/null | wc -l)
        echo "   - Kraken2 reports: $REPORTS"
        find results_short/Taxonomy -name "*.kraken2.report.txt" -exec echo "     - {}" \;
    fi
    
    if [ -d "results_short/Assembly" ]; then
        echo ""
        echo "📊 Assembly results:"
        find results_short/Assembly -name "*.fa.gz" -exec ls -lh {} \;
    fi
else
    echo "❌ Short reads analysis failed!"
    echo "Exit code: $EXIT_CODE"
fi

echo ""
echo "End time: $(date)"
echo "=========================================="
