#!/bin/bash
#SBATCH --job-name=Short_Read_Metagenome
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Short_Read_Metagenome_%j.out
#SBATCH --error=Short_Read_Metagenome_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

echo "=========================================="
echo "🧬  Short-Read Metagenome Workflow"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate nextflow_env

# Verify tools
echo "🧪 2. Verifying tools..."
NEXTFLOW_PATH=$(which nextflow)
echo "✅ Nextflow: $NEXTFLOW_PATH"

if command -v apptainer &> /dev/null; then
    echo "✅ Apptainer: $(which apptainer)"
elif command -v singularity &> /dev/null; then
    echo "✅ Singularity: $(which singularity)"
else
    echo "❌ Apptainer/Singularity not found"
    exit 1
fi

echo ""
echo "ℹ️  Workflow: Short-read only (MEGAHIT + SPAdes + Bowtie2 + Kraken2)"
echo ""

# Set database paths
KRAKEN2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref"

# Set short-read samplesheet path
if [ -f "samplesheet_short.csv" ]; then
    SAMPLESHEET_SHORT="samplesheet_short.csv"
elif [ -f "/scratch/sp96859/Meta-genome-data-analysis/Apptainer/yitiaolong/data/reads/samplesheet_short.csv" ]; then
    SAMPLESHEET_SHORT="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/yitiaolong/data/reads/samplesheet_short.csv"
else
    echo "❌ Short-read samplesheet not found!"
    echo "   Checked: ./samplesheet_short.csv"
    echo "   Checked: /scratch/sp96859/Meta-genome-data-analysis/Apptainer/yitiaolong/data/reads/samplesheet_short.csv"
    exit 1
fi

# Verify databases
echo "🗄️ 3. Verifying databases..."
if [ -d "$KRAKEN2_DB" ]; then
    echo "✅ Kraken2 database: $KRAKEN2_DB"
else
    echo "❌ Kraken2 database not found: $KRAKEN2_DB"
    exit 1
fi
echo ""

# Verify input files
echo "📁 4. Verifying input files..."
if [ -f "$SAMPLESHEET_SHORT" ]; then
    echo "✅ Short-read samplesheet: $SAMPLESHEET_SHORT"
    SHORT_SAMPLES=$(tail -n +2 "$SAMPLESHEET_SHORT" | wc -l)
    echo "   📊 Found $SHORT_SAMPLES short-read samples"
else
    echo "❌ Samplesheet not found"
    exit 1
fi
echo ""

# Clean previous results
echo "🧹 5. Cleaning previous results..."
if [ -d "results_short" ]; then
    echo "Removing previous short-read results..."
    rm -rf results_short
fi
echo ""

# Set Singularity bind paths
export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases:/databases"

# Run workflow
echo "🚀 6. Running short-read workflow..."
echo ""
echo "📝 Workflow steps:"
echo "   1. fastp quality control"
echo "   2. MEGAHIT and metaSPAdes assembly"
echo "   3. Bowtie2 mapping and abundance (RPM/RPKM)"
echo "   4. Kraken2 classification"
echo "   5. Comprehensive report generation"
echo ""

nextflow run metagenome_hybrid_workflow.nf \
    -c metagenome_hybrid_workflow.config \
    --input_short "$SAMPLESHEET_SHORT" \
    --outdir_short results_short \
    --kraken2_db "$KRAKEN2_DB"

WORKFLOW_EXIT_CODE=$?

# Check results
echo ""
echo "=========================================="
echo "🎯 Workflow Results"
echo "=========================================="

if [ $WORKFLOW_EXIT_CODE -eq 0 ]; then
    echo "✅ Workflow completed successfully!"
    
    if [ -d "results_short" ]; then
        echo ""
        echo "📊 Short-read results (results_short/):"
        
        if [ -d "results_short/fastp" ]; then
            FASTP_COUNT=$(find results_short/fastp -name "*.html" 2>/dev/null | wc -l)
            echo "  ✅ fastp reports: $FASTP_COUNT files"
        fi
        
        if [ -d "results_short/abundance_megahit" ]; then
            ABUND_M=$(find results_short/abundance_megahit -name "*_abundance.txt" 2>/dev/null | wc -l)
            echo "  ✅ MEGAHIT abundance: $ABUND_M files"
        fi
        
        if [ -d "results_short/abundance_spades" ]; then
            ABUND_S=$(find results_short/abundance_spades -name "*_abundance.txt" 2>/dev/null | wc -l)
            echo "  ✅ SPAdes abundance: $ABUND_S files"
        fi
        
        if [ -d "results_short/kraken2_megahit" ]; then
            echo "  ✅ Kraken2 MEGAHIT classification"
        fi
        
        if [ -d "results_short/kraken2_spades" ]; then
            echo "  ✅ Kraken2 SPAdes classification"
        fi
        
        if [ -d "results_short/merged_reports" ]; then
            MERGED=$(find results_short/merged_reports -name "*_merged_report.csv" 2>/dev/null | wc -l)
            CONSENSUS=$(find results_short/merged_reports -name "*_virus_consensus.txt" 2>/dev/null | wc -l)
            echo "  ✅ Merged reports: $MERGED files"
            if [ "$CONSENSUS" -gt 0 ]; then
                echo "  ✅ Virus consensus analysis: $CONSENSUS files ⭐"
            fi
        fi
        
        echo ""
        echo "📋 Summary:"
        echo "  Total files: $(find results_short -type f 2>/dev/null | wc -l)"
        echo ""
        echo "📁 Results location: results_short/"
    else
        echo "❌ Results directory not found"
    fi
    
else
    echo "❌ Workflow failed with exit code: $WORKFLOW_EXIT_CODE"
    echo "🔍 Check the error log for details"
fi

echo ""
echo "End time: $(date)"
echo "=========================================="
