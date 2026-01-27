#!/bin/bash
#SBATCH --job-name=Long_Read_Metagenome
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Long_Read_Metagenome_%j.out
#SBATCH --error=Long_Read_Metagenome_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

echo "=========================================="
echo "🧬  Long-Read Metagenome Workflow"
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

# Verify viralFlye_env exists
if conda env list | grep -q "viralFlye_env"; then
    echo "✅ viralFlye_env: Found"
else
    echo "❌ viralFlye_env not found!"
    echo "   Please create viralFlye_env and install viralFlye:"
    echo "   conda create -n viralFlye_env python=3.8"
    echo "   conda activate viralFlye_env"
    echo "   git clone https://github.com/Dmitry-Antipov/viralFlye.git"
    echo "   cd viralFlye && pip install -r requirements.txt"
    exit 1
fi

echo ""
echo "ℹ️  Workflow execution environment:"
echo "   - Assembly (metaFlye): Apptainer"
echo "   - Mapping (Minimap2, Samtools): Conda + auto symlink"
echo "   - Abundance calculation: Conda + auto symlink"
echo "   - Classification (Kraken2): Conda"
echo ""

# Set database paths
KRAKEN2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref"

# Set long-read samplesheet path
if [ -f "samplesheet_long.csv" ]; then
    SAMPLESHEET_LONG="samplesheet_long.csv"
elif [ -f "/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/data/samplesheet_long.csv" ]; then
    SAMPLESHEET_LONG="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/data/samplesheet_long.csv"
else
    echo "❌ Long-read samplesheet not found!"
    echo "   Checked: ./samplesheet_long.csv"
    echo "   Checked: /scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/data/samplesheet_long.csv"
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
if [ -f "$SAMPLESHEET_LONG" ]; then
    echo "✅ Long-read samplesheet: $SAMPLESHEET_LONG"
    LONG_SAMPLES=$(tail -n +2 "$SAMPLESHEET_LONG" | wc -l)
    echo "   📊 Found $LONG_SAMPLES long-read samples"
else
    echo "❌ Samplesheet not found"
    exit 1
fi
echo ""

# Clean previous results
echo "🧹 5. Cleaning previous results..."
if [ -d "results_long" ]; then
    echo "Removing previous long-read results..."
    rm -rf results_long
fi
echo ""

# Set Singularity bind paths
export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases:/databases"

# Run workflow
echo "🚀 6. Running long-read workflow..."
echo ""
echo "📝 Workflow steps:"
echo "   1. metaFlye assembly (general metagenome)"
echo "   2. viralFlye identify viral contigs (from metaFlye results)"
echo "      - Linear viral contigs"
echo "      - Circular viral contigs"  
echo "   3. Minimap2 mapping and abundance (RPM/RPKM)"
echo "      - metaFlye all contigs"
echo "      - viralFlye linear viral contigs"
echo "      - viralFlye circular viral contigs"
echo "   4. Kraken2 classification for all contig sets"
echo ""

nextflow run metagenome_hybrid_workflow.nf \
    -c metagenome_hybrid_workflow.config \
    --input_long "$SAMPLESHEET_LONG" \
    --outdir_long results_long \
    --kraken2_db "$KRAKEN2_DB"

WORKFLOW_EXIT_CODE=$?

# Check results
echo ""
echo "=========================================="
echo "🎯 Workflow Results"
echo "=========================================="

if [ $WORKFLOW_EXIT_CODE -eq 0 ]; then
    echo "✅ Workflow completed successfully!"
    
    if [ -d "results_long" ]; then
        echo ""
        echo "📊 Long-read results (results_long/):"
        
        if [ -d "results_long/abundance_flye" ]; then
            ABUND_F=$(find results_long/abundance_flye -name "*_abundance.txt" 2>/dev/null | wc -l)
            echo "  ✅ metaFlye abundance: $ABUND_F files"
        fi
        
        if [ -d "results_long/viralflye" ]; then
            echo "  ✅ viralFlye identified viral contigs"
        fi
        
        if [ -d "results_long/abundance_viralflye_linear" ]; then
            ABUND_VFL=$(find results_long/abundance_viralflye_linear -name "*_abundance.txt" 2>/dev/null | wc -l)
            echo "  ✅ viralFlye linear viral abundance: $ABUND_VFL files"
        fi
        
        if [ -d "results_long/abundance_viralflye_circular" ]; then
            ABUND_VFC=$(find results_long/abundance_viralflye_circular -name "*_abundance.txt" 2>/dev/null | wc -l)
            echo "  ✅ viralFlye circular viral abundance: $ABUND_VFC files"
        fi
        
        if [ -d "results_long/kraken2_flye" ]; then
            KRAKEN_F=$(find results_long/kraken2_flye -name "*_report.txt" 2>/dev/null | wc -l)
            echo "  ✅ Kraken2 metaFlye classification: $KRAKEN_F files"
        fi
        
        if [ -d "results_long/kraken2_viralflye_linear" ]; then
            KRAKEN_VFL=$(find results_long/kraken2_viralflye_linear -name "*_report.txt" 2>/dev/null | wc -l)
            echo "  ✅ Kraken2 linear viral classification: $KRAKEN_VFL files"
        fi
        
        if [ -d "results_long/kraken2_viralflye_circular" ]; then
            KRAKEN_VFC=$(find results_long/kraken2_viralflye_circular -name "*_report.txt" 2>/dev/null | wc -l)
            echo "  ✅ Kraken2 circular viral classification: $KRAKEN_VFC files"
        fi
        
        echo ""
        echo "📋 Summary:"
        echo "  Total files: $(find results_long -type f 2>/dev/null | wc -l)"
        echo ""
        echo "📁 Results location: results_long/"
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
