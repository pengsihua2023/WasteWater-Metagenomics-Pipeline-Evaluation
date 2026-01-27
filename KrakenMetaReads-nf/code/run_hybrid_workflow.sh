#!/bin/bash
#SBATCH --job-name=Hybrid_Metagenome
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Hybrid_Metagenome_%j.out
#SBATCH --error=Hybrid_Metagenome_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

echo "=========================================="
echo "🧬  Hybrid Metagenome Assembly Workflow"
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

# Verify viralFlye_env if processing long reads
if [ -f "samplesheet_long.csv" ] || [ -f "/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/data/samplesheet_long.csv" ]; then
    if conda env list | grep -q "viralFlye_env"; then
        echo "✅ viralFlye_env: Found (required for long-read analysis)"
    else
        echo "⚠️  viralFlye_env not found (needed for viralFlye viral identification)"
        echo "   Long-read workflow will fail if viralFlye is enabled"
        echo "   To disable viralFlye, edit metagenome_hybrid_workflow.config: run_viralflye = false"
    fi
fi

echo ""
echo "ℹ️  Workflow execution environment:"
echo "   Short reads:"
echo "     - Quality control (fastp): Conda"
echo "     - Assembly (MEGAHIT, SPAdes): Apptainer"
echo "     - Mapping (Bowtie2): Apptainer"
echo "     - Abundance calculation: Conda + auto symlink"
echo "     - Classification (Kraken2): Conda"
echo "   Long reads:"
echo "     - Assembly (metaFlye): Apptainer"
echo "     - Mapping (Minimap2, Samtools): Conda + auto symlink"
echo "     - Abundance calculation: Conda + auto symlink"
echo "     - Classification (Kraken2): Conda"
echo ""

# Set database paths
KRAKEN2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref"

# Set samplesheet paths (check current directory first, then absolute paths)
if [ -f "samplesheet_short.csv" ]; then
    SAMPLESHEET_SHORT="samplesheet_short.csv"
elif [ -f "/scratch/sp96859/Meta-genome-data-analysis/Apptainer/yitiaolong/data/reads/samplesheet_short.csv" ]; then
    SAMPLESHEET_SHORT="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/yitiaolong/data/reads/samplesheet_short.csv"
else
    SAMPLESHEET_SHORT=""
fi

if [ -f "samplesheet_long.csv" ]; then
    SAMPLESHEET_LONG="samplesheet_long.csv"
elif [ -f "/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/data/samplesheet_long.csv" ]; then
    SAMPLESHEET_LONG="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/data/samplesheet_long.csv"
else
    SAMPLESHEET_LONG=""
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

# Check short-read samplesheet
if [ -n "$SAMPLESHEET_SHORT" ] && [ -f "$SAMPLESHEET_SHORT" ]; then
    echo "✅ Short-read samplesheet: $SAMPLESHEET_SHORT"
    SHORT_SAMPLES=$(tail -n +2 "$SAMPLESHEET_SHORT" | wc -l)
    echo "   📊 Found $SHORT_SAMPLES short-read samples"
    PROCESS_SHORT=true
else
    echo "⚠️  Short-read samplesheet not found, skipping short-read analysis"
    echo "   Checked: ./samplesheet_short.csv"
    echo "   Checked: /scratch/sp96859/Meta-genome-data-analysis/Apptainer/yitiaolong/data/reads/samplesheet_short.csv"
    PROCESS_SHORT=false
fi

# Check long-read samplesheet
if [ -n "$SAMPLESHEET_LONG" ] && [ -f "$SAMPLESHEET_LONG" ]; then
    echo "✅ Long-read samplesheet: $SAMPLESHEET_LONG"
    LONG_SAMPLES=$(tail -n +2 "$SAMPLESHEET_LONG" | wc -l)
    echo "   📊 Found $LONG_SAMPLES long-read samples"
    PROCESS_LONG=true
else
    echo "⚠️  Long-read samplesheet not found, skipping long-read analysis"
    echo "   Checked: ./samplesheet_long.csv"
    echo "   Checked: /scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/data/samplesheet_long.csv"
    PROCESS_LONG=false
fi

if [ "$PROCESS_SHORT" = false ] && [ "$PROCESS_LONG" = false ]; then
    echo ""
    echo "❌ No valid samplesheets found!"
    echo ""
    echo "💡 Please create at least one of the following files:"
    echo "   - samplesheet_short.csv (in current directory)"
    echo "   - samplesheet_long.csv (in current directory)"
    echo ""
    exit 1
fi
echo ""

# Clean previous results
echo "🧹 5. Cleaning previous results..."
if [ "$PROCESS_SHORT" = true ] && [ -d "results_short" ]; then
    echo "Removing previous short-read results..."
    rm -rf results_short
fi
if [ "$PROCESS_LONG" = true ] && [ -d "results_long" ]; then
    echo "Removing previous long-read results..."
    rm -rf results_long
fi
echo ""

# Set Singularity bind paths
export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases:/databases"

# Build command
echo "🚀 6. Running hybrid metagenome workflow..."
echo ""

CMD="nextflow run metagenome_hybrid_workflow.nf \
    -c metagenome_hybrid_workflow.config \
    --kraken2_db \"$KRAKEN2_DB\""

if [ "$PROCESS_SHORT" = true ]; then
    CMD="$CMD --input_short \"$SAMPLESHEET_SHORT\" --outdir_short results_short"
fi

if [ "$PROCESS_LONG" = true ]; then
    CMD="$CMD --input_long \"$SAMPLESHEET_LONG\" --outdir_long results_long"
fi

echo "📝 Workflow steps:"
if [ "$PROCESS_SHORT" = true ]; then
    echo "   Short reads:"
    echo "     1. fastp quality control"
    echo "     2. MEGAHIT and metaSPAdes assembly"
    echo "     3. Bowtie2 mapping and abundance (RPM/RPKM)"
    echo "     4. Kraken2 classification"
    echo "     5. Comprehensive report generation"
fi

if [ "$PROCESS_LONG" = true ]; then
    echo "   Long reads:"
    echo "     1. metaFlye assembly"
    echo "     2. viralFlye viral identification (linear + circular)"
    echo "     3. Minimap2 mapping and abundance (RPM/RPKM)"
    echo "     4. Kraken2 classification"
fi
echo ""

echo "Command: $CMD"
echo ""

eval $CMD
WORKFLOW_EXIT_CODE=$?

# Check results
echo ""
echo "=========================================="
echo "🎯 Workflow Results"
echo "=========================================="

if [ $WORKFLOW_EXIT_CODE -eq 0 ]; then
    echo "✅ Workflow completed successfully!"
    
    # Check short-read results
    if [ "$PROCESS_SHORT" = true ] && [ -d "results_short" ]; then
        echo ""
        echo "📊 Short-read results (results_short/):"
        
        if [ -d "results_short/fastp" ]; then
            FASTP_COUNT=$(find results_short/fastp -name "*.html" | wc -l)
            echo "  ✅ fastp reports: $FASTP_COUNT files"
        fi
        
        if [ -d "results_short/abundance_megahit" ]; then
            ABUND_M=$(find results_short/abundance_megahit -name "*_abundance.txt" | wc -l)
            echo "  ✅ MEGAHIT abundance: $ABUND_M files"
        fi
        
        if [ -d "results_short/abundance_spades" ]; then
            ABUND_S=$(find results_short/abundance_spades -name "*_abundance.txt" | wc -l)
            echo "  ✅ SPAdes abundance: $ABUND_S files"
        fi
        
        if [ -d "results_short/kraken2_megahit" ]; then
            echo "  ✅ Kraken2 MEGAHIT classification"
        fi
        
        if [ -d "results_short/kraken2_spades" ]; then
            echo "  ✅ Kraken2 SPAdes classification"
        fi
        
        if [ -d "results_short/merged_reports" ]; then
            MERGED=$(find results_short/merged_reports -name "*_merged_report.csv" | wc -l)
            CONSENSUS=$(find results_short/merged_reports -name "*_virus_consensus.txt" | wc -l)
            echo "  ✅ Merged reports: $MERGED files"
            if [ "$CONSENSUS" -gt 0 ]; then
                echo "  ✅ Virus consensus analysis: $CONSENSUS files ⭐"
            fi
        fi
    fi
    
    # Check long-read results
    if [ "$PROCESS_LONG" = true ] && [ -d "results_long" ]; then
        echo ""
        echo "📊 Long-read results (results_long/):"
        
        if [ -d "results_long/abundance_flye" ]; then
            ABUND_F=$(find results_long/abundance_flye -name "*_abundance.txt" 2>/dev/null | wc -l)
            echo "  ✅ metaFlye abundance: $ABUND_F files"
        fi
        
        if [ -d "results_long/viralflye" ]; then
            echo "  ✅ viralFlye viral identification"
        fi
        
        if [ -d "results_long/abundance_viralflye_linear" ]; then
            echo "  ✅ viralFlye linear viral abundance"
        fi
        
        if [ -d "results_long/abundance_viralflye_circular" ]; then
            echo "  ✅ viralFlye circular viral abundance"
        fi
        
        if [ -d "results_long/kraken2_flye" ]; then
            echo "  ✅ Kraken2 metaFlye classification"
        fi
        
        if [ -d "results_long/kraken2_viralflye_linear" ]; then
            echo "  ✅ Kraken2 linear viral classification"
        fi
        
        if [ -d "results_long/kraken2_viralflye_circular" ]; then
            echo "  ✅ Kraken2 circular viral classification"
        fi
    fi
    
    echo ""
    echo "📋 Summary:"
    echo "  Short-read files: $(find results_short -type f 2>/dev/null | wc -l)"
    echo "  Long-read files: $(find results_long -type f 2>/dev/null | wc -l)"
    
else
    echo "❌ Workflow failed with exit code: $WORKFLOW_EXIT_CODE"
    echo "🔍 Check the error log for details"
fi

echo ""
echo "End time: $(date)"
echo "=========================================="
