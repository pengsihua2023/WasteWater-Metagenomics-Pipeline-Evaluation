#!/bin/bash
#SBATCH --job-name=Viral_Classification_ShortRead
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Viral_Classification_ShortRead_%j.out
#SBATCH --error=Viral_Classification_ShortRead_%j.err

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "🦠  Metagenome Viral Classification Workflow (Short-Read Mode)"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0

# User's conda environment path
USER_CONDA_ENV="/home/sp96859/.conda/envs/nextflow_env"

# Get conda base path
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo "⚠️  Warning: conda info --base failed, trying alternative method..."
    CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
fi

echo "   Conda base: $CONDA_BASE"
echo "   Target env: $USER_CONDA_ENV"

# Initialize conda
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "❌ Cannot find conda.sh at $CONDA_BASE/etc/profile.d/conda.sh"
    exit 1
fi

# Check if user environment exists
if [ ! -d "$USER_CONDA_ENV" ]; then
    echo "❌ Conda environment not found: $USER_CONDA_ENV"
    echo "   Available environments:"
    conda env list
    exit 1
fi

# Activate environment using absolute path
conda activate "$USER_CONDA_ENV"

# Force update PATH
export PATH="$USER_CONDA_ENV/bin:$PATH"

# Set conda-related environment variables
export CONDA_PREFIX="$USER_CONDA_ENV"
export CONDA_DEFAULT_ENV="nextflow_env"

# Get Python version and set PYTHONPATH
PYTHON_VERSION=$("$USER_CONDA_ENV/bin/python" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "3.9")
export PYTHONPATH="$USER_CONDA_ENV/lib/python${PYTHON_VERSION}/site-packages:${PYTHONPATH:-}"

# Verify environment activation
PYTHON_PATH=$(which python)
echo "   After PATH update:"
echo "   - Python path: $PYTHON_PATH"
echo "   - CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV"

if [[ "$PYTHON_PATH" == *"$USER_CONDA_ENV"* ]]; then
    echo "✅ Conda environment activated successfully!"
    echo "   Environment: nextflow_env"
    echo "   Python: $PYTHON_PATH"
else
    echo "❌ Failed to activate user conda environment!"
    echo "   Expected Python in: $USER_CONDA_ENV"
    echo "   Actual Python: $PYTHON_PATH"
    exit 1
fi

# Verify tools
echo "🧪 2. Verifying tools..."
echo "✅ Nextflow: $(which nextflow)"

# Check for Apptainer/Singularity
if command -v apptainer &> /dev/null; then
    echo "✅ Apptainer: $(which apptainer)"
elif command -v singularity &> /dev/null; then
    echo "✅ Singularity: $(which singularity)"
else
    echo "❌ Apptainer/Singularity not found (required for containers)"
    exit 1
fi

echo ""
echo "ℹ️  Note: Short-read workflow execution environment"
echo "   - Quality control (fastp): Conda environment"
echo "   - Assembly tools (MEGAHIT, metaSPAdes): Apptainer containers"
echo "   - Viral identification tools:"
echo "     * VirSorter2: Pre-installed in nextflow_env ✅"
echo "     * DeepVirFinder: Pre-installed in dvf environment ✅"
echo ""

# Set database paths
VIRSORTER2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/virsorter2/db"
DEEPVIRFINDER_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/DeepVirFinder"

# Verify databases and tools
echo "🗄️ 3. Verifying databases and tools..."

# VirSorter2 database
if [ -d "$VIRSORTER2_DB" ]; then
    echo "✅ VirSorter2 database: $VIRSORTER2_DB"
    echo "   Database size: $(du -sh $VIRSORTER2_DB | cut -f1)"
else
    echo "❌ VirSorter2 database not found: $VIRSORTER2_DB"
    exit 1
fi

# DeepVirFinder installation
if [ -d "$DEEPVIRFINDER_DIR" ] && [ -f "$DEEPVIRFINDER_DIR/dvf.py" ]; then
    echo "✅ DeepVirFinder: $DEEPVIRFINDER_DIR"
else
    echo "❌ DeepVirFinder not found: $DEEPVIRFINDER_DIR"
    exit 1
fi

echo ""

# Verify input files
echo "📁 4. Verifying input files..."
if [ -f "samplesheet_short.csv" ]; then
    SAMPLESHEET="samplesheet_short.csv"
    echo "✅ Samplesheet: samplesheet_short.csv"
    echo "📊 Found $(tail -n +2 samplesheet_short.csv | grep -v '^$' | wc -l) samples"
elif [ -f "samplesheet.csv" ]; then
    SAMPLESHEET="samplesheet.csv"
    echo "✅ Samplesheet: samplesheet.csv"
    echo "📊 Found $(tail -n +2 samplesheet.csv | grep -v '^$' | wc -l) samples"
else
    echo "❌ Samplesheet not found: samplesheet_short.csv or samplesheet.csv"
    echo "   Please create samplesheet_short.csv with format:"
    echo "   sample,fastq_1,fastq_2"
    echo "   sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz"
    exit 1
fi

# Clean previous results and Nextflow cache
echo "🧹 5. Cleaning previous results and Nextflow cache..."
if [ -d "results_short" ]; then
    echo "Removing previous results directory..."
    rm -rf results_short
fi
# Clean Nextflow work directory to force fresh execution (optional)
if [ -d ".nextflow" ]; then
    echo "Clearing Nextflow cache..."
    rm -rf .nextflow/cache
fi

# Set Singularity/Apptainer bind paths and ensure it's enabled
export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases:/databases"
# Ensure Apptainer/Singularity is enabled for Nextflow (don't disable it)
unset NXF_DISABLE_SINGULARITY

# Build nextflow command
NEXTFLOW_CMD="nextflow run metagenome_assembly_classification_workflow.nf \
    -c metagenome_assembly_classification.config \
    --input $SAMPLESHEET \
    --outdir results_short \
    --virsorter2_db \"$VIRSORTER2_DB\" \
    --deepvirfinder_dir \"$DEEPVIRFINDER_DIR\" \
    --longread false"

# Run workflow
echo "🚀 6. Running Metagenome Viral Classification workflow (Short-Read Mode)..."
echo "Command: $NEXTFLOW_CMD"
echo ""
echo "📝 Workflow steps:"
echo "   1. fastp quality control (auto adapter removal, low-quality read filtering)"
echo "   2. MEGAHIT and metaSPAdes parallel assembly"
echo "   3. VirSorter2 viral sequence identification"
echo "   4. DeepVirFinder viral prediction"
echo "   5. Tool result merging (VirSorter2 + DeepVirFinder per assembler)"
echo "   6. Assembler comparison (MEGAHIT vs SPAdes) → Final consensus viral list ⭐"
echo "   7. Viral abundance calculation (RPM and RPKM) for all identified viruses ⭐⭐⭐"
echo ""
echo "✅ Note: Using short-read mode ( Illumina paired-end reads )"
echo "✅ Note: Full run from scratch (no resume mode)"
echo ""

eval $NEXTFLOW_CMD

# Check results
echo ""
echo "=========================================="
echo "🎯 Workflow Results"
echo "=========================================="

if [ $? -eq 0 ]; then
    echo "✅ Workflow completed successfully!"
    
    if [ -d "results_short" ]; then
        echo "📁 Results directory created: results_short/"
        echo "📊 Generated results:"
        
        # Check fastp results
        if [ -d "results_short/fastp" ]; then
            echo "  ✅ fastp quality reports: results_short/fastp/"
            FASTP_HTML=$(find results_short/fastp -name "*.html" | wc -l)
            echo "     - Generated $FASTP_HTML HTML quality reports"
        fi
        
        # Check clean reads
        if [ -d "results_short/clean_reads" ]; then
            echo "  ✅ Clean reads: results_short/clean_reads/"
            CLEAN_READS=$(find results_short/clean_reads -name "*_clean_R*.fastq.gz" | wc -l)
            echo "     - Generated $CLEAN_READS clean read files"
        fi
        
        # Check assembly results
        if [ -d "results_short/assembly_megahit" ]; then
            echo "  ✅ MEGAHIT assembly: results_short/assembly_megahit/"
            MEGAHIT_CONTIGS=$(find results_short/assembly_megahit -name "*_megahit_contigs.fa" | wc -l)
            echo "     - Generated $MEGAHIT_CONTIGS contig files"
        fi
        
        if [ -d "results_short/assembly_spades" ]; then
            echo "  ✅ SPAdes assembly: results_short/assembly_spades/"
            SPADES_CONTIGS=$(find results_short/assembly_spades -name "*_spades_contigs.fa" | wc -l)
            echo "     - Generated $SPADES_CONTIGS contig files"
        fi
        
        # Check VirSorter2 results
        if [ -d "results_short/virsorter2_megahit" ]; then
            echo "  ✅ VirSorter2 MEGAHIT results: results_short/virsorter2_megahit/"
            VS2_MEGAHIT=$(find results_short/virsorter2_megahit -name "*_vs2_final-viral-score.tsv" | wc -l)
            echo "     - Generated $VS2_MEGAHIT viral identification reports"
        fi
        
        if [ -d "results_short/virsorter2_spades" ]; then
            echo "  ✅ VirSorter2 SPAdes results: results_short/virsorter2_spades/"
            VS2_SPADES=$(find results_short/virsorter2_spades -name "*_vs2_final-viral-score.tsv" | wc -l)
            echo "     - Generated $VS2_SPADES viral identification reports"
        fi
        
        # Check DeepVirFinder results
        if [ -d "results_short/deepvirfinder_megahit" ]; then
            echo "  ✅ DeepVirFinder MEGAHIT results: results_short/deepvirfinder_megahit/"
            DVF_MEGAHIT=$(find results_short/deepvirfinder_megahit -name "*_dvf_output.txt" | wc -l)
            echo "     - Generated $DVF_MEGAHIT viral prediction reports"
        fi
        
        if [ -d "results_short/deepvirfinder_spades" ]; then
            echo "  ✅ DeepVirFinder SPAdes results: results_short/deepvirfinder_spades/"
            DVF_SPADES=$(find results_short/deepvirfinder_spades -name "*_dvf_output.txt" | wc -l)
            echo "     - Generated $DVF_SPADES viral prediction reports"
        fi
        
        # Check merged viral reports
        if [ -d "results_short/merged_viral_reports_megahit" ]; then
            echo "  ✅ Merged viral reports (MEGAHIT): results_short/merged_viral_reports_megahit/"
            MERGED_VIRAL_MEGAHIT=$(find results_short/merged_viral_reports_megahit -name "*_viral_merged_report.txt" | wc -l)
            echo "     - Generated $MERGED_VIRAL_MEGAHIT comprehensive viral analysis reports"
        fi
        
        if [ -d "results_short/merged_viral_reports_spades" ]; then
            echo "  ✅ Merged viral reports (SPAdes): results_short/merged_viral_reports_spades/"
            MERGED_VIRAL_SPADES=$(find results_short/merged_viral_reports_spades -name "*_viral_merged_report.txt" | wc -l)
            echo "     - Generated $MERGED_VIRAL_SPADES comprehensive viral analysis reports"
        fi
        
        # Check assembler comparison
        if [ -d "results_short/assembler_comparison" ]; then
            echo "  ✅ Assembler comparison (MEGAHIT vs SPAdes): results_short/assembler_comparison/"
            ASSEMBLER_COMP=$(find results_short/assembler_comparison -name "*_assembler_comparison.txt" | wc -l)
            echo "     - Generated $ASSEMBLER_COMP assembler comparison reports"
            CONSENSUS_SEQS=$(find results_short/assembler_comparison -name "*_consensus_viral_sequences.txt" | wc -l)
            echo "     - Generated $CONSENSUS_SEQS final consensus viral sequence lists"
        fi
        
        # Check viral abundance results
        if [ -d "results_short/abundance" ]; then
            echo "  ⭐⭐⭐ Viral abundance analysis: results_short/abundance/"
            
            if [ -d "results_short/abundance/megahit" ]; then
                ABUNDANCE_MEGAHIT=$(find results_short/abundance/megahit -name "*_abundance.csv" | wc -l)
                echo "     - MEGAHIT viral abundance: $ABUNDANCE_MEGAHIT CSV files (RPM & RPKM)"
            fi
            
            if [ -d "results_short/abundance/spades" ]; then
                ABUNDANCE_SPADES=$(find results_short/abundance/spades -name "*_abundance.csv" | wc -l)
                echo "     - SPAdes viral abundance: $ABUNDANCE_SPADES CSV files (RPM & RPKM)"
            fi
        fi
        
        echo ""
        echo "📋 Summary of key viral identification files:"
        echo "  VirSorter2 viral scores:"
        find results_short/virsorter2_* -name "*_vs2_final-viral-score.tsv" 2>/dev/null | head -10
        echo ""
        echo "  DeepVirFinder predictions:"
        find results_short/deepvirfinder_* -name "*_dvf_output.txt" 2>/dev/null | head -10
        echo ""
        echo "  Merged viral reports:"
        find results_short/merged_viral_reports_* -name "*.txt" -o -name "*.csv" 2>/dev/null | head -10
        echo ""
        echo "  Assembler comparison (MEGAHIT vs SPAdes):"
        find results_short/assembler_comparison -name "*_assembler_comparison.txt" 2>/dev/null | head -10
        echo ""
        echo "  ⭐ Final consensus viral sequences (recommended for downstream analysis):"
        find results_short/assembler_comparison -name "*_consensus_viral_sequences.txt" 2>/dev/null | head -10
        echo ""
        echo "  ⭐⭐⭐ Viral abundance analysis (RPM & RPKM):"
        find results_short/abundance -name "*_abundance.csv" 2>/dev/null | head -10
        echo ""
        echo "  📊 Abundance summary reports:"
        find results_short/abundance -name "*_abundance_summary.txt" 2>/dev/null | head -10
        echo ""
        echo "Total files: $(find results_short -type f | wc -l)"
        
    else
        echo "❌ Results directory not found"
    fi
    
else
    echo "❌ Workflow failed with exit code: $?"
    echo "🔍 Check the error log for details"
fi

echo ""
echo "End time: $(date)"
echo "=========================================="
