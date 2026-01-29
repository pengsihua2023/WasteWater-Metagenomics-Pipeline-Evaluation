#!/bin/bash
#SBATCH --job-name=Metagenome_Assembly_Classification
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Metagenome_Assembly_Classification_%j.out
#SBATCH --error=Metagenome_Assembly_Classification_%j.err

# Usage:
#   sbatch run_metagenome_assembly_classification_en.sh [read_type] [samplesheet]
#
# Arguments:
#   read_type    - 'short' for short reads (Illumina paired-end) or 'long' for long reads (Nanopore/PacBio)
#                  Default: 'short'
#   samplesheet  - Path to samplesheet CSV file (optional)
#                  If not provided, automatically uses 'samplesheet_short.csv' for short reads
#                  or 'samplesheet_long.csv' for long reads
#
# Examples:
#   # Short reads (uses samplesheet_short.csv automatically)
#   sbatch run_metagenome_assembly_classification_en.sh short
#
#   # Long reads (uses samplesheet_long.csv automatically)
#   sbatch run_metagenome_assembly_classification_en.sh long
#
#   # Custom samplesheet
#   sbatch run_metagenome_assembly_classification_en.sh short my_samplesheet.csv

cd "$SLURM_SUBMIT_DIR" || exit 1

echo "=========================================="
echo "🧬 Metagenome Assembly and Classification Workflow"
echo "   Supports both short (Illumina) and long (Nanopore/PacBio) reads"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Determine read type based on command line argument
READ_TYPE="${1:-short}"  # Default to 'short' if not specified

# Auto-select samplesheet based on read type if not provided
if [ -z "$2" ]; then
    # No samplesheet provided, auto-select based on read type
    if [ "$READ_TYPE" == "long" ]; then
        SAMPLESHEET="samplesheet_long.csv"
    else
        SAMPLESHEET="samplesheet_short.csv"
    fi
else
    # Use provided samplesheet
    SAMPLESHEET="$2"
fi

# Set output directory based on read type
if [ "$READ_TYPE" == "long" ]; then
    OUTDIR="results_long"
else
    OUTDIR="results_short"
fi

echo "ℹ️  Detected read type: $READ_TYPE"
echo "ℹ️  Using samplesheet: $SAMPLESHEET"
echo "ℹ️  Output directory: $OUTDIR"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow_env

# Verify tools
echo "🧪 2. Verifying tools..."
echo "✅ Nextflow: $(which nextflow)"

# Check for Apptainer/Singularity (optional, not required since using Conda for all tools)
if command -v apptainer &> /dev/null; then
    echo "✅ Apptainer: $(which apptainer) (optional, not required)"
elif command -v singularity &> /dev/null; then
    echo "✅ Singularity: $(which singularity) (optional, not required)"
else
    echo "ℹ️  Apptainer/Singularity not found (not required, all tools use Conda)"
fi

echo ""
echo "ℹ️  Note: Workflow execution environment"
if [ "$READ_TYPE" == "long" ]; then
    echo "   - Long read assembly (MetaFlye): Conda environment (auto-created)"
    echo "   - Refinement (viralFlye, optional): Conda environment (auto-created)"
else
    echo "   - Quality control (fastp): Conda environment (auto-created)"
    echo "   - Assembly tools (MEGAHIT, metaSPAdes): nextflow_env (pre-installed)"
fi
echo "   - Gene prediction (Prodigal): Conda environment (auto-created)"
echo "   - Classification tool (Diamond): Conda environment (auto-created)"
echo "   Note: MEGAHIT and SPAdes use nextflow_env (pre-installed), other tools auto-create Conda environments"
echo ""

# Set database paths
# Diamond database for classification (RVDB - Reference Viral DataBase)
DIAMOND_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB/RVDB_prot_ref.dmnd"

# Pfam-A HMM database for viralFlye (required for long reads)
PFAM_HMM="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/Pfam/Pfam-A.hmm"

# Verify databases
echo "🗄️ 3. Verifying databases..."
if [ -f "$DIAMOND_DB" ]; then
    echo "✅ Diamond database: $DIAMOND_DB"
    echo "   Database size: $(du -h $DIAMOND_DB | cut -f1)"
else
    echo "❌ Diamond database not found: $DIAMOND_DB"
    echo ""
    echo "📝 Note: RVDB database information:"
    echo "   Database: Reference Viral DataBase (RVDB)"
    echo "   Expected file: RVDB_prot_ref.dmnd"
    echo "   Location: /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB/"
    echo ""
    echo "   To obtain RVDB database:"
    echo "   1. Download from: https://rvdb-prot.pasteur.fr/"
    echo "   2. Build Diamond index: diamond makedb --in RVDB.fasta -d RVDB_prot_ref"
    exit 1
fi

# Verify Pfam-A HMM database (required for long reads  )
if [ "$READ_TYPE" == "long" ]; then
    if [ -f "$PFAM_HMM" ]; then
        echo "✅ Pfam-A HMM database: $PFAM_HMM"
        echo "   Database size: $(du -h $PFAM_HMM | cut -f1)"
    else
        echo "⚠️  Pfam-A HMM database not found: $PFAM_HMM"
        echo ""
        echo "📝 Note: This database is required for viralFlye refinement"
        echo "   To download:"
        echo "   mkdir -p /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/Pfam"
        echo "   cd /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/Pfam"
        echo "   wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
        echo "   gunzip Pfam-A.hmm.gz"
        echo ""
        echo "   If you want to skip viralFlye, set skip_viralflye = true in config file"
        exit 1
    fi
    
    # Check viralFlye_env environment
    if conda env list | grep -q "viralFlye_env"; then
        echo "✅ viralFlye_env conda environment found"
    else
        echo "⚠️  viralFlye_env environment not found"
        echo ""
        echo "📝 Note: viralFlye requires a separate conda environment"
        echo "   To create it:"
        echo "   conda create -n viralFlye_env python=3.7"
        echo "   conda activate viralFlye_env"
        echo "   conda install -c bioconda viralflye"
        echo ""
        echo "   Or skip viralFlye by setting skip_viralflye = true in config file"
        exit 1
    fi
    
    # Check NCBI Taxonomy database (for taxonomic resolution)
    TAXONOMY_NAMES="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB/names.dmp"
    TAXONOMY_NODES="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB/nodes.dmp"
    
    if [ -f "$TAXONOMY_NAMES" ] && [ -f "$TAXONOMY_NODES" ]; then
        echo "✅ NCBI Taxonomy database found"
        echo "   names.dmp size: $(du -h $TAXONOMY_NAMES | cut -f1)"
        echo "   nodes.dmp size: $(du -h $TAXONOMY_NODES | cut -f1)"
    else
        echo "⚠️  NCBI Taxonomy database not found"
        echo ""
        echo "📝 Note: Taxonomy database is used to resolve virus names and lineage"
        echo "   Missing files:"
        [ ! -f "$TAXONOMY_NAMES" ] && echo "   - $TAXONOMY_NAMES"
        [ ! -f "$TAXONOMY_NODES" ] && echo "   - $TAXONOMY_NODES"
        echo ""
        echo "   Download from: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
        exit 1
    fi
fi
echo ""

# Verify input files
echo "📁 4. Verifying input files..."
if [ -f "$SAMPLESHEET" ]; then
    echo "✅ Samplesheet: $SAMPLESHEET"
    echo "📊 Found $(wc -l < $SAMPLESHEET) samples"
else
    echo "❌ Samplesheet not found: $SAMPLESHEET"
    exit 1
fi

# Clean previous results
echo "🧹 5. Cleaning previous results..."
if [ -d "$OUTDIR" ]; then
    echo "Removing previous results directory: $OUTDIR"
    rm -rf "$OUTDIR"
fi

# Container binding configuration is no longer needed (all tools now use Conda environments)
# export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases:/databases,/tmp:/lscratch"

# Run workflow
echo "🚀 6. Running Metagenome Assembly and Classification workflow..."
echo "Command: nextflow run metagenome_assembly_classification_workflow_en.nf -c metagenome_assembly_classification_en.config --input $SAMPLESHEET --outdir $OUTDIR --diamond_db $DIAMOND_DB --read_type $READ_TYPE"
echo ""
echo "📝 Workflow steps:"
if [ "$READ_TYPE" == "long" ]; then
    echo "   1. MetaFlye long read assembly"
    echo "   2. viralFlye refinement (optional)"
    echo "   3. Prodigal gene prediction (metagenome mode)"
    echo "   4. Diamond BLASTP classification against protein database"
else
    echo "   1. fastp quality control (auto adapter removal, low-quality read filtering)"
    echo "   2. MEGAHIT and metaSPAdes parallel assembly"
    echo "   3. Prodigal gene prediction (metagenome mode)"
    echo "   4. Diamond BLASTP classification against protein database"
    echo "   5. Comprehensive report generation"
    echo "   6. Viral abundance calculation (RPM and RPKM, if enabled)"
fi
echo ""

nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --diamond_db "$DIAMOND_DB" \
    --read_type "$READ_TYPE"

# Check results
echo ""
echo "=========================================="
echo "🎯 Workflow Results"
echo "=========================================="

if [ $? -eq 0 ]; then
    echo "✅ Workflow completed successfully!"
    
    if [ -d "$OUTDIR" ]; then
        echo "📁 Results directory created: $OUTDIR/"
        echo "📊 Generated results:"
        
        # Check results based on read type
        if [ "$READ_TYPE" == "long" ]; then
            # Long read results
            if [ -d "$OUTDIR/abundance_metaflye" ]; then
                echo "  ✅ Viral abundance (MetaFlye): $OUTDIR/abundance_metaflye/"
                METAFLYE_ABUNDANCE=$(find $OUTDIR/abundance_metaflye -name "*.txt" -o -name "*.csv" | wc -l)
                echo "     - Generated $METAFLYE_ABUNDANCE abundance reports (RPM and RPKM)"
            fi
            
            if [ -d "$OUTDIR/abundance_viralflye" ]; then
                echo "  ✅ Viral abundance (viralFlye): $OUTDIR/abundance_viralflye/"
                VIRALFLYE_ABUNDANCE=$(find $OUTDIR/abundance_viralflye -name "*.txt" -o -name "*.csv" | wc -l)
                echo "     - Generated $VIRALFLYE_ABUNDANCE abundance reports (RPM and RPKM)"
            fi
            if [ -d "$OUTDIR/assembly_metaflye" ]; then
                echo "  ✅ MetaFlye assembly contigs: $OUTDIR/assembly_metaflye/"
                METAFLYE_CONTIGS=$(find $OUTDIR/assembly_metaflye -name "*_contigs.fa" | wc -l)
                echo "     - Generated $METAFLYE_CONTIGS contig files"
            fi
            
            if [ -d "$OUTDIR/assembly_viralflye" ]; then
                echo "  ✅ viralFlye refined contigs: $OUTDIR/assembly_viralflye/"
                VIRALFLYE_CONTIGS=$(find $OUTDIR/assembly_viralflye -name "*_contigs.fa" | wc -l)
                echo "     - Generated $VIRALFLYE_CONTIGS refined contig files"
            fi
            
            if [ -d "$OUTDIR/prodigal_metaflye" ]; then
                echo "  ✅ Prodigal MetaFlye gene predictions (all sequences): $OUTDIR/prodigal_metaflye/"
                METAFLYE_GENES=$(find $OUTDIR/prodigal_metaflye -name "*.faa" | wc -l)
                echo "     - Generated $METAFLYE_GENES protein sequence files"
            fi
            
            if [ -d "$OUTDIR/prodigal_viralflye" ]; then
                echo "  ✅ Prodigal viralFlye gene predictions (viral only): $OUTDIR/prodigal_viralflye/"
                VIRALFLYE_GENES=$(find $OUTDIR/prodigal_viralflye -name "*.faa" | wc -l)
                echo "     - Generated $VIRALFLYE_GENES viral protein sequence files"
            fi
            
            if [ -d "$OUTDIR/diamond_metaflye" ]; then
                echo "  ✅ Diamond MetaFlye results (all sequences): $OUTDIR/diamond_metaflye/"
            fi
            
            if [ -d "$OUTDIR/diamond_viralflye" ]; then
                echo "  ✅ Diamond viralFlye results (viral only): $OUTDIR/diamond_viralflye/"
            fi
            
            if [ -d "$OUTDIR/taxonomy_metaflye" ]; then
                echo "  ✅ Taxonomy analysis (all sequences): $OUTDIR/taxonomy_metaflye/"
                TAXONOMY_META=$(find $OUTDIR/taxonomy_metaflye -name "*_with_taxonomy.txt" | wc -l)
                echo "     - Generated $TAXONOMY_META reports with Kingdom/Class/Family/Species"
            fi
            
            if [ -d "$OUTDIR/taxonomy_viralflye" ]; then
                echo "  ✅ Taxonomy analysis (viral only): $OUTDIR/taxonomy_viralflye/"
                TAXONOMY_VIRAL=$(find $OUTDIR/taxonomy_viralflye -name "*_with_taxonomy.txt" | wc -l)
                echo "     - Generated $TAXONOMY_VIRAL viral-specific reports"
            fi
            
            if [ -d "$OUTDIR/consensus_analysis" ]; then
                echo "  ✅ Dual-track consensus analysis: $OUTDIR/consensus_analysis/"
                CONSENSUS_REPORTS=$(find $OUTDIR/consensus_analysis -name "*_consensus_viruses.txt" | wc -l)
                echo "     - Generated $CONSENSUS_REPORTS consensus virus reports"
                echo "     - Consensus viruses: highest confidence (both methods confirm)"
            fi
        else
            # Short read results
            if [ -d "$OUTDIR/fastp" ]; then
                echo "  ✅ fastp quality reports: $OUTDIR/fastp/"
                FASTP_HTML=$(find $OUTDIR/fastp -name "*.html" | wc -l)
                echo "     - Generated $FASTP_HTML HTML quality reports"
            fi
            
            if [ -d "$OUTDIR/assembly_megahit" ]; then
                echo "  ✅ MEGAHIT assembly contigs: $OUTDIR/assembly_megahit/"
                MEGAHIT_CONTIGS=$(find $OUTDIR/assembly_megahit -name "*_contigs.fa" | wc -l)
                echo "     - Generated $MEGAHIT_CONTIGS contig files"
            fi
            
            if [ -d "$OUTDIR/assembly_spades" ]; then
                echo "  ✅ SPAdes assembly contigs: $OUTDIR/assembly_spades/"
                SPADES_CONTIGS=$(find $OUTDIR/assembly_spades -name "*_contigs.fa" | wc -l)
                echo "     - Generated $SPADES_CONTIGS contig files"
            fi
            
            if [ -d "$OUTDIR/prodigal_megahit" ]; then
                echo "  ✅ Prodigal MEGAHIT gene predictions: $OUTDIR/prodigal_megahit/"
                MEGAHIT_GENES=$(find $OUTDIR/prodigal_megahit -name "*.faa" | wc -l)
                echo "     - Generated $MEGAHIT_GENES protein sequence files"
            fi
            
            if [ -d "$OUTDIR/prodigal_spades" ]; then
                echo "  ✅ Prodigal SPAdes gene predictions: $OUTDIR/prodigal_spades/"
                SPADES_GENES=$(find $OUTDIR/prodigal_spades -name "*.faa" | wc -l)
                echo "     - Generated $SPADES_GENES protein sequence files"
            fi
            
            if [ -d "$OUTDIR/diamond_megahit" ]; then
                echo "  ✅ Diamond MEGAHIT results: $OUTDIR/diamond_megahit/"
            fi
            
            if [ -d "$OUTDIR/diamond_spades" ]; then
                echo "  ✅ Diamond SPAdes results: $OUTDIR/diamond_spades/"
            fi
            
            if [ -d "$OUTDIR/merged_reports" ]; then
                echo "  ✅ Comprehensive analysis reports: $OUTDIR/merged_reports/"
                MERGED_REPORTS=$(find $OUTDIR/merged_reports -name "*.txt" | wc -l)
                echo "     - Generated $MERGED_REPORTS comprehensive reports"
            fi
            
            if [ -d "$OUTDIR/abundance_megahit" ]; then
                echo "  ✅ Viral abundance (MEGAHIT): $OUTDIR/abundance_megahit/"
                MEGAHIT_ABUNDANCE=$(find $OUTDIR/abundance_megahit -name "*.txt" -o -name "*.csv" | wc -l)
                echo "     - Generated $MEGAHIT_ABUNDANCE abundance reports (RPM and RPKM)"
            fi
            
            if [ -d "$OUTDIR/abundance_spades" ]; then
                echo "  ✅ Viral abundance (SPAdes): $OUTDIR/abundance_spades/"
                SPADES_ABUNDANCE=$(find $OUTDIR/abundance_spades -name "*.txt" -o -name "*.csv" | wc -l)
                echo "     - Generated $SPADES_ABUNDANCE abundance reports (RPM and RPKM)"
            fi
        fi
        
        echo ""
        echo "📋 Summary of generated files:"
        echo "  fastp reports:"
        find $OUTDIR/fastp -name "*.html" -o -name "*.json" 2>/dev/null | head -10
        echo ""
        echo "  Assembly contigs:"
        find $OUTDIR/assembly_* -name "*_contigs.fa" 2>/dev/null | head -10
        echo ""
        echo "  Gene prediction results:"
        find $OUTDIR/prodigal_* -name "*.faa" 2>/dev/null | head -10
        echo ""
        echo "  Diamond classification results:"
        find $OUTDIR/diamond_* -name "*.txt" 2>/dev/null | head -10
        echo ""
        echo "  Merged reports:"
        find $OUTDIR/merged_reports -name "*.txt" -o -name "*.csv" 2>/dev/null | head -10
        echo ""
        echo "  Abundance reports:"
        find $OUTDIR/abundance_* -name "*.txt" -o -name "*.csv" 2>/dev/null | head -10
        echo ""
        echo "Total files: $(find $OUTDIR -type f | wc -l)"
        
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
