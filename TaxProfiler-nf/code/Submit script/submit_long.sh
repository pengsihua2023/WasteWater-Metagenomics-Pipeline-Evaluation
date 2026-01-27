#!/bin/bash
#SBATCH --job-name=Viral_Long
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Viral_Long_%j.out
#SBATCH --error=Viral_Long_%j.err

# === 1. Activate nextflow_env (critical! contains Java 17 + apptainer) ===
source activate nextflow_env

# === 2. Load Nextflow (if not in environment) ===
module load nextflow 2>/dev/null || true

# === 3. Set Apptainer cache ===
export APPTAINER_CACHEDIR=$HOME/.apptainer_cache
export NXF_SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR

# === 4. Verify environment (optional, for debugging) ===
echo "=== Environment Verification ==="
echo "Conda env: $CONDA_DEFAULT_ENV"
echo "Java:"
java -version
echo "Apptainer:"
apptainer --version
echo "Nextflow:"
nextflow -v

# === 5. Run nf-core/taxprofiler workflow (long-read data) ===
echo ""
echo "========================================"
echo "▶️  Starting TaxProfiler Analysis"
echo "========================================"

nextflow run nf-core/taxprofiler \
  -r 1.2.0 \
  -profile apptainer \
  -c nextflow_long.config \
  -resume

# Check if workflow completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "✅ TaxProfiler Analysis Completed"
    echo "========================================"
    
    # === 6. Automatic viral abundance calculation (RPM & RPKM) ===
    echo ""
    echo "========================================"
    echo "🧮 Auto-calculating Viral Abundance (RPM & RPKM)"
    echo "========================================"
    
    # Get output directory (read from config file)
    OUTDIR=$(grep "outdir" nextflow_long.config | head -1 | sed "s/.*= *['\"]//;s/['\"].*//")
    
    if [ -z "$OUTDIR" ]; then
        OUTDIR="results_viral_long"
    fi
    
    echo "Output directory: $OUTDIR"
    
    # Run abundance calculation script (long-read specific)
    # Note: Long-read data doesn't use Bracken, calculates abundance directly from Kraken2 report
    if [ -f "batch_calculate_abundance_longread_en.sh" ]; then
        bash batch_calculate_abundance_longread_en.sh "$OUTDIR"
        
        if [ $? -eq 0 ]; then
            echo ""
            echo "========================================"
            echo "🎉 All Tasks Completed!"
            echo "========================================"
            echo "📊 Classification results: $OUTDIR/kraken2/"
            echo "📈 Abundance results: $OUTDIR/abundance/"
            echo "📋 Comprehensive report: $OUTDIR/multiqc/multiqc_report.html"
            echo "========================================"
        else
            echo "⚠️  Warning: Abundance calculation failed, but main analysis completed"
        fi
    else
        echo "⚠️  Warning: batch_calculate_abundance_longread_en.sh script not found"
        echo "    You can run manually: bash batch_calculate_abundance_longread_en.sh $OUTDIR"
    fi
else
    echo ""
    echo "========================================"
    echo "❌ TaxProfiler Analysis Failed"
    echo "========================================"
    echo "Please check log files for detailed information"
    exit 1
fi
