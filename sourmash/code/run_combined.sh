#!/bin/bash

#SBATCH --job-name=Sourmash
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=512G
#SBATCH --time=5-04:00:00
#SBATCH --output=Combined_Metagenome_%j.out
#SBATCH --error=Combined_Metagenome_%j.err

# Switch to the submission directory
cd "$SLURM_SUBMIT_DIR" || exit 1

# Print job info
echo "=========================================="
echo "  SLURM job info - short + long reads"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node(s): $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: ${SLURM_MEM_PER_NODE:-N/A} MB"
echo "Workdir: $SLURM_SUBMIT_DIR"
echo "Start time: $(date)"
echo "=========================================="
echo ""

# Initialize conda
echo "Initializing conda..."

# Method 1: Try loading a conda module (module systems)
if command -v module &> /dev/null; then
    module load Miniforge3/24.11.3-0 2>/dev/null || true
fi

# Method 2: Try common conda locations
CONDA_INITIALIZED=false

# Try sourcing conda.sh from common locations
for conda_path in \
    "/apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "/opt/conda/etc/profile.d/conda.sh"; do
    if [ -f "$conda_path" ]; then
        echo "Loading conda from: $conda_path"
        # shellcheck disable=SC1090
        source "$conda_path"
        CONDA_INITIALIZED=true
        break
    fi
done

# Method 3: If conda is already on PATH, use the shell hook
if [ "$CONDA_INITIALIZED" = false ]; then
    if command -v conda &> /dev/null; then
        eval "$(conda shell.bash hook)" 2>/dev/null && CONDA_INITIALIZED=true
    fi
fi

# Method 4: Derive conda base via `conda info --base`
if [ "$CONDA_INITIALIZED" = false ]; then
    if command -v conda &> /dev/null; then
        CONDA_BASE=$(conda info --base 2>/dev/null)
        if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
            # shellcheck disable=SC1090
            source "$CONDA_BASE/etc/profile.d/conda.sh"
            CONDA_INITIALIZED=true
        fi
    fi
fi

if [ "$CONDA_INITIALIZED" = false ]; then
    echo "[ERROR] Failed to initialize conda."
    echo "Please check that conda is installed or that the correct module is loaded."
    exit 1
fi

echo "[OK] Conda initialized."

# Activate the nextflow_env environment
echo "Activating conda env: nextflow_env"
# Activate using an absolute env path if present
CONDA_ENVS_PATH="/scratch/sp96859/conda/envs"
if [ -d "$CONDA_ENVS_PATH/nextflow_env" ]; then
    conda activate "$CONDA_ENVS_PATH/nextflow_env" || conda activate nextflow_env
else
    conda activate nextflow_env
fi

if [ "$CONDA_DEFAULT_ENV" != "nextflow_env" ] && [ "$CONDA_DEFAULT_ENV" != "$CONDA_ENVS_PATH/nextflow_env" ]; then
    echo "[ERROR] Failed to activate conda env: nextflow_env"
    exit 1
fi

echo "[OK] Active conda env: $CONDA_DEFAULT_ENV"
echo ""

# Note: The Nextflow driver itself can use small resources.
# The heavy compute tasks are submitted as separate Slurm jobs via the slurm executor (resources are in nextflow.config).

# Put Nextflow home on scratch to avoid many small files under $HOME
NXF_HOME_DIR="/scratch/${USER}/nextflow_home"
mkdir -p "$NXF_HOME_DIR"
export NXF_HOME="$NXF_HOME_DIR"
echo "[OK] NXF_HOME: $NXF_HOME"

# Add Nextflow to PATH (if needed)
NEXTFLOW_DIR="/scratch/sp96859/bin"
if [ -d "$NEXTFLOW_DIR" ] && [[ ":$PATH:" != *":$NEXTFLOW_DIR:"* ]]; then
    export PATH="$NEXTFLOW_DIR:$PATH"
    echo "[OK] Added to PATH: $NEXTFLOW_DIR"
fi

# Set Apptainer cache directory
APPTAINER_CACHE="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/cache"
export APPTAINER_CACHEDIR="$APPTAINER_CACHE"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHE"
echo "[OK] Apptainer cache: $APPTAINER_CACHE"

# Validate Nextflow
echo "Validating Nextflow..."
NEXTFLOW_BIN="$(command -v nextflow 2>/dev/null || true)"
if [ -z "$NEXTFLOW_BIN" ]; then
    # Try an explicit path
    if [ -f "$NEXTFLOW_DIR/nextflow" ]; then
        echo "[WARN] Nextflow not on PATH, but found: $NEXTFLOW_DIR/nextflow"
        NEXTFLOW_BIN="$NEXTFLOW_DIR/nextflow"
        export PATH="$NEXTFLOW_DIR:$PATH"
        echo "[OK] Nextflow path configured."
    else
        echo "[ERROR] Nextflow not found."
        echo "Please check that Nextflow is installed under: $NEXTFLOW_DIR"
        exit 1
    fi
fi

echo "Nextflow version:"
"$NEXTFLOW_BIN" -version
echo ""

# This script runs "short + long reads" combined analysis.
# If you ONLY have short reads, use `run_short.sh` (it will not trigger the metaflye branch).

# Validate samplesheets
SAMPLESHEET_SHORT="samplesheet_short.csv"
SAMPLESHEET_LONG="samplesheet_long.csv"

if [ ! -f "$SAMPLESHEET_SHORT" ]; then
    echo "[ERROR] samplesheet_short.csv not found."
    echo "Please make sure samplesheet_short.csv is in the current directory."
    exit 1
fi

echo "[OK] Found short-read samplesheet: $SAMPLESHEET_SHORT"

# Long reads are optional. If you don't have long reads, do NOT pass --samplesheet_long (avoids metaflye branch).
ENABLE_LONG_READS=${ENABLE_LONG_READS:-true}
if [ "$ENABLE_LONG_READS" = "true" ]; then
    if [ ! -f "$SAMPLESHEET_LONG" ]; then
        echo "[ERROR] ENABLE_LONG_READS=true but samplesheet_long.csv not found."
        echo "Either provide samplesheet_long.csv, or run:"
        echo "  ENABLE_LONG_READS=false sbatch run_combined.sh"
        exit 1
    fi
    echo "[OK] Found long-read samplesheet: $SAMPLESHEET_LONG"
else
    echo "[WARN] Long reads disabled: samplesheet_long.csv will not be used and metaflye will not run."
fi
echo ""

# Nextflow parameters (short + long reads)
OUTDIR="results_combined"
WORKDIR="work_combined"
SLURM_PARTITION="bahl_p"

# MEGAHIT resources (critical for Slurm/cgroups; prevents SIGKILL/OOM)
MEGAHIT_TASK_MEMORY="240 GB"
MEGAHIT_TASK_CPUS="16"
MEGAHIT_TASK_TIME="144h"

# METAFLYE resources (also set a generous walltime)
METAFLYE_TASK_MEMORY="128 GB"
METAFLYE_TASK_CPUS="16"
METAFLYE_TASK_TIME="144h"

# Run the workflow
echo "=========================================="
echo "  Running Nextflow workflow - short + long reads"
echo "=========================================="
echo ""
echo "Parameters:"
echo "  samplesheet_short: $SAMPLESHEET_SHORT"
echo "  samplesheet_long: $SAMPLESHEET_LONG"
echo "  outdir: $OUTDIR"
echo "  workdir: $WORKDIR"
if [ "$ENABLE_LONG_READS" = "true" ]; then
    echo "  data type: short + long reads"
else
    echo "  data type: short reads only"
fi
echo "  skip short-read assembly: false (will assemble)"
echo "  short-read assembler: megahit"
if [ "$ENABLE_LONG_READS" = "true" ]; then
    echo "  skip long-read assembly: false (will assemble)"
    echo "  long-read assembler: metaflye"
    echo "  use local Flye: true"
else
    echo "  skip long-read assembly: true (long-read branch will not run)"
fi
echo "  Slurm partition: $SLURM_PARTITION"
echo "  MEGAHIT resources: memory=$MEGAHIT_TASK_MEMORY cpus=$MEGAHIT_TASK_CPUS time=$MEGAHIT_TASK_TIME"
if [ "$ENABLE_LONG_READS" = "true" ]; then
    echo "  METAFLYE resources: memory=$METAFLYE_TASK_MEMORY cpus=$METAFLYE_TASK_CPUS time=$METAFLYE_TASK_TIME"
fi
echo ""

NF_ARGS=()
NF_ARGS+=(--samplesheet_short "$SAMPLESHEET_SHORT")
if [ "$ENABLE_LONG_READS" = "true" ]; then
    NF_ARGS+=(--samplesheet_long "$SAMPLESHEET_LONG")
    NF_ARGS+=(--skip_long_assembly false)
    NF_ARGS+=(--long_assembler metaflye)
    NF_ARGS+=(--use_local_flye true)
    NF_ARGS+=(--metaflye_task_memory "$METAFLYE_TASK_MEMORY")
    NF_ARGS+=(--metaflye_task_cpus "$METAFLYE_TASK_CPUS")
    NF_ARGS+=(--metaflye_task_time "$METAFLYE_TASK_TIME")
else
    # Even if a long-read samplesheet exists, force-skip the long-read branch
    NF_ARGS+=(--skip_long_assembly true)
fi

"$NEXTFLOW_BIN" run sourmash_workflow.nf \
    "${NF_ARGS[@]}" \
    --skip_assembly false \
    --assembler megahit \
    --use_local_sourmash true \
    --use_local_megahit true \
    --use_local_spades true \
    --outdir "$OUTDIR" \
    --work_dir "$WORKDIR" \
    --slurm_partition "$SLURM_PARTITION" \
    --megahit_task_memory "$MEGAHIT_TASK_MEMORY" \
    --megahit_task_cpus "$MEGAHIT_TASK_CPUS" \
    --megahit_task_time "$MEGAHIT_TASK_TIME" \
    -with-report "${OUTDIR}/nextflow_report.html" \
    -with-trace "${OUTDIR}/nextflow_trace.txt" \
    -with-timeline "${OUTDIR}/nextflow_timeline.html"

# Exit status
EXIT_CODE=$?

echo ""
echo "=========================================="
echo "  Job finished"
echo "=========================================="
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "[OK] Workflow completed successfully."
    echo ""
    echo "Results directory: $OUTDIR"
    echo "To inspect results:"
    echo "  ls -lh $OUTDIR"
else
    echo "[ERROR] Workflow failed."
    echo "Please check the log file:"
    echo "  tail -f ${OUTDIR}/.nextflow.log"
    exit $EXIT_CODE
fi

