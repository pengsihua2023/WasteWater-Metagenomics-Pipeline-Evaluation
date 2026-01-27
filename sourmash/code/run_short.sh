#!/bin/bash

#SBATCH --job-name=Sourmash_short
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=512G
#SBATCH --time=5-04:00:00
#SBATCH --output=Short_Metagenome_%j.out
#SBATCH --error=Short_Metagenome_%j.err

# Short-read only run (will not trigger the metaflye/long-read branch)

cd "$SLURM_SUBMIT_DIR" || exit 1

echo "=========================================="
echo "  SLURM job info - short reads only"
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

# Initialize conda (kept consistent with run_combined.sh)
echo "Initializing conda..."
if command -v module &> /dev/null; then
    module load Miniforge3/24.11.3-0 2>/dev/null || true
fi

CONDA_INITIALIZED=false
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

if [ "$CONDA_INITIALIZED" = false ]; then
    if command -v conda &> /dev/null; then
        eval "$(conda shell.bash hook)" 2>/dev/null && CONDA_INITIALIZED=true
    fi
fi

if [ "$CONDA_INITIALIZED" = false ]; then
    echo "[ERROR] Failed to initialize conda."
    exit 1
fi
echo "[OK] Conda initialized."

echo "Activating conda env: nextflow_env"
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

NXF_HOME_DIR="/scratch/${USER}/nextflow_home"
mkdir -p "$NXF_HOME_DIR"
export NXF_HOME="$NXF_HOME_DIR"
echo "[OK] NXF_HOME: $NXF_HOME"

NEXTFLOW_DIR="/scratch/sp96859/bin"
if [ -d "$NEXTFLOW_DIR" ] && [[ ":$PATH:" != *":$NEXTFLOW_DIR:"* ]]; then
    export PATH="$NEXTFLOW_DIR:$PATH"
fi

APPTAINER_CACHE="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/cache"
export APPTAINER_CACHEDIR="$APPTAINER_CACHE"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHE"

NEXTFLOW_BIN="$(command -v nextflow 2>/dev/null || true)"
if [ -z "$NEXTFLOW_BIN" ]; then
    if [ -f "$NEXTFLOW_DIR/nextflow" ]; then
        NEXTFLOW_BIN="$NEXTFLOW_DIR/nextflow"
        export PATH="$NEXTFLOW_DIR:$PATH"
    else
        echo "[ERROR] Nextflow not found."
        exit 1
    fi
fi

SAMPLESHEET_SHORT="samplesheet_short.csv"
if [ ! -f "$SAMPLESHEET_SHORT" ]; then
    echo "[ERROR] samplesheet_short.csv not found."
    exit 1
fi
echo "[OK] Found short-read samplesheet: $SAMPLESHEET_SHORT"
echo ""

OUTDIR="results_short"
WORKDIR="work_short"
SLURM_PARTITION="bahl_p"

MEGAHIT_TASK_MEMORY="240 GB"
MEGAHIT_TASK_CPUS="16"
MEGAHIT_TASK_TIME="144h"

echo "=========================================="
echo "  Running Nextflow workflow - short reads only"
echo "=========================================="
echo ""

"$NEXTFLOW_BIN" run sourmash_workflow.nf \
    --samplesheet_short "$SAMPLESHEET_SHORT" \
    --skip_assembly false \
    --assembler megahit \
    --skip_long_assembly true \
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

EXIT_CODE=$?
echo ""
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
exit $EXIT_CODE

