#!/bin/bash
#SBATCH --job-name=CLARK_Build_Virus_DB
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=CLARK_Build_Virus_DB_%j.out
#SBATCH --error=CLARK_Build_Virus_DB_%j.err

cd "$SLURM_SUBMIT_DIR" || exit 1

# Set database directory and CLARK software directory
# CLARK_DIR: CLARK software installation directory, should contain set_targets.sh, classify_metagenome.sh, etc.
export CLARK_DB_DIR=${CLARK_DB_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/CLARK_db"}
export CLARK_DIR=${CLARK_DIR:-"/scratch/sp96859/Meta-genome-data-analysis/Apptainer/CLARK/CLARK"}

echo "=========================================="
echo "Start building CLARK virus database"
echo "Time: $(date)"
echo "Database directory: $CLARK_DB_DIR"
echo "CLARK directory: $CLARK_DIR"
echo "=========================================="

# Check if CLARK is installed
if [ ! -f "$CLARK_DIR/set_targets.sh" ]; then
    echo "Error: Cannot find set_targets.sh script!"
    echo "Please set CLARK_DIR environment variable or ensure running this script in CLARK directory"
    exit 1
fi

# Create database directory
mkdir -p "$CLARK_DB_DIR"

# Enter CLARK directory
cd "$CLARK_DIR" || exit 1

# Build virus database
echo "Start downloading and building virus database..."
./set_targets.sh "$CLARK_DB_DIR" viruses

if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "Virus database built successfully!"
    echo "Completion time: $(date)"
    echo "Database location: $CLARK_DB_DIR/Viruses"
    echo "=========================================="
else
    echo "=========================================="
    echo "Error: Database build failed!"
    echo "Failure time: $(date)"
    echo "=========================================="
    exit 1
fi
