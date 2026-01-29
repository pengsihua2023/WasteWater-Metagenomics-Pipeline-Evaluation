# Sourmash Metagenome Analysis Pipeline

This repository contains SLURM scripts for running sourmash sketch and gather analysis on metagenomic contig files from different assembly methods.

[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue.svg)](https://www.python.org/)
---

## Overview

Sourmash is a tool for quickly searching, comparing, and analyzing genomic and metagenomic sequence data using MinHash sketches. This pipeline includes scripts for:

1. **Sketch**: Creating signature files from contig assemblies
2. **Gather**: Searching signatures against reference databases (e.g., NCBI viruses)

## Files

### Sketch Scripts (Create Signatures)

- `sourmash_sketch_longread.slurm` - Create signatures from long-read assembly (MetaFlye)
- `sourmash_sketch_spades.slurm` - Create signatures from short-read SPAdes assembly
- `sourmash_sketch_MEGAHIT.slurm` - Create signatures from short-read MEGAHIT assembly

### Gather Scripts (Database Search)

- `sourmash_gather_longread.slurm` - Search long-read signatures against database
- `sourmash_gather_spades.slurm` - Search SPAdes signatures against database
- `sourmash_gather_MEGAHIT.slurm` - Search MEGAHIT signatures against database

## Input Files

### Long-read Assembly
- **Directory**: `/scratch/sp96859/Meta-genome-data-analysis/Apptainer/sourmash/no-nf/data/long_contig`
- **File**: `llnl_66d1047e_metaflye_contigs.fa`
- **Output signature**: `llnl_66d1047e_metaflye_contigs.sig`

### SPAdes Assembly
- **Directory**: `/scratch/sp96859/Meta-genome-data-analysis/Apptainer/sourmash/no-nf/data/short_contig_spades`
- **File**: `llnl_66ce4dde_spades_contigs.fa`
- **Output signature**: `llnl_66ce4dde_spades_contigs.sig`

### MEGAHIT Assembly
- **Directory**: `/scratch/sp96859/Meta-genome-data-analysis/Apptainer/sourmash/no-nf/data/short_contig_MEGAHIT`
- **File**: `llnl_66ce4dde_megahit_contigs.fa`
- **Output signature**: `llnl_66ce4dde_megahit_contigs.sig`

## Database

- **Database directory**: `/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/sourmash_db`
- **Database file**: `ncbi-viruses-2025.01.k31.zip`
- **Database type**: NCBI Viruses (k=31)

## Usage

### Step 1: Create Signatures

First, create signature files from your contig assemblies:

```bash
# For long-read assembly
sbatch sourmash_sketch_longread.slurm

# For SPAdes assembly
sbatch sourmash_sketch_spades.slurm

# For MEGAHIT assembly
sbatch sourmash_sketch_MEGAHIT.slurm
```

### Step 2: Search Database

After the sketch job completes, search against the database:

```bash
# For long-read signatures
sbatch sourmash_gather_longread.slurm

# For SPAdes signatures
sbatch sourmash_gather_spades.slurm

# For MEGAHIT signatures
sbatch sourmash_gather_MEGAHIT.slurm
```

### Monitor Jobs

```bash
# Check job status
squeue -u $USER

# View output log
tail -f sourmash_sketch_longread_<JOB_ID>.out

# View error log
tail -f sourmash_sketch_longread_<JOB_ID>.err
```

## Parameters

### Sketch Parameters

- **k-mer sizes**: k=21, k=31, k=51
- **Scaled**: 1000 (1 in 1000 k-mers retained)
- **Abundance tracking**: Enabled (`abund`)
- **Molecule type**: DNA

### Gather Parameters

- **k-mer size**: k=31 (matches database)
- **Threshold**: 300 bp minimum overlap
- **Database**: NCBI Viruses 2025.01

## Output Files

### Signature Files (.sig)

- Contains k-mer sketches for k=21, 31, and 51
- Includes abundance information
- Format: JSON-based signature file

### Gather Results (.csv)

CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `intersect_bp` | Overlapping base pairs |
| `f_orig_query` | Fraction in original query |
| `f_match` | Fraction in matched sequence |
| `average_abund` | Average k-mer abundance |
| `median_abund` | Median k-mer abundance |
| `std_abund` | Standard deviation of abundance |
| `name` | Matched sequence name |
| `gather_result_rank` | Match ranking (0 = best match) |
| `remaining_bp` | Remaining unmatched base pairs |
| `query_abundance` | Whether query has abundance info |
| `query_containment_ani` | Query containment ANI |
| `match_containment_ani` | Match containment ANI |
| `average_containment_ani` | Average containment ANI |

## Resource Requirements

All scripts use the following SLURM resources:

- **Partition**: `bahl_p`
- **CPUs**: 32 cores
- **Memory**: 256 GB
- **Time limit**: 72 hours

## Output File Names

### Sketch Outputs
- Long-read: `llnl_66d1047e_metaflye_contigs.sig`
- SPAdes: `llnl_66ce4dde_spades_contigs.sig`
- MEGAHIT: `llnl_66ce4dde_megahit_contigs.sig`

### Gather Outputs
- Long-read: `gather_results_ncbi_viruses.csv`
- SPAdes: `gather_results_spades_ncbi_viruses.csv`
- MEGAHIT: `gather_results_MEGAHIT_ncbi_viruses.csv`

## Interpreting Results

### Abundance Information

The results include abundance information because signatures were created with the `abund` parameter:

- **average_abund**: Average k-mer occurrence count
- **median_abund**: Median k-mer occurrence count
- **std_abund**: Standard deviation of k-mer abundances

### Match Quality

- **gather_result_rank**: Lower rank = better match (0 is best)
- **intersect_bp**: Larger overlap = more confident match
- **containment_ani**: Higher ANI = more similar sequences
- **remaining_bp**: Shows how much of the query remains unmatched

### Example Interpretation

A match with:
- `gather_result_rank: 0`
- `intersect_bp: 17000`
- `average_containment_ani: 0.76`
- `average_abund: 1.06`

Indicates:
- Best match found
- 17,000 bp overlap
- 76% average containment ANI
- Low abundance (k-mers appear ~1 time on average)

## Notes

1. **Threshold**: The gather threshold is set to 300 bp. Lower thresholds detect more matches but may include false positives.

2. **k-mer Matching**: Only k=31 is used for gather (matching the database). The sketch files contain k=21, 31, and 51 for flexibility.

3. **Abundance Tracking**: Enabled by default. This allows analysis of k-mer abundance in addition to presence/absence.

4. **File Locations**: Ensure input files exist at the specified paths before submitting jobs.

5. **Job Dependencies**: Gather jobs require sketch jobs to complete first.

## Contact

For questions or issues, please check the sourmash documentation:
- https://sourmash.readthedocs.io/


