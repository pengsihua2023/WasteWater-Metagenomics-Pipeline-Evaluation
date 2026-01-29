# GOTTCHA2 Metagenome Viral Classification

This repository contains SLURM job submission scripts for performing taxonomic profiling and viral classification using GOTTCHA2 on metagenomic sequencing data.

[![Python](https://img.shields.io/badge/python-%3D%203.10-blue.svg)](https://www.python.org/)

---

## Overview

GOTTCHA2 (Genomic Origin Through Taxonomic CHAllenge) is a taxonomic profiling tool for metagenomic data. This workflow supports both short-read (Illumina) and long-read (Oxford Nanopore, PacBio) sequencing data.

## Files

- `run_metagenome_assembly_classification_shortread.sh` - Script for short-read paired-end data analysis
- `run_metagenome_assembly_classification_longread.sh` - Script for long-read single-end data analysis

## Requirements

### Software Dependencies
- SLURM workload manager
- Miniforge3/24.11.3-0 (or compatible version)
- GOTTCHA2 (installed in conda environment)

### Conda Environment
- Environment name: `gottcha2_env`
- Location: `/home/sp96859/.conda/envs/gottcha2_env`

### Database
- GOTTCHA2 database location: `/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/gottcha2_db`
- Database file: `gottcha_db.species.fna.mmi`

## Usage

### Short-Read Analysis (Illumina)

**Input Data:**
- Paired-end FASTQ files: `llnl_66ce4dde_R1.fastq.gz` and `llnl_66ce4dde_R2.fastq.gz`

**Submit Job:**
```bash
sbatch run_metagenome_assembly_classification_shortread.sh
```

**Output Directory:** `result_shortread/`

### Long-Read Analysis (Oxford Nanopore / PacBio)

**Input Data:**
- Single FASTQ file: `llnl_66d1047e.fastq.gz`

**Submit Job:**
```bash
sbatch run_metagenome_assembly_classification_longread.sh
```

**Output Directory:** `result_longread/`

**Platform-Specific Settings:**
- For Oxford Nanopore data: use `--nanopore` parameter (default)
- For PacBio data: change to `--pacbio` parameter in the script

## Resource Requirements

Both scripts are configured with the following SLURM resources:

| Resource | Value |
|----------|-------|
| Partition | bahl_p |
| CPUs | 32 |
| Memory | 256 GB |
| Time Limit | 72 hours |

## Output Files

Each analysis produces the following output files:

- `*.tsv` - Classification results table
- `*.summary.tsv` - Classification summary
- `*.full.tsv` - Full classification information

## Customization

### Modifying Input Files

**For short-read analysis:**
Edit line 39 in `run_metagenome_assembly_classification_shortread.sh`:
```bash
  -i your_R1.fastq.gz your_R2.fastq.gz \
```

**For long-read analysis:**
Edit line 45 in `run_metagenome_assembly_classification_longread.sh`:
```bash
  -i your_longread_data.fastq.gz \
```

### Changing Thread Count

Modify the `-t` parameter (default: 32) to match your available CPU cores:
```bash
  -t 32 \
```

### Adjusting Memory and Time Limits

Edit the SLURM directives at the top of each script:
```bash
#SBATCH --mem=256G
#SBATCH --time=72:00:00
```



## Supported File Formats

- `.fastq.gz` (compressed FASTQ, recommended)
- `.fq.gz` (compressed FASTQ)
- `.fastq` (uncompressed FASTQ)
- `.fq` (uncompressed FASTQ)


## References

- GOTTCHA2: [https://github.com/poeli/GOTTCHA2](https://github.com/poeli/GOTTCHA2)
- Publication: Freitas et al. (2018) "Accurate read-based metagenome characterization using a hierarchical suite of unique signatures"
