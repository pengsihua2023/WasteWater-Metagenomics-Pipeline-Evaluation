## Metagenome Viral Classification Pipeline (Short-Read Mode)

## Overview

This pipeline performs **taxonomic profiling of viral sequences** from short-read metagenomic data. It takes paired-end FASTQ files as input and produces abundance tables and summary reports for viral taxa using GOTTCHA2.

[![Python](https://img.shields.io/badge/python-%3D%203.10-blue.svg)](https://www.python.org/)

---

## Workflow Diagram

```
Raw/Preprocessed FASTQ (R1, R2)
           │
           ▼
   ┌───────────────────┐
   │ 1. Environment    │  Miniforge3, conda env (gottcha2_env)
   │    setup           │
   └─────────┬─────────┘
             │
             ▼
   ┌───────────────────┐
   │ 2. Taxonomic      │  GOTTCHA2 profile
   │    profiling       │  Database: gottcha_db.species.fna
   └─────────┬─────────┘
             │
             ▼
   Result directory (result_shortread/)
   • Summary TSV
   • Abundance tables
   • Per-sample reports
```

---

## Software Used

### Execution environment
| Software / Component | Role |
|----------------------|------|
| **SLURM** | Job scheduling and resource allocation (partition, CPUs, memory, wall time) |
| **Miniforge3** (v24.11.3-0) | Conda distribution (module-loaded) |
| **Conda** | Environment manager; runs GOTTCHA2 from `gottcha2_env` |

### Data preprocessing
The current script **does not** run preprocessing inside the pipeline. Input is expected to be:
- **Paired-end FASTQ** (e.g. `*_1.fastq` / `*_2.fastq` or `*_R1.fastq.gz` / `*_R2.fastq.gz`).

If preprocessing is done **before** this pipeline, typical tools include:
- **fastp** or **Trimmomatic**: adapter trimming, quality filtering
- **BBduk** (BBTools): contaminant removal, quality trimming
- **PRINSEQ**: quality and length filtering

So: **preprocessing is optional and external**; the pipeline starts from FASTQ and runs taxonomic profiling.

### Abundance / taxonomic profiling
| Software | Version / note | Role |
|----------|----------------|------|
| **GOTTCHA2** | Run via `gottcha2 profile` | Reference-based taxonomic profiling; assigns reads to taxa and produces **abundance** (read counts and fractions) for viruses and other taxa. |

- **Database**: GOTTCHA2 database built from reference sequences (e.g. `gottcha_db.species.fna` with `.mmi` index under `DB_DIR`).
- **Algorithm**: Minimap2-based alignment of reads to the database; assignment and aggregation to produce per-taxon abundances.

No other abundance or viral-detection tools are invoked in the current script (e.g. no Kraken2, MetaPhlAn, or MetaTaxProfiler in this single script).

---

## Pipeline Steps (as in the script)

1. **Environment setup**
   - Change to `$SLURM_SUBMIT_DIR`.
   - Load Miniforge3 module.
   - Activate conda environment `gottcha2_env`.

2. **GOTTCHA2 taxonomic profiling**
   - **Input**: Paired-end FASTQ (e.g. `SRR35987572_1.fastq`, `SRR35987572_2.fastq`).
   - **Database**: `DB_PREFIX` = `$DB_DIR/gottcha_db.species.fna` (with `.mmi` in `DB_DIR`).
   - **Parameters**: `-t 32` (threads), output directory `-o result_shortread`.
   - **Command**: `gottcha2 profile -i <R1> <R2> -d <DB_PREFIX> -t 32 -o result_shortread`.

3. **Logging**
   - Start/end time and duration are printed to the SLURM output file (via a shell `trap` on EXIT).

---

## Result Files

All profiling results are written under the **output directory** specified by `-o` (default: **`result_shortread/`**).

### Main result files (typical GOTTCHA2 output)

| File / directory | Description |
|------------------|-------------|
| **`result_shortread/`** | Top-level result directory for the run. |
| **`result_shortread/summary.tsv`** | **Main summary table**: taxon IDs (e.g. TaxID), names, read counts, and **abundance fractions** (relative abundance) per taxon. Primary file for downstream analysis and plotting. |
| **`result_shortread/*.tsv`** (other TSV) | Additional tables (e.g. per-level or filtered summaries) depending on GOTTCHA2 version and options. |
| **`result_shortread/*.txt`** (if any) | Logs or text summaries produced by GOTTCHA2. |

### Job logs (SLURM)

| File | Description |
|------|-------------|
| **`Viral_Classification_ShortRead_<jobid>.out`** | Standard output: step messages, start/end time, duration. |
| **`Viral_Classification_ShortRead_<jobid>.err`** | Standard error: warnings and errors from GOTTCHA2 or the environment. |

---

## Summary

- **Input**: Paired-end FASTQ (optionally preprocessed elsewhere).
- **Preprocessing**: Not performed in this script; done externally if at all.
- **Abundance/taxonomic profiling**: **GOTTCHA2** only, using a species-level GOTTCHA2 database.
- **Final result files**: **`result_shortread/`** (especially **`result_shortread/summary.tsv`**) plus SLURM `.out` and `.err` logs.

For virus-focused analysis, downstream steps typically use the viral taxa and their abundance fractions from **`summary.tsv`** (e.g. filtering to viruses, comparing with other tools, or building plots).
