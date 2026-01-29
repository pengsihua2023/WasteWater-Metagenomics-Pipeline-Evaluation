# MetaTaxProfiler Pipeline Workflow

This document describes the end-to-end workflow of MetaTaxProfiler: software used at each stage (including data preprocessing and abundance analysis) and the final result files produced.

[![Apptainer](https://img.shields.io/badge/Apptainer-%E2%89%A51.3-42b983)](https://apptainer.org/)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue.svg)](https://www.python.org/)

---

## 1. Overview

MetaTaxProfiler is built on **nf-core/taxprofiler** (v1.2.0) and adds automated viral abundance calculation. The pipeline supports two data types:

| Data type   | Entry script      | Output directory       |
|------------|-------------------|------------------------|
| Short-read | `submit_short.sh` | `results_viral_short/` |
| Long-read  | `submit_long.sh`  | `results_viral_long/`  |

Execution is containerized with **Apptainer** and orchestrated by **Nextflow**.

---

## 2. Workflow Diagram

### 2.1 Short-read (Illumina)

```
Raw FASTQ (e.g. *_1.fastq.gz, *_2.fastq.gz)
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  DATA PREPROCESSING                                             │
│  • FastQC          – quality metrics per sample                 │
│  • fastp            – adapter trimming, QC, paired/single       │
│  • Run merging      – merge multiple runs per sample (optional) │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  TAXONOMIC CLASSIFICATION                                       │
│  • Kraken2         – k-mer–based taxonomic assignment           │
│  • Bracken         – species-level read count re-estimation     │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌───────────────────────────────────────────────────────────────────────┐
│  ABUNDANCE ANALYSIS (MetaTaxProfiler)                                 │
│  • batch_calculate_abundance_en.sh  – batch driver                    │
│  • calculate_abundance_en.py        – Bracken → RPM (or fallback      │
│  • calculate_abundance_longread_en.py – Kraken2 → RPM if no Bracken)  │
└───────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  REPORTING & STANDARDISATION                                    │
│  • Taxpasta        – standardised profile tables (e.g. MOTU)    │
│  • MultiQC         – aggregated QC and pipeline report          │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
   Final results (see Section 4)
```

### 2.2 Long-read (Nanopore / PacBio)

```
Raw FASTQ (long reads, e.g. *.fastq.gz)
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  DATA PREPROCESSING                                             │
│  • NanoPlot / NanoStat  – long-read QC metrics                  │
│  • Porechop (or Porechop_ABI) – adapter removal                 │
│  • Run merging          – merge runs per sample (optional)      │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  TAXONOMIC CLASSIFICATION                                       │
│  • Kraken2         – taxonomic assignment (Bracken not used)    │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  ABUNDANCE ANALYSIS (MetaTaxProfiler)                               │
│  • batch_calculate_abundance_longread_en.sh – batch driver          │
│  • calculate_abundance_longread_en.py       – Kraken2 report → RPM  │
└─────────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  REPORTING & STANDARDISATION                                    │
│  • Taxpasta        – standardised profile tables                │
│  • MultiQC         – aggregated report                          │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
   Final results (see Section 4)
```

---

## 3. Software Used

### 3.1 Workflow engine and containers

| Software    | Role |
|------------|------|
| **Nextflow** | Workflow orchestration (DSL2). |
| **Apptainer** | Runs pipeline tools in containers (replaces Singularity). |

### 3.2 Data preprocessing

| Software    | Data type   | Role |
|------------|-------------|------|
| **FastQC** | Short-read  | Per-sample quality control (e.g. base quality, GC content, adapter content). |
| **fastp**  | Short-read  | Adapter trimming, quality filtering, optional merging; outputs analysis-ready FASTQs. |
| **NanoPlot / NanoStat** | Long-read | Long-read quality and length statistics. |
| **Porechop / Porechop_ABI** | Long-read | Adapter detection and removal for Nanopore (and compatible) data. |

### 3.3 Taxonomic classification

| Software    | Data type   | Role |
|------------|-------------|------|
| **Kraken2** | Short & long | K-mer–based taxonomic classification against a reference database (e.g. viral). |
| **Bracken** | Short-read only | Re-estimates species-level abundances from Kraken2 output; used when Bracken database is configured. |

### 3.4 Abundance analysis (MetaTaxProfiler)

| Component | Role |
|----------|------|
| **batch_calculate_abundance_en.sh** | Short-read: finds Bracken/Kraken2 outputs, dispatches to Python; falls back to long-read method if no Bracken. |
| **batch_calculate_abundance_longread_en.sh** | Long-read (or short-read fallback): finds Kraken2 reports, runs Python per sample and builds summary tables. |
| **calculate_abundance_en.py** | Reads Bracken (+ optional Kraken2) output; computes **RPM** (reads per million); writes per-sample TSV. Uses **Python 3**, **pandas**, **numpy**. |
| **calculate_abundance_longread_en.py** | Reads Kraken2 report; computes **RPM**; writes per-sample TSV. Uses **Python 3**, **pandas**, **numpy**. |

### 3.5 Reporting and standardisation

| Software    | Role |
|------------|------|
| **Taxpasta** | Converts profiler outputs to standardised tables (e.g. MOTU, BIOM) when `run_profile_standardisation` is enabled. |
| **MultiQC**  | Aggregates FastQC, fastp, Kraken2/Bracken, and other tool outputs into a single HTML report. |

---

## 4. Final result files

All paths are relative to the run output directory: **`results_viral_short/`** or **`results_viral_long/`** (set by `outdir` in the corresponding Nextflow config).

### 4.1 Classification and abundance (main results)

| Path / file | Description |
|-------------|-------------|
| **kraken2/** | Kraken2 classification results (e.g. `*.kraken2.report.txt`, report format). One report per sample/database. |
| **bracken/** | Bracken output (short-read only, if enabled). Species-level re-estimated read counts (e.g. `*_bracken*.tsv`, `*.bracken_species.tsv`). |
| **abundance/** | **MetaTaxProfiler abundance outputs:** |
| **abundance/<sample>_abundance.tsv** | Per-sample abundance table: taxonomy, assigned reads, fraction, **RPM**. |
| **abundance/all_samples_abundance_summary.tsv** | Combined table: all samples, all taxa, with RPM. |
| **abundance/top_viruses_summary.tsv** | Subset of taxa with RPM ≥ 10 (or equivalent threshold used in the script). |

### 4.2 Preprocessing and QC

| Path / file | Description |
|-------------|-------------|
| **fastqc/** (short-read) | FastQC HTML and JSON per sample. |
| **nanoq/** (long-read)   | Long-read QC statistics (NanoPlot/NanoQ pipeline outputs). |
| Preprocessed/analysis-ready FASTQs | If `save_preprocessed_reads` / `save_analysis_ready_fastqs` are set, these appear in the Nextflow `work` output or published directories as configured. |

### 4.3 Aggregated report

| Path / file | Description |
|-------------|-------------|
| **multiqc/multiqc_report.html** | MultiQC report: QC, classification, and pipeline summary in one HTML file. |

### 4.4 Standardised profiles (optional)

| Path / file | Description |
|-------------|-------------|
| **taxpasta/** (or equivalent) | Standardised profile tables (e.g. MOTU, BIOM) when `run_profile_standardisation` and related options are enabled. |

---

## 5. Abundance table columns (RPM-only)

After removal of RPKM, the main abundance TSVs typically contain:

| Column          | Description |
|-----------------|-------------|
| Sample          | Sample ID (in summary tables). |
| Species         | Taxonomic name (e.g. species or lowest assigned rank). |
| Taxonomy_ID     | NCBI taxonomy ID. |
| Assigned_Reads  | Read count assigned to that taxon. |
| Fraction        | Fraction of total reads (0–1). |
| **RPM**         | Reads per million: relative abundance normalized to 1 million reads. |

---

## 6. Quick reference

- **Short-read run:** `sbatch submit_short.sh` → preprocessing (FastQC, fastp) → Kraken2 → Bracken (if configured) → abundance scripts → MultiQC.  
- **Long-read run:** `sbatch submit_long.sh` → preprocessing (NanoPlot, Porechop) → Kraken2 → abundance scripts → MultiQC.  
- **Containers:** All pipeline tools run via Apptainer; config uses `apptainer.runOptions = '--no-mount /lscratch'` where needed.  
- **Main deliverables:** Kraken2 reports, Bracken tables (short-read), **abundance/*.tsv** (RPM), and **multiqc/multiqc_report.html**.

For sample sheets, database setup, and options, see **README.md**, **SHORTREAD_ABUNDANCE_GUIDE.md**, and **LONGREAD_GUIDE.md**.



