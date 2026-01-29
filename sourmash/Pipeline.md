# sourmash Pipeline Overview

This repository contains a Nextflow DSL2 pipeline for metagenomic screening using **sourmash** MinHash signatures. It supports **short reads** and **long reads**, with optional **assembly** prior to sketching.

[![Apptainer](https://img.shields.io/badge/Apptainer-%3D%201.3-42b983)](https://apptainer.org/)
[![Nextflow](https://img.shields.io/badge/nextflow-%3D%2024.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%3D%203.10-blue.svg)](https://www.python.org/)

---

## Overview

At a high level, the pipeline:

1. **Reads inputs** (short/long reads) from a samplesheet or glob patterns.
2. Optionally **assembles** reads into contigs (short reads: MEGAHIT and/or metaSPAdes; long reads: metaFlye).
3. **Creates sourmash DNA sketches** (MinHash signatures) from reads and/or contigs.
4. Runs **sourmash gather** against a reference database (viral by default; bacterial optional).
5. Writes a lightweight **HTML report** and Nextflow run reports (if enabled).

## Inputs

The pipeline accepts any of the following:

- **Short-read samplesheet**: `--samplesheet_short <csv>`
- **Long-read samplesheet**: `--samplesheet_long <csv>`
- **Legacy combined samplesheet**: `--samplesheet <csv>`
- **Globs**:
  - `--reads "path/*_{R1,R2}.fastq.gz"` (paired-end)
  - `--long_reads "path/*.fastq.gz"` (single-end)

### Samplesheet columns

The pipeline is permissive about column names and tries multiple aliases:

- **Short reads (paired-end)**:
  - R1: `fastq_1` or `R1` or `read1` or `fastq_r1`
  - R2: `fastq_2` or `R2` or `read2` or `fastq_r2`
  - Optional ID: `sample_id` or `sample` or `id` (otherwise inferred from filename)

- **Long reads (single-end)**:
  - long read: `fastq_long` or `long_read` or `long_reads`
  - Optional ID: `sample_id` or `sample` or `id` (otherwise inferred from filename)

## Workflow steps

### 1) Input handling & minimal preprocessing

- For **compressed inputs**, steps that require plain FASTQ will decompress using `zcat` into temporary files.
- For `SKETCH_READS` (no assembly), paired-end reads are concatenated into a single temporary FASTQ file.
- For short-read assembly, the pipeline includes a small robustness feature: if a samplesheet path ends with `.fastq` but the actual file is `.fastq.gz` (or vice versa), it tries to auto-correct the suffix before failing.

**Note:** The pipeline does **not** perform quality control or trimming (no `fastp`, `Trimmomatic`, etc.). If you need QC/trimming, run it upstream and point the samplesheet to the cleaned FASTQs.

### 2) Assembly (optional)

#### Short reads

Controlled by:

- `--skip_assembly` (default: `false`)
- `--assembler` (`megahit`, `metaspades`, or `both`)
- `--run_both_assemblers` (alternative flag to run both)

Assemblers:

- **MEGAHIT** (`MEGAHIT_ASSEMBLY`)
- **metaSPAdes** (`META_SPADES`)

#### Long reads

Controlled by:

- `--skip_long_assembly` (default: `false`)
- `--long_assembler` (currently supports `metaflye`)

Assembler:

- **metaFlye** (`METAFLYE_ASSEMBLY`, uses Flye `--meta`)

### 3) Sketching (MinHash signatures)

The pipeline creates **sourmash DNA signatures** using:

- `sourmash sketch dna -p k=<k_value>,scaled=<scaled>`

Sketching modes:

- **From contigs**: `SKETCH_CONTIGS`
  - Output naming: `<sample_id>_<assembler>_contigs.sig`
- **From short reads (no assembly)**: `SKETCH_READS`
  - Output naming: `<sample_id>_short_reads.sig`
- **From long reads (no assembly)**: `SKETCH_LONG_READS`
  - Output naming: `<sample_id>_long_reads.sig`

### 4) Screening against reference databases (gather)

The pipeline runs:

- **Viral screening** (always): `VIRAL_GATHER`
- **Bacterial screening** (optional): `BACTERIAL_GATHER` (enabled when `--bacterial_db` is provided)

Command used (conceptually):

- `sourmash gather <query.sig> <db.zip> --threshold-bp <threshold_bp> -o <out.csv> --save-matches <matches.sig>`

The pipeline also writes a small per-sample text summary file that includes either:

- the first lines of the CSV (if produced), or
- a note that no matches above threshold were found.

### 5) Report generation

`GENERATE_REPORT` produces a simple HTML page (`analysis_report.html`) containing run parameters and a pointer to the CSV results.

## Software used

- **Nextflow** (workflow engine; DSL2)
- **Slurm** (executor)
- **Conda** (used by run scripts and by the workflow when `use_local_*` is enabled)
- **Assemblers**:
  - **MEGAHIT**
  - **metaSPAdes** (SPAdes)
  - **Flye** in metagenome mode (**metaFlye**)
- **sourmash**:
  - `sourmash sketch dna`
  - `sourmash gather`
- Common Unix tools: `zcat`, `cat`, `awk`, `find`, `head`, `ls`

### Containerization (Apptainer/Singularity)

The configuration enables Apptainer (`apptainer.enabled = true`), and the workflow defines container images for some steps. However, many steps can be forced to run locally by setting:

- `--use_local_sourmash true`
- `--use_local_megahit true`
- `--use_local_spades true`
- `--use_local_flye true`

When these are `true`, the pipeline uses `conda run -n sourmash_env ...` and does **not** execute the container image for that step.

## Abundance / quantification (important note)

This pipeline performs **screening based on MinHash signature overlap**, not read-mapping-based quantification.

- The produced `*_viral.csv` / `*_bacterial.csv` files contain **gather metrics** (e.g., signature overlap and related statistics) generated by sourmash.
- The workflow uses `sourmash sketch dna` without explicit abundance tracking flags, so it is **not** generating abundance-tracking signatures by default.

If you need true abundance/relative abundance estimates, you typically need either:

- abundance-tracking sketches (sourmash `--track-abundance` workflow), and/or
- read mapping / coverage-based approaches (outside the scope of this pipeline).

## Outputs (results directory)

All outputs are written under `--outdir` (default: `results`, commonly `results_short`).

### Assemblies

- `assembly/megahit/<sample_id>/contigs.fasta`
- `assembly/metaspades/<sample_id>/contigs.fasta`
- `assembly/metaflye/<sample_id>/contigs.fasta`

Which ones appear depends on your chosen inputs and assembler settings.

### Signatures

In `signatures/`:

- `<sample_id>_<assembler>_contigs.sig`
- `<sample_id>_short_reads.sig` (when `--skip_assembly true` for short reads)
- `<sample_id>_long_reads.sig` (when `--skip_long_assembly true` for long reads)

### Viral screening results

In `viral_detection/`:

- `<sample_id>_viral.csv`
- `<sample_id>_viral_summary.txt`

### Bacterial screening results (optional)

In `bacterial_detection/` (only if `--bacterial_db` is set):

- `<sample_id>_bacterial.csv`
- `<sample_id>_bacterial_summary.txt`

### Report

In the output root:

- `analysis_report.html`

### Nextflow run reports (if enabled)

If you run with `-with-report`, `-with-trace`, and `-with-timeline` (as in `run_short.sh`), you will also get:

- `nextflow_report.html`
- `nextflow_trace.txt`
- `nextflow_timeline.html`

## Example: short reads only (MEGAHIT)

The provided `run_short.sh` launches a short-read-only run using:

- MEGAHIT assembly
- sourmash sketch on contigs
- sourmash gather against the configured viral database

