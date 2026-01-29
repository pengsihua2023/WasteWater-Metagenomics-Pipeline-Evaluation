# Workflows of nf-core/mag

This document summarizes the **workflows provided in this repository** for metagenome analysis, including:

- **nf-core/mag** runs (short-read and hybrid patterns, using **Nextflow + Apptainer**)
- additional **long-read utilities** (Flye + Kraken2), implemented as standalone scripts

It also lists **major tools used**, including data pre-processing and **taxonomy/abundance-style** reporting, and highlights the **key output files**.

[![Apptainer](https://img.shields.io/badge/Apptainer-%E2%89%A51.3-42b983)](https://apptainer.org/)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue.svg)](https://www.python.org/)

---

## What is implemented here (scripts)

### nf-core/mag (Nextflow) runs

- **Short reads, co-assembly, viral-focused**: `run_short_reads.sh`  
  - Output: `results_short/`
  - Post-step: Kraken2 classification of assembled contigs (MEGAHIT/SPAdes)

- **Short reads only (simpler wrapper)**: `run_nfcore_mag_v3_short_only.sh`  
  - Output: `results_mag_v3/`

- **Mixed samplesheet (short + long rows in one CSV)**: `run_nfcore_mag_v3.sh`  
  - Output: `results_mag_v3/`
  - Note: still configured as “viral-focused” and skips binning/QC/GTDB-Tk in the script.

- **Older stable nf-core/mag run (v2.5.4)**: `run_nfcore_mag_viral.sh`  
  - Output: `results_mag/`

### Long-read utility runs (standalone, not nf-core/mag)

- **Flye assembly + Kraken2 on contigs and reads**: `run_long_reads.sh`  
  - Output: `results_long/`

- **Kraken2 on long reads only**: `run_long_reads_kraken2_only.sh`  
  - Output: `results_long_kraken2/`

### Experimental / workaround scripts

- **Attempt nf-core/mag dev for long-reads-only**: `run_long_reads_nfcore_dev.sh`
- **Workaround: dummy short reads + real long reads** (to bypass strict validation): `run_long_reads_with_dummy.sh`

## Common runtime settings

- **Workflow engine**: Nextflow
- **Primary pipeline**: `nf-core/mag` (commonly `3.1.0` in this repo)
- **Container runtime**: **Apptainer** (`-profile apptainer`)

## High-level nf-core/mag workflow (generic)

nf-core/mag is a modular metagenome pipeline. Depending on parameters, it can run:

1) **Read pre-processing / QC** (short reads and/or long reads)
2) **Assembly** (short-read, long-read, or hybrid)
3) **Mapping / coverage estimation** (used for binning/abundance/coverage tables when enabled)
4) **Binning** (MAG reconstruction)
5) **Bin quality control** (BUSCO / CheckM)
6) **Taxonomic classification** (e.g., Kraken2/CAT/GTDB-Tk depending on settings)
7) **Reporting** (MultiQC + Nextflow reports)

In this repository, most nf-core/mag runs are configured to **focus on assembly + viral classification** and to **skip binning and MAG-oriented steps**.

## What actually ran in your successful short-read run (`run_short_reads.sh`)

### 1) Input & samplesheet

- Input is a Nextflow CSV samplesheet (`--input samplesheet_short.csv`) with paired-end short reads.
- Group-level co-assembly is enabled (`--coassemble_group`).

### 2) Data pre-processing (QC + filtering)

- **FastQC**: raw read quality reports
- **fastp**: trimming/filtering (e.g., `--fastp_qualified_quality 20`)
- **Bowtie2 PhiX removal**: read mapping to PhiX and removal of PhiX-like reads (enabled in your log)

### 3) Assembly

Assemblers executed in your run:

- **MEGAHIT**
- **SPAdes (meta)**

### 4) Gene prediction

- **Prodigal** (shown as `...DIGAL...` in Nextflow process names)

> Genome annotation with **Prokka** was skipped in the script.

### 5) Taxonomy / abundance-style reporting

- **Kraken2** (read classification) produces a report with **percentages and counts per taxon** (commonly used as an abundance proxy).
- **Krona** generates an interactive HTML summary.

Repo-specific detail:

- `--kraken2_db` is passed a placeholder file to satisfy schema validation in v3.x.
- The real Kraken2 DB is injected at process level via `nfcore_mag_v3_no_validation.config`:
  - `process.withName: 'NFCORE_MAG:MAG:KRAKEN2' { ext.args = "--db /scratch/.../kraken2_Viral_ref" }`

### 6) Reporting

- **MultiQC** aggregates QC and tool reports
- Nextflow execution artifacts were enabled: `-with-report`, `-with-timeline`

### 7) Extra post-processing (custom; outside nf-core/mag)

`run_short_reads.sh` additionally runs **Kraken2 on assembled contigs** (MEGAHIT and SPAdes) using Apptainer and writes separate contig classification outputs.

## Tools/software used (across these workflows)

### Workflow & containers

- **Nextflow**
- **nf-core/mag**
- **Apptainer**

### Short-read pre-processing / QC

- **FastQC**
- **fastp**
- **Bowtie2** (PhiX removal)

### Long-read pre-processing / QC (available in nf-core/mag; may be enabled in other modes)

- **NanoPlot**
- **porechop / porechop_abi**
- **Filtlong**

### Assembly

- **MEGAHIT** (short reads)
- **SPAdes** / **SPAdesHybrid** (short / hybrid)
- **Flye** (long-read utility script; and may be used in some nf-core/mag long-read contexts depending on version/profile)

### Gene prediction

- **Prodigal**

### Taxonomy / abundance-style profiling

- **Kraken2**
- **Krona**

### MAG-oriented steps (available in nf-core/mag but typically skipped in this repo)

- **Binning**: MetaBAT2, MaxBin2, CONCOCT (optionally DAS Tool refinement)
- **Bin QC**: BUSCO or CheckM
- **GTDB-Tk** classification (taxonomy of MAGs)

### Reporting

- **MultiQC**

## Key output files (successful short-read run: `results_short/`)

All paths below are relative to `results_short/`.

### Read-level taxonomy (Kraken2)

- `Taxonomy/kraken2/<SAMPLE>/<SAMPLE>.kraken2_report.txt`  
  Example:
  - `Taxonomy/kraken2/llnl_66ce4dde/llnl_66ce4dde.kraken2_report.txt`

### Assemblies

- `Assembly/MEGAHIT/MEGAHIT-<GROUP>.contigs.fa.gz`
- `Assembly/SPAdes/SPAdes-<GROUP>_contigs.fasta.gz`
- (optional) `Assembly/SPAdes/SPAdes-<GROUP>_scaffolds.fasta.gz`

Examples from your run summary:

- `Assembly/MEGAHIT/MEGAHIT-group-viral_group1.contigs.fa.gz`
- `Assembly/SPAdes/SPAdes-group-viral_group1_contigs.fasta.gz`
- `Assembly/SPAdes/SPAdes-group-viral_group1_scaffolds.fasta.gz`

### Contig-level taxonomy (custom post-step)

- `Taxonomy/contigs/llnl_66ce4dde_MEGAHIT_classification.txt`
- `Taxonomy/contigs/llnl_66ce4dde_MEGAHIT_kraken2_report.txt`
- `Taxonomy/contigs/llnl_66ce4dde_SPAdes_classification.txt`
- `Taxonomy/contigs/llnl_66ce4dde_SPAdes_kraken2_report.txt`

### MultiQC

The exact folder naming can vary by pipeline version/config. Search under `results_short/` for:

- `*multiqc*html`

### Nextflow execution reports

- `execution_report.html`
- `execution_timeline.html`

## How to interpret “abundance” in these runs

- **Read-level abundance proxy**: the Kraken2 report (`*.kraken2_report.txt`) provides per-taxon read counts and percentages.
- **Contig-level abundance proxy**: contig Kraken2 reports show the taxonomic distribution of assembled sequences.
- **MAG abundance / coverage tables**: only produced if binning + mapping/coverage steps are enabled; in most scripts here, they are intentionally skipped.

## Notes about skips (viral-focused configuration)

Many scripts in this repository explicitly disable time-consuming MAG-oriented steps:

- `--skip_binning`
- `--skip_binqc`
- `--skip_gtdbtk`
- `--skip_prokka`
- `--skip_quast`

This is intended when the goal is **assembly + viral taxonomy profiling** rather than full MAG reconstruction.



