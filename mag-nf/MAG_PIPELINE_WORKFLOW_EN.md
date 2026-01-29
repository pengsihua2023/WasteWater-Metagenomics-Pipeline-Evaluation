## nf-core/mag (v3.1.0) Short-read workflow 

This document describes the **actual workflow executed by `run_short_reads.sh`** in this project and summarizes **tools**, **major steps**, and **result files** (under `results_short/`).

### Overview

- **Workflow engine**: Nextflow
- **Pipeline**: `nf-core/mag` (release `3.1.0`)
- **Container runtime**: **Apptainer** (`-profile apptainer`)
- **Input type**: paired-end short reads (Illumina-like)
- **Co-assembly**: enabled (`--coassemble_group`) — samples within the same group are co-assembled
- **Binning / MAG QC / GTDB-Tk**: **skipped** by design (方案 A)
  - `skip_binning = true`, `skip_binqc = true`, `skip_gtdbtk = true` in `nfcore_mag_v3_no_validation.config`

### High-level flow (what happens, in order)

#### 1) Input & samplesheet

- The run uses a Nextflow CSV samplesheet (`--input samplesheet_short.csv`) containing:
  - sample id, group id
  - `short_reads_1`, `short_reads_2`
  - `long_reads` left empty for this run

#### 2) Read pre-processing (QC + filtering)

Goal: improve read quality before assembly and downstream profiling.

Main steps (typical for nf-core/mag short-read runs, and visible in your log):

- **FastQC** (raw reads)  
  Produces per-sample quality reports for raw input reads.

- **fastp** (adapter/quality trimming + filtering)  
  Performs read trimming and filtering (e.g. qualified quality threshold).  
  In this project run: `--fastp_qualified_quality 20` and `--fastp_save_trimmed_fail`.

- **Bowtie2 PhiX removal** (optional but enabled in the pipeline run you showed)  
  Maps reads to a PhiX reference and removes PhiX-like reads.

#### 3) Assembly (co-assembly per group)

Goal: build contigs from reads.

Assemblers used in this run:

- **MEGAHIT** (`MEGAHIT-group-viral_group1...`)
- **SPAdes (meta)** (`SPAdes-group-viral_group1...`)

Both assemblies were executed for the group and produced gzipped contig FASTA files.

#### 4) Gene prediction (for assemblies)

Goal: predict genes on assembled contigs to support downstream reporting.

- **Prodigal** (process shows as `...DIGAL...` in your log)

> Note: genome annotation with **Prokka** was explicitly skipped (`--skip_prokka`), so you get gene prediction outputs but not Prokka-annotated genomes.

#### 5) Taxonomic profiling (abundance-style reporting from reads)

Goal: estimate **taxonomic composition** of reads and provide **relative abundance** style summaries.

- **Kraken2** (read-level classification)  
  Produces a Kraken2 classification report (`*.kraken2_report.txt`) which contains per-taxon counts and percentages (commonly used as a taxonomic abundance proxy).

- **Krona** (interactive visualization)  
  Converts taxonomic reports into an interactive HTML view.

Implementation detail in this repo:

- The pipeline parameter `--kraken2_db` is passed a placeholder file to satisfy schema validation.
- The actual Kraken2 DB path is injected at process level in `nfcore_mag_v3_no_validation.config`:
  - `process.withName: 'NFCORE_MAG:MAG:KRAKEN2' { ext.args = "--db /scratch/.../kraken2_Viral_ref" }`

#### 6) Reporting

- **MultiQC** aggregates QC metrics and reports from multiple tools into a single report.
- Nextflow also generates execution reports/timelines as requested by `-with-report` and `-with-timeline`.

#### 7) Additional (custom) analysis: contig classification (outside nf-core/mag)

After the pipeline finishes, `run_short_reads.sh` performs an extra step:

- Run **Kraken2** on the **assembled contigs** (MEGAHIT and SPAdes) using Apptainer
- Writes separate classification and report files for each assembly

This is **not** part of nf-core/mag itself; it is an added post-processing step in your script.

### Tools/software used (this run)

#### Workflow & containers

- **Nextflow**
- **nf-core/mag** `3.1.0`
- **Apptainer** (container execution)

#### Read pre-processing / QC

- **FastQC** (raw read QC)
- **fastp** (short-read trimming/filtering)
- **Bowtie2** (PhiX removal mapping)

#### Assembly

- **MEGAHIT**
- **SPAdes** (meta assembly mode)

#### Gene prediction

- **Prodigal**

#### Taxonomy / “abundance” profiling

- **Kraken2** (read classification; report includes per-taxon % and counts)
- **Krona** (interactive taxonomy visualization)

#### Reporting

- **MultiQC**

### What is intentionally skipped (important for interpretation)

Because 方案 A is enabled in config/CLI, the following are **not** produced:

- **Metagenome binning** (MetaBAT2 / MaxBin2 / CONCOCT / DAS Tool)
- **MAG QC** (BUSCO / CheckM)
- **GTDB-Tk classification**
- **Genome annotation** with Prokka
- **Assembly QC** with QUAST

Therefore, **MAG-level abundance tables (bin depths/coverage summaries)** will not exist in this run. Your “abundance-style” outputs are primarily **taxonomic composition** from Kraken2 reports (reads and contigs).

### Key output files (this run)

All paths below are relative to the output directory: `results_short/`.

#### 1) Read-level taxonomy (Kraken2)

- `Taxonomy/kraken2/<SAMPLE>/<SAMPLE>.kraken2_report.txt`  
  Example from your script:
  - `Taxonomy/kraken2/llnl_66ce4dde/llnl_66ce4dde.kraken2_report.txt`

> If you also need the full per-read assignment output, look for `*.kraken2_output.txt` or similarly named files in the same folder (exact names can vary by pipeline version/config).

#### 2) Assembly outputs

- `Assembly/MEGAHIT/MEGAHIT-<GROUP>.contigs.fa.gz`
- `Assembly/SPAdes/SPAdes-<GROUP>_contigs.fasta.gz`
- (optional) `Assembly/SPAdes/SPAdes-<GROUP>_scaffolds.fasta.gz`

Examples shown in your run summary:

- `Assembly/MEGAHIT/MEGAHIT-group-viral_group1.contigs.fa.gz`
- `Assembly/SPAdes/SPAdes-group-viral_group1_contigs.fasta.gz`
- `Assembly/SPAdes/SPAdes-group-viral_group1_scaffolds.fasta.gz`

#### 3) Contig-level taxonomy (custom post-step)

Created by the post-processing section of `run_short_reads.sh`:

- `Taxonomy/contigs/llnl_66ce4dde_MEGAHIT_classification.txt`
- `Taxonomy/contigs/llnl_66ce4dde_MEGAHIT_kraken2_report.txt`
- `Taxonomy/contigs/llnl_66ce4dde_SPAdes_classification.txt`
- `Taxonomy/contigs/llnl_66ce4dde_SPAdes_kraken2_report.txt`

#### 4) MultiQC

- `MultiQC/multiqc_report.html` (typical location; exact subfolder name can vary)

If you do not see this exact path, search under `results_short/` for:

- `*multiqc*html`

#### 5) Nextflow execution reports

Requested explicitly in the script:

- `execution_report.html`
- `execution_timeline.html`

### How to interpret “abundance” in this run

- **Read-level abundance (taxonomic composition)**: use `*.kraken2_report.txt` to get per-taxon read counts and percentages.
- **Contig-level abundance proxy**: use the contig Kraken2 reports to see how assembled sequences distribute across taxa.
- **MAG abundance / depth (coverage) tables**: **not available** because binning and bin QC were skipped.

### Reproducibility notes

- This run uses Apptainer containers pulled from `quay.io/biocontainers/...` images.
- The exact command-line arguments are printed by the script (debug block) and can be copied into the log for auditing.



