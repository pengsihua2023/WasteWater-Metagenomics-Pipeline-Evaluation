## rvdb-viral-metagenome-nf Pipeline Overview

This repository provides a **Nextflow DSL2** pipeline for metagenome assembly and protein-level classification against a **DIAMOND** database (e.g., RVDB). It supports:

- **Short reads** (Illumina paired-end)
- **Long reads** (Nanopore / PacBio)

The pipeline runs on an HPC scheduler (e.g., **SLURM**) and primarily uses **Conda** environments (some assemblers can be provided via pre-installed environments added to `PATH`). **Apptainer/Singularity is optional and not enabled by default** in the provided config.

[![Nextflow](https://img.shields.io/badge/nextflow-%3D%2024.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%3D%203.9-blue.svg)](https://www.python.org/)

---

## Key software used

### Workflow engine
- **Nextflow** (DSL2)

### Short-read track (Illumina)
- **fastp**: read QC and adapter/quality trimming (optional)
- **MEGAHIT**: metagenome assembly
- **metaSPAdes** (`metaspades.py`): metagenome assembly (resource intensive)
- **Prodigal**: gene prediction (`-p meta`)
- **DIAMOND**: `blastp` classification against a protein database (`.dmnd`)
- **Python + pandas**: report generation (merge + taxonomy enrichment)
- **BWA + samtools**: read mapping back to contigs for abundance (RPM/RPKM)

### Long-read track (Nanopore / PacBio)
- **Flye** (MetaFlye mode): long-read metagenome assembly (`flye --meta`)
- **viralFlye**: optional refinement to focus on viral contigs (requires Pfam-A HMM)
- **Prodigal**: gene prediction (`-p meta`)
- **DIAMOND**: `blastp` classification
- **Python + pandas**: taxonomy enrichment and summary statistics
- **minimap2**: long-read mapping back to contigs for abundance (RPM/RPKM)

### Databases / reference files
- **DIAMOND database**: `--diamond_db /path/to/*.dmnd`
- **NCBI taxonomy dumps** (for lineage resolution):
  - `--taxonomy_names /path/to/names.dmp`
  - `--taxonomy_nodes /path/to/nodes.dmp`
- **Pfam-A HMM** (only for viralFlye): `--pfam_hmm /path/to/Pfam-A.hmm`

---

## Inputs

### 1) Short reads samplesheet (paired-end)
CSV columns:

```text
sample,fastq_1,fastq_2
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz
```

### 2) Long reads samplesheet
CSV columns:

```text
sample,fastq_long
SAMPLE1,/path/to/SAMPLE1.fastq.gz
```

---

## How the pipeline works

The pipeline chooses the branch based on `--read_type`:

- `--read_type short` → Illumina paired-end workflow
- `--read_type long`  → Nanopore/PacBio workflow

### Short-read workflow (Illumina paired-end)

#### Stage 0 — QC (optional)
1. **FASTP**
   - Input: `R1/R2`
   - Output: cleaned reads `*_clean_R1.fastq.gz`, `*_clean_R2.fastq.gz` + `*_fastp.html/json`

#### Stage 1 — Assembly (parallel assemblers)
2. **MEGAHIT_ASSEMBLY**
   - Input: cleaned `R1/R2`
   - Output: `${sample}_megahit_contigs.fa`

3. **SPADES_ASSEMBLY (metaSPAdes)**
   - Input: cleaned `R1/R2`
   - Output: `${sample}_spades_contigs.fa`
   - Note: metaSPAdes can be extremely memory-hungry on very large datasets.
     - In this pipeline, SPAdes failures can be handled by emitting an empty contigs file and continuing downstream (so MEGAHIT-based results can still complete).

#### Stage 2 — Gene prediction
4. **PRODIGAL_MEGAHIT**
   - Input: MEGAHIT contigs
   - Output: `${sample}_megahit_proteins.faa` + `${sample}_megahit_genes.fna`

5. **PRODIGAL_SPADES**
   - Input: SPAdes contigs
   - Output: `${sample}_spades_proteins.faa` + `${sample}_spades_genes.fna`
   - If contigs are empty, empty Prodigal outputs are created to keep the workflow consistent.

#### Stage 3 — Protein classification
6. **DIAMOND_CLASSIFICATION_MEGAHIT**
   - Input: MEGAHIT proteins + `--diamond_db`
   - Output: `${sample}_megahit_diamond.txt`

7. **DIAMOND_CLASSIFICATION_SPADES**
   - Input: SPAdes proteins + `--diamond_db`
   - Output: `${sample}_spades_diamond.txt`
   - If proteins are empty, an empty Diamond output is created.

#### Stage 4 — Comparative reporting + taxonomy (optional but enabled by default)
8. **MERGE_DIAMOND_REPORTS**
   - Inputs:
     - MEGAHIT Diamond report
     - SPAdes Diamond report
     - `names.dmp` + `nodes.dmp`
   - What it does:
     - Parses Diamond tabular output
     - Adds NCBI taxonomy lineage (superkingdom → species) per hit
     - Produces a merged comparison report and per-assembler enhanced tables

#### Stage 5 — Viral abundance estimation (optional but enabled by default)
9. **CALCULATE_ABUNDANCE_MEGAHIT_SHORT**
10. **CALCULATE_ABUNDANCE_SPADES_SHORT**

Abundance logic (short reads):
- Uses **Diamond results** to determine which contigs are considered “viral candidates”
  (by extracting contig IDs from `qseqid`).
- Maps paired-end reads to contigs using:
  - `bwa index`
  - `bwa mem`
  - `samtools view/sort/index/idxstats`
- Computes:
  - **RPM**: reads per million mapped reads
  - **RPKM**: reads per kilobase per million mapped reads

---

### Long-read workflow (Nanopore / PacBio)

#### Stage 1 — Assembly
1. **METAFLYE_ASSEMBLY**
   - Uses **Flye** in metagenome mode (`flye --meta`)
   - Output: `${sample}_metaflye_contigs.fa`

2. **VIRALFLYE_REFINEMENT** (optional, controlled by `--skip_viralflye`)
   - Refines / filters viral contigs using **viralFlye**
   - Requires `--pfam_hmm` (Pfam-A HMM database)
   - Output: `${sample}_viralflye_contigs.fa`
   - If viralFlye fails, the pipeline can fall back to using MetaFlye contigs as output.

#### Stage 2 — Gene prediction
3. **PRODIGAL_METAFLYE**
   - Output: `${sample}_metaflye_proteins.faa` + `${sample}_metaflye_genes.fna`

4. **PRODIGAL_VIRALFLYE** (if viralFlye enabled)
   - Output: `${sample}_viralflye_proteins.faa` + `${sample}_viralflye_genes.fna`

#### Stage 3 — Protein classification
5. **DIAMOND_CLASSIFICATION_METAFLYE**
   - Output: `${sample}_metaflye_diamond.txt`

6. **DIAMOND_CLASSIFICATION_VIRALFLYE** (if viralFlye enabled)
   - Output: `${sample}_viralflye_diamond.txt`

#### Stage 4 — Taxonomy enrichment
7. **ADD_TAXONOMY_METAFLYE**
   - Outputs:
     - `${sample}_metaflye_diamond_with_taxonomy.txt`
     - `${sample}_metaflye_taxonomy_summary.txt`

8. **ADD_TAXONOMY_VIRALFLYE** (if viralFlye enabled)
   - Outputs:
     - `${sample}_viralflye_diamond_with_taxonomy.txt`
     - `${sample}_viralflye_taxonomy_summary.txt`

#### Stage 5 — Dual-track comparison (optional, when viralFlye is enabled)
9. **COMPARE_DUAL_TRACKS**
   - Compares MetaFlye (all contigs) vs viralFlye (viral-focused contigs)
   - Produces “consensus” and “track-specific” virus sets

#### Stage 6 — Viral abundance estimation (optional)
10. **CALCULATE_ABUNDANCE_METAFLYE**
11. **CALCULATE_ABUNDANCE_VIRALFLYE** (if viralFlye enabled)

Abundance logic (long reads):
- Uses **Diamond results** to identify candidate viral contigs.
- Maps long reads to contigs using **minimap2** (preset `map-ont` by default, configurable via `long_read_preset`).
- Counts mapped reads per contig from the SAM alignment and computes **RPM/RPKM**.

---

## Output structure (results directory)

The pipeline writes outputs under `--outdir` with per-stage subfolders.

### Short-read outputs (`--read_type short`)
- `fastp/`
  - `${sample}_fastp.html`
  - `${sample}_fastp.json`
- `assembly_megahit/`
  - `${sample}_megahit_contigs.fa`
- `assembly_spades/`
  - `${sample}_spades_contigs.fa`
- `prodigal_megahit/`
  - `${sample}_megahit_proteins.faa`
  - `${sample}_megahit_genes.fna`
- `prodigal_spades/`
  - `${sample}_spades_proteins.faa`
  - `${sample}_spades_genes.fna`
- `diamond_megahit/`
  - `${sample}_megahit_diamond.txt`
- `diamond_spades/`
  - `${sample}_spades_diamond.txt`
- `merged_reports/` (if enabled)
  - `${sample}_merged_report.txt`
  - `${sample}_merged_report.csv`
  - `${sample}_megahit_with_taxonomy.txt`
  - `${sample}_spades_with_taxonomy.txt`
- `abundance_megahit/` (if enabled)
  - `${sample}_megahit_abundance.txt`
  - `${sample}_megahit_abundance.csv`
- `abundance_spades/` (if enabled)
  - `${sample}_spades_abundance.txt`
  - `${sample}_spades_abundance.csv`

### Long-read outputs (`--read_type long`)
- `assembly_metaflye/`
  - `${sample}_metaflye_contigs.fa`
- `assembly_viralflye/` (if enabled)
  - `${sample}_viralflye_contigs.fa`
- `prodigal_metaflye/`
  - `${sample}_metaflye_proteins.faa`
  - `${sample}_metaflye_genes.fna`
- `prodigal_viralflye/` (if enabled)
  - `${sample}_viralflye_proteins.faa`
  - `${sample}_viralflye_genes.fna`
- `diamond_metaflye/`
  - `${sample}_metaflye_diamond.txt`
- `diamond_viralflye/` (if enabled)
  - `${sample}_viralflye_diamond.txt`
- `taxonomy_metaflye/`
  - `${sample}_metaflye_diamond_with_taxonomy.txt`
  - `${sample}_metaflye_taxonomy_summary.txt`
- `taxonomy_viralflye/` (if enabled)
  - `${sample}_viralflye_diamond_with_taxonomy.txt`
  - `${sample}_viralflye_taxonomy_summary.txt`
- `consensus_analysis/` (if enabled)
  - `${sample}_consensus_viruses.txt`
  - `${sample}_metaflye_only_viruses.txt`
  - `${sample}_viralflye_only_viruses.txt`
  - `${sample}_dual_track_comparison.txt`
- `abundance_metaflye/` (if enabled)
  - `${sample}_metaflye_abundance.txt`
  - `${sample}_metaflye_abundance.csv`
- `abundance_viralflye/` (if enabled)
  - `${sample}_viralflye_abundance.txt`
  - `${sample}_viralflye_abundance.csv`

---

## Notes on execution environment

- The provided configuration typically uses:
  - Nextflow `conda { enabled = true }` for tools like fastp/Prodigal/Diamond/pandas.
  - Pre-installed assembler environments can be injected into task `PATH` via `beforeScript`.
- **Apptainer/Singularity** is not enabled by default in the config (the block is commented out).

