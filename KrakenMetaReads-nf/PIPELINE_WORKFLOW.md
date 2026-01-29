# Hybrid Metagenome Assembly and Kraken2 Taxonomic Classification Pipeline

**Version:** 3.0.0  
**Workflow engine:** Nextflow (DSL2)  
**Container runtime:** Apptainer/Singularity  
**Environment:** Conda (for selected processes)

---

## 1. Pipeline Overview

This pipeline performs metagenomic assembly, read mapping, abundance estimation, and taxonomic classification for both **short-read** (Illumina) and **long-read** (Nanopore/PacBio) sequencing data. It supports running short-read only, long-read only, or both in a single execution.

---

## 2. Workflow Diagram

### 2.1 Short-Read Workflow

```
Raw paired-end reads (FASTQ)
        │
        ▼
┌───────────────────┐
│  Data Preprocessing │  fastp (QC, adapter trimming, filtering)
└─────────┬─────────┘
          │
          ▼
┌───────────────────┐
│     Assembly      │  MEGAHIT + metaSPAdes (parallel)
└─────────┬─────────┘
          │
          ├──────────────────────────────┐
          ▼                              ▼
┌───────────────────┐          ┌───────────────────┐
│  Bowtie2 index    │          │  Bowtie2 index    │
│  (MEGAHIT contigs)│          │  (SPAdes contigs)  │
└─────────┬─────────┘          └─────────┬───────────┘
          │                              │
          ▼                              ▼
┌───────────────────┐          ┌───────────────────┐
│  Bowtie2 align    │          │  Bowtie2 align     │
│  (clean reads →   │          │  (clean reads →    │
│   MEGAHIT)        │          │   SPAdes)          │
└─────────┬─────────┘          └─────────┬───────────┘
          │                              │
          ▼                              ▼
┌───────────────────┐          ┌───────────────────┐
│  Abundance        │          │  Abundance        │
│  (RPM, RPKM)       │          │  (RPM, RPKM)      │
│  MEGAHIT           │          │  SPAdes           │
└─────────┬─────────┘          └─────────┬───────────┘
          │                              │
          ▼                              ▼
┌───────────────────┐          ┌───────────────────┐
│  Kraken2          │          │  Kraken2           │
│  (MEGAHIT contigs) │          │  (SPAdes contigs) │
└─────────┬─────────┘          └─────────┬───────────┘
          │                              │
          └──────────────┬───────────────┘
                         ▼
              ┌───────────────────┐
              │  Merge Kraken2     │
              │  reports + virus   │
              │  consensus         │
              └───────────────────┘
```

### 2.2 Long-Read Workflow

```
Raw long reads (FASTQ)
        │
        ▼
┌───────────────────┐
│  metaFlye Assembly │  (general metagenome)
└─────────┬─────────┘
          │
          ├──────────────────────────────────────────┐
          ▼                                            │
┌───────────────────┐                                  │
│  Minimap2 align   │                                  │
│  (reads → Flye    │                                  │
│   contigs)        │                                  │
└─────────┬─────────┘                                  │
          │                                            │
          ▼                                            │
┌───────────────────┐                                  │
│  Abundance (Flye) │                                  │
│  RPM, RPKM         │                                  │
└─────────┬─────────┘                                  │
          │                                            │
          ▼                                            │
┌───────────────────┐                                  │
│  Kraken2 (Flye    │                                  │
│  contigs)         │                                  │
└───────────────────┘                                  │
                                                       │
          [Optional: viralFlye enabled]                │
                       │                               │
                       ▼                               │
              ┌───────────────────┐                    │
              │  viralFlye         │  (viral contig     │
              │  (linear +         │   identification   │
              │   circular)        │   from Flye)       │
              └─────────┬─────────┘                    │
                       │                               │
          ├────────────┼────────────┐                  │
          ▼            ▼            ▼                  │
    Minimap2      Minimap2      Minimap2               │
    (linear)      (circular)    (→ abundance +         │
          │            │         Kraken2 each)         │
          └────────────┴───────────────────────────────┘
```

---

## 3. Software Used

### 3.1 Data Preprocessing

| Software | Version | Purpose | Mode |
|----------|---------|---------|------|
| **fastp** | 0.23.4 | Quality control, adapter detection/trimming, length filtering, quality filtering | Conda |

*Parameters (configurable):* qualified quality Phred ≥20, unqualified percent limit 40%, min length 50 bp, auto-detect adapters for PE.

### 3.2 Assembly

| Software | Version | Purpose | Mode |
|----------|---------|---------|------|
| **MEGAHIT** | 1.2.9 | Short-read metagenomic assembly (k-mer based) | Apptainer |
| **metaSPAdes** | 3.15.5 | Short-read metagenomic assembly (multi-k) | Apptainer |
| **metaFlye** | 2.9.2 | Long-read metagenomic assembly | Apptainer |
| **viralFlye** | (conda env) | Viral contig identification from metaFlye output (linear + circular) | Conda (viralFlye_env) |

### 3.3 Read Mapping

| Software | Version | Purpose | Mode |
|----------|---------|---------|------|
| **Bowtie2** | 2.5.1 | Build index and align short reads to MEGAHIT/SPAdes contigs | Apptainer |
| **samtools** | (in container) | Sort and index BAM | Apptainer |
| **Minimap2** | (in container) | Align long reads to Flye/viralFlye contigs | Apptainer |
| **samtools** | 1.17 | Used in abundance scripts (idxstats) | Conda |

### 3.4 Abundance Analysis

| Component | Purpose | Mode |
|-----------|---------|------|
| **Custom Python script** (Biopython, samtools) | Compute **RPM** (reads per million) and **RPKM** (reads per kilobase per million) per contig from BAM + FASTA | Conda |

*Formulas:*
- RPM = (mapped reads to contig / total mapped reads) × 10⁶  
- RPKM = reads / (length/1000) / (total_mapped_reads/10⁶)

### 3.5 Taxonomic Classification

| Software | Version | Purpose | Mode |
|----------|---------|---------|------|
| **Kraken2** | 2.1.3 | Taxonomic classification of contigs using a user-provided database | Conda |

### 3.6 Report Merging and Virus Consensus (Short-Read Only)

| Component | Purpose | Mode |
|-----------|---------|------|
| **Custom Python script** (pandas, numpy) | Merge Kraken2 reports (MEGAHIT vs SPAdes); extract virus consensus (both assemblers) and detection category (Consensus/SPAdes only/MEGAHIT only) | Conda |

---

## 4. Final Output Files

Outputs are written under `results_short/` (short-read) and `results_long/` (long-read) by default.

### 4.1 Short-Read Outputs (`results_short/`)

| Subdirectory | Contents |
|--------------|----------|
| **fastp/** | `*_fastp.html`, `*_fastp.json` — QC reports and cleaned read stats |
| **abundance_megahit/** | `*_megahit_abundance.txt` (Contig_ID, Length, Mapped_Reads, RPM, RPKM), `*_megahit_abundance_summary.txt` |
| **abundance_spades/** | `*_spades_abundance.txt`, `*_spades_abundance_summary.txt` |
| **kraken2_megahit/** | `*_megahit_classification.txt`, `*_megahit_report.txt` |
| **kraken2_spades/** | `*_spades_classification.txt`, `*_spades_report.txt` |
| **merged_reports/** | `*_merged_report.txt`, `*_merged_report.csv` (merged Kraken2 taxonomy); `*_virus_consensus.txt`, `*_virus_consensus.csv` (virus-focused consensus) |

*Note:* Assembled contigs and BAMs are produced in the pipeline but are not necessarily published to these folders by default; the main delivered results are QC, abundance tables, Kraken2 reports, and merged/virus consensus reports.

### 4.2 Long-Read Outputs (`results_long/`)

| Subdirectory | Contents |
|--------------|----------|
| **flye_assembly/** | `*_flye_assembly/` — metaFlye assembly directory (e.g. `assembly.fasta`) |
| **viralflye/** | viralFlye outputs: linear/circular viral contigs and components (when `run_viralflye=true`) |
| **abundance_flye/** | `*_flye_abundance.txt`, `*_flye_abundance_summary.txt` |
| **abundance_viralflye_linear/** | `*_viralflye_linear_abundance.txt`, `*_viralflye_linear_abundance_summary.txt` |
| **abundance_viralflye_circular/** | `*_viralflye_circular_abundance.txt`, `*_viralflye_circular_abundance_summary.txt` |
| **kraken2_flye/** | `*_flye_classification.txt`, `*_flye_report.txt` |
| **kraken2_viralflye_linear/** | `*_viralflye_linear_*.txt` (classification + report) |
| **kraken2_viralflye_circular/** | `*_viralflye_circular_*.txt` (classification + report) |

---

## 5. Summary Table of Key Result Files

| Category | Short-Read | Long-Read |
|----------|------------|-----------|
| **Preprocessing** | `*_fastp.html`, `*_fastp.json` | — |
| **Abundance** | `*_megahit_abundance.txt`, `*_spades_abundance.txt` (+ summaries) | `*_flye_abundance.txt` (+ viralflye linear/circular if enabled) (+ summaries) |
| **Taxonomy** | `*_megahit_report.txt`, `*_spades_report.txt` | `*_flye_*.txt` (+ viralflye linear/circular reports if enabled) |
| **Merged / virus** | `*_merged_report.csv`, `*_virus_consensus.csv` (and `.txt`) | — |

---

## 6. Input Requirements

- **Short-read:** CSV with columns `sample`, `fastq_1`, `fastq_2` (paired-end FASTQ paths).
- **Long-read:** CSV with columns `sample`, `fastq_long` (long-read FASTQ path).
- **Kraken2:** Path to a pre-built Kraken2 database (required).
- **viralFlye (optional):** Path to Pfam HMM file (e.g. Pfam-A.hmm) when running long-read viral identification.

---

## 7. Running the Pipeline

- **Short-read only:** e.g. `run_short_only.sh` or Nextflow with `--input_short samplesheet_short.csv` and `--kraken2_db /path/to/db`.
- **Long-read only:** e.g. `run_long_only.sh` or Nextflow with `--input_long samplesheet_long.csv` and `--kraken2_db /path/to/db`.
- **Both:** e.g. `run_hybrid_workflow.sh` or Nextflow with both `--input_short` and `--input_long` and `--kraken2_db`.

Configuration (resources, paths, options) is set in `metagenome_hybrid_workflow.config`.
