# Summary of 8 Viral Metagenomics Pipelines (Overview)

This document summarizes the 8 README markdown files in this folder (each describing one pipeline/toolchain). It highlights each pipeline’s goal, typical inputs/outputs, key steps, major tools and databases, and recommended use cases for quick comparison and selection.

## Comparison Table (At a Glance)

| Pipeline | Primary goal | Supported data types | Assembly | Core identification/classification method | Key databases | Main outputs |
|---|---|---|---|---|---|---|
| **Sourmash** | Fast similarity search via MinHash (contigs vs reference DB) | Assembled contigs (short/long-read assemblies) | Uses external assemblies (MetaFlye/SPAdes/MEGAHIT) | `sourmash sketch` + `sourmash gather` | NCBI viruses (k=31, zip) | `.sig` signatures; `gather_results*.csv` hits table |
| **rvdb-viral-metagenome-nf** | Viral protein similarity-based identification + (short-read dual-assembler comparison / long-read dual-track comparison) | Illumina (short) + Nanopore/PacBio (long) | ✅ (short: MEGAHIT+SPAdes; long: MetaFlye) | Prodigal gene calling + Diamond BLASTP (vs RVDB) + full lineage comparison | RVDB protein DB (Diamond index) + NCBI taxonomy; long-read also uses Pfam-A (viralFlye) | 7-level taxonomy comparison reports; dual-track consensus/difference lists; optional RPM/RPKM |
| **mag-nf (nf-core/mag approach)** | Assembly + Kraken2 classification (reads and contigs) | Short reads (nf-core/mag); long reads (direct containers) | ✅ (short: MEGAHIT+SPAdes; long: Flye) | Kraken2 (reads + contigs) | Kraken2 viral DB | `results_short/` and `results_long/` (Kraken2 reports, contig classification, etc.) |
| **TaxProfiler-nf (MetaTaxProfiler)** | Standardized taxonomic profiling + automatic viral abundance (RPM/RPKM) | Illumina + Nanopore/PacBio | ❌ (classification-focused) | Kraken2; optional Bracken correction for short reads; automatic summaries/reports | Kraken2 DB (Bracken DB recommended for short reads) | `results_viral_short/abundance/` or `results_viral_long/abundance/` (summary tables, top viruses) |
| **MLMVD-nf (ML-enhanced discovery)** | “Novel/distant virus” discovery with multi-tool parallel calls + cross-validation | Short + long reads | ✅ (short: MEGAHIT+SPAdes; long: metaFlye) | VirSorter2 (hybrid ML) + DeepVirFinder (deep learning) + viralFlye (Pfam domain validation; long-read) + intersection stratification | VirSorter2 DB; DeepVirFinder models; Pfam-A (viralFlye) | 3-tool comparison report (1/2/3-tool consensus tiers); high-confidence lists; optional RPM/RPKM |
| **KrakenMetaReads-nf** | TaxProfiler-centered Kraken2/Bracken profiling + (optional multi-assembly) + abundance | Short + long reads | ✅ (short: MEGAHIT+SPAdes; long: Flye+ViralFlye) | Kraken2; Bracken for short reads; batch processing + RPM/RPKM; BIOM output (as described) | Kraken2 DB + (short) Bracken DB | `results_viral_short/` and `results_viral_long/` (abundance_* folders, merged reports, etc.) |
| **GOTTCHA2** | Metagenome profiling / viral classification with GOTTCHA2 | Short + long reads | ❌ (mostly read-based) | GOTTCHA2 (k-mer signature profiling) | GOTTCHA2 DB (`.mmi`) | `*.tsv` / `*.summary.tsv` / `*.full.tsv` |
| **CLARK** | Fast classification using discriminative k-mers (incl. virus DB building + abundance estimation) | Single-end or paired-end (read-centric) | ❌ (tool itself does not assemble) | CLARK / CLARK-l / CLARK-S; built-in abundance estimation or custom scripts | NCBI RefSeq (can build viruses-only targets) | Classification CSV; optional abundance tables via scripts |

## Key Points by Pipeline

### 1) `sourmash README.md` (Sourmash Metagenome Analysis Pipeline)
- **Positioning**: Convert contigs into MinHash signatures (sketch), then perform containment-style searches against a reference DB (gather). Good for fast screening/similarity lookups.
- **Input**: Contigs from different assemblers (e.g., MetaFlye/SPAdes/MEGAHIT).
- **Key parameters**:
  - sketch: k=21/31/51, scaled=1000, abundance tracking enabled (`abund`).
  - gather: k=31 (matches DB), threshold 300 bp.
- **Outputs**:
  - `.sig`: signature files.
  - `gather_results*.csv`: overlap (`intersect_bp`), containment ANI, abundance statistics, etc.
- **Best for**: You already have assemblies and want a quick indication of “what known viral references these contigs resemble,” plus fast comparative screening.

### 2) `rvdb-viral-metagenome-nf README.md` (RVDB + Dual-Strategy Viral Metagenomics)
- **Positioning**: Protein-similarity-driven viral identification using RVDB + Diamond. Short reads use **dual assembler comparison**; long reads use **viralFlye feature filtering + Diamond dual-track comparison** for complementarity and coverage.
- **Short-read flow**: fastp → (MEGAHIT ∥ SPAdes) → Prodigal → Diamond vs RVDB → 7-level taxonomy comparison → (optional) BWA-based RPM/RPKM.
- **Long-read flow**: MetaFlye → viralFlye (Pfam + viralVerify) → dual tracks:
  - Track 1: all contigs → Prodigal → Diamond → taxonomy
  - Track 2: viralFlye-filtered viral contigs → Prodigal → Diamond → taxonomy
  - Outputs consensus/differences (e.g., “Diamond-positive but viralFlye-filtered” distant candidates).
- **Output highlights**: Short-read side-by-side taxonomy comparisons (Kingdom→Species); long-read consensus virus lists and MetaFlye-only distant-candidate tiers.
- **Best for**: You want assembler-to-assembler comparisons (short reads) and/or feature-based vs similarity-based complementarity (long reads) with fully interpretable taxonomy statistics.

### 3) `mag-nf README.md` (nf-core/mag + Custom Contig Classification; Long-read Direct Run)
- **Positioning**: Treat nf-core/mag (short reads) as a stable, engineered pipeline, and add the missing piece: **contig-level Kraken2 classification**. For long reads, use a simpler Flye + Kraken2 containerized approach.
- **Short reads**: fastp, PhiX removal (Bowtie2), MEGAHIT+SPAdes, Kraken2 on reads (nf-core/mag) + Kraken2 on contigs (custom), MultiQC.
- **Long reads**: Flye assembly, then Kraken2 on contigs and raw reads separately.
- **Best for**: You want nf-core standardization for short reads but still need contig-level classification; for long reads you want a straightforward deployable workflow.

### 4) `TaxProfiler-nf README.md` (MetaTaxProfiler: From Classification to Standardized Abundance Outputs)
- **Positioning**: A one-command workflow packaging “Kraken2 classification + (optional Bracken for short reads) + automatic RPM/RPKM + reports” to produce standardized abundance tables.
- **Input**: Samplesheet (short/long formats), database config via `databases.csv`.
- **Outputs**:
  - Under `abundance/`: per-sample detailed tables, all-samples summaries, and a “Top viruses” table (default RPM≥10).
- **Best for**: Your primary deliverable is a clean viral abundance matrix ready for downstream statistics/plots (rather than deep novel-virus discovery or complex cross-validation).

### 5) `MLMVD-nf README.md` (Machine Learning-Enhanced Novel Viral Discovery + Multi-Tool Consensus)
- **Positioning**: Designed for discovering novel/distant viruses using **VirSorter2 + DeepVirFinder** (short+long) plus **viralFlye (Pfam domain validation; long-read)**, then stratifying results by tool agreement (confidence tiers).
- **Long-read (key highlight)**: metaFlye → (VS2 ∥ DVF ∥ viralFlye) in parallel → 3-tool comparison:
  - 3-tool consensus: highest confidence
  - 2-tool consensus: medium confidence
  - single-tool calls: exploratory set
- **Short-read**: MEGAHIT+SPAdes in parallel; run VS2 and DVF on each assembly and compare across assemblers.
- **Best for**: Maximizing candidate discovery while still keeping a clear, multi-evidence confidence framework—especially useful for viral “dark matter” exploration in new environments.

### 6) `KrakenMetaReads-nf README.md` (TaxProfiler-Centered Kraken2/Bracken + Abundance)
- **Positioning**: Standardize short/long read profiling around nf-core/taxprofiler, add RPM/RPKM abundance calculation and batch processing, and (as described) generate BIOM outputs for downstream analysis.
- **Short reads**: Kraken2 (optionally with Bracken) on MEGAHIT and SPAdes assemblies → `abundance_megahit/` and `abundance_spades/`.
- **Long reads**: Flye + ViralFlye (circular/linear separation) followed by Kraken2 → corresponding abundance directories.
- **Best for**: You want standardized Kraken2/Bracken-based profiling and abundance outputs across both short and long reads, managed as an engineered Nextflow workflow with batch support.

### 7) `GOTTCHA2 README.md` (GOTTCHA2 Profiling Scripts)
- **Positioning**: Read-based profiling/viral classification using GOTTCHA2, wrapped as SLURM submission scripts; relatively lightweight for quick profiling.
- **Input**: Paired-end short-read FASTQ or single-end long-read FASTQ.
- **Outputs**: `.tsv`, `.summary.tsv`, `.full.tsv`.
- **Best for**: A simpler profiling path when you already have GOTTCHA2 DB and operational habits.

### 8) `CLARK README.md` (CLARK Installation & Virus Classification Guide)
- **Positioning**: CLARK family (including CLARK-l/CLARK-S) performs fast classification using discriminative k-mers. The README focuses on installation, building viruses-only targets, classification, and abundance estimation.
- **Key notes**:
  - Building **viruses-only** targets can reduce DB size and cost.
  - k-mer size impacts sensitivity: k=20/21 is often used for higher sensitivity in virus detection; k=31 is a balanced default.
  - Abundance can be computed using `estimate_abundance.sh` or custom scripts.
- **Best for**: High-speed CLARK-based classification, and as a reference for building a RefSeq-derived viruses-only database workflow.

## Selection Guide (How to Choose)

- **If you mainly need standardized viral abundance tables (RPM/RPKM)**: start with `TaxProfiler-nf` or `KrakenMetaReads-nf` (Kraken2/Bracken ecosystem, stats-friendly outputs).
- **If you need short-read dual-assembler comparisons and/or long-read dual-track complementarity (feature-based vs similarity-based)**: use `rvdb-viral-metagenome-nf`.
- **If you want to maximize discovery of novel/distant viruses with confidence tiers**: use `MLMVD-nf` (VS2 + DVF + viralFlye multi-evidence stratification).
- **If you already have contigs and want fast similarity screening**: use `Sourmash`.
- **If you prefer non-Kraken2 k-mer profiling/classification**: consider `GOTTCHA2` (signature profiling) or `CLARK` (discriminative k-mers + viruses-only DB building).

