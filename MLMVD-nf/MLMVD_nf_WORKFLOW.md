# Metagenome Viral Classification Workflow (Nextflow DSL2)

This repository provides a **metagenome viral classification pipeline** implemented in **Nextflow DSL2** (`metagenome_assembly_classification_workflow.nf`).

It supports two execution modes:
- **Short-read mode (Illumina paired-end)**: QC → assembly (MEGAHIT + metaSPAdes) → viral identification (VirSorter2 + DeepVirFinder) → merge reports → assembler comparison → abundance (RPM/RPKM)
- **Long-read mode (Nanopore/PacBio)**: assembly (metaFlye) → viral identification (VirSorter2 + DeepVirFinder + optional viralFlye) → optional three-tool comparison → abundance (RPM/RPKM)

---

## Inputs

### Short-read samplesheet (paired-end)
CSV with header:
- `sample`: sample ID
- `fastq_1`: R1 FASTQ/FASTQ.GZ path
- `fastq_2`: R2 FASTQ/FASTQ.GZ path

Example:

```csv
sample,fastq_1,fastq_2
S1,/path/to/S1_1.fastq.gz,/path/to/S1_2.fastq.gz
```

### Long-read samplesheet (single-end)
CSV with header:
- `sample`: sample ID
- `fastq_long`: long-read FASTQ/FASTQ.GZ path

Example:

```csv
sample,fastq_long
S1,/path/to/S1_nanopore.fastq.gz
```

---

## Software used (by step)

### Core workflow engine
- **Nextflow** (DSL2)

### Short-read mode
- **fastp**: read QC / filtering (Conda: `bioconda::fastp=0.23.4`)
- **MEGAHIT**: metagenome assembly (container: `docker://quay.io/biocontainers/megahit:1.2.9--h2e03b76_1`)
- **metaSPAdes (SPAdes)**: metagenome assembly (container: `docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1`)
- **VirSorter2**: viral sequence identification (Conda env path configured in the pipeline, e.g. `/home/sp96859/.conda/envs/VirSorter2_env`)
- **DeepVirFinder**: viral prediction (runs `dvf.py` from `--deepvirfinder_dir`, activates an existing DVF conda environment)
- **Python + pandas/numpy**: report merging and comparisons
- **bowtie2 + samtools**: read mapping to viral contigs and abundance calculation (RPM/RPKM)

### Long-read mode
- **(Optional placeholder) Long-read QC**: standardizes output naming
- **metaFlye (Flye --meta)**: long-read metagenome assembly (Conda: `bioconda::flye=2.9`)
- **VirSorter2**: viral sequence identification (same as above)
- **DeepVirFinder**: viral prediction (same as above)
- **viralFlye** (optional): viral identification/refinement with Pfam validation (uses a pre-existing conda env path `--viralflye_env` and Pfam database `--pfam_db`)
- **minimap2 + samtools**: long-read mapping to viral contigs and abundance calculation (RPM/RPKM)

### Container runtime
- **Apptainer** is enabled via `metagenome_assembly_classification.config` and is used when processes declare a `container` image (e.g., MEGAHIT and SPAdes).

---

## Workflow overview

### A) Short-read mode (Illumina paired-end)

1. **FASTP (QC)**
   - Input: `fastq_1`, `fastq_2`
   - Output: cleaned reads (`*_clean_R1.fastq.gz`, `*_clean_R2.fastq.gz`) + QC reports (`*.html`, `*.json`)

2. **MEGAHIT_ASSEMBLY**
   - Input: cleaned reads
   - Output: assembled contigs (`*_megahit_contigs.fa`)

3. **SPADES_ASSEMBLY (metaSPAdes)**
   - Input: cleaned reads
   - Output: assembled contigs (`*_spades_contigs.fa`)

4. **VIRSORTER2 (per assembler)**
   - Runs on both MEGAHIT and SPAdes contigs
   - Output (per assembler): viral score table (`*_vs2_final-viral-score.tsv`) and optional viral contigs FASTA (`*_vs2_final-viral-combined.fa`)

5. **DEEPVIRFINDER (per assembler)**
   - Runs on both MEGAHIT and SPAdes contigs
   - Output: `*_dvf_output.txt` (scores and p-values)

6. **MERGE_VIRAL_REPORTS (per assembler)**
   - Integrates VirSorter2 + DeepVirFinder results
   - Output:
     - `*_viral_merged_report.txt` (human-readable summary)
     - `*_viral_merged_report.csv` (machine-readable table)
     - `*_viral_consensus.txt` (IDs identified by both tools for that assembler; “high confidence” within that assembler)

7. **COMPARE_ASSEMBLERS (MEGAHIT vs SPAdes)**
   - Compares the two merged reports (MEGAHIT vs SPAdes)
   - Output:
     - `*_assembler_comparison.txt`
     - `*_assembler_comparison.csv`
     - `*_consensus_viral_sequences.txt` (final high-confidence consensus across both assemblers)

8. **CALCULATE_ABUNDANCE (per assembler)**
   - Maps reads to viral contigs and computes **RPM** and **RPKM**
   - Output (per assembler):
     - `*_{assembler}_abundance.csv`
     - `*_{assembler}_abundance_summary.txt`

### B) Long-read mode (Nanopore / PacBio)

1. **METAFLYE_ASSEMBLY (Flye --meta)**
   - Input: long reads
   - Output:
     - `*_metaflye_contigs.fa`
     - `*_flye_output/` (complete Flye output directory, published for downstream use)

2. **Viral identification on metaFlye output**
   - **VirSorter2** → `*_metaflye_vs2_final-viral-score.tsv`, `*_metaflye_vs2_final-viral-combined.fa` (optional)
   - **DeepVirFinder** → `*_metaflye_dvf_output.txt`
   - **viralFlye** (optional) → `*_viralflye_contigs.fa`, `*_viralflye_summary.csv` (+ full output dir)

3. **COMPARE_THREE_VIRAL_TOOLS** (optional; when viralFlye + VirSorter2 + DeepVirFinder are all enabled)
   - Produces a three-tool comparison and a prioritized high-confidence list:
     - `*_three_tools_comparison.txt`
     - `*_three_tools_comparison.csv`
     - `*_high_confidence_viruses.txt`
     - `*_high_confidence_viruses.fa` (optional)

4. **CALCULATE_ABUNDANCE_LONGREAD**
   - Maps long reads to viral contigs and computes **RPM** and **RPKM**
   - Output (per source, e.g. `metaflye`, `viralflye`):
     - `*_{source}_abundance.csv`
     - `*_{source}_abundance_summary.txt`

---

## Output directory structure (by default under `--outdir`)

### Short-read mode outputs
- `fastp/`
  - `*_fastp.html`
  - `*_fastp.json`
- `clean_reads/` (if `save_clean_reads=true`)
  - `*_clean_R1.fastq.gz`
  - `*_clean_R2.fastq.gz`
- `assembly_megahit/`
  - `*_megahit_contigs.fa`
- `assembly_spades/`
  - `*_spades_contigs.fa`
- `virsorter2_megahit/`
  - `*_megahit_vs2_final-viral-score.tsv`
  - `*_megahit_vs2_final-viral-combined.fa` (optional)
- `virsorter2_spades/`
  - `*_spades_vs2_final-viral-score.tsv`
  - `*_spades_vs2_final-viral-combined.fa` (optional)
- `deepvirfinder_megahit/`
  - `*_megahit_dvf_output.txt`
- `deepvirfinder_spades/`
  - `*_spades_dvf_output.txt`
- `merged_viral_reports_megahit/`
  - `*_megahit_viral_merged_report.txt`
  - `*_megahit_viral_merged_report.csv`
  - `*_megahit_viral_consensus.txt`
- `merged_viral_reports_spades/`
  - `*_spades_viral_merged_report.txt`
  - `*_spades_viral_merged_report.csv`
  - `*_spades_viral_consensus.txt`
- `assembler_comparison/`
  - `*_assembler_comparison.txt`
  - `*_assembler_comparison.csv`
  - `*_consensus_viral_sequences.txt`  **(final consensus)**
- `abundance/megahit/`
  - `*_megahit_abundance.csv`
  - `*_megahit_abundance_summary.txt`
- `abundance/spades/`
  - `*_spades_abundance.csv`
  - `*_spades_abundance_summary.txt`

### Long-read mode outputs
- `assembly_metaflye/`
  - `*_metaflye_contigs.fa`
- `metaflye_full_output/`
  - `*_flye_output/` (complete Flye output)
- `virsorter2_metaflye/`
  - `*_metaflye_vs2_final-viral-score.tsv`
  - `*_metaflye_vs2_final-viral-combined.fa` (optional)
- `deepvirfinder_metaflye/`
  - `*_metaflye_dvf_output.txt`
- `viralflye_results/` (if enabled)
  - `*_viralflye_contigs.fa`
  - `*_viralflye_summary.csv`
- `viralflye_full_output/` (if enabled)
  - `*_viralflye_output/`
- `three_tools_comparison/` (only when all three tools are enabled)
  - `*_three_tools_comparison.txt`
  - `*_three_tools_comparison.csv`
  - `*_high_confidence_viruses.txt`
  - `*_high_confidence_viruses.fa` (optional)
- `merged_viral_reports_metaflye/` (used when three-tool comparison is not active)
  - `*_metaflye_viral_merged_report.txt`
  - `*_metaflye_viral_merged_report.csv`
  - `*_metaflye_viral_consensus.txt`
- `abundance/metaflye/`
  - `*_metaflye_abundance.csv`
  - `*_metaflye_abundance_summary.txt`
- `abundance/viralflye/` (if enabled)
  - `*_viralflye_abundance.csv`
  - `*_viralflye_abundance_summary.txt`

---

## Key “final” result files to focus on

### Short-read mode
- `assembler_comparison/*_consensus_viral_sequences.txt`  
  Final high-confidence viral ID list after combining **(VirSorter2 + DeepVirFinder)** and comparing **MEGAHIT vs SPAdes**.
- `abundance/*/*_abundance.csv`  
  Abundance table for each contig (RPM and RPKM).
- `abundance/*/*_abundance_summary.txt`  
  Top abundant contigs summary for quick review.

### Long-read mode
- `three_tools_comparison/*_high_confidence_viruses.txt` (if enabled)  
  Prioritized high-confidence virus list across **VirSorter2 + DeepVirFinder + viralFlye**.
- `abundance/*/*_abundance.csv` and `*_abundance_summary.txt`  
  Abundance tables and summaries.

---

## Notes
- Apptainer is enabled in the config and will be used for any process that declares a `container` image.
- VirSorter2 and DeepVirFinder steps assume **pre-existing environments** on the system (paths are referenced inside the pipeline and/or provided via parameters).
- Abundance metrics:
  - **RPM** = reads per million total reads
  - **RPKM** = reads per kilobase per million (normalized by contig length and total reads)

