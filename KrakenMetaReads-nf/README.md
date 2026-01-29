# KrakenMetaReads-nf

A comprehensive Nextflow workflow for viral metagenomic classification and abundance analysis using nf-core/taxprofiler. This pipeline supports both short-read (Illumina) and long-read (Nanopore/PacBio) sequencing data, providing automated taxonomic classification with Kraken2 and abundance quantification using RPM (Reads Per Million) and RPKM (Reads Per Kilobase Million) metrics.

[![Apptainer](https://img.shields.io/badge/Apptainer-%3D%201.3-42b983)](https://apptainer.org/)
[![Nextflow](https://img.shields.io/badge/nextflow-%3D%2024.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%3D%203.9-blue.svg)](https://www.python.org/)
---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output Structure](#output-structure)
- [Abundance Calculation](#abundance-calculation)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## Features

- **Dual Sequencing Support**: Optimized workflows for both short-read (Illumina) and long-read (Nanopore/PacBio) data
- **Assembly-Based Analysis**: 
  - Short-read: MEGAHIT and SPAdes for metagenomic assembly
  - Long-read: Flye and ViralFlye for long-read assembly (with circular/linear distinction)
- **Automated Classification**: Uses Kraken2 for fast and accurate taxonomic classification on assembled contigs
- **Abundance Quantification**: 
  - Short-read: Kraken2 + Bracken statistical correction (if available)
  - Long-read: Direct extraction from Kraken2 reports (Bracken not needed)
- **Multiple Metrics**: Calculates both RPM and RPKM for comprehensive abundance analysis
- **Batch Processing**: Automated batch processing of multiple samples
- **Containerized**: Uses Apptainer/Singularity for reproducible analysis
- **Quality Control**: Integrated QC steps (Fastp for short-read, Nanoplot for long-read)
- **Standardized Output**: Generates BIOM format files for downstream analysis

## Requirements

### Software Dependencies

- **Nextflow** (>= 22.10.0)
- **Java** (>= 17)
- **Apptainer/Singularity** (>= 1.1.0)
- **Python 3** (>= 3.7) with pandas
- **nf-core/taxprofiler** (>= 1.2.0)

### System Requirements

- **CPU**: 32 cores recommended
- **Memory**: 256 GB RAM recommended
- **Storage**: Sufficient space for databases and results (varies by dataset size)
- **SLURM**: For cluster execution (scripts include SLURM directives)

### Database Requirements

- Kraken2 viral reference database (specified in `databases.csv`)
- Bracken database (for short-read analysis, built from Kraken2 database)

## Installation

### 1. Clone the Repository

```bash
git clone <repository-url>
cd KrakenMetaReads-nf
```

### 2. Set Up Nextflow Environment

Create a conda environment with required dependencies:

```bash
conda create -n nextflow_env python=3.9 java=17
conda activate nextflow_env
conda install -c bioconda nextflow apptainer
```

### 3. Prepare Databases

Ensure your Kraken2 viral database is built and accessible. Update the database path in `databases.csv`:

```csv
tool,db_name,db_params,db_path
kraken2,Viral_ref,"",/path/to/kraken2_Viral_ref
bracken,Viral_ref,";-r 150",/path/to/kraken2_Viral_ref
```

### 4. Install Python Dependencies

```bash
pip install pandas
```

## Quick Start

### For Short-Read Data (Illumina)

1. **Prepare samplesheet**: Edit `samplesheet_short.csv` with your sample information
2. **Submit job**: 
   ```bash
   sbatch submit_short.sh
   ```

### For Long-Read Data (Nanopore/PacBio)

1. **Prepare samplesheet**: Edit `samplesheet_long.csv` with your sample information
2. **Submit job**:
   ```bash
   sbatch submit_long.sh
   ```

## Configuration

### Samplesheet Format

#### Short-Read Samplesheet (`samplesheet_short.csv`)

```csv
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
sample1,run1,ILLUMINA,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,
sample2,run2,ILLUMINA,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,
```

#### Long-Read Samplesheet (`samplesheet_long.csv`)

```csv
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
sample1,run1,OXFORD_NANOPORE,/path/to/sample1.fastq.gz,,
sample2,run2,OXFORD_NANOPORE,/path/to/sample2.fastq.gz,,
```

### Configuration Files

#### `nextflow_short.config` (Short-Read Configuration)

Key parameters:
- `input`: Path to short-read samplesheet
- `outdir`: Output directory (default: `results_viral_short`)
- `databases`: Path to databases CSV file
- `run_kraken2`: Enable Kraken2 classification (default: `true`)
- `run_bracken`: Enable Bracken abundance estimation (default: `true`)
- `bracken_precision`: Taxonomic level for Bracken (default: `'S'` for species)
- `bracken_readlen`: Read length for Bracken (default: `150`)

#### `nextflow_long.config` (Long-Read Configuration)

Key parameters:
- `input`: Path to long-read samplesheet
- `outdir`: Output directory (default: `results_viral_long`)
- `run_kraken2`: Enable Kraken2 classification (default: `true`)
- `run_bracken`: Disabled for long-read data (default: `false`)
- `perform_longread_qc`: Enable long-read QC (default: `true`)

### Database Configuration (`databases.csv`)

```csv
tool,db_name,db_params,db_path
kraken2,Viral_ref,"",/path/to/kraken2_Viral_ref
bracken,Viral_ref,";-r 150",/path/to/kraken2_Viral_ref
```

## Usage

### Automated Workflow (Recommended)

The submission scripts (`submit_short.sh` and `submit_long.sh`) automatically:
1. Run the nf-core/taxprofiler workflow
2. Calculate abundance metrics (RPM & RPKM) for all samples
3. Generate summary tables

### Manual Execution

#### Step 1: Run TaxProfiler

```bash
# Short-read
nextflow run nf-core/taxprofiler \
  -r 1.2.0 \
  -profile apptainer \
  -c nextflow_short.config \
  -resume

# Long-read
nextflow run nf-core/taxprofiler \
  -r 1.2.0 \
  -profile apptainer \
  -c nextflow_long.config \
  -resume
```

#### Step 2: Calculate Abundance

```bash
# Short-read (with Bracken)
bash batch_calculate_abundance_en.sh results_viral_short

# Long-read (Kraken2 only)
bash batch_calculate_abundance_longread_en.sh results_viral_long
```

### Single Sample Abundance Calculation

#### Short-Read Data

**For MEGAHIT assemblies:**
```bash
# Note: Short-read data may use Bracken if available, or calculate directly from Kraken2
# For MEGAHIT assembly results
python3 calculate_abundance_en.py \
  --bracken results_viral_short/bracken/sample1_bracken.tsv \
  --kraken results_viral_short/kraken2_megahit/sample1.report \
  --output results_viral_short/abundance_megahit/sample1_abundance.tsv

# If Bracken is not available, use long-read script method:
python3 calculate_abundance_longread_en.py \
  --kraken results_viral_short/kraken2_megahit/sample1.report \
  --output results_viral_short/abundance_megahit/sample1_abundance.tsv \
  --level S
```

**For SPAdes assemblies:**
```bash
python3 calculate_abundance_en.py \
  --bracken results_viral_short/bracken/sample1_bracken.tsv \
  --kraken results_viral_short/kraken2_spades/sample1.report \
  --output results_viral_short/abundance_spades/sample1_abundance.tsv
```

#### Long-Read Data

**For Flye assemblies:**
```bash
# Note: Long-read data does NOT use Bracken
python3 calculate_abundance_longread_en.py \
  --kraken results_viral_long/kraken2_flye/sample1.report \
  --output results_viral_long/abundance_flye/sample1_abundance.tsv \
  --level S
```

**For ViralFlye circular contigs:**
```bash
python3 calculate_abundance_longread_en.py \
  --kraken results_viral_long/kraken2_viralflye_circular/sample1.report \
  --output results_viral_long/abundance_viralflye_circular/sample1_abundance.tsv \
  --level S
```

**For ViralFlye linear contigs:**
```bash
python3 calculate_abundance_longread_en.py \
  --kraken results_viral_long/kraken2_viralflye_linear/sample1.report \
  --output results_viral_long/abundance_viralflye_linear/sample1_abundance.tsv \
  --level S
```

## Output Structure

### Short-Read Data Structure

```
results_viral_short/
├── fastp/                            # Quality control and preprocessing (Fastp)
├── kraken2_megahit/                  # Kraken2 classification on MEGAHIT assemblies
│   └── [sample reports and classification files]
├── kraken2_spades/                   # Kraken2 classification on SPAdes assemblies
│   └── [sample reports and classification files]
├── abundance_megahit/                # Abundance metrics for MEGAHIT assemblies
│   ├── sample1_abundance.tsv
│   ├── sample2_abundance.tsv
│   ├── all_samples_abundance_summary.tsv
│   └── top_viruses_summary.tsv
├── abundance_spades/                 # Abundance metrics for SPAdes assemblies
│   ├── sample1_abundance.tsv
│   ├── sample2_abundance.tsv
│   ├── all_samples_abundance_summary.tsv
│   └── top_viruses_summary.tsv
└── merged_reports/                   # Merged classification reports
```

**Note**: 
- Short-read data uses multiple assembly tools (MEGAHIT and SPAdes)
- Kraken2 classification is performed on each assembly separately
- Abundance calculations are generated for each assembly method
- The `merged_reports/` directory contains consolidated results

### Long-Read Data Structure

```
results_viral_long/
├── flye_assembly/                    # Flye assembly results
├── viralflye/                        # ViralFlye assembly results
├── kraken2_flye/                     # Kraken2 classification on Flye assemblies
│   └── [sample reports and classification files]
├── kraken2_viralflye_circular/       # Kraken2 classification on ViralFlye circular contigs
│   └── [sample reports and classification files]
├── kraken2_viralflye_linear/         # Kraken2 classification on ViralFlye linear contigs
│   └── [sample reports and classification files]
├── abundance_flye/                   # Abundance metrics for Flye assemblies
│   ├── sample1_abundance.tsv
│   ├── sample2_abundance.tsv
│   ├── all_samples_abundance_summary.tsv
│   └── top_viruses_summary.tsv
├── abundance_viralflye_circular/     # Abundance metrics for ViralFlye circular contigs
│   ├── sample1_abundance.tsv
│   ├── sample2_abundance.tsv
│   ├── all_samples_abundance_summary.tsv
│   └── top_viruses_summary.tsv
└── abundance_viralflye_linear/       # Abundance metrics for ViralFlye linear contigs
    ├── sample1_abundance.tsv
    ├── sample2_abundance.tsv
    ├── all_samples_abundance_summary.tsv
    └── top_viruses_summary.tsv
```

**Note**: 
- Long-read data uses Flye and ViralFlye for assembly
- ViralFlye distinguishes between circular and linear viral contigs
- Kraken2 classification is performed on each assembly type separately
- Abundance calculations are generated for each assembly method and contig type
- **Bracken directory is NOT present for long-read data** (Bracken is designed for short reads only)

### Abundance Output Format

Each sample abundance file (`*_abundance.tsv`) contains:

| Column | Description |
|--------|-------------|
| Species | Taxonomic species name |
| Taxonomy_ID | NCBI taxonomy ID |
| Assigned_Reads | Number of reads assigned to this species |
| Fraction | Fraction of total reads (0-1) |
| RPM | Reads Per Million |
| Genome_Length_bp | Genome length in base pairs (if available) |
| RPKM | Reads Per Kilobase Million (if genome length available) |

## Abundance Calculation

### Metrics Explained

#### RPM (Reads Per Million)
```
RPM = (assigned_reads / total_reads) × 1,000,000
```
- Normalizes read counts by total sequencing depth
- Useful for comparing abundance across samples with different sequencing depths

#### RPKM (Reads Per Kilobase Million)
```
RPKM = assigned_reads / (genome_length_kb × total_reads_million)
```
- Normalizes by both sequencing depth and genome length
- Accounts for genome size differences between viruses
- More accurate for comparing abundance across different viral species

### Short-Read vs Long-Read

**Short-Read (Illumina)**:
- **Assembly tools**: Uses MEGAHIT and SPAdes for metagenomic assembly
- Uses Kraken2 for classification on assembled contigs
- Applies Bracken statistical correction to improve abundance estimates (if available)
- Bracken is designed for short reads (50-300 bp)
- **Output structure**: 
  - `kraken2_megahit/` and `kraken2_spades/` for classification results
  - `abundance_megahit/` and `abundance_spades/` for abundance metrics
  - `fastp/` for quality control and preprocessing
  - `merged_reports/` for consolidated results
- **Abundance calculation**: Uses both Kraken2 reports and Bracken output files (or Kraken2 only if Bracken unavailable)
- **QC tools**: Fastp for quality control

**Long-Read (Nanopore/PacBio)**:
- **Assembly tools**: Uses Flye and ViralFlye for long-read assembly
- ViralFlye distinguishes between circular and linear viral contigs
- Uses Kraken2 for classification on assembled contigs
- No Bracken correction needed (long reads contain more information)
- Direct extraction from Kraken2 reports is sufficient
- **Output structure**: 
  - `kraken2_flye/`, `kraken2_viralflye_circular/`, `kraken2_viralflye_linear/` for classification results
  - `abundance_flye/`, `abundance_viralflye_circular/`, `abundance_viralflye_linear/` for abundance metrics
  - `flye_assembly/` and `viralflye/` for assembly results
- **Abundance calculation**: Uses only Kraken2 report files
- **QC tools**: Nanoplot for long-read quality control

### Custom Genome Length Database

You can provide a custom genome length database:

```bash
python3 calculate_abundance_en.py \
  --bracken sample1_bracken.tsv \
  --kraken sample1.report \
  --output sample1_abundance.tsv \
  --genome-db custom_genome_lengths.tsv
```

Format of `custom_genome_lengths.tsv`:
```
Species Name<TAB>Genome_Length_bp
Severe acute respiratory syndrome coronavirus 2<TAB>29903
Influenza A virus<TAB>13588
```

## Troubleshooting

### Common Issues

#### 1. Bracken Results Not Found (Short-Read)

**Symptom**: Warning message about missing Bracken results

**Solution**: 
- Check that `run_bracken = true` in configuration
- Verify Bracken database is properly configured in `databases.csv`
- The script will automatically fall back to Kraken2-only method if Bracken fails

#### 2. Container Mount Errors

**Symptom**: Apptainer mount point errors

**Solution**: The configuration includes `--no-mount /lscratch` option. If you encounter other mount issues, update `apptainer.runOptions` in the config file.

#### 3. Memory Issues

**Symptom**: Out of memory errors

**Solution**: 
- Increase `max_memory` in configuration file
- Reduce `max_cpus` to limit parallel processes
- Process samples in smaller batches

#### 4. Database Path Errors

**Symptom**: Database not found errors

**Solution**:
- Verify database paths in `databases.csv` are absolute paths
- Ensure database directories are accessible
- Check that Kraken2 database is properly built

#### 5. Python Dependencies Missing

**Symptom**: `ModuleNotFoundError: No module named 'pandas'`

**Solution**:
```bash
pip install pandas
# or
conda install pandas
```

### Getting Help

1. Check Nextflow logs: `work/` directory contains detailed execution logs
2. Review SLURM output files: `*_%j.out` and `*_%j.err`
3. Check QC results: `fastp/` directory for short-read data, assembly directories for long-read data

## Citation

If you use this workflow in your research, please cite:

1. **nf-core/taxprofiler**: 
   - Ewels, P. A., et al. (2023). The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology*, 41(2), 178-183.

2. **Kraken2**:
   - Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257.

3. **Bracken**:
   - Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.

4. **Nextflow**:
   - Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4), 316-319.

## License

MIT License.

## Contact

- Email: sihua.peng@uga.edu, Workflow code programmer  
- Email: justin.bahl@uga.edu, Project supervisor  
- Suggestion: [Click here!](https://github.com/pengsihua2023/rvdb-viral-metagenome-nf/issues/new)

---







