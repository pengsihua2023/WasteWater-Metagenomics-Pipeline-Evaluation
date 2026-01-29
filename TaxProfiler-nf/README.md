
# TaxProfiler - Viral Metagenomics Analysis Workflow

MetaTaxProfiler is a viral metagenomics analysis tool built upon the TaxProfiler-nf workflow. It provides an automated, containerized pipeline that converts raw metagenomic sequencing data into standardized viral abundance outputs. Supporting both Illumina and Nanopore/PacBio platforms, it performs quality control, taxonomic classification (Kraken2), optional statistical correction (Bracken), and automatic RPM/RPKM calculation—all with a single command.

With its containerized design, MetaTaxProfiler ensures high reproducibility, delivers publication-ready results, and intelligently adapts to available computational resources. It is an ideal solution for researchers seeking efficient and standardized viral metagenomic analysis.

[![Apptainer](https://img.shields.io/badge/Apptainer-%E2%89%A51.3-42b983)](https://apptainer.org/)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue.svg)](https://www.python.org/)

---

## Quick Start

### Short-read Data (Illumina)

```bash
sbatch submit_short.sh
```

**Automatic pipeline**: QC → Kraken2 classification → (Bracken correction) → Auto RPM/RPKM calculation → Report generation

**Output**: `results_viral_short/abundance/`

### Long-read Data (Nanopore/PacBio)

```bash
sbatch submit_long.sh
```

**Automatic pipeline**: QC → Kraken2 classification → Auto RPM/RPKM calculation → Report generation

**Output**: `results_viral_long/abundance/`

---

## Project Files

### Core Runtime Files (10)

```
├── Configuration Files (5)
│   ├── databases.csv              # Database configuration
│   ├── samplesheet_short.csv      # Short-read sample list
│   ├── samplesheet_long.csv       # Long-read sample list
│   ├── nextflow_short.config      # Short-read configuration
│   └── nextflow_long.config       # Long-read configuration
│
├── Submission Scripts (2)
│   ├── submit_short.sh            # Short-read analysis entry point
│   └── submit_long.sh             # Long-read analysis entry point
│
└── Abundance Calculation (4)
    ├── batch_calculate_abundance_en.sh              # Short-read abundance (smart dispatcher)
    ├── batch_calculate_abundance_longread_en.sh     # Universal abundance (Kraken2)
    ├── calculate_abundance_en.py                    # Bracken → RPM/RPKM
    └── calculate_abundance_longread_en.py           # Kraken2 → RPM/RPKM
```

---

## Sample Sheet Format

### Short-read (samplesheet_short.csv)

```csv
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
sample1,run1,ILLUMINA,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,
```

### Long-read (samplesheet_long.csv)

```csv
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
sample1,run1,OXFORD_NANOPORE,/path/to/reads.fastq.gz,,
```

---

## Database Configuration (databases.csv)

```csv
tool,db_name,db_params,db_path
kraken2,Viral_ref,"",/path/to/kraken2_viral_database
bracken,Viral_ref,";-r 150",/path/to/kraken2_viral_database
```

**Requirements**:
- **Short-read**: Kraken2 database + Bracken database (recommended for improved accuracy)
  - Bracken database files should be in the same directory as Kraken2 database
  - Required files: `database150mers.kmer_distrib`, `database150mers.kraken`
- **Long-read**: Kraken2 database only (Bracken not needed)

---

## Output Results

### Directory Structure

```
results_viral_short/  (or results_viral_long/)
├── kraken2/                   # Kraken2 classification results
├── bracken/                   # Bracken abundance estimation (short-read only, if available)
├── abundance/                 # ⭐ Abundance calculation results (auto-generated)
│   ├── sample_abundance.tsv              # Individual sample detailed table
│   ├── all_samples_abundance_summary.tsv # All samples summary
│   └── top_viruses_summary.tsv           # High abundance viruses (RPM≥10)
├── fastqc/                    # QC reports
└── multiqc/                   # Comprehensive report
    └── multiqc_report.html    # HTML report (recommended)
```

### Abundance Table Format

| Column | Description |
|--------|-------------|
| Species | Viral species name |
| Taxonomy_ID | NCBI taxonomy ID |
| Assigned_Reads | Number of assigned reads |
| Fraction | Relative abundance (0-1) |
| **RPM** | Normalized abundance per million reads ⭐ |
| Genome_Length_bp | Genome length (for known viruses) |
| **RPKM** | Normalized abundance accounting for genome length ⭐ |

---

## Abundance Metrics

### RPM (Reads Per Million) - Relative Abundance
- **Formula**: `RPM = (viral reads / total reads) × 1,000,000`
- **Purpose**: Compare the same virus across different samples
- **Advantage**: Simple and intuitive, no genome length required
- **Type**: Relative abundance metric (normalized to total reads)

### RPKM (Reads Per Kilobase Million) - Relative Abundance (Genome-length Normalized)
- **Formula**: `RPKM = (viral reads) / (genome length kb × total reads million)`
- **Purpose**: Compare viral loads between different viruses (accounts for genome size)
- **Type**: Relative abundance metric (normalized to total reads and genome length)
- **Limitation**: Requires known genome length (unknown viruses show NA)

**Note**: Both RPM and RPKM are **relative abundance** metrics, not absolute abundance. For absolute abundance, additional calibration methods (e.g., spike-in controls) are required.

---

## Workflow

### Short-read Pipeline

```
1. Submit job: sbatch submit_short.sh
2. TaxProfiler runs:
   ├─ FastQC/fastp QC
   ├─ Kraken2 classification
   └─ Bracken abundance estimation (if database available)
3. Auto abundance calculation:
   ├─ With Bracken → Bracken-based calculation (more accurate)
   └─ Without Bracken → Kraken2-based calculation (still valid)
4. Generate results: results_viral_short/abundance/
```

### Long-read Pipeline

```
1. Submit job: sbatch submit_long.sh
2. TaxProfiler runs:
   ├─ NanoPlot/Porechop QC
   └─ Kraken2 classification (higher accuracy for long reads)
3. Auto abundance calculation:
   └─ Kraken2 report-based (no Bracken needed)
4. Generate results: results_viral_long/abundance/
```

---

## Short-read vs Long-read Comparison

| Feature | Short-read (Illumina) | Long-read (Nanopore/PacBio) |
|---------|---------------------|---------------------------|
| **Read length** | 100-300 bp | 1,000-100,000+ bp |
| **Kraken2 classification** | May reach genus/family level | Usually reaches species level ✅ |
| **Uses Bracken** | Recommended (improves accuracy) | Not needed |
| **Database requirements** | Kraken2 (+optional Bracken) | Kraken2 only |
| **Abundance source** | Bracken or Kraken2 | Kraken2 |
| **Classification accuracy** | 85-98% | >95% |
| **Setup time** | 0-6 hours (if building Bracken) | Ready immediately |

---

## Configuration

### Key Parameters (nextflow_short.config)

```groovy
params {
    input       = 'samplesheet_short.csv'
    outdir      = 'results_viral_short'
    databases   = 'databases.csv'
    
    run_kraken2 = true
    run_bracken = true   // Recommended if Bracken database available
    
    max_cpus    = 32
    max_memory  = '256.GB'
}
```

### Key Parameters (nextflow_long.config)

```groovy
params {
    input       = 'samplesheet_long.csv'
    outdir      = 'results_viral_long'
    databases   = 'databases.csv'
    
    run_kraken2 = true
    run_bracken = false  // Long reads don't need Bracken
    
    max_cpus    = 32
    max_memory  = '256.GB'
}
```

---

## Usage Tips

### View Abundance Results

```bash
# View individual sample details
head -30 results_viral_short/abundance/sample_abundance.tsv

# View all samples summary
head -50 results_viral_short/abundance/all_samples_abundance_summary.tsv

# View high abundance viruses (RPM ≥ 10)
cat results_viral_short/abundance/top_viruses_summary.tsv

# Count detected viruses
tail -n +2 results_viral_short/abundance/all_samples_abundance_summary.tsv | wc -l
```

### Extract Specific Viruses

```bash
# Search for specific virus
grep -i "coronavirus" results_viral_short/abundance/all_samples_abundance_summary.tsv

# Extract TOP 50 viruses
head -51 results_viral_short/abundance/all_samples_abundance_summary.tsv
```

---

## Environment Requirements

- **Java**: 17+ (OpenJDK 17.0.3)
- **Nextflow**: 25.04+ (current 25.04.7)
- **Apptainer**: 1.3+ (current 1.3.6)
  - **Configuration**: `--no-mount /lscratch` (handles missing mount points)
- **Conda**: nextflow_env environment
- **Python**: 3.9+ (for abundance calculation)
- **Python packages**: pandas (for abundance calculation)

---

## Data Interpretation

### Typical Abundance Ranges

**Human viral infection samples**:
- High load: RPM > 10,000
- Medium load: RPM 100-10,000
- Low load: RPM < 100

**Environmental/metagenomic samples**:
- Major viruses: RPM 100-1,000
- Minor viruses: RPM 10-100
- Rare viruses: RPM < 10

### About RPKM Showing NA

**Normal occurrence**:
- Environmental phages are usually not in the genome length database
- Only RPM available for these viruses
- **RPM is sufficient** for most analyses

---

## Troubleshooting

### Issue: No abundance directory generated

**Solution**: Check if TaxProfiler completed successfully

```bash
# View logs
tail -100 Viral_Short_*.out

# Check if Kraken2 results exist
ls -lh results_viral_short/kraken2/
```

### Issue: top_viruses_summary.tsv is empty

**Cause**: All viruses have RPM < 10

**Solution**: Manually generate list with lower threshold

```bash
# Generate list with RPM >= 5
awk -F'\t' 'NR==1 || $6 >= 5' results_viral_short/abundance/all_samples_abundance_summary.tsv | \
  sort -t$'\t' -k6,6nr > results_viral_short/abundance/top_viruses_rpm5.tsv
```

### Issue: Bracken didn't run

**This doesn't affect abundance calculation!** Scripts will automatically use Kraken2 results.

To use Bracken (optional, more accurate):
1. Check if database has `database150mers.kmer_distrib` file
2. Ensure `run_bracken = true` in config file
3. Re-run analysis

---

## Related Documentation

- **ABUNDANCE_USAGE.md** - Detailed abundance calculation guide
- **LONGREAD_GUIDE.md** - Long-read data specific guide
- **SHORTREAD_ABUNDANCE_GUIDE.md** - Short-read abundance guide

---

## Technical Support

- **Database path**: `/scratch/sp96859/Meta-genome-data-analysis/`
- **Working environment**: `nextflow_env` (conda)
- **Compute partition**: `bahl_p`
- **nf-core/taxprofiler**: https://nf-co.re/taxprofiler

---

### Contact

- Email: sihua.peng@uga.edu, Workflow code programmer  
- Email: justin.bahl@uga.edu, Project supervisor  
- Suggestion: [Click here!](https://github.com/pengsihua2023/rvdb-viral-metagenome-nf/issues/new)
