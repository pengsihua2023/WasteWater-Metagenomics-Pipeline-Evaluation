# CLARK Software Installation and User Guide

## Table of Contents
- [CLARK Introduction](#clark-introduction)
- [System Requirements](#system-requirements)
- [Installation Steps](#installation-steps)
- [Database Preparation](#database-preparation)
  - [NCBI RefSeq Database Download](#ncbi-refseq-database-download)
- [Usage](#usage)
- [Common Command Examples](#common-command-examples)
  - [Virus Detection and Classification Examples](#example-7-virus-detection-and-classification-specifically-for-virus-analysis)
- [Result Interpretation](#result-interpretation)
  - [Abundance Calculation](#method-2-manual-abundance-calculation-python-script)
- [Complete Virus Detection and Classification Pipeline](#complete-virus-detection-and-classification-pipeline)
- [Common Issues](#common-issues)

## CLARK Introduction

CLARK (Classification of metagenomic and genomic sequences) is an efficient metagenomic sequence classification tool that can quickly and accurately classify metagenomic data. CLARK supports multiple classification modes, including species classification, genus classification, etc., and is suitable for large-scale metagenomic data analysis.


[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue.svg)](https://www.python.org/)
[![Apptainer](https://img.shields.io/badge/Apptainer-%E2%89%A51.3-42b983)](https://apptainer.org/)

---

## Quick Start: Virus Detection and Classification

If you need to quickly perform virus detection and classification, you can follow these steps:

### 1. Install CLARK

```bash
# Clone repository from GitHub
git clone https://github.com/rouni001/CLARK.git
cd CLARK

# Run installation script
./install.sh

# Or compile manually
make
```

### 2. Download and Build Virus Database

```bash
# Set database directory
export CLARK_DB_DIR=/path/to/CLARK_DB
mkdir -p $CLARK_DB_DIR

# Download and build virus database (this may take some time)
./set_targets.sh $CLARK_DB_DIR viruses
```

### 3. Run Virus Detection

```bash
# Single-end sequencing data (using classification script, recommended)
./classify_metagenome.sh -O sample.fastq -R virus_results -m 0 -n 8

# Or directly use CLARK command
./CLARK -T ./targets_addresses.txt -D $CLARK_DB_DIR -O sample.fastq -R virus_results -m 0 -n 8

# Or use the provided batch processing script
chmod +x virus_detection_pipeline.sh
./virus_detection_pipeline.sh
```

### 4. Calculate Virus Abundance

```bash
# Use Python script to calculate abundance
python calculate_virus_abundance.py \
    results/virus_results.csv \
    results/virus_abundance.csv
```

**Note**: On Linux/macOS systems, ensure scripts have execute permissions: `chmod +x *.sh`

## System Requirements

### Hardware Requirements

CLARK memory requirements depend on the version used and database:

- **CLARK (default version)**:
  - Loading bacteria database: 58 GB RAM
  - Building bacteria database: 156 GB RAM
  
- **CLARK-l** (lightweight version):
  - Suitable for 4 GB RAM laptops
  - Suitable for small metagenomic analysis
  
- **CLARK-S** (sparse k-mer version):
  - Maximum 101 GB RAM required (for spaced k-mer classification)

- **Storage Space**: 
  - At least 50GB available space (for storing databases and results)
  - Virus database approximately 5-10 GB
  - Complete database approximately 200-300 GB

- **CPU**: Multi-core processor (CLARK supports multithreading using OpenMP)

### Software Dependencies
- **Operating System**: 64-bit Linux or macOS
- **Compiler**: GNU GCC 4.4 or higher (64-bit support)
- **OpenMP**: For multithreading operations
- **Python**: Python 2.7 or Python 3.x (for auxiliary scripts)
- **Make**: For compilation

## Installation Steps

### Method 1: Install from GitHub (Recommended)

```bash
# 1. Clone CLARK repository
git clone https://github.com/rouni001/CLARK.git
cd CLARK

# 2. Run installation script (recommended)
./install.sh

# Or compile manually
make

# 3. Verify installation
./CLARK -v
```

### Method 2: Download from Official Website

```bash
# 1. Download latest version from CLARK official website
# Official website: http://clark.cs.ucr.edu
# Or download from GitHub Releases: https://github.com/rouni001/CLARK/releases

# 2. Extract
tar -xvf CLARKV1.3.0.tar.gz
cd CLARKV1.3.0

# 3. Run installation script
./install.sh

# 4. Add to PATH (optional)
export PATH=$PATH:$(pwd)
```

### Verify Installation

```bash
# Check if CLARK is working properly
./CLARK -v

# View help information
./CLARK -h
```

## Database Preparation

CLARK requires reference databases for classification. CLARK supports various databases, including NCBI RefSeq, NCBI nt, etc.

### NCBI RefSeq Database Download

#### Method 1: Use CLARK Auto-Download (Recommended)

CLARK provides scripts for automatic downloading and building NCBI RefSeq databases:

```bash
# 1. Set database directory
export CLARK_DB_DIR=/path/to/CLARK_DB
mkdir -p $CLARK_DB_DIR

# 2. Download only virus database (for virus detection and classification)
./set_targets.sh $CLARK_DB_DIR viruses

# 3. Download virus and other databases
./set_targets.sh $CLARK_DB_DIR viruses bacteria

# Note: CLARK v1.3.0 set_targets.sh syntax:
# ./set_targets.sh <DIR_DB/> <target1> <target2> ...
# No need to use --db parameter
```

#### Method 2: Manually Download NCBI RefSeq Virus Database

If you need to manually download or use mirror sites, you can follow these steps:

**Step 1: Access NCBI RefSeq Database**

- **NCBI RefSeq Homepage**: https://www.ncbi.nlm.nih.gov/refseq/
- **Virus RefSeq Database**: https://www.ncbi.nlm.nih.gov/refseq/about/viral/

**Step 2: Download Virus Genome Sequences**

```bash
# Create database directory
export CLARK_DB_DIR=/path/to/CLARK_DB/viruses
mkdir -p $CLARK_DB_DIR

# Method A: Direct download from NCBI FTP (recommended)
# FTP address: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
cd $CLARK_DB_DIR

# Download all virus genome FASTA files
wget -r -np -nH --cut-dirs=3 \
     ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.*.genomic.fna.gz

# Or use rsync (faster)
rsync -av --include="*.fna.gz" --exclude="*" \
      rsync://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ \
      $CLARK_DB_DIR/

# Decompress all files
gunzip *.fna.gz

# Merge all FASTA files
cat *.fna > viral_genomes.fasta
```

**Step 3: Download Classification Information Files**

```bash
# Download NCBI classification information
cd $CLARK_DB_DIR

# Download taxonomy database
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz

# Download accession to taxid mapping file
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
```

**Step 4: Build Database Using CLARK**

```bash
# If using custom downloaded database
./set_targets.sh $CLARK_DB_DIR custom

# Or use download script
./download_RefSeqDB.sh $CLARK_DB_DIR viruses
```

#### Method 3: Use NCBI Data Download Tools

```bash
# Use NCBI datasets tool (need to install first)
# Installation: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

# Download all virus genomes
datasets download genome taxon viruses --refseq --filename viral_genomes.zip

# Extract
unzip viral_genomes.zip
```

### Database Build Options

CLARK supports the following database types:
- **bacteria**: Bacteria database
- **viruses**: Virus database (**for virus detection and classification**)
- **fungi**: Fungi database
- **human**: Human genome
- **custom**: Custom database

### Database Build Examples

```bash
# Build only virus database (recommended for virus detection)
./set_targets.sh $CLARK_DB_DIR viruses

# Build virus and bacteria databases
./set_targets.sh $CLARK_DB_DIR viruses bacteria

# Build complete database (including all types)
./set_targets.sh $CLARK_DB_DIR bacteria viruses fungi human

# Use download script (if available)
./download_RefSeqDB.sh $CLARK_DB_DIR viruses
```

### Database Size Reference

- **Virus database**: Approximately 5-10 GB (compressed)
- **Bacteria database**: Approximately 50-100 GB (compressed)
- **Complete database**: Approximately 200-300 GB (compressed)

**Note**: Virus database is relatively small, suitable for dedicated virus detection and classification.

## Usage

### Basic Command Format

CLARK v1.3.0 provides two usage methods:

**Method 1: Directly use CLARK command**
```bash
CLARK [options] -T <target file> -D <database directory> -O <input file> -R <result file>
```

**Method 2: Use classification script (recommended)**
```bash
./classify_metagenome.sh -O <input file> -R <result file prefix> [other options]
```

### How to Choose K-mer Length

Choosing the appropriate k-mer length is important for classification results:

- **k = 20 or 21**: High sensitivity, suitable for detecting low abundance species (recommended for virus detection)
- **k = 31**: Balanced speed, accuracy and RAM usage (default value, recommended for most cases)
- **k = 27**: Default for CLARK-l

**Recommendations**:
- For virus detection, if pursuing high sensitivity, use `-k 20` or `-k 21`
- For general analysis, use the default `-k 31`
- If memory is limited, you can use CLARK-l (default k=27)

### Main Parameter Description

| Parameter | Description |
|-----------|-------------|
| `-T <file>` | Target definition file (required, generated by set_targets.sh) |
| `-D <directory>` | Database directory path (required) |
| `-O <file>` | Input sequence file (FASTA/FASTQ format, required) |
| `-R <file>` | Result output file (required) |
| `-P <file1> <file2>` | Paired-end sequencing files |
| `-m <mode>` | Execution mode: 0=full, 1=default, 2=express, 3=spectrum |
| `-k <integer>` | k-mer size (default 31, recommended 20-21 for high sensitivity) |
| `-n <integer>` | Number of threads (default 1) |
| `-s <factor>` | Sampling factor (default 2) |
| `--gzipped` | Input file is gzip compressed format |
| `--light` | Use CLARK-l mode (low memory) |
| `--spaced` | Use CLARK-S mode (spaced k-mer) |
| `--long` | For ultra-long sequences (e.g., long contigs, Nanopore/Pacbio reads) |
| `--extended` | Extended output (full mode only) |

## Common Command Examples

### Example 1: Using Classification Script (Recommended Method)

```bash
# Basic classification command (using classify_metagenome.sh script)
./classify_metagenome.sh -O sample.fa -R result

# Use 20-mers for higher sensitivity
./classify_metagenome.sh -O sample.fa -R result -k 20

# Full mode, 8 threads
./classify_metagenome.sh -O sample.fa -R result -m 0 -n 8

# Process gzip compressed files
./classify_metagenome.sh -O sample.fa.gz -R result -m 0 -n 8 --gzipped
```

### Example 1b: Directly Use CLARK Command

```bash
# Basic classification command (need to run set_targets.sh first to generate targets file)
./CLARK -T ./targets_addresses.txt \
        -D $CLARK_DB_DIR \
        -O sample.fa \
        -R result \
        -k 31 \
        -n 8
```

### Example 2: Paired-end Sequencing Data Classification

```bash
# Use classification script to process paired-end data
./classify_metagenome.sh -O sample1.fastq sample2.fastq -R result

# Use CLARK command to process paired-end data
./CLARK -T ./targets_addresses.txt \
        -D $CLARK_DB_DIR \
        -P sample1.fastq sample2.fastq \
        -R paired.results \
        -k 20 \
        -n 8
```

### Example 3: Compressed File Classification

```bash
# Use classification script to process gzip compressed files
./classify_metagenome.sh -O sample.fa.gz -R result -m 0 -n 8 --gzipped

# Use CLARK-S mode to process compressed files
./classify_metagenome.sh -O sample.fa.gz -R result -m 2 -n 8 --spaced --gzipped
```

### Example 4: FASTA File Classification

```bash
# Use classification script to process FASTA files
./classify_metagenome.sh -O sequences.fasta -R sequences_results -m 0 -k 31 -n 8

# Full mode, extended output
./classify_metagenome.sh -O sequences.fasta -R sequences_results -m 0 --extended
```

### Example 5: Using CLARK-l (Low Memory Mode)

```bash
# Use CLARK-l to process paired-end data (suitable for 4GB RAM laptops)
./classify_metagenome.sh -P sample1.fastq sample2.fastq -R result --light

# CLARK-l uses 27-mers (default)
```

### Example 6: Using CLARK-S (Spaced k-mer Mode)

```bash
# Use CLARK-S to process paired-end data
./classify_metagenome.sh -P sample1.fastq sample2.fastq -R result --spaced

# CLARK-S with full mode, 8 threads
./classify_metagenome.sh -O sample.fa -R result -m 0 -n 8 --spaced

# Reduce CLARK-S RAM usage
./classify_metagenome.sh -O sample.fa -R result --spaced -s 2
```

### Example 7: Virus Detection and Classification (Specifically for Virus Analysis)

```bash
# Method 1: Use classification script (recommended)
./classify_metagenome.sh -O sample.fastq -R virus_classification -m 0 -k 31 -n 8

# Method 2: Directly use CLARK command
./CLARK -T ./targets_addresses.txt \
        -D $CLARK_DB_DIR \
        -O sample.fastq \
        -R virus_classification \
        -m 0 \
        -k 31 \
        -n 8

# Notes:
# -D points to virus database directory (database built using ./set_targets.sh $CLARK_DB_DIR viruses)
# -m 0 means full mode (complete mode)
# -k 31 is k-mer size, using 20-21 can achieve higher sensitivity
```

### Example 8: Virus Classification (Paired-end Sequencing Data)

```bash
# Use classification script to process paired-end virus data
./classify_metagenome.sh -O sample_R1.fastq sample_R2.fastq -R virus_classification -m 0 -n 8

# Use CLARK command to process paired-end virus data
./CLARK -T ./targets_addresses.txt \
        -D $CLARK_DB_DIR \
        -P sample_R1.fastq sample_R2.fastq \
        -R virus_classification \
        -m 0 \
        -n 8
```

## Result Interpretation

### Output File Description

CLARK output format depends on the mode used:

#### Full Mode Output Format (-m 0)

```csv
<Object_ID>,<hit count in target 1>,...,<hit count in target N>,<Length of object>,<Gamma>,<first assignment>,<hit count of first>,<second assignment>,<hit count of second>,<confidence score>
```

#### Default/Express Mode Output Format (-m 1 or -m 2)

```csv
<Object_ID>,<Length of object>,<1st_assignment>
```

### Result File Format Examples

**Default Mode Example:**
```csv
read_1,150,Escherichia coli
read_2,150,Escherichia coli
read_3,150,Staphylococcus aureus
```

**Full Mode Example (Extended Output):**
```csv
read_1,95,5,150,0.95,Escherichia coli,95,Staphylococcus aureus,5,0.95
```

### Result Statistics and Abundance Calculation

#### Method 1: Use CLARK Built-in Abundance Estimation Tool (Recommended)

CLARK v1.3.0 provides dedicated abundance estimation script:

```bash
# Basic abundance estimation
./estimate_abundance.sh -F ./result.csv -D $CLARK_DB_DIR

# High confidence assignments
./estimate_abundance.sh -F ./result.csv -D $CLARK_DB_DIR --highconfidence

# Filter by confidence score (e.g., 0.8)
./estimate_abundance.sh -F ./result.csv -D $CLARK_DB_DIR -c 0.80

# Filter by gamma score (e.g., 0.03)
./estimate_abundance.sh -F ./result.csv -D $CLARK_DB_DIR -g 0.03

# Output MetaPhlAn format
./estimate_abundance.sh -F ./result.csv -D $CLARK_DB_DIR -g 0.03 --mpa

# Filter by abundance (e.g., >2%)
./estimate_abundance.sh -F ./result.csv -D $CLARK_DB_DIR -a 2
```

#### Method 2: Manual Abundance Calculation (Python Script)

Create a Python script to calculate virus abundance:

```bash
# Create abundance calculation script
cat > calculate_virus_abundance.py << 'EOF'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLARK result abundance calculation script
For calculating relative abundance of viruses or other taxa
"""

import pandas as pd
import sys
from collections import Counter

def calculate_abundance(input_file, output_file):
    """
    Calculate abundance of taxa
    
    Parameters:
        input_file: CLARK output CSV file
        output_file: Abundance report output file
    """
    # Read CLARK result file
    df = pd.read_csv(input_file)
    
    # Count sequences for each taxon
    classification_counts = Counter(df['Classification_Name'])
    
    # Calculate total number of sequences
    total_sequences = len(df)
    
    # Calculate relative abundance
    abundance_data = []
    for classification, count in classification_counts.items():
        relative_abundance = (count / total_sequences) * 100
        abundance_data.append({
            'Classification': classification,
            'Read_Count': count,
            'Relative_Abundance_%': round(relative_abundance, 4),
            'Total_Reads': total_sequences
        })
    
    # Create DataFrame and sort
    abundance_df = pd.DataFrame(abundance_data)
    abundance_df = abundance_df.sort_values('Relative_Abundance_%', ascending=False)
    
    # Save results
    abundance_df.to_csv(output_file, index=False)
    
    # Print summary
    print(f"Total sequences: {total_sequences}")
    print(f"Number of taxa detected: {len(abundance_df)}")
    print(f"\nTop 10 most abundant taxa:")
    print(abundance_df.head(10).to_string(index=False))
    
    return abundance_df

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python calculate_virus_abundance.py <input CSV file> <output file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    calculate_abundance(input_file, output_file)
EOF

# Use script to calculate abundance
python calculate_virus_abundance.py \
    virus_classification.csv \
    virus_abundance_report.csv
```

#### Method 3: Use R for Abundance Analysis and Visualization

```r
# Create R script for abundance analysis
cat > analyze_abundance.R << 'EOF'
#!/usr/bin/env Rscript
# CLARK result abundance analysis and visualization

library(ggplot2)
library(dplyr)

# Read CLARK results
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript analyze_abundance.R <CLARK result CSV file>")
}

input_file <- args[1]
df <- read.csv(input_file)

# Calculate abundance
abundance <- df %>%
  count(Classification_Name, name = "Count") %>%
  mutate(Relative_Abundance = (Count / sum(Count)) * 100) %>%
  arrange(desc(Relative_Abundance))

# Save abundance table
write.csv(abundance, "abundance_report.csv", row.names = FALSE)

# Plot top 20 most abundant viruses
top20 <- head(abundance, 20)

p <- ggplot(top20, aes(x = reorder(Classification_Name, Relative_Abundance), 
                       y = Relative_Abundance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Virus Relative Abundance",
       x = "Virus Classification",
       y = "Relative Abundance (%)") +
  theme_minimal()

ggsave("virus_abundance_plot.png", p, width = 12, height = 8, dpi = 300)

print("Analysis completed!")
print(paste("Number of virus species detected:", nrow(abundance)))
print(head(abundance, 10))
EOF

# Run R script
Rscript analyze_abundance.R virus_classification.csv
```

#### Method 4: Simple Shell Script Abundance Calculation

```bash
# Create simple abundance calculation script
cat > calculate_abundance.sh << 'EOF'
#!/bin/bash
# Simple abundance calculation script

INPUT=$1
OUTPUT=$2

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: ./calculate_abundance.sh <input CSV> <output file>"
    exit 1
fi

# Extract classification name column (assume column 3 is classification name)
# Count occurrences of each classification
awk -F',' 'NR>1 {print $3}' $INPUT | sort | uniq -c | sort -rn | \
awk 'BEGIN {
    total=0
    print "Classification,Read_Count,Relative_Abundance_%"
}
{
    count[$2] = $1
    total += $1
}
END {
    for (class in count) {
        abundance = (count[class] / total) * 100
        printf "%s,%d,%.4f\n", class, count[class], abundance
    }
}' > $OUTPUT

echo "Abundance calculation completed! Results saved in $OUTPUT"
EOF

chmod +x calculate_abundance.sh

# Use script
./calculate_abundance.sh virus_classification.csv virus_abundance.csv
```


## Complete Virus Detection and Classification Pipeline

### Complete Workflow Example

```bash
#!/bin/bash
# Complete virus detection and classification pipeline script

# 1. Set environment variables
export CLARK_DB_DIR=/path/to/CLARK_DB
export INPUT_DIR=/path/to/input_data
export OUTPUT_DIR=/path/to/output

# 2. Ensure virus database is built
if [ ! -d "$CLARK_DB_DIR/viruses" ]; then
    echo "Building virus database..."
    ./set_targets.sh $CLARK_DB_DIR viruses
fi

# 3. Virus classification for each sample
for sample in $INPUT_DIR/*.fastq; do
    sample_name=$(basename $sample .fastq)
    
    echo "Processing sample: $sample_name"
    
    ./CLARK -T $sample \
            -O $OUTPUT_DIR \
            -R ${sample_name}_virus \
            -D $CLARK_DB_DIR \
            -m 0 \
            -k 31 \
            -t 8 \
            --threshold 0.7
done

# 4. Calculate virus abundance for each sample
for result in $OUTPUT_DIR/*_virus.csv; do
    sample_name=$(basename $result _virus.csv)
    python calculate_virus_abundance.py \
        $result \
        ${OUTPUT_DIR}/${sample_name}_abundance.csv
done

echo "Virus detection and classification completed!"
```

## Advanced Usage

### Batch Processing Multiple Samples

```bash
# Create batch processing script
for sample in sample1 sample2 sample3; do
    ./CLARK -T ${sample}.fastq \
            -O ./results \
            -R ${sample}_results \
            -D $CLARK_DB_DIR \
            -m 0 \
            -t 8
done
```

### Result Merging and Analysis

```bash
# Merge multiple result files
cat sample1_results.csv sample2_results.csv > combined_results.csv

# Use Python or other tools for subsequent analysis
python analyze_clark_results.py combined_results.csv
```

### Virus Detection Dedicated Batch Processing Script

```bash
#!/bin/bash
# Virus detection batch processing script

export CLARK_DB_DIR=/path/to/CLARK_DB/viruses
export THREADS=16

# Process all FASTQ files
for fastq_file in *.fastq; do
    sample_name=$(basename $fastq_file .fastq)
    
    echo "=== Processing sample: $sample_name ==="
    
    # Run CLARK classification (using classification script)
    ./classify_metagenome.sh -O $fastq_file \
            -R ./virus_results/${sample_name}_virus \
            -m 0 \
            -k 31 \
            -n $THREADS
    
    # Calculate abundance
    if [ -f "./virus_results/${sample_name}_virus.csv" ]; then
        python calculate_virus_abundance.py \
            ./virus_results/${sample_name}_virus.csv \
            ./virus_results/${sample_name}_abundance.csv
    fi
done

echo "All samples processed!"
```

## References

- **Official GitHub**: https://github.com/rouni001/CLARK
- **Official Website**: http://clark.cs.ucr.edu
- **User Community**: https://groups.google.com/g/clarkusers
- **Latest Version**: v1.3.0.0 (Released May 17, 2024)
- **Paper**: Ounit R, Wanamaker S, Close TJ, Lonardi S. CLARK: fast and accurate classification of metagenomic and genomic sequences using discriminative k-mers. BMC Genomics. 2015;16:236.

## Version Information

This guide is based on CLARK v1.3.0.0. Specific parameters may vary slightly between versions, please refer to the documentation for the corresponding version.

### Version History

- **v1.3.0.0 (May 17, 2024)**: Latest version
- **v1.2.6.1**: Previous stable version
- **v1.1.3 (June 3, 2015)**: Extended features and simplified output
- **v1.1.2 (April 22, 2015)**: Added abundance estimation script and gzip support
- **v1.1.1 (April 15, 2015)**: Improved database loading efficiency
- **v1.1 (February 20, 2015)**: Improved storage efficiency and multithreaded database loading
- **v1.0 (September 1, 2014)**: Initial version

### Version Compatibility

- Databases created with v1.0 are incompatible with v1.1 or newer versions, and need to be rebuilt using v1.1 or newer versions.

## License

CLARK software follows its official license. Please read the license terms carefully before use.

---

**Note**: This guide provides basic installation and usage methods for CLARK. For specific research needs, it is recommended to refer to the official documentation or contact the CLARK development team for support.

## Provided Auxiliary Scripts

This guide includes the following auxiliary scripts that can be used directly:

1. **`calculate_virus_abundance.py`**: Python script for calculating virus abundance
2. **`analyze_abundance.R`**: R script for abundance analysis and visualization
3. **`calculate_abundance.sh`**: Shell script using awk to calculate abundance (no Python required)
4. **`virus_detection_pipeline.sh`**: Complete virus detection pipeline script

### Using Auxiliary Scripts

```bash
# Python abundance calculation (recommended)
python calculate_virus_abundance.py input.csv output.csv

# R abundance analysis and visualization (requires ggplot2 and dplyr)
Rscript analyze_abundance.R input.csv

# Shell script abundance calculation
./calculate_abundance.sh input.csv output.csv

# Complete virus detection pipeline
export CLARK_DB_DIR=/path/to/CLARK_DB
export INPUT_DIR=/path/to/fastq_files
./virus_detection_pipeline.sh
```

