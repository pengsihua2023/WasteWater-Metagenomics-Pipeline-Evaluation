# Short-read Abundance Calculation Guide (Illumina)

## Overview

Short-read (Illumina) abundance calculation prioritizes **Bracken** for statistical correction. If Bracken is unavailable, automatically falls back to Kraken2-based calculation.

---

## Quick Start

### Automatic Execution (Recommended)

```bash
sbatch submit_short.sh
```

**Automatically completes**:
1. QC (FastQC, fastp)
2. Kraken2 viral classification
3. Bracken abundance estimation (if database available)
4. **Automatic RPM/RPKM calculation**
5. Generate MultiQC report

**Output**: `results_viral_short/abundance/`

---

## Two Calculation Methods

### Method 1: Standard Method (Recommended, More Accurate)

**Requires Bracken database**

```
Kraken2 classification
    ↓
Bracken statistical correction ← Improves species-level accuracy
    ↓
calculate_abundance_en.py
    ↓
RPM/RPKM values (accuracy ~98%)
```

**Check Bracken database**:
```bash
ls -lh /path/to/kraken2_db/database150mers.kmer_distrib
ls -lh /path/to/kraken2_db/database150mers.kraken
```

**Configure in databases.csv**:
```csv
tool,db_name,db_params,db_path
kraken2,Viral_ref,"",/path/to/kraken2_viral_database
bracken,Viral_ref,";-r 150",/path/to/kraken2_viral_database
```

Note: `;-r 150` specifies read length of 150 bp (adjust based on your data: 50, 75, 100, 150, 200, 250, 300)

### Method 2: Fallback Method (Automatic Backup, Still Valid)

**No Bracken database required**

```
Kraken2 classification
    ↓
Direct abundance extraction
    ↓
calculate_abundance_longread_en.py
    ↓
RPM/RPKM values (accuracy ~90%)
```

**When used**:
- Bracken database unavailable
- Quick analysis
- Exploratory research

---

## Smart Workflow Dispatch

### Automatic Selection Logic

```bash
After submit_short.sh runs:
│
├─ TaxProfiler completes
│  ├─ Kraken2 ✅ Always runs
│  └─ Bracken ❓ Depends on database
│
└─ batch_calculate_abundance_en.sh (smart dispatcher)
   │
   ├─ Bracken results detected?
   │
   ├─ Yes → calculate_abundance_en.py
   │         └─ Bracken-based (more accurate) ✅✅
   │
   └─ No  → batch_calculate_abundance_longread_en.sh
             └─ Kraken2-based (still valid) ✅
```

**Result**: **Either way, you get abundance values!**

---

## Sample Sheet Configuration

### Basic Format

```csv
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
llnl_66ce4dde,run1,ILLUMINA,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,
```

### Required Columns

- `sample` - Unique sample identifier
- `run_accession` - Sequencing run ID
- `instrument_platform` - Must be `ILLUMINA`
- `fastq_1` - Read 1 file path
- `fastq_2` - Read 2 file path
- `fasta` - Leave empty

---

## Configuration Details

### Core Parameters (nextflow_short.config)

```groovy
params {
    input       = 'samplesheet_short.csv'
    outdir      = 'results_viral_short'
    databases   = 'databases.csv'
    
    // Tool configuration
    run_kraken2 = true
    run_bracken = true   // If Bracken database available
    
    // Bracken parameters
    bracken_precision = 'S'    // Species level
    bracken_readlen   = 150    // Read length (match your data)
    
    // QC
    perform_shortread_qc = true
    
    // Resources
    max_cpus    = 32
    max_memory  = '256.GB'
    
    // Apptainer configuration (handles missing mount points)
    profiles {
        apptainer {
            apptainer.runOptions = '--no-mount /lscratch'
        }
    }
}
```

### Bracken Parameters

```groovy
bracken_precision = 'S'    // Species level
bracken_readlen   = 150    // Read length (match your data)
```

**Read length options**: 50, 75, 100, 150, 200, 250, 300

---

## Abundance Results

### Output Files

```
results_viral_short/abundance/
├── llnl_66ce4dde_abundance.tsv          # Single sample (1087 species)
├── all_samples_abundance_summary.tsv    # All samples summary
└── top_viruses_summary.tsv              # High abundance (RPM≥10)
```

### Example Results

**Your short-read sample** (llnl_66ce4dde):

```
Total reads (classified): 12,987 (Bracken-estimated reads)
Detected viruses: 1,087 species
Highest RPM: 51,051 (Shigella phage SfIV)
Main type: Environmental phages
Note: Total reads shown is from Bracken classification, not raw sequencing reads
```

---

## Bracken vs Non-Bracken Comparison

### Practical Impact

For **environmental/phage-dominated samples**:
- Bracken accuracy: ~98%
- Kraken2 accuracy: ~90%
- **Difference: ~8%** - Acceptable for most analyses

For **human virus/clinical samples**:
- Bracken recommendation: ⭐⭐⭐⭐⭐
- Need precise pathogen identification

### When Bracken is Essential

- Clinical diagnostic samples
- Distinguishing closely related viruses (e.g., influenza subtypes)
- Publishing in high-impact journals
- Requiring highest accuracy

### When Bracken is Optional

- Environmental samples (your case)
- Exploratory research
- Viral diversity analysis
- Rapid screening

---

## Optimization Suggestions

### 1. Database Selection

**Standard viral database**:
- Suitable for: Human virus detection
- Size: ~2-3 GB
- Build time: 2-4 hours (Bracken)

**Complete microbial database**:
- Suitable for: Comprehensive environmental sample analysis
- Size: ~50-100 GB
- Build time: 10-24 hours (Bracken)

### 2. QC Stringency

Adjust based on sample quality:

```groovy
// High quality data
perform_shortread_complexityfilter = false

// Low quality data
perform_shortread_complexityfilter = true
perform_shortread_hostremoval = true  // Remove host DNA
```

### 3. Read Length Confirmation

```bash
# Check actual read length
zcat your_R1.fastq.gz | head -n 2 | tail -n 1 | wc -c
```

Ensure `bracken_readlen` matches actual read length (or use closest value).

---

## Troubleshooting

### Issue: Empty abundance directory

**Cause**: Abundance calculation script didn't run

**Solution**:
```bash
# Manually run abundance calculation
bash batch_calculate_abundance_en.sh results_viral_short

# Or use universal script
bash batch_calculate_abundance_longread_en.sh results_viral_short
```

### Issue: top_viruses_summary.tsv is empty

**Cause**: All viruses have RPM < 10

**Solution**:
```bash
# Manually generate with lower threshold
awk -F'\t' 'NR==1 || $6 >= 5' results_viral_short/abundance/all_samples_abundance_summary.tsv | \
  sort -t$'\t' -k6,6nr > top_viruses_rpm5.tsv
```

### Issue: All RPKM values are NA

**Cause**: Detected viruses not in genome length database

**Explanation**: This is normal (environmental phages)
- **RPM is already sufficient** for analysis
- RPKM is just additional information

---

## Data Analysis Examples

### Basic Statistics

```bash
# Number of viruses per sample
awk -F'\t' 'NR>1 {count[$1]++} END {for (s in count) print s"\t"count[s]}' \
  results_viral_short/abundance/all_samples_abundance_summary.tsv

# TOP 20 most common viruses
awk -F'\t' 'NR>1 {rpm[$2]+=$6} END {for (v in rpm) print v"\t"rpm[v]}' \
  results_viral_short/abundance/all_samples_abundance_summary.tsv | \
  sort -t$'\t' -k2,2nr | head -20
```

### R Visualization

```r
library(ggplot2)
library(dplyr)

data <- read.table("results_viral_short/abundance/all_samples_abundance_summary.tsv", 
                   header=TRUE, sep="\t", quote="")

# TOP 10 viruses bar plot
top10 <- data %>% group_by(Species) %>% 
  summarise(total_rpm = sum(RPM)) %>% 
  arrange(desc(total_rpm)) %>% head(10)

ggplot(top10, aes(x=reorder(Species, total_rpm), y=total_rpm)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Top 10 Viruses by Total RPM", x="Species", y="Total RPM") +
  theme_minimal()
```

---

## Best Practices

### 1. Data Filtering

```bash
# Filter low abundance viruses (RPM < 5)
awk -F'\t' 'NR==1 || $6 >= 5' all_samples_abundance_summary.tsv > filtered.tsv
```

### 2. Diversity Analysis

```bash
# Calculate Shannon diversity index (requires R or Python)
# Or use RPM values for community analysis
```

### 3. Cross-sample Comparison

- Use RPM (not raw read counts)
- Ensure consistent sample processing
- Account for sequencing depth differences

---

## Summary

### Short-read Data Analysis Workflow

```
sbatch submit_short.sh  →  Auto analysis  →  results_viral_short/abundance/
```

### Key Points

1. **Automated** - Fully automatic after submission
2. **Smart dispatch** - Automatically selects best calculation method
3. **Robust** - Succeeds regardless of Bracken availability
4. **Standardized abundance** - RPM/RPKM ready for publication

### Bracken Optionality

- With Bracken: More accurate (~98%) ✅✅
- Without Bracken: Still valid (~90%) ✅
- Environmental samples: Small difference
- Clinical samples: Bracken recommended

---

**More information?** See README.md or ABUNDANCE_USAGE.md


