# RVDB-Based Dual-Method Viral Metagenome Analysis Pipeline

## Project Overview

This is a Nextflow-based viral metagenome analysis pipeline using **dual-method strategies** to maximize virus identification coverage. Features **dual-assembler comparison** for short reads (MEGAHIT + SPAdes) and **dual-track analysis** for long reads (viralFlye feature-based + Diamond similarity-based). All viral classification is performed using **RVDB (Reference Viral Database)** via Diamond BLASTP. Supports both short-read (**Illumina**) and long-read (**Nanopore/PacBio**) data types.

**Version**: 4.2.0  
**Author**: Sihua Peng  
**Last Updated**: 2025-11-17

---

## Key Features

### Short-Read Workflow (Illumina Paired-end) - Dual-Assembler Comparison
1. **Quality Control**: fastp automatic adapter removal and quality filtering
2. **Dual Assembler Parallel**: MEGAHIT + metaSPAdes
3. **Gene Prediction**: Prodigal (metagenome mode)
4. **Viral Protein Classification**: Diamond BLASTP vs RVDB
5. **Comprehensive Comparison**: Merge and compare results from both assemblers
   - Complete 7-level taxonomic comparison: Kingdom → Phylum → Class → Order → Family → Genus → Species
   - Side-by-side comparison of MEGAHIT vs SPAdes at each taxonomic level
   - Identify assembler-specific and shared discoveries
6. **Viral Abundance Calculation** (if enabled): Calculate RPM and RPKM for identified viral contigs
   - BWA alignment of reads to viral contigs (short reads)
   - Read counting and normalization
   - Generate abundance reports in both text and CSV formats

### Long-Read Workflow (Nanopore/PacBio) - Dual-Track Complementary Analysis
1. **Metagenome Assembly**: MetaFlye (--meta mode)
2. **Viral Identification and Filtering**: viralFlye analysis
   - Use Pfam-A HMM and viralVerify to identify viral features
   - Filter viral sequences from all contigs (typical filtering rate ~18%)
   - Generate lists of linear and circular viral contigs
3. **Dual-Track Parallel Analysis**:
   - **Track 1 (All Sequences)**: No filtering, analyze all ~1,200 MetaFlye contigs
     - Includes 222 viruses identified by viralFlye
     - Includes 991 sequences filtered by viralFlye (bacteria/eukaryotes/unknown/potential new viruses)
   - **Track 2 (Viral Sequences)**: Analyze only ~222 high-confidence viral contigs identified by viralFlye
4. **Gene Prediction**: Prodigal predicts genes for both sequence sets separately
5. **Viral Protein Classification**: Diamond BLASTP vs RVDB (independent classification for both sets)
6. **Taxonomic Resolution**: Both result sets provide complete taxonomic lineage (Kingdom/Phylum/Family/Species)
7. **Comparative Analysis**: Compare results from both tracks for complementary virus identification
   - Missed by viralFlye but found by Diamond → distant viruses
   - Not found by Diamond but identified by viralFlye → new viral features without similar sequences
8. **Viral Abundance Calculation** (if enabled): Calculate RPM and RPKM for identified viral contigs
   - minimap2 alignment of long reads to viral contigs
   - Read counting and normalization
   - Generate abundance reports in both text and CSV formats

---

## System Requirements

### Required Software
- **Nextflow** >= 24.04
- **Conda/Mamba** (environment management)
- **Apptainer/Singularity** (container support)
- **SLURM** cluster scheduler

### Computational Resources
- **Short reads**:
  - CPU: 32 cores
  - Memory: 512 GB (SPAdes requires large memory)
  - Time: 48-72 hours
  
- **Long reads**:
  - CPU: 32 cores
  - Memory: 128 GB
  - Time: 24-48 hours

---

## Installation and Configuration

### 1. Create Nextflow Environment

```bash
# Load Conda
module load Miniforge3/24.11.3-0

# Create environment
conda create -n nextflow_env -c bioconda nextflow diamond-aligner prodigal flye fastp pandas -y

# Activate environment
conda activate nextflow_env
```

### 2. Create viralFlye Environment (only required for long-read viral analysis)

```bash
# Create separate environment
conda create -n viralFlye_env python=3.7 -y
conda activate viralFlye_env

# Install viralFlye and dependencies
conda install -c bioconda viralflye seqtk biopython -y

# Verify installation
viralFlye.py --help
seqtk
python -c "from Bio import SeqIO; print('Biopython OK')"
```

**Important**:
- `seqtk`: Used to extract viral sequences from assembly results
- `biopython`: Backup sequence extraction tool
- If these tools are missing, viralFlye will fall back to using all MetaFlye contigs

### 3. Download Databases

#### RVDB Viral Protein Database (Required)
```bash
# Create directory
mkdir -p /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB

cd /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB

# Download RVDB
wget https://rvdb-prot.pasteur.fr/files/RVDB_prot.fasta.xz

# Decompress and build Diamond index
xz -d RVDB_prot.fasta.xz
diamond makedb --in RVDB_prot.fasta -d RVDB_prot_ref
```

#### Pfam-A HMM Database (required for viralFlye)
```bash
mkdir -p /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/Pfam

cd /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/Pfam

# Download Pfam-A
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
```

#### NCBI Taxonomy Database (required for taxonomic resolution)
```bash
# Usually included in RVDB directory
# If separate download is needed:
cd /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
```

---

## Usage

### Prepare Input Data

#### Short-Read Data (samplesheet_short.csv)
```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

#### Long-Read Data (samplesheet_long.csv)
```csv
sample,fastq_long
sample1,/path/to/sample1.fastq.gz
sample2,/path/to/sample2.fastq.gz
```

### Run Analysis

#### Method 1: Using SLURM Script (Recommended)

```bash
# Short-read analysis
sbatch run_metagenome_assembly_classification_en.sh short

# Long-read viral analysis
sbatch run_metagenome_assembly_classification_en.sh long

# Using custom samplesheet
sbatch run_metagenome_assembly_classification_en.sh long my_samples.csv
```

#### Method 2: Run Nextflow Directly

```bash
# Short reads
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet_short.csv \
    --outdir results_short \
    --diamond_db /path/to/RVDB_prot_ref.dmnd \
    --read_type short

# Long reads
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet_long.csv \
    --outdir results_long \
    --diamond_db /path/to/RVDB_prot_ref.dmnd \
    --read_type long
```

### Resume Interrupted Analysis

```bash
# Add -resume parameter
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet_long.csv \
    --outdir results_long \
    --diamond_db /path/to/RVDB_prot_ref.dmnd \
    --read_type long \
    -resume
```

---

## Output Results

### Short-Read Output (results_short/)

```
results_short/
├── fastp/                      # Quality control reports
│   ├── *_fastp.html           # HTML quality reports
│   └── *_fastp.json           # JSON quality data
├── assembly_megahit/           # MEGAHIT assembly results
│   └── *_megahit_contigs.fa   # Assembled contigs
├── assembly_spades/            # SPAdes assembly results
│   └── *_spades_contigs.fa    # Assembled contigs
├── prodigal_megahit/           # MEGAHIT gene prediction
│   └── *_megahit_proteins.faa # Predicted protein sequences
├── prodigal_spades/            # SPAdes gene prediction
│   └── *_spades_proteins.faa  # Predicted protein sequences
├── diamond_megahit/            # MEGAHIT Diamond classification
│   └── *_megahit_diamond.txt  # Classification results
├── diamond_spades/             # SPAdes Diamond classification
│   └── *_spades_diamond.txt   # Classification results
├── merged_reports/             # Comprehensive analysis reports
│   ├── *_merged_report.txt    # Text format comprehensive report
│   ├── *_merged_report.csv    # CSV format detailed comparison
│   ├── *_megahit_with_taxonomy.txt  # With taxonomic information
│   └── *_spades_with_taxonomy.txt   # With taxonomic information
├── abundance_megahit/          # Viral abundance for MEGAHIT (if enabled)
│   ├── *_megahit_abundance.txt # Abundance report with RPM and RPKM (text format)
│   └── *_megahit_abundance.csv # Abundance data with RPM and RPKM (CSV format)
└── abundance_spades/           # Viral abundance for SPAdes (if enabled)
    ├── *_spades_abundance.txt  # Abundance report with RPM and RPKM (text format)
    └── *_spades_abundance.csv  # Abundance data with RPM and RPKM (CSV format)
```

### Long-Read Output (results_long/) - Dual-Track Analysis

```
results_long/
├── assembly_metaflye/                    # MetaFlye assembly results (all sequences)
│   └── *_metaflye_contigs.fa            # Metagenome contigs
├── assembly_viralflye/                   # viralFlye viral sequence refinement (viral only)
│   └── *_viralflye_contigs.fa           # Viral contigs (filtered and optimized)
│
├── [Track 1: All Sequences Analysis - Comprehensive coverage, avoid omissions]
├── prodigal_metaflye/                    # MetaFlye gene prediction (all sequences)
│   ├── *_metaflye_proteins.faa          # Predicted protein sequences
│   └── *_metaflye_genes.fna             # Predicted gene sequences
├── diamond_metaflye/                     # Diamond classification (all sequences)
│   └── *_metaflye_diamond.txt           # BLAST format classification results
├── taxonomy_metaflye/                    # Taxonomic resolution (all sequences)
│   ├── *_metaflye_diamond_with_taxonomy.txt  # Complete taxonomic info (22 columns)
│   └── *_metaflye_taxonomy_summary.txt       # Statistical summary
│
├── [Track 2: Viral Sequences Analysis - High-purity viral dataset]
├── prodigal_viralflye/                   # viralFlye gene prediction (viral only)
│   ├── *_viralflye_proteins.faa         # Predicted viral protein sequences
│   └── *_viralflye_genes.fna            # Predicted viral gene sequences
├── diamond_viralflye/                    # Diamond classification (viral only)
│   └── *_viralflye_diamond.txt          # BLAST format classification results
├── taxonomy_viralflye/                   # Taxonomic resolution (viral only)
│   ├── *_viralflye_diamond_with_taxonomy.txt  # Complete taxonomic info (22 columns)
│   └── *_viralflye_taxonomy_summary.txt       # Statistical summary
│
├── consensus_analysis/                   # Dual-track comparison and consensus analysis ⭐ Recommended!
│   ├── *_consensus_viruses.txt           # ★★★ Consensus viruses (confirmed by both methods, highest confidence)
│   │                                     #     Typical: ~180 contigs, ~850 viral proteins
│   ├── *_metaflye_only_viruses.txt       # ★ MetaFlye-only (distant viral candidates, needs verification)
│   │                                     #     Typical: ~500 contigs, ~950 viral proteins
│   ├── *_viralflye_only_viruses.txt      # ★★ viralFlye-only (feature-based but sequence-unique)
│   │                                     #     Typical: 0-10 contigs (rare)
│   └── *_dual_track_comparison.txt       # Detailed comparison report (with complete stats and taxonomic distribution)
│                                          #     Includes: Kingdom/Phylum/Family/Genus/Species
├── abundance_metaflye/                   # Viral abundance for MetaFlye (if enabled)
│   ├── *_metaflye_abundance.txt          # Abundance report with RPM and RPKM (text format)
│   └── *_metaflye_abundance.csv         # Abundance data with RPM and RPKM (CSV format)
└── abundance_viralflye/                  # Viral abundance for viralFlye (if enabled)
    ├── *_viralflye_abundance.txt         # Abundance report with RPM and RPKM (text format)
    └── *_viralflye_abundance.csv         # Abundance data with RPM and RPKM (CSV format)
```

---

## Short-Read Comprehensive Report Format

The `*_merged_report.txt` in `results_short/merged_reports/` contains complete taxonomic comparison across **7 levels**:

```
================================================================================
Diamond Comprehensive Analysis Report - MEGAHIT vs SPAdes Assembly Comparison
================================================================================

[Overall Statistics]
--------------------------------------------------------------------------------
SPAdes total hits:          51,835
MEGAHIT total hits:         16,226

SPAdes unique queries:      51,835
MEGAHIT unique queries:     16,226

SPAdes average identity:    41.56%
MEGAHIT average identity:   34.40%

SPAdes average length:      96.2 aa
MEGAHIT average length:     212.8 aa


[Kingdom Level Comparison (Top 15)]
--------------------------------------------------------------------------------
Kingdom                        SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
Bamfordvirae                   35680           11248           46928          
Orthornavirae                  530             219             749            
...

[Phylum Level Comparison (Top 15)]
--------------------------------------------------------------------------------
Phylum                         SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
Nucleocytoviricota             35680           11248           46928          
Peploviricota                  530             219             749            
...

[Class Level Comparison (Top 15)]
--------------------------------------------------------------------------------
Class                          SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
Megaviricetes                  25000           8000            33000          
...

[Order Level Comparison (Top 15)]
--------------------------------------------------------------------------------
Order                          SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
Algavirales                    15000           5000            20000          
...

[Family Level Comparison (Top 15)]
--------------------------------------------------------------------------------
Family                         SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
Mimiviridae                    18133           5445            23578          
Phycodnaviridae                6092            1936            8028           
Marseilleviridae               3652            1267            4919           
...

[Genus Level Comparison (Top 15)]
--------------------------------------------------------------------------------
Genus                          SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
Mimivirus                      12000           4000            16000          
Cafeteriavirus                 5000            1500            6500           
...

[Species Level Comparison (Top 20)]
--------------------------------------------------------------------------------
Species                                  SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
Acanthamoeba polyphaga mimivirus        2914            780             3694           
Fadolivirus algeromassiliense           2096            703             2799           
Cafeteria roenbergensis virus           2059            703             2762           
...

[Taxonomic ID Comparison (Top 20)]
[Unique Findings]
```

**Key Points**:
- All numbers represent **protein match counts**
- Sorted by total count (SPAdes + MEGAHIT) in descending order
- Provides complete taxonomic comparison from Kingdom to Species level

---

## Consensus Analysis Result Interpretation (Long-Read Dual-Track)

### Dual-Track Comparison Report Example

```
[Contig-Level Comparison]
--------------------------------------------------------------------------------
MetaFlye Viral Contigs:          709  ← Viral contigs identified by Diamond
viralFlye Viral Contigs:         181  ← Viral contigs identified by viralFlye
Consensus Viral Contigs (Both):  181 ★★★ Highest confidence
MetaFlye-Only Viral Contigs:     528 → Distant viral candidates
viralFlye-Only Viral Contigs:    0   → Feature-based but sequence-unique
```

### Meaning of Numbers:

**Important**: Distinguish between **contig counts** and **protein match counts**

| Statistic | Contig Count | Protein Matches | Description |
|-----------|--------------|-----------------|-------------|
| MetaFlye Total | 1,213 | 2,434 | All sequences |
| MetaFlye Viral | 709 | 1,809 | Viruses identified by Diamond |
| viralFlye Total | 222 | 1,129 | Sequences selected by viralFlye |
| viralFlye Viral | 181 | 844 | Viruses confirmed by both methods |
| **Consensus Viruses** | **181** | **844** | ★★★ Highest confidence |
| MetaFlye-Only | 528 | 965 | Distant viral candidates |
| viralFlye-Only | 0 | 0 | - |

### Key Findings:

1. **100% Consensus Rate**:
   - 181 viral contigs identified by viralFlye
   - Diamond **matched all of them** to viral proteins
   - Indicates viralFlye identification is very accurate!

2. **MetaFlye Found 528 Additional Viral Contigs**:
   - Matched to viral proteins by Diamond
   - But filtered by viralFlye (features don't meet standards)
   - Possibly:
     - Low gene density viruses
     - Distant viruses
     - Partial viral genome fragments
     - False positives (need verification)

3. **Gene Density Differences**:
   ```
   Consensus viruses: 844 proteins ÷ 181 contigs = 4.7 proteins/contig
   MetaFlye-only: 965 proteins ÷ 528 contigs = 1.8 proteins/contig
   ```
   - Consensus virus gene density is **2.6 times** that of MetaFlye-only
   - Proves viralFlye indeed selects high-quality, complete viral genomes

### Species Statistics Interpretation

```
Species Distribution (Top 15):
  Mimiviridae sp. ChoanoV1                     57  ← 57 protein matches
  Fadolivirus algeromassiliense                49
  ...
```

**The second column numbers** (e.g., 57) = **protein match count** for that viral species

**Meaning**:
- **Abundance indicator**: Larger numbers indicate higher viral abundance in the sample
- **Completeness indicator**: Larger numbers may indicate more complete genome assembly
- **Confidence indicator**: Multiple protein matches increase identification confidence

**Example Analysis**:
```
Mimiviridae sp. ChoanoV1: 57 proteins
  → Likely from 5-10 different contigs
  → Each contig has 5-10 gene matches
  → This virus has high abundance or good assembly quality in the sample
```

### Usage Recommendations:

1. **High Confidence Research**: Use `consensus_viruses.txt` (181 contigs, 844 proteins)
   - Confirmed by both methods
   - Lowest false positive rate
   
2. **Distant Virus Exploration**: Analyze `metaflye_only_viruses.txt` (528 contigs, 965 proteins)
   - Matched by Diamond but not identified by viralFlye
   - May contain distant viruses with atypical features
   - Requires further verification

3. **Abundance Analysis**: Focus on viruses with high protein counts in Species statistics
   - Top 3 viruses (57, 49, 45 proteins) are dominant viruses in the sample
   - Prioritize for in-depth study

---

## Diamond Output Format Explanation

### Raw Diamond Output (13 columns)

| Column | Name | Description | Example |
|--------|------|-------------|---------|
| 1 | qseqid | Query sequence ID | contig_1001_2 |
| 2 | sseqid | Matched viral protein ID | YP_009449147.1 |
| 3 | pident | Amino acid identity (%) | 36.2 |
| 4 | length | Alignment length (amino acids) | 348 |
| 5 | mismatch | Number of mismatches | 201 |
| 6 | gapopen | Number of gap openings | 8 |
| 7 | qstart | Query sequence start | 28 |
| 8 | qend | Query sequence end | 363 |
| 9 | sstart | Target sequence start | 14 |
| 10 | send | Target sequence end | 352 |
| 11 | evalue | E-value (smaller is better) | 6.93e-56 |
| 12 | bitscore | Bit score (higher is better) | 192 |
| 13 | staxids | NCBI taxonomy ID | 2023057 |

### Enhanced Output (22 columns, including taxonomic information)

Original 13 columns + the following taxonomic columns:

| Column | Name | Description | Example |
|--------|------|-------------|---------|
| 14 | organism_name | Virus name | Tomato bushy stunt virus |
| 15 | superkingdom | Superkingdom | N/A (viruses lack this level in new system) |
| 16 | kingdom | Kingdom | Bamfordvirae |
| 17 | phylum | Phylum | Nucleocytoviricota |
| 18 | class | Class | N/A |
| 19 | order | Order | Pimascovirales |
| 20 | family | Family | Mimiviridae |
| 21 | genus | Genus | Mimivirus |
| 22 | species | Species | Acanthamoeba polyphaga mimivirus |

---

## Quality Assessment Standards

### 🌟🌟🌟 High Confidence Viral Protein Matches
- **E-value** < 1e-20
- **Bit score** > 100
- **Alignment length** > 150 aa
- **Identity** > 30%

### 🌟🌟 Medium Confidence
- **E-value** < 1e-10
- **Bit score** > 60
- **Alignment length** > 80 aa

### ⚠️ Requires Further Verification
- **E-value** > 1e-5
- **Identity** < 25%
- **Alignment length** < 50 aa

---

## Workflow Parameters

### General Parameters
```bash
--input                    # Input samplesheet (CSV format)
--outdir                   # Output directory
--read_type                # 'short' or 'long'
--diamond_db              # Diamond database path
--help                    # Display help information
```

### Short-Read Parameters
```bash
--skip_fastp              # Skip quality control (default: false)
--fastp_qualified_quality # Minimum quality value (default: 20)
--fastp_min_length        # Minimum read length (default: 50)
--megahit_min_contig_len  # MEGAHIT minimum contig length (default: 1000)
--skip_merge_reports      # Skip comprehensive reports (default: false)
--skip_abundance          # Skip viral abundance calculation (default: false)
```

### Long-Read Parameters
```bash
--metaflye_genome_size    # Genome size estimate (e.g., '10m', optional)
--skip_viralflye          # Skip viralFlye (default: false)
--pfam_hmm                # Pfam-A HMM database path
--viralflye_min_length    # Minimum viral length (default: 5000)
--viralflye_completeness  # Completeness threshold (default: 0.5)
--skip_abundance          # Skip viral abundance calculation (default: false)
--long_read_preset        # minimap2 preset for long reads (default: 'map-ont' for Nanopore, 'map-pb' for PacBio)
```

### Taxonomy Parameters
```bash
--taxonomy_names          # NCBI taxonomy names.dmp path
--taxonomy_nodes          # NCBI taxonomy nodes.dmp path
```

### Diamond Parameters
```bash
--diamond_evalue          # E-value threshold (default: 1e-5)
--diamond_max_target_seqs # Maximum target sequences (default: 1)
```

---

## Usage Examples

### Example 1: Short-Read Viral Detection

```bash
# Prepare samplesheet
cat > samplesheet_short.csv << EOF
sample,fastq_1,fastq_2
virus_sample1,/data/sample1_R1.fastq.gz,/data/sample1_R2.fastq.gz
EOF

# Run analysis
sbatch run_metagenome_assembly_classification_en.sh short
```

**Output**: `results_short/` directory contains dual assembler results and comprehensive reports

### Example 2: Long-Read Viral Detection (with viralFlye dual-track analysis)

```bash
# Prepare samplesheet
cat > samplesheet_long.csv << EOF
sample,fastq_long
virus_nanopore,/data/nanopore_reads.fastq.gz
EOF

# Run analysis
sbatch run_metagenome_assembly_classification_en.sh long
```

**Output**: `results_long/` directory contains dual-track analysis results:
- **Track 1 (MetaFlye)**: Comprehensive analysis of all ~1,200 contigs
  - Based on sequence similarity via Diamond alignment
  - May identify ~700 viral contigs, ~1,800 viral proteins
- **Track 2 (viralFlye)**: ~200 high-purity viral contigs
  - viralFlye based on feature identification
  - May identify ~180 viral contigs, ~800 viral proteins
- **Consensus Analysis**: Automatically compare both result sets
  - Consensus viruses: ~180 contigs (highest confidence)
  - MetaFlye-only: ~500 contigs (distant viral candidates)
  - Complete taxonomic statistics (Kingdom → Species)

### Example 3: Custom Parameters

```bash
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input my_samples.csv \
    --outdir my_results \
    --diamond_db /path/to/RVDB.dmnd \
    --read_type long \
    --metaflye_genome_size 20m \
    --viralflye_min_length 3000
```

---

## Result Interpretation

### Long-Read Consensus Virus Analysis Result Example

```
================================================================================
Dual-Track Analysis Comparison Report - Viral Identification Consensus Analysis
================================================================================

Sample Name: llnl_66d1047e

[Overall Statistics]
--------------------------------------------------------------------------------
MetaFlye Total Hits:              2,434  ← Protein matches
MetaFlye Viral Hits:              1,809
viralFlye Total Hits:             1,129
viralFlye Viral Hits:             844

[Contig-Level Comparison]
--------------------------------------------------------------------------------
MetaFlye Viral Contigs:           709   ← Contig count
viralFlye Viral Contigs:          181
Consensus Viral Contigs (Both):   181 ★★★ Highest confidence
MetaFlye-Only Viral Contigs:      528 → Distant viral candidates
viralFlye-Only Viral Contigs:     0

[Taxonomic Distribution of Consensus Viruses]  ← Based on consensus viruses (highest confidence)
--------------------------------------------------------------------------------

Kingdom Distribution (Top 10):
  Bamfordvirae                                761  ← 761 viral proteins belong to this kingdom
  Orthornavirae                                44
  Heunggongvirae                               36

Phylum Distribution (Top 10):
  Nucleocytoviricota                          761  ← Nucleocytoplasmic large DNA viruses phylum
  Peploviricota                                29
  Lenarviricota                                22

Family Distribution (Top 10):
  Mimiviridae                                 353  ← Mimiviridae (giant viruses)
  Phycodnaviridae                             128  ← Algal DNA viruses
  Marseilleviridae                             93  ← Marseilleviruses
  Poxviridae                                   41

Genus Distribution (Top 15):
  Mimivirus                                    xx
  Marseillevirus                               xx
  ...

Species Distribution (Top 15):
  Mimiviridae sp. ChoanoV1                     57  ← 57 proteins matched to this species
  Fadolivirus algeromassiliense                49  ← Second most abundant virus in sample
  Megaviridae environmental sample             45
  ...

[MetaFlye-Only Viruses (Distant Viral Candidates)]
--------------------------------------------------------------------------------
Count: 528 contigs, 965 protein matches
Average Identity: 34.18%  ← Similar identity to consensus viruses, but filtered by viralFlye

Family Distribution (Top 5):
  Mimiviridae                                 445
  Phycodnaviridae                             143
  ...
```

**Note**: All numbers in distribution statistics are **protein match counts**, not contig counts!

### Viral Classification System Notes

**Note**: NCBI reformed the viral classification system in 2018-2020

- **Old system**: Superkingdom: Viruses → Family → Genus → Species
- **New system**: Kingdom (e.g., Bamfordvirae) → Phylum → Class → Order → Family → Genus → Species

Therefore:
- ✅ **Superkingdom showing N/A is normal** (viruses lack this level in the new system)
- ✅ **Kingdom shows viral kingdoms** (e.g., Bamfordvirae, Orthornavirae)
- ✅ **Phylum/Family/Species provide detailed classification**

---

## Dual-Track Analysis: Complementary Virus Identification Strategy

### 🔬 Methodological Principles

**Important Note**: Dual-track analysis **cannot discover completely unknown new viruses**, but rather **maximizes viral identification coverage** through two complementary methods.

**Limitations**:
- ❌ Diamond: Can only find viruses with protein sequence similarity to RVDB
- ❌ viralFlye: Can only find viruses conforming to known viral genome features
- ❌ Completely unknown viruses (different in both sequence and features): Neither method can identify

**Value of Dual-Track Analysis**:
- ✅ Complementary identification: Expand virus discovery range
- ✅ Distant viruses: Identifiable despite low similarity
- ✅ Feature-atypical viruses: Discoverable through sequence similarity
- ✅ High-purity dataset: viralFlye provides verified viruses

The workflow simultaneously generates two result sets:
- **Track 1 (MetaFlye + Diamond)**: Contains all sequences (1,213 contigs)
  - No filtering, Diamond alignment on all sequences
  - Identifies viruses through **sequence similarity**
  - **Can find**: Sequences with protein similarity to RVDB viruses (even 20-30% similarity)
  - **Cannot find**: Viruses with completely different protein sequences
  
- **Track 2 (viralFlye + Diamond)**: Only viruses identified by viralFlye (222 contigs)
  - Identifies based on **viral genome features** (gene density, Pfam domains, viralVerify scores)
  - **Can find**: Sequences conforming to known viral genome features
  - **Cannot find**: Feature-atypical viruses

**Complementary Value**:
- Potential distant viruses filtered by viralFlye but found by Diamond (feature-atypical but with similar proteins)
- Potential viruses identified by viralFlye but not found by Diamond (features conform but proteins dissimilar)
- Combined provides most comprehensive viral identification

### Analysis Strategy 1: Identify Sequences Filtered by viralFlye but Matched by Diamond

**Purpose**: Discover distant viruses with atypical features but similar proteins

**Principle**:
- viralFlye filters based on features (e.g., gene density, Pfam domains)
- Diamond matches based on sequence similarity
- Filtered by viralFlye but viral matches found by Diamond → distant viral candidates

**Actual Data** (example):
- MetaFlye total contigs: 1,213
- viralFlye viral contigs: 222 (18.3%)
- **Filtered sequences: 991** (81.7%)
  - Most are bacteria, eukaryotes, low-quality sequences
  - A small portion may be feature-atypical but protein-similar distant viruses

```bash
# Extract viralFlye contig IDs
grep ">" results_long/assembly_viralflye/*_viralflye_contigs.fa | \
  sed 's/>//' | cut -d' ' -f1 > viralflye_contig_ids.txt

# Extract MetaFlye contig IDs
grep ">" results_long/assembly_metaflye/*_metaflye_contigs.fa | \
  sed 's/>//' | cut -d' ' -f1 > metaflye_contig_ids.txt

# Find contigs in MetaFlye but not in viralFlye (filtered ~991)
comm -23 <(sort metaflye_contig_ids.txt) <(sort viralflye_contig_ids.txt) > viralflye_filtered_contigs.txt

echo "Number of contigs filtered by viralFlye:"
wc -l viralflye_filtered_contigs.txt

# View Diamond classification results for these filtered contigs
while read contig; do
  grep "^${contig}_" results_long/taxonomy_metaflye/*_with_taxonomy.txt
done < viralflye_filtered_contigs.txt > viralflye_filtered_analysis.txt

# Count viral matches in filtered sequences
grep -c "virae" viralflye_filtered_analysis.txt
# If count > 0, these sequences:
# - Have viral protein similarity (found by Diamond)
# - But don't meet viralFlye viral feature standards
# - Possibly distant viruses or feature-atypical viruses
```

### Analysis Strategy 2: Identify Low-Similarity Viral Proteins (Distant Viruses)

```bash
# Extract matches with identity < 30% but E-value < 1e-10
awk '$3 < 30 && $11 < 1e-10' \
  results_long/taxonomy_metaflye/*_with_taxonomy.txt \
  > low_similarity_high_confidence.txt

# These sequences may be:
# - Distant viruses (low similarity to known viral proteins)
# - New viral genera/species (within known viral families/orders)
# - Highly variable viral proteins
# 
# Note: Still have some similarity to known viruses, not completely unknown
```

### Analysis Strategy 3: Identify Viral Sequences with Incomplete Classification

```bash
# Extract sequences with viral Kingdom but Species = N/A
awk -F'\t' '$16 ~ /virae/ && $22 == "N/A"' \
  results_long/taxonomy_metaflye/*_with_taxonomy.txt \
  > unclassified_viral_kingdom.txt

# These may be:
# - Known viral kingdom but unnamed species
# - New viral species (within known viral kingdoms/phyla)
# - Viruses with incomplete annotation in NCBI database
```

### Analysis Strategy 4: Compare Both Result Sets to Find Complementary Discoveries

```bash
# Count contigs in both sets
echo "MetaFlye total contigs:"
grep -c ">" results_long/assembly_metaflye/*_contigs.fa
# Expected: ~1,200+

echo "viralFlye viral contigs:"
grep -c ">" results_long/assembly_viralflye/*_contigs.fa  
# Expected: ~200+ (if counts are the same, sequence extraction failed, need to install seqtk)

# Count Diamond matches in both sets
echo "MetaFlye Diamond matches:"
wc -l results_long/diamond_metaflye/*_diamond.txt

echo "viralFlye Diamond matches:"
wc -l results_long/diamond_viralflye/*_diamond.txt

# Compare viral proportions
echo "Viral proportion in MetaFlye:"
awk -F'\t' '$16 ~ /virae/ || $17 ~ /viricota/' \
  results_long/taxonomy_metaflye/*_with_taxonomy.txt | wc -l

echo "Viral proportion in viralFlye:"
awk -F'\t' '$16 ~ /virae/ || $17 ~ /viricota/' \
  results_long/taxonomy_viralflye/*_with_taxonomy.txt | wc -l
```

**Expected Result Interpretation**:
- If viralFlye contigs significantly fewer than MetaFlye (e.g., 222 vs 1,213)
  → ✅ viralFlye successfully filtered, difference (~991) includes filtered sequences
- If counts are the same
  → ⚠️ Sequence extraction failed, need to install seqtk or biopython

**Actual Run Data Example**:
```
MetaFlye: 1,213 contigs → 2,434 protein matches → 709 viral contigs (1,809 viral proteins)
viralFlye: 222 contigs → 1,129 protein matches → 181 viral contigs (844 viral proteins)
Consensus viruses: 181 contigs (844 proteins) - confirmed by both methods
MetaFlye-only: 528 contigs (965 proteins) - distant viral candidates
```

**Key Insights**:
- viralFlye retained 18% of contigs but 46% of protein matches
- Indicates viralFlye selects sequences with **high gene density**
- Consensus virus gene density: 4.7 proteins/contig
- MetaFlye-only gene density: 1.8 proteins/contig

### Distant Virus/New Viral Species Determination Criteria

#### 🌟 Strong Candidates (Distant Viruses or New Viral Genera):
- ✅ Detected in MetaFlye results (Track 1)
- ✅ **Filtered by viralFlye** (not in 222 virus list)
- ✅ But matched to viral proteins by Diamond (proves it's viral)
- ✅ Low similarity match (< 30%, indicates distant or new virus)
- ✅ Significant E-value (< 1e-10, statistically reliable)
- ✅ Kingdom belongs to viral kingdom but Species is N/A

**Explanation**:
- viralFlye cannot identify feature-atypical viruses
- But if the virus has proteins similar to known viruses
- Diamond may still discover through low similarity matching
- In this case, sequences filtered by viralFlye but matched by Diamond merit in-depth study

**Important Note**:
- Completely unknown viruses (different in both protein sequence and features) cannot be identified by either method
- Dual-track analysis maximizes identification range but cannot guarantee finding all new viruses
- Needs combination with other methods (e.g., genome structure analysis, host signals)

#### ⭐ Medium Candidates (New Viral Species):
- ✅ Detected in both result sets
- ✅ All protein similarities < 40%
- ✅ Belongs to known viral family but species name is N/A

#### ⚠️ Requires Further Verification:
- ⚠️ Only detected in MetaFlye
- ⚠️ Filtered by viralFlye
- ⚠️ No Diamond match or E-value > 1e-5
- ⚠️ Possibly non-viral sequence or low-quality assembly

### Follow-up Steps for Validating New Viral Candidates

1. **Sequence Quality Check**: Use CheckV to assess genome completeness
2. **Genome Structure Analysis**:
   - Genome size, gene density
   - Genome terminal structures (ITR, DTR, etc.)
   - Presence of virus-specific gene arrangement patterns
3. **Phylogenetic Analysis**: Construct phylogenetic tree to confirm taxonomic position
4. **Host Signal Analysis**:
   - CRISPR spacer matching
   - Integration site analysis
   - Similarity to host genomes
5. **Functional Annotation**:
   - Complete annotation of all proteins
   - Identify conserved viral replication/packaging proteins
6. **Experimental Validation**: PCR validation, cell culture, metatranscriptome validation, etc.

**Note**: For sequences not found by either Diamond or viralFlye, more evidence is needed to determine if they are viral.

---

## Workflow Architecture

### Short-Read Workflow Diagram

```
Illumina Paired-end Reads
         ↓
    fastp Quality Control
         ↓
    ┌────┴────┐
    ↓         ↓
 MEGAHIT   SPAdes
    ↓         ↓
Prodigal  Prodigal
    ↓         ↓
 Diamond   Diamond
    └────┬────┘
         ↓
    Comprehensive Analysis Report
         ↓
  Complete Taxonomic Information
```

### Long-Read Workflow Diagram (Dual-Track Analysis)

```
Nanopore/PacBio Long Reads
         ↓
    MetaFlye Assembly
    (--meta mode, all sequences)
         ↓
    ┌────┴─────────────────┐
    ↓                      ↓
[Track 1: All Sequences] [Track 2: Viral Only]
    ↓                      ↓
Prodigal               viralFlye
(MetaFlye)          (Viral Sequence ID)
All genes                  ↓
    ↓                  Prodigal
Diamond             (viralFlye)
(All sequence       Viral genes
classification)         ↓
    ↓                  Diamond
Taxonomy            (Viral
resolution          classification)
    ↓                      ↓
Comprehensive           Taxonomy
analysis                resolution
No omissions               ↓
                    High confidence
                    viral sequences
```

**Notes**:
- Dual-track analysis provides two complementary identification methods (feature-based vs sequence similarity)
- Both depend on known viral information, can only identify sequences related to known viruses
- Completely unknown new viruses (different in both sequence and features) require other identification methods

---

## Key Configuration Files

### metagenome_assembly_classification_en.config

Main configuration items:
- **executor**: SLURM
- **queue**: bahl_p
- **Resource allocation**: CPU/memory/time configuration for different processes
- **Database paths**: Diamond, Pfam-A, Taxonomy
- **Parameter settings**: Assembly, classification, filtering parameters

### Custom Configuration

Edit `metagenome_assembly_classification_en.config` file to adjust:
- Computational resource allocation
- Database paths
- Analysis parameters

---

## Troubleshooting

### Issue 1: viralFlye Environment Activation Failed

**Error**: `viralFlye.py not found`

**Solution**: Workflow is configured to use absolute paths, no manual intervention needed. For debugging:

```bash
# Verify viralFlye installation
conda activate viralFlye_env
which viralFlye.py
viralFlye.py --help
```

### Issue 1b: viralFlye Identified Viruses but Extraction Failed

**Symptom**: `llnl_66d1047e_viralflye_contigs.fa` identical to `llnl_66d1047e_metaflye_contigs.fa`

**Cause**: Missing `seqtk` or `biopython` tools for sequence extraction

**Solution**:
```bash
# Install sequence extraction tools in viralFlye_env
conda activate viralFlye_env
conda install -c bioconda seqtk
# Or
conda install -c conda-forge biopython

# Re-run
rm -rf work/ .nextflow/
sbatch run_metagenome_assembly_classification_en.sh long
```

**Verify Fix**:
```bash
# Compare contig counts in both files
grep -c ">" results_long/assembly_metaflye/*_contigs.fa    # Should be ~1,200
grep -c ">" results_long/assembly_viralflye/*_contigs.fa   # Should be ~200

# If counts differ significantly, extraction was successful
```

### Issue 2: SPAdes Out of Memory

**Error**: `SPAdes killed - out of memory`

**Solution**: Increase memory limit or use MEGAHIT results

```groovy
// In config file
withName: 'SPADES_ASSEMBLY' {
    memory = '1024 GB'  // Increase memory
}
```

### Issue 3: Taxonomy Database Missing

**Error**: `NCBI Taxonomy database not found`

**Solution**: Download taxonomy database (see Installation and Configuration section)

### Issue 4: Classification Results Show Many N/A

**Cause**:
- Some sequences in RVDB database lack taxonomy IDs
- Normal phenomenon, typically 20-30% of sequences may lack complete classification information

**Solution**: Focus on classified results, N/A may be new viruses or incompletely annotated sequences

---

## Performance Optimization Suggestions

### Short-Read Optimization
1. **Skip Quality Control** (if data already filtered): `--skip_fastp true`
2. **Use Only One Assembler**: Comment out unneeded assembler in workflow
3. **Increase MEGAHIT Memory**: Improve assembly quality

### Long-Read Optimization
1. **Set Genome Size**: `--metaflye_genome_size 10m` speeds up assembly
2. **Dual-Track Analysis Recommendations**:
   - For new virus discovery: Keep dual-track analysis (default)
   - Only analyzing known viruses: Can skip MetaFlye track (requires workflow modification)
   - Don't need virus-specific: `--skip_viralflye true` (MetaFlye track only)
3. **Adjust Minimum Length**: `--viralflye_min_length 3000` retains more short viruses
4. **Parallelization**: Prodigal and Diamond in dual-track automatically run in parallel

---

## Frequently Asked Questions (FAQ)

### Q1: Is viralFlye required?
**A**: No. viralFlye is mainly used for viral research, providing additional viral sequence optimization. If not specifically studying viruses, this step can be skipped.

### Q2: Why are two assemblers needed (MEGAHIT and SPAdes)?
**A**: Different assemblers have different advantages, dual assemblers can:
- Obtain more comprehensive results
- Mutually verify assembly quality
- Discover more viral sequences

### Q3: Can short-read and long-read data be analyzed simultaneously?
**A**: Cannot analyze both in a single run. Need to run separately:
```bash
sbatch run_metagenome_assembly_classification_en.sh short
sbatch run_metagenome_assembly_classification_en.sh long
```

### Q4: Can Diamond database be replaced with other databases?
**A**: Yes. Any Diamond-format protein database can be used, such as:
- NCBI nr
- UniProt
- Custom viral databases

### Q5: Will output directories be automatically distinguished?
**A**: Yes. Workflow automatically sets:
- Short reads → `results_short/`
- Long reads → `results_long/`

### Q6: How much computational time does dual-track analysis add?
**A**: About 40-60% increase, but viralFlye track is faster (fewer sequences):
- MetaFlye assembly: 1 time (shared)
- viralFlye viral identification: ~1 time (~3 minutes)
- Prodigal: 2 times (MetaFlye ~3 minutes + viralFlye ~30 seconds)
- Diamond: 2 times (MetaFlye ~1 minute + viralFlye ~10 seconds)
- Taxonomy: 2 times (each ~20 seconds)
- Consensus analysis: 1 time (~20 seconds)
- **Total additional time**: about 6-10 minutes

**Example Run Time** (1 sample, ~1,200 contigs):
- Single-track (MetaFlye only): ~10 minutes
- Dual-track analysis + consensus: ~16 minutes

**Computational Resource Consumption** (dual-track):
- MetaFlye track: 2,434 protein alignments
- viralFlye track: 1,129 protein alignments (54% less)
- Total computation increases by about 46%

### Q7: Will viralFlye filter out new viruses?
**A**: Yes. Based on actual run data:
- **MetaFlye assembly**: ~1,213 contigs (all sequences)
- **viralFlye identification**: ~222 viral contigs (18.3%)
- **Filtered**: ~991 contigs (81.7%)

Filtered sequences may include:
- ✅ Bacteria, eukaryotes, host sequences (normal filtering)
- ⚠️ Completely new viruses (don't conform to known viral features)
- ⚠️ Distant viruses (atypical features)

**This is why we use dual-track analysis (complementary identification)**:

| Method | Identification Principle | Can Find | Cannot Find |
|--------|--------------------------|----------|-------------|
| **Diamond** | Sequence similarity | Viruses with similar proteins to known viruses | Viruses with completely different protein sequences |
| **viralFlye** | Viral features | Viruses conforming to viral genome features | Feature-atypical viruses |

**Complementarity Examples**:
- **Case A**: Distant virus (proteins have low similarity but features atypical)
  - Diamond ✅ Found (low similarity match)
  - viralFlye ❌ Filtered (features don't conform)
  - **Conclusion**: Discovered in MetaFlye track

- **Case B**: New virus (features conform but sequences dissimilar)
  - Diamond ❌ Cannot find (no similar sequences)
  - viralFlye ✅ Found (features conform)
  - **Conclusion**: Discovered in viralFlye track

- **Case C**: Completely unknown new virus (different in both features and sequences)
  - Diamond ❌ Cannot find
  - viralFlye ❌ Cannot find
  - **Conclusion**: Neither method can identify

**Dual-track analysis provides maximum viral discovery range, but still has limitations.**

### Q8: What do the numbers in classification statistics mean?
**A**: They are **protein match counts**, not contig counts.

**Example**:
```
Species Distribution:
  Mimiviridae sp. ChoanoV1     57  ← 57 protein matches
```

**Meaning**:
- 57 proteins were matched by Diamond to this viral species
- These proteins may come from 5-10 different contigs
- **Larger numbers** → Higher viral abundance or more complete assembly

**Relationship to Contigs**:
```
Consensus viruses: 181 contigs contain 844 viral proteins
  ↓
Average 4.7 protein matches per contig
  ↓
Grouped by Species:
  - Mimiviridae sp. ChoanoV1: 57 proteins (possibly from ~10 contigs)
  - Fadolivirus algeromassiliense: 49 proteins (possibly from ~10 contigs)
  - ...
```

### Q9: How to determine if a sequence is a distant virus or new viral species?
**A**: Comprehensive consideration of multiple indicators:

1. **Filtered by viralFlye but matched by Diamond**
   - Has viral protein similarity but atypical features
   - Possibly distant virus or new viral genus

2. **Diamond low similarity (20-30%) but significant E-value (< 1e-10)**
   - Some similarity to known viral proteins but large differences
   - Possibly new genus/species within known viral family/order

3. **Kingdom is viral kingdom but Species is N/A**
   - Belongs to known viral classification system but species unnamed
   - Possibly new viral species

4. **Identified by viralFlye but no Diamond match or low similarity**
   - Has viral features but proteins dissimilar
   - Possibly new viral genus

**Important**:
- These are all "new" relative to known viruses (new genus, new species, distant)
- Completely unknown viruses (no similarity, no features) require other identification methods
- Ultimately requires phylogenetic analysis and experimental validation

---

## Citations and Acknowledgments

### Tools Used
- **Nextflow**: Di Tommaso et al. (2017) Nature Biotechnology
- **MEGAHIT**: Li et al. (2015) Bioinformatics
- **SPAdes**: Bankevich et al. (2012) Journal of Computational Biology
- **MetaFlye**: Kolmogorov et al. (2020) Nature Methods
- **viralFlye**: Antipov et al. (2020)
- **Prodigal**: Hyatt et al. (2010) BMC Bioinformatics
- **Diamond**: Buchfink et al. (2021) Nature Methods
- **fastp**: Chen et al. (2018) Bioinformatics

### Databases
- **RVDB**: Reference Viral DataBase - https://rvdb-prot.pasteur.fr/
- **NCBI Taxonomy**: https://www.ncbi.nlm.nih.gov/taxonomy
- **Pfam**: https://pfam.xfam.org/

---

## Contact

- Email: sihua.peng@uga.edu, Workflow code programmer  
- Email: justin.bahl@uga.edu, Project supervisor  
- Suggestion: [Click here!](https://github.com/pengsihua2023/rvdb-viral-metagenome-nf/issues/new)

---


## License

MIT License

---

## Appendix

### File List

```
.
├── metagenome_assembly_classification_workflow_en.nf  # Main workflow file
├── metagenome_assembly_classification_en.config       # Configuration file
├── run_metagenome_assembly_classification_en.sh       # SLURM submission script
├── samplesheet_short.csv                              # Short-read sample sheet
└── samplesheet_long.csv                               # Long-read sample sheet
```

### Recommended Directory Structure

```
project_root/
├── data/
│   ├── short_reads/
│   └── long_reads/
├── databases/
│   ├── RVDB/
│   │   ├── RVDB_prot_ref.dmnd
│   │   ├── names.dmp
│   │   └── nodes.dmp
│   └── Pfam/
│       └── Pfam-A.hmm
├── results_short/
├── results_long/
└── work/  (Nextflow working directory)
```









