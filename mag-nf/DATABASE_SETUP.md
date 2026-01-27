# Viral Metagenome Database Setup Guide

## Recommended Viral Classification Databases

### 1. Kraken2 Viral Database Options

#### Option A: Standard Viral Database (Recommended for virus-specific studies)

```bash
# Download and build standard viral database (~465 MB, fast)

DB_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_viral"
mkdir -p $DB_DIR

# Download viral sequence library
kraken2-build --download-library viral --db $DB_DIR --threads 32

# Download NCBI taxonomy
kraken2-build --download-taxonomy --db $DB_DIR

# Build database
kraken2-build --build --db $DB_DIR --threads 32

# Clean intermediate files
kraken2-build --clean --db $DB_DIR
```

**Advantages:**
- ✅ Small size, fast build (~10-20 minutes)
- ✅ Focused on viral classification with high accuracy
- ✅ Contains NCBI RefSeq viral genomes

**Disadvantages:**
- ⚠️ Viruses only, cannot identify host or bacterial contamination


#### Option B: PlusPFP Database (Recommended for comprehensive analysis)

```bash
# Download and build PlusPFP database (~60 GB, more comprehensive)

DB_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/k2_pluspfp"
mkdir -p $DB_DIR

# Download multiple taxonomic groups (viral, bacterial, archaeal, eukaryotic, plasmid)
kraken2-build --download-library viral --db $DB_DIR --threads 32
kraken2-build --download-library bacteria --db $DB_DIR --threads 32
kraken2-build --download-library archaea --db $DB_DIR --threads 32
kraken2-build --download-library plasmid --db $DB_DIR --threads 32
kraken2-build --download-library human --db $DB_DIR --threads 32
kraken2-build --download-library fungi --db $DB_DIR --threads 32
kraken2-build --download-library protozoa --db $DB_DIR --threads 32

# Download NCBI taxonomy
kraken2-build --download-taxonomy --db $DB_DIR

# Build database (requires ~150 GB RAM)
kraken2-build --build --db $DB_DIR --threads 32 --max-db-size 60000000000

# Clean intermediate files
kraken2-build --clean --db $DB_DIR
```

**Advantages:**
- ✅ Comprehensive coverage of viruses, bacteria, archaea, eukaryotes
- ✅ Can identify hosts, contamination, symbiotic microbes
- ✅ Suitable for complex sample analysis

**Disadvantages:**
- ⚠️ Large size (~60 GB)
- ⚠️ Time-consuming build (2-4 hours)
- ⚠️ Requires high memory (~150 GB)


#### Option C: Custom Viral Database (Recommended for specific viral research)

```bash
# Download specific viral genomes from IMG/VR or RVDB

DB_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_custom_viral"
mkdir -p $DB_DIR/library/added

# 1. Download RVDB (Reference Viral DataBase)
wget https://hive.biochemistry.gwu.edu/dna.cgi?cmd=objFile&ids=1928374&filename=U-RVDBv27.0-prot.fasta.gz
gunzip U-RVDBv27.0-prot.fasta.gz

# 2. Add sequences to Kraken2 database
kraken2-build --add-to-library U-RVDBv27.0-prot.fasta --db $DB_DIR

# 3. Or download specific viral families from NCBI (e.g., Herpesviridae)
# Visit NCBI Genome: https://www.ncbi.nlm.nih.gov/genome/viruses/
# Download target viral genome FASTA files
# kraken2-build --add-to-library herpesvirus_genomes.fasta --db $DB_DIR

# 4. Download taxonomy
kraken2-build --download-taxonomy --db $DB_DIR

# 5. Build database
kraken2-build --build --db $DB_DIR --threads 32

# 6. Clean
kraken2-build --clean --db $DB_DIR
```

**Advantages:**
- ✅ Highly customizable for specific research goals
- ✅ Can include latest viral genomes
- ✅ Controllable size

**Disadvantages:**
- ⚠️ Requires manual sequence collection
- ⚠️ Coverage depends on input data


### 2. Other Recommended Viral Databases

#### IMG/VR (Integrated Microbial Genomes/Viral Resources)
```bash
# IMG/VR: Largest viral metagenome database

# Access: https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html
# Download complete database or environment-specific viral sequences
```

**Features:**
- Over 10 million viral sequences
- Contains uncultivated viral genomes
- Suitable for discovering novel viruses


#### RVDB (Reference Viral DataBase)
```bash
# RVDB: Curated reference viral database

wget https://hive.biochemistry.gwu.edu/dna.cgi?cmd=objFile&ids=1928374&filename=U-RVDBv27.0-prot.fasta.gz
```

**Features:**
- High-quality viral protein sequences
- Regular updates (2-3 times per year)
- Suitable for functional annotation


#### NCBI RefSeq Viral
```bash
# NCBI RefSeq: Standard reference viral genomes

# Already included in Kraken2 --download-library viral
```

**Features:**
- Manually curated reference genomes
- High-quality annotation
- Monthly updates


### 3. Special Considerations for Long Reads

```bash
# For Nanopore/PacBio long reads, recommend:

# 1. Use larger k-mer (k=35)

# 2. Lower confidence threshold (0.01-0.05)

# 3. Consider using specialized long-read viral classifiers:
#    - Centrifuge
#    - Kaiju
#    - MMseqs2
```


## Database Maintenance

### Regular Updates

```bash
# Update Kraken2 database every 3-6 months

DB_DIR="/path/to/kraken2/db"

# Backup old database
cp -r $DB_DIR ${DB_DIR}_backup_$(date +%Y%m%d)

# Re-download and rebuild
kraken2-build --download-library viral --db $DB_DIR --threads 32
kraken2-build --download-taxonomy --db $DB_DIR --no-masking
kraken2-build --build --db $DB_DIR --threads 32
kraken2-build --clean --db $DB_DIR
```


## Recommended Configuration

### Virus-Specific Studies

```bash
# Use Standard Viral database + custom sequences

DB="/scratch/sp96859/.../kraken2_viral"
```

### Environmental Sample Analysis

```bash
# Use PlusPFP database (comprehensive coverage)

DB="/scratch/sp96859/.../k2_pluspfp"
```

### Clinical Viral Diagnostics

```bash
# Use Standard Viral + human genome (remove host sequences)

DB="/scratch/sp96859/.../kraken2_viral_clinical"
```


## Database Size and Performance Comparison

| Database | Size | Build Time | RAM | Speed | Use Case |
|----------|------|------------|-----|-------|----------|
| Viral | ~500 MB | 10-20 min | ~8 GB | Very Fast | Viral-specific |
| PlusPFP | ~60 GB | 2-4 hours | ~150 GB | Fast | Comprehensive |
| Custom | Variable | 30-60 min | ~8-32 GB | Fast | Targeted |
| IMG/VR | ~200 GB | - | ~200 GB | Medium | Viral discovery |


## Validate Database

```bash
# Test if database works correctly

DB="/path/to/kraken2/db"
TEST_READS="test_reads.fastq"

kraken2 --db $DB \
    --threads 8 \
    --report test_report.txt \
    --output test_classification.txt \
    $TEST_READS

# Check report
head -20 test_report.txt
```


## Resource Links

- **Kraken2 Official Documentation**: https://github.com/DerrickWood/kraken2/wiki
- **IMG/VR**: https://img.jgi.doe.gov/vr/
- **RVDB**: https://hive.biochemistry.gwu.edu/rvdb/
- **NCBI Virus**: https://www.ncbi.nlm.nih.gov/labs/virus/
- **nf-core/mag**: https://nf-co.re/mag/


## Quick Start Examples

### Example 1: Build Standard Viral Database (Fastest)

```bash
#!/bin/bash
# Quick viral database setup (10-20 minutes)

DB_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_viral"

mkdir -p $DB_DIR
kraken2-build --download-library viral --db $DB_DIR --threads 32
kraken2-build --download-taxonomy --db $DB_DIR
kraken2-build --build --db $DB_DIR --threads 32
kraken2-build --clean --db $DB_DIR

echo "✅ Viral database ready: $DB_DIR"
```

### Example 2: Build PlusPFP for Comprehensive Analysis

```bash
#!/bin/bash
# Comprehensive database with viruses, bacteria, archaea (2-4 hours)

DB_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/k2_pluspfp"

mkdir -p $DB_DIR

# Download all libraries
for LIB in viral bacteria archaea plasmid human fungi protozoa; do
    echo "Downloading $LIB library..."
    kraken2-build --download-library $LIB --db $DB_DIR --threads 32
done

# Build database
kraken2-build --download-taxonomy --db $DB_DIR
kraken2-build --build --db $DB_DIR --threads 32 --max-db-size 60000000000
kraken2-build --clean --db $DB_DIR

echo "✅ PlusPFP database ready: $DB_DIR"
```

### Example 3: Add Custom Viral Genomes

```bash
#!/bin/bash
# Add custom viral sequences to existing database

DB_DIR="/path/to/kraken2/db"
CUSTOM_FASTA="my_viral_genomes.fasta"

# Add custom sequences
kraken2-build --add-to-library $CUSTOM_FASTA --db $DB_DIR

# Rebuild database
kraken2-build --build --db $DB_DIR --threads 32
kraken2-build --clean --db $DB_DIR

echo "✅ Custom sequences added to database"
```
