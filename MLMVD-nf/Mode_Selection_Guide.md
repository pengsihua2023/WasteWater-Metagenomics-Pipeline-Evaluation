# Workflow Mode Selection Guide

## Quick Selection

**What data do you have?**

```
Illumina paired-end sequencing (150-300 bp)
    ↓
Use【Short-read Mode】

Nanopore / PacBio long-read (>1 kb)
    ↓
Use【Long-read Mode】+ Enable viralFlye ⭐
```

---

## Detailed Comparison of Two Modes

### Short-read Mode (Illumina)

#### Design Strategy: Dual Assemblers in Parallel

```
Illumina Reads
     ↓
  ┌──┴──┐
  ↓     ↓
MEGAHIT SPAdes
  ↓     ↓
  ├─VS2 ├─VS2
  └─DVF └─DVF
  ↓     ↓
  └──┬──┘
     ↓
Assembler Comparison
 Consensus Viruses
```

#### Core Features

| Feature | Content |
|---------|---------|
| **Assemblers** | MEGAHIT + metaSPAdes (dual assemblers in parallel) |
| **Viral Identification** | VirSorter2 + DeepVirFinder (run on each assembler) |
| **Comprehensive Comparison** | Assembler comparison (MEGAHIT vs SPAdes) |
| **Consensus Strategy** | Viruses identified by both assemblers = High confidence |
| **viralFlye** | ❌ **Not supported** (viralFlye only supports long-read) |

#### Advantages

- ✅ Dual-assembler cross-validation
- ✅ Identify assembler consensus viruses
- ✅ Compare MEGAHIT and SPAdes performance
- ✅ Suitable for Illumina data

#### Output Example

```
results/
├── assembly_megahit/             # MEGAHIT assembly
├── assembly_spades/              # SPAdes assembly
├── virsorter2_megahit/           # VS2 (MEGAHIT)
├── virsorter2_spades/            # VS2 (SPAdes)
├── deepvirfinder_megahit/        # DVF (MEGAHIT)
├── deepvirfinder_spades/         # DVF (SPAdes)
├── merged_viral_reports_megahit/ # Merged reports (MEGAHIT)
├── merged_viral_reports_spades/  # Merged reports (SPAdes)
├── assembler_comparison/         # Assembler comparison ⭐
│   ├── sample_assembler_comparison.txt
│   └── sample_consensus_viral_sequences.txt  # Dual-assembler consensus
└── abundance/                    # Viral abundance analysis (RPM & RPKM) ⭐⭐⭐
    ├── megahit/
    │   ├── sample_megahit_abundance.csv
    │   └── sample_megahit_abundance_summary.txt
    └── spades/
        ├── sample_spades_abundance.csv
        └── sample_spades_abundance_summary.txt
```

---

### Long-read Mode (Nanopore/PacBio) ⭐

#### Design Strategy: Three-Tool Parallel Analysis

```
Nanopore/PacBio Reads
         ↓
    metaFlye Assembly
         ↓
    All contigs
         ↓
   ┌─────┼─────┐
   ↓     ↓     ↓
  VS2   DVF  viralFlye
 Viral  ML   Pfam
Feature Learn Validation
   |     |     |
   └─────┼─────┘
         ↓
  Three-Tool Comparison
   - 3-tool consensus ⭐⭐⭐
   - 2-tool consensus ⭐⭐
   - Single-tool ⭐
```

#### Core Features

| Feature | Content |
|---------|---------|
| **Assembler** | metaFlye (single assembler, but saves complete output) |
| **Viral Identification** | VirSorter2 + DeepVirFinder + **viralFlye** (three tools parallel) ⭐ |
| **Comprehensive Comparison** | **Three-tool comparison** (VS2 vs DVF vs viralFlye) ⭐⭐⭐ |
| **Consensus Strategy** | 1-3 tool consensus stratification (fine-grained confidence levels) |
| **viralFlye** | ✅ **Supported** (Pfam protein validation) ⭐ |

#### Advantages

- ✅ **Three independent methods** validation
- ✅ **Pfam protein domain** validation (viralFlye exclusive)
- ✅ **Confidence stratification** (1/2/3 tool consensus)
- ✅ Can identify **novel viruses** (viralFlye function-based, not sequence similarity-dependent)
- ✅ Assess viral **genome completeness** (viralComplete)
- ✅ Save **complete metaFlye output** (assembly_info.txt, assembly_graph.gfa)

#### Output Example

```
results_long/
├── assembly_metaflye/            # metaFlye contigs
├── metaflye_full_output/         # metaFlye complete output ⭐
│   └── sample_flye_output/
│       ├── assembly_info.txt     # Assembly statistics
│       └── assembly_graph.gfa    # Assembly graph
├── virsorter2_metaflye/          # VS2 identification
├── deepvirfinder_metaflye/       # DVF identification
├── viralflye_results/            # viralFlye identification ⭐
│   ├── sample_viralflye_contigs.fa
│   └── sample_viralflye_summary.csv
├── three_tools_comparison/       # Three-tool comparison ⭐⭐⭐
│   ├── sample_three_tools_comparison.txt
│   ├── sample_three_tools_comparison.csv
│   └── sample_high_confidence_viruses.txt
└── abundance/                    # Viral abundance analysis (RPM & RPKM) ⭐⭐⭐
    ├── metaflye/
    │   ├── sample_metaflye_abundance.csv
    │   └── sample_metaflye_abundance_summary.txt
    └── viralflye/
        ├── sample_viralflye_abundance.csv
        └── sample_viralflye_abundance_summary.txt
```

---

## Detailed Comparison Table

| Comparison Item | Short-read Mode | Long-read Mode |
|----------------|----------------|----------------|
| **Data Type** | Illumina paired-end | Nanopore/PacBio single-end |
| **Read Length** | 150-300 bp | 1-100 kb |
| **Number of Assemblers** | 2 (MEGAHIT + SPAdes) | 1 (metaFlye) |
| **Viral ID Tools** | 2 (VS2 + DVF) | 3 (VS2 + DVF + viralFlye) ⭐ |
| **Signature Tool** | - | **viralFlye** (Pfam validation) ⭐ |
| **Comparison Strategy** | Assembler comparison | Tool comparison |
| **Confidence Levels** | 2 levels (dual-tool, single-tool) | 4 levels (3-tool, 2-tool, 1-tool, viralFlye-unique) ⭐ |
| **Pfam Database** | ❌ Not required | ✅ Required |
| **Computation Time** | 24-48h (medium data) | 48-72h (medium data) |
| **Memory Requirements** | 64-512 GB (SPAdes) | 128-256 GB |
| **Virus Quality** | Medium-High | High (viralFlye validated) ⭐ |
| **Novel Virus Discovery** | Medium | Strong (viralFlye function-based) ⭐ |
| **Viral Abundance Analysis** | ✅ RPM & RPKM | ✅ RPM & RPKM |

---

## How to Choose a Mode?

### Scenario 1: Only Illumina Data

```
Sequencing Data: Illumina HiSeq/NovaSeq
Read Type: Paired-end, 150 bp

Choice: Short-read mode
Command: sbatch run_metagenome_assembly_classification_shortread.sh

Obtain:
- Dual-assembler results
- Dual-tool viral identification
- Assembler consensus viruses
```

### Scenario 2: Only Long-read Data

```
Sequencing Data: Nanopore MinION/PromethION
Read Type: Single-end, average 5-10 kb

Choice: Long-read mode + Enable viralFlye ⭐
Command: sbatch run_metagenome_assembly_classification_longread.sh
      (Ensure ENABLE_VIRALFLYE="true")

Obtain:
- metaFlye assembly
- Three-tool viral identification
- Pfam protein validation
- Three-tool comprehensive comparison ⭐⭐⭐
```

### Scenario 3: Both Short-read and Long-read

```
Sequencing Data: Illumina + Nanopore (hybrid sequencing)

Choice: Run twice
1. Short-read mode (Illumina)
2. Long-read mode (Nanopore + viralFlye)

Then manually compare results from both modes
```

---

## Why Does viralFlye Only Support Long-read?

### Technical Reasons

1. **Requires Flye Output Structure**
   - viralFlye is designed to process Flye/metaFlye output directory
   - Requires assembly_graph.gfa (assembly graph)
   - Requires assembly_info.txt (assembly statistics)

2. **Flye Only Supports Long-read**
   - Flye is a long-read-specific assembler
   - Does not support short-read data
   - MEGAHIT/SPAdes output structure is completely different

3. **Long-read Advantages**
   - Easier to obtain complete viral genomes
   - viralFlye more effective (requires multiple protein-coding genes)
   - Pfam validation needs sufficiently long sequences

### Alternative for Short-read Data

**Short-read mode already has excellent validation strategies**:

1. **Dual-assembler Consensus** ⭐
   - Viruses identified by both MEGAHIT and SPAdes
   - Very reliable

2. **Dual-tool Consensus**
   - Identified by both VirSorter2 and DeepVirFinder
   - High confidence

3. **Tool Combination**
   - VS2 + DVF combination is already powerful
   - Covers different types of viruses

**Conclusion**: Although short-read mode lacks viralFlye, the combination of dual-assemblers + dual-tools already provides reliable viral identification.

---

## Recommended Configurations

### Best Configuration for Short-read Mode

```groovy
// Dual assemblers + Dual tools
params {
    skip_fastp = false              // Enable QC
    
    // VirSorter2
    virsorter2_min_score = 0.5
    virsorter2_min_length = 1000
    
    // DeepVirFinder
    deepvirfinder_pvalue = 0.05     // ⭐ Recommended 0.05 (higher sensitivity for novel viruses)
    deepvirfinder_min_length = 1000
}

Focus on:
- assembler_comparison/: Assembler comparison ⭐⭐⭐
- consensus_viral_sequences.txt: Dual-assembler consensus viruses
```

### Best Configuration for Long-read Mode (Three Tools) ⭐

```groovy
// Three-tool parallel analysis
params {
    longread = true
    longread_platform = 'nano'      // or 'pacbio'
    enable_viralflye = true         // ⭐ Enable third tool
    
    // viralFlye (Pfam validation) ⭐⭐⭐
    viralflye_min_length = 500          // ⭐ Optimized: identify shorter viruses
    viralflye_completeness = 0.3        // ⭐⭐⭐ Completeness threshold (KEY parameter!)
    // Explanation:
    //   0.5 (default) = Only highly complete viruses (2-5)
    //   0.3 (recommended) = Balance quality and quantity (10-20) ⭐⭐
    //   0.2 (aggressive) = Identify more viruses (20-40)
    pfam_db = '/path/to/Pfam-A.hmm'
    
    // VirSorter2
    virsorter2_min_score = 0.5
    virsorter2_min_length = 1000
    
    // DeepVirFinder
    deepvirfinder_pvalue = 0.05     // ⭐ High sensitivity for novel viruses
    deepvirfinder_min_length = 1000
}

Focus on:
- three_tools_comparison/: Three-tool comprehensive comparison ⭐⭐⭐
- high_confidence_viruses.txt: Virus list stratified by consensus
- viralflye_results/: Pfam-validated viruses (significantly increased) ⭐
```

---

## Feature Comparison

### Short-read Mode Features

✅ **Available Features**:
- fastp QC
- MEGAHIT assembly
- metaSPAdes assembly
- VirSorter2 viral identification (MEGAHIT)
- VirSorter2 viral identification (SPAdes)
- DeepVirFinder identification (MEGAHIT)
- DeepVirFinder identification (SPAdes)
- Dual-assembler result integration
- Assembler comparison and consensus

❌ **Unavailable Features**:
- viralFlye viral identification
- Pfam protein validation
- Three-tool comprehensive comparison
- viralComplete completeness assessment

### Long-read Mode Features

✅ **Available Features**:
- metaFlye assembly
- Save complete metaFlye output (assembly_info, assembly_graph) ⭐
- VirSorter2 viral identification
- DeepVirFinder identification
- **viralFlye viral identification** (Pfam validation) ⭐
- **Three-tool comprehensive comparison** (VS2 vs DVF vs viralFlye) ⭐⭐⭐
- **Confidence stratification** (1-3 tool consensus)
- viralComplete completeness assessment
- viralVerify virus/bacteria/plasmid classification

❌ **Unavailable Features**:
- Dual-assembler comparison (only one assembler)

---

## Result Quality Comparison

### Short-read Mode

**Validation Strategy**: Dual-assemblers + Dual-tools

```
Dual-assembler Consensus Viruses:
- MEGAHIT identified ∩ SPAdes identified
- Confidence: High ⭐⭐
- Quantity: Medium

Dual-tool Consensus Viruses:
- VirSorter2 ∩ DeepVirFinder
- Confidence: High ⭐⭐
- Quantity: For each assembler

Highest Confidence:
- Dual-assembler ∩ Dual-tool
- Confidence: Very High ⭐⭐⭐
```

### Long-read Mode (viralFlye Enabled)

**Validation Strategy**: Three-tool parallel

```
Three-tool Consensus Viruses:
- VS2 ∩ DVF ∩ viralFlye
- Confidence: Highest ⭐⭐⭐
- Quantity: 5-15

Two-tool Consensus Viruses:
- VS2 ∩ DVF, VS2 ∩ viralFlye, DVF ∩ viralFlye
- Confidence: High ⭐⭐
- Quantity: 20-50

viralFlye Unique:
- Pfam protein validated
- Confidence: High (high specificity) ⭐⭐
- Quantity: May have viralFlye-exclusive discoveries

Single-tool Identification:
- Only one tool identified
- Confidence: Medium-Low ⭐
- Quantity: More (exploratory)
```

---

## Unique Value of viralFlye (Long-read Exclusive)

### Why is Long-read Mode More Powerful?

**Because of viralFlye!** ⭐

#### 1. Pfam Protein Validation

```
viralFlye uses Pfam database (~20,000 protein families)
    ↓
Identifies protein domains
    ↓
Judges based on domain combination patterns:
- Viral feature combination → Virus ✅
- Bacterial feature combination → Chromosome
- Plasmid feature combination → Plasmid
```

**Advantages**:
- Based on protein **function**, not just sequence similarity
- Can discover **completely novel viruses**
- Accurately distinguishes viruses, bacteria, plasmids

#### 2. Completeness Assessment

viralFlye includes viralComplete component:
- Assesses viral genome completeness
- Marks complete/partial viral genomes
- Provides confidence scores

#### 3. Circular/Linear Identification

- Identifies circular viral genomes
- Identifies linear viral genomes
- Many phages are circular ⭐

### Real Case Example

**Sample**: llnl_66d1047e

**2 viruses identified by viralFlye** (using default parameters):

```
contig_1085: 39,632 bp
- Type: Circular virus
- Pfam features: Phage_Mu_F, Portal_Mu, Phage_tail_terminator
- Score: 33.68
- Classification: Mu phage

contig_1192: 17,488 bp
- Type: Circular virus
- Pfam features: TerL_ATPase, Phage_capsid, Phage_portal
- Score: 25.93
- Classification: Long-tailed phage (Siphoviridae)
```

**Both are high-quality complete phage genomes!**

**After using optimized parameters** ⭐⭐⭐:
```groovy
viralflye_min_length = 500          // ⭐ Identify shorter viruses
viralflye_completeness = 0.3        // ⭐⭐⭐ Lower completeness threshold (KEY!)
```

**Expected Effect**:
- viralFlye identification: **10-20 viruses** (increased from 2 to 10-20) ⭐⭐
- Including:
  - Complete circular viruses (original 2)
  - Relatively complete linear viruses (new)
  - Large viral fragments (new)
- Quality: Still high (Pfam validated)

**Key Understanding**:
- `min_length = 500`: Relax length restriction
- `completeness = 0.3`: **This is the KEY parameter!** ⭐⭐⭐
  - Default 0.5 requires 50% completeness (very strict)
  - Lowering to 0.3 requires 30% completeness (balanced)
  - This parameter has the biggest impact on virus count

---

## Mode Configuration

### Short-read Mode Configuration

```bash
# samplesheet.csv
sample,fastq_1,fastq_2
s1,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz

# Run
sbatch run_metagenome_assembly_classification_shortread.sh

# Or command line
nextflow run metagenome_assembly_classification_workflow.nf \
    --input samplesheet.csv \
    --outdir results \
    --virsorter2_db /path/to/db \
    --deepvirfinder_dir /path/to/DVF
    # No need for --longread or --enable_viralflye
```

### Long-read Mode Configuration (Basic)

```bash
# samplesheet_long.csv
sample,fastq_long
s1,/path/to/nanopore.fastq.gz

# Run (without viralFlye)
nextflow run metagenome_assembly_classification_workflow.nf \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/db \
    --deepvirfinder_dir /path/to/DVF \
    --longread true \
    --longread_platform nano
    # Only VS2 + DVF, no viralFlye
```

### Long-read Mode Configuration (Complete Three Tools) ⭐

```bash
# samplesheet_long.csv
sample,fastq_long
s1,/path/to/nanopore.fastq.gz

# Run (with viralFlye enabled)
sbatch run_metagenome_assembly_classification_longread.sh
# Or
nextflow run metagenome_assembly_classification_workflow.nf \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/db \
    --deepvirfinder_dir /path/to/DVF \
    --longread true \
    --longread_platform nano \
    --enable_viralflye true \               # ⭐ Enable third tool
    --pfam_db /path/to/Pfam-A.hmm \
    --viralflye_env /path/to/viralFlye_env \
    --viralflye_min_length 500 \            # ⭐ Optimized: identify shorter viruses
    --viralflye_completeness 0.3            # ⭐⭐⭐ Completeness threshold (KEY!)
```

---

## Recommended Choice

### If You Have Illumina Data

```
✅ Use short-read mode
✅ Dual assemblers in parallel (MEGAHIT + SPAdes)
✅ Dual-tool validation (VirSorter2 + DeepVirFinder)
✅ Focus on assembler consensus viruses

Run: sbatch run_metagenome_assembly_classification_shortread.sh
```

### If You Have Nanopore/PacBio Data (Recommended) ⭐⭐⭐

```
✅ Use long-read mode
✅ Enable viralFlye (obtain Pfam validation) ⭐
✅ Three-tool parallel analysis
✅ Focus on three-tool consensus and viralFlye-unique viruses

Run: sbatch run_metagenome_assembly_classification_longread.sh
      Ensure ENABLE_VIRALFLYE="true"
```

---

## viralFlye Parameter Optimization Important Tips ⭐⭐⭐

### Problem: viralFlye Only Identifies 2 Viruses?

If your viralFlye output has very few viruses (e.g., only 2), **don't just lower `min_length`**!

### Solution

**The KEY parameter is `viralflye_completeness`** (completeness threshold):

```groovy
// ❌ Insufficient optimization
viralflye_min_length = 500          // Only changing this has limited effect

// ✅ Correct optimization
viralflye_min_length = 500          // ⭐ Identify shorter viruses
viralflye_completeness = 0.3        // ⭐⭐⭐ This is the KEY!
```

### Why?

viralFlye has multiple filtering layers:
1. **Length filtering** (min_length) - Can pass after relaxing
2. Pfam classification
3. **Completeness assessment** (completeness) - **This is the main bottleneck!** ⭐⭐⭐
4. viralComplete validation

Default `completeness = 0.5` (requires 50% completeness) is too strict!

### Recommended Configuration

| completeness | Virus Count | Quality | Use Case |
|--------------|------------|---------|----------|
| 0.5 (default) | 2-5 | Very High | Publication-grade |
| **0.3 (recommended)** | **10-20** | High ⭐⭐ | **Routine analysis** ⭐⭐⭐ |
| 0.2 | 20-40 | Medium-High | Exploratory |

### Detailed Guide

See:
- **`viralFlye_parameter_optimization_guide.md`** - Complete optimization guide (English)
- **`viralFlye_quick_reference.md`** - Quick reference card (English)

---

## Summary

### Core Differences

| | Short-read Mode | Long-read Mode |
|-|----------------|----------------|
| **Strategy** | Dual-assembler strategy | Three-tool strategy ⭐ |
| **Unique Feature** | Assembler comparison | viralFlye + Pfam validation ⭐ |
| **Confidence** | Dual validation | Triple validation ⭐⭐⭐ |

### Key Understanding

1. **viralFlye is the unique advantage of long-read mode** ⭐
   - Provides Pfam protein validation
   - Can discover novel viruses
   - Assesses genome completeness

2. **Short-read mode has its own advantages**
   - Dual-assembler cross-validation
   - Lower cost
   - Mature technology

3. **Both modes have dual-tool viral identification**
   - VirSorter2 + DeepVirFinder
   - This is the common foundation of both modes

4. **Both modes support viral abundance analysis** ⭐⭐⭐
   - Calculates RPM (Reads Per Million) and RPKM (Reads Per Kilobase per Million)
   - Identifies the most abundant viral species
   - Enables cross-sample abundance comparison

---

**Choose the appropriate mode and start your viral metagenomic analysis!** 🚀

**Related Documents**:
- **`README.md`** - Complete user guide (English)
- **`Mode_Selection_Guide.md`** - Mode selection guide (this document)
- **`Plan_A_three_tool_parallel_analysis.md`** - Long-read mode detailed design (English)
- **`parameter_optimization_guide.md`** - Comprehensive parameter optimization guide (English)
- **`viralFlye_parameter_optimization_guide.md`** ⭐⭐⭐ - viralFlye specialized optimization (Important! English)
- **`viralFlye_quick_reference.md`** ⭐⭐ - viralFlye quick reference card (English)



