# Viral Abundance Calculation Feature

## Overview

This workflow now includes automatic viral abundance calculation for all identified viral contigs. The abundance analysis provides two key metrics: RPM (Reads Per Million) and RPKM (Reads Per Kilobase per Million).

## Output Location

All abundance results are saved in the `results/abundance/` directory, organized by assembler/method:

### Short-Read Mode
- `results/abundance/megahit/` - Abundance for MEGAHIT-assembled viral contigs
- `results/abundance/spades/` - Abundance for SPAdes-assembled viral contigs

### Long-Read Mode
- `results/abundance/metaflye/` - Abundance for metaFlye-assembled viral contigs
- `results/abundance/viralflye/` - Abundance for viralFlye-identified viral contigs (if enabled)

## Output Files

For each assembler/method, two files are generated:

1. **`{sample}_{assembler}_abundance.csv`** - Detailed abundance data
   - Columns:
     - `contig_id`: Viral contig identifier
     - `contig_length`: Length of the viral contig (bp)
     - `mapped_reads`: Number of reads mapped to this contig
     - `total_reads`: Total number of reads in the sample
     - `rpm`: Reads Per Million
     - `rpkm`: Reads Per Kilobase per Million
   - Sorted by RPKM (descending)

2. **`{sample}_{assembler}_abundance_summary.txt`** - Summary report
   - Total reads in sample
   - Total viral contigs identified
   - Total reads mapped to viral contigs
   - Mapping rate
   - Top 10 most abundant viral contigs

## Metrics Explanation

### RPM (Reads Per Million)
```
RPM = (Mapped Reads / Total Reads) × 1,000,000
```
- Normalizes by total read count
- Useful for comparing relative abundance across samples
- Independent of contig length

### RPKM (Reads Per Kilobase per Million)
```
RPKM = (Mapped Reads × 1,000,000,000) / (Contig Length × Total Reads)
```
- Normalizes by both contig length and total read count
- Best metric for comparing abundance of contigs with different lengths
- Accounts for the fact that longer contigs will naturally receive more reads
- **Recommended metric for identifying most abundant viral species**

## Methodology

### Short-Read Mode (Illumina)
1. Reads are mapped to viral contigs using **bowtie2**
2. Alignment results are processed with **samtools**
3. Read counts are calculated per contig
4. RPM and RPKM are computed

### Long-Read Mode (PacBio/Nanopore)
1. Reads are mapped to viral contigs using **minimap2** (with appropriate preset: map-ont or map-pb)
2. Alignment results are processed with **samtools**
3. Read counts are calculated per contig
4. RPM and RPKM are computed

## Usage

The abundance calculation is automatically enabled when VirSorter2 is run (i.e., when `--skip_virsorter2` is not set to true). No additional parameters are required.

### Example Command (Short-Read)
```bash
nextflow run metagenome_assembly_classification_workflow.nf \
    --input samplesheet.csv \
    --outdir results \
    --virsorter2_db /path/to/virsorter2/db
```

### Example Command (Long-Read with viralFlye)
```bash
nextflow run metagenome_assembly_classification_workflow.nf \
    --input samplesheet_long.csv \
    --outdir results_long \
    --virsorter2_db /path/to/virsorter2/db \
    --longread true \
    --longread_platform nano \
    --enable_viralflye true
```

## Interpreting Results

1. **Top Abundant Viruses**: Check the summary file for the top 10 most abundant viruses (by RPKM)

2. **Comparing Across Assemblers**: 
   - In short-read mode, compare MEGAHIT and SPAdes abundance results
   - High-confidence viruses should show consistent abundance patterns across assemblers

3. **Cross-Sample Comparison**:
   - Use RPM for quick cross-sample comparisons
   - Use RPKM for more accurate abundance estimates, especially when contig lengths vary

4. **Low Abundance Viruses**:
   - Viruses with RPKM < 1 are typically low abundance
   - Viruses with RPKM > 100 are typically high abundance
   - These thresholds are approximate and depend on sequencing depth

## Dependencies

The following tools are automatically installed via conda/bioconda:

### Short-Read Mode
- bowtie2 (v2.5.1 or later)
- samtools (v1.18 or later)

### Long-Read Mode
- minimap2 (v2.28 or later)
- samtools (v1.18 or later)

## Notes

- If no viral contigs are identified for a sample, empty output files are created
- The abundance calculation uses only primary alignments (for long reads) or uniquely mapped reads
- Total read count includes ALL reads in the input, not just those mapped to viruses
- Abundance calculation is performed independently for each assembler/method

## Troubleshooting

### No abundance files generated
- Check that VirSorter2 identified at least one viral contig
- Verify that viral contig FASTA files are not empty
- Check the process logs for mapping errors

### Very low mapping rates
- May indicate that identified viral contigs are not truly viral
- May indicate issues with read quality
- Compare results across different assemblers for validation

### Different abundance values across assemblers
- Expected due to different assembly algorithms
- Focus on high-confidence viruses (identified by multiple tools)
- Use consensus results for most reliable abundance estimates

## Citation

If you use this abundance calculation feature in your research, please cite:
- Bowtie2: Langmead & Salzberg (2012) Nature Methods
- Minimap2: Li (2018) Bioinformatics
- SAMtools: Li et al. (2009) Bioinformatics
