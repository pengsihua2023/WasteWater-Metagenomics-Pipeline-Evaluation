#!/usr/bin/env nextflow

/*
 * Hybrid Metagenome Assembly and Kraken2 Taxonomic Classification Workflow
 * 
 * This workflow supports both short-read and long-read data:
 * 1. Short reads (Illumina): fastp QC → MEGAHIT/SPAdes assembly → abundance → Kraken2
 * 2. Long reads (Nanopore/PacBio): metaFlye assembly → abundance → Kraken2
 * 
 * Author: Assistant
 * Version: 3.0.0
 */

nextflow.enable.dsl = 2

// Workflow parameters
params.input_short = null
params.input_long = null
params.outdir_short = './results_short'
params.outdir_long = './results_long'
params.help = false

// Short-read parameters
params.skip_fastp = false
params.fastp_qualified_quality = 20
params.fastp_unqualified_percent = 40
params.fastp_min_length = 50

params.megahit_memory = 0.8
params.megahit_min_contig_len = 1000
params.spades_meta = true

// Long-read parameters
params.flye_genome_size = '5m'  // Estimated metagenome size
params.flye_min_overlap = 3000
params.long_read_type = 'nanopore'  // Options: 'nanopore', 'pacbio', 'pacbio-hifi'
params.run_viralflye = false  // Set to true and provide --viralflye_hmm to enable viralFlye
params.viralflye_hmm = null  // REQUIRED for viralFlye: Path to Pfam-A.hmm.gz
params.viralflye_min_length = 5000  // Minimum viral contig length
params.viralflye_completeness = 0.5  // Completeness cutoff for viralComplete
params.viralflye_threads = 10  // Threads for viralFlye

// Kraken2 parameters
params.kraken2_db = null

// Merge reports
params.skip_merge_reports = false

// Resource parameters
params.max_cpus = 32
params.max_memory = '256.GB'
params.max_time = '72.h'

// Print help information
if (params.help) {
    log.info """
    ==========================================
    Hybrid Metagenome Assembly and Kraken2 Workflow
    ==========================================
    
    Usage:
    nextflow run metagenome_hybrid_workflow.nf \\
        --input_short samplesheet_short.csv \\
        --input_long samplesheet_long.csv \\
        --outdir_short results_short \\
        --outdir_long results_long \\
        --kraken2_db /path/to/db
    
    Parameters:
    --input_short              Short-read input samplesheet (optional)
    --input_long               Long-read input samplesheet (optional)
    --outdir_short             Short-read output directory (default: results_short)
    --outdir_long              Long-read output directory (default: results_long)
    --kraken2_db              Kraken2 database path (required)
    
    At least one of --input_short or --input_long must be provided.
    """
    exit 0
}

// Validate parameters
if (!params.input_short && !params.input_long) {
    error "At least one input (--input_short or --input_long) is required."
}

if (!params.kraken2_db) {
    error "Kraken2 database path is required. Use --kraken2_db parameter."
}

// Print workflow information
log.info """
==========================================
🧬 Hybrid Metagenome Assembly and Kraken2 Workflow
==========================================
Workflow version: 3.0.0
Short-read input: ${params.input_short ?: 'Not provided'}
Long-read input: ${params.input_long ?: 'Not provided'}
Short-read output: ${params.outdir_short}
Long-read output: ${params.outdir_long}
Kraken2 database: ${params.kraken2_db}
==========================================
"""

// Create input channels
if (params.input_short) {
    Channel
        .fromPath(params.input_short)
        .splitCsv(header: true)
        .map { row -> 
            def sample = row.sample
            def read1 = file(row.fastq_1)
            def read2 = file(row.fastq_2)
            return tuple(sample, [read1, read2])
        }
        .set { ch_short_reads }
} else {
    ch_short_reads = Channel.empty()
}

if (params.input_long) {
    Channel
        .fromPath(params.input_long)
        .splitCsv(header: true)
        .map { row -> 
            def sample = row.sample
            def fastq = file(row.fastq_long)
            return tuple(sample, fastq)
        }
        .set { ch_long_reads }
} else {
    ch_long_reads = Channel.empty()
}

// Define workflow
workflow {
    // ========================================
    // SHORT-READ WORKFLOW
    // ========================================
    if (params.input_short) {
        // Stage 0: Quality Control
        if (!params.skip_fastp) {
            FASTP(ch_short_reads)
            ch_clean_reads = FASTP.out.clean_reads
        } else {
            ch_clean_reads = ch_short_reads
        }
        
        // Stage 1: Assembly
        MEGAHIT_ASSEMBLY(ch_clean_reads)
        SPADES_ASSEMBLY(ch_clean_reads)
        
        // Stage 2: Abundance calculation
        BOWTIE2_BUILD_MEGAHIT(MEGAHIT_ASSEMBLY.out.contigs)
        BOWTIE2_BUILD_SPADES(SPADES_ASSEMBLY.out.contigs)
        
        BOWTIE2_ALIGN_MEGAHIT(ch_clean_reads.join(BOWTIE2_BUILD_MEGAHIT.out.index))
        BOWTIE2_ALIGN_SPADES(ch_clean_reads.join(BOWTIE2_BUILD_SPADES.out.index))
        
        CALCULATE_ABUNDANCE_MEGAHIT(
            BOWTIE2_ALIGN_MEGAHIT.out.bam.join(MEGAHIT_ASSEMBLY.out.contigs)
        )
        CALCULATE_ABUNDANCE_SPADES(
            BOWTIE2_ALIGN_SPADES.out.bam.join(SPADES_ASSEMBLY.out.contigs)
        )
        
        // Stage 3: Kraken2 Classification
        KRAKEN2_CLASSIFICATION_MEGAHIT(
            MEGAHIT_ASSEMBLY.out.contigs,
            params.kraken2_db
        )
        KRAKEN2_CLASSIFICATION_SPADES(
            SPADES_ASSEMBLY.out.contigs,
            params.kraken2_db
        )
        
        // Stage 4: Merge Reports
        if (!params.skip_merge_reports) {
            KRAKEN2_CLASSIFICATION_MEGAHIT.out.kraken2_megahit
                .join(KRAKEN2_CLASSIFICATION_SPADES.out.kraken2_spades)
                .set { ch_short_reports }
            
            MERGE_KRAKEN2_REPORTS_SHORT(ch_short_reports)
        }
    }
    
    // ========================================
    // LONG-READ WORKFLOW
    // ========================================
    if (params.input_long) {
        // Stage 1: Assembly with metaFlye
        FLYE_ASSEMBLY(ch_long_reads)
        
        // Stage 2: Abundance calculation for metaFlye
        MINIMAP2_ALIGN_FLYE(
            ch_long_reads.join(FLYE_ASSEMBLY.out.contigs)
        )
        
        CALCULATE_ABUNDANCE_FLYE(
            MINIMAP2_ALIGN_FLYE.out.bam.join(FLYE_ASSEMBLY.out.contigs)
        )
        
        // Stage 3: Kraken2 Classification for metaFlye
        KRAKEN2_CLASSIFICATION_FLYE(
            FLYE_ASSEMBLY.out.contigs,
            params.kraken2_db
        )
        
        // Stage 4: viralFlye - Identify viral contigs from metaFlye results
        if (params.run_viralflye) {
            VIRALFLYE_IDENTIFY(
                ch_long_reads.join(FLYE_ASSEMBLY.out.assembly_dir)
            )
            
            // Calculate abundance for viral contigs (linear + circular)
            MINIMAP2_ALIGN_VIRALFLYE_LINEAR(
                ch_long_reads.join(VIRALFLYE_IDENTIFY.out.linear_contigs)
            )
            
            MINIMAP2_ALIGN_VIRALFLYE_CIRCULAR(
                ch_long_reads.join(VIRALFLYE_IDENTIFY.out.circular_contigs)
            )
            
            CALCULATE_ABUNDANCE_VIRALFLYE_LINEAR(
                MINIMAP2_ALIGN_VIRALFLYE_LINEAR.out.bam.join(VIRALFLYE_IDENTIFY.out.linear_contigs)
            )
            
            CALCULATE_ABUNDANCE_VIRALFLYE_CIRCULAR(
                MINIMAP2_ALIGN_VIRALFLYE_CIRCULAR.out.bam.join(VIRALFLYE_IDENTIFY.out.circular_contigs)
            )
            
            // Kraken2 classification for viral contigs
            KRAKEN2_CLASSIFICATION_VIRALFLYE_LINEAR(
                VIRALFLYE_IDENTIFY.out.linear_contigs,
                params.kraken2_db
            )
            
            KRAKEN2_CLASSIFICATION_VIRALFLYE_CIRCULAR(
                VIRALFLYE_IDENTIFY.out.circular_contigs,
                params.kraken2_db
            )
        }
    }
}

// ================================================================================
// Process Definitions - SHORT-READ PROCESSES
// ================================================================================

// Process: fastp Quality Control
process FASTP {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.outdir_short}/fastp", mode: 'copy', pattern: "*.{html,json}"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_clean_R{1,2}.fastq.gz"), emit: clean_reads
    path("${sample}_fastp.html"), emit: html
    path("${sample}_fastp.json"), emit: json
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    echo "=== fastp Quality Control: ${sample} ==="
    
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample}_clean_R1.fastq.gz \\
        -O ${sample}_clean_R2.fastq.gz \\
        --thread ${task.cpus} \\
        --qualified_quality_phred ${params.fastp_qualified_quality} \\
        --unqualified_percent_limit ${params.fastp_unqualified_percent} \\
        --length_required ${params.fastp_min_length} \\
        --detect_adapter_for_pe \\
        --compression 6 \\
        --html ${sample}_fastp.html \\
        --json ${sample}_fastp.json
    
    echo "fastp: Quality control completed for ${sample}"
    """
}

// Process: MEGAHIT Assembly
process MEGAHIT_ASSEMBLY {
    tag "${sample}_MEGAHIT"
    label 'process_high'
    container 'docker://quay.io/biocontainers/megahit:1.2.9--h2e03b76_1'
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("megahit_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== MEGAHIT Assembly: ${sample} ==="
    
    megahit \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o megahit_output \\
        -t ${task.cpus} \\
        --memory ${params.megahit_memory} \\
        --min-contig-len ${params.megahit_min_contig_len}
    
    cp megahit_output/final.contigs.fa megahit_contigs.fa
    
    echo "MEGAHIT: Generated \$(grep -c ">" megahit_contigs.fa) contigs"
    """
}

// Process: SPAdes Assembly
process SPADES_ASSEMBLY {
    tag "${sample}_SPAdes"
    label 'process_high'
    container 'docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("spades_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== metaSPAdes Assembly: ${sample} ==="
    
    metaspades.py \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o spades_output \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()} \\
        --only-assembler
    
    cp spades_output/contigs.fasta spades_contigs.fa
    
    echo "metaSPAdes: Generated \$(grep -c ">" spades_contigs.fa) contigs"
    """
}

// Process: Build Bowtie2 index for MEGAHIT contigs
process BOWTIE2_BUILD_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_medium'
    container 'docker://quay.io/biocontainers/bowtie2:2.5.1--py310h8d7afc0_0'
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("megahit_index*"), emit: index
    
    script:
    """
    bowtie2-build --threads ${task.cpus} ${contigs} megahit_index
    """
}

// Process: Build Bowtie2 index for SPAdes contigs
process BOWTIE2_BUILD_SPADES {
    tag "${sample}_SPAdes"
    label 'process_medium'
    container 'docker://quay.io/biocontainers/bowtie2:2.5.1--py310h8d7afc0_0'
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("spades_index*"), emit: index
    
    script:
    """
    bowtie2-build --threads ${task.cpus} ${contigs} spades_index
    """
}

// Process: Align reads to MEGAHIT contigs
process BOWTIE2_ALIGN_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_high'
    container 'docker://quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0'
    
    input:
    tuple val(sample), path(reads), path(index)
    
    output:
    tuple val(sample), path("${sample}_megahit.sorted.bam"), path("${sample}_megahit.sorted.bam.bai"), emit: bam
    
    script:
    """
    bowtie2 \\
        -x megahit_index \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --threads ${task.cpus} \\
        --no-unal \\
        | samtools view -bS - \\
        | samtools sort -@ ${task.cpus} -o ${sample}_megahit.sorted.bam
    
    samtools index ${sample}_megahit.sorted.bam
    """
}

// Process: Align reads to SPAdes contigs
process BOWTIE2_ALIGN_SPADES {
    tag "${sample}_SPAdes"
    label 'process_high'
    container 'docker://quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0'
    
    input:
    tuple val(sample), path(reads), path(index)
    
    output:
    tuple val(sample), path("${sample}_spades.sorted.bam"), path("${sample}_spades.sorted.bam.bai"), emit: bam
    
    script:
    """
    bowtie2 \\
        -x spades_index \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --threads ${task.cpus} \\
        --no-unal \\
        | samtools view -bS - \\
        | samtools sort -@ ${task.cpus} -o ${sample}_spades.sorted.bam
    
    samtools index ${sample}_spades.sorted.bam
    """
}

// Process: Calculate RPM and RPKM for MEGAHIT contigs
process CALCULATE_ABUNDANCE_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_low'
    conda 'conda-forge::python=3.10 bioconda::samtools=1.17 conda-forge::biopython=1.81'
    publishDir "${params.outdir_short}/abundance_megahit", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(bam), path(bai), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_megahit_abundance.txt"), emit: abundance
    path("${sample}_megahit_abundance_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env bash
    
    # Setup symbolic link to resolve libbz2.so.1.0 dependency
    mkdir -p \$HOME/.local/lib_tmp
    for libdir in /usr/lib /usr/lib64 /lib /lib64 /usr/lib/x86_64-linux-gnu; do
        if [ -f "\$libdir/libbz2.so.1" ]; then
            ln -sf "\$libdir/libbz2.so.1" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        elif [ -f "\$libdir/libbz2.so" ]; then
            ln -sf "\$libdir/libbz2.so" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        fi
    done
    export LD_LIBRARY_PATH=\$HOME/.local/lib_tmp:\$LD_LIBRARY_PATH
    
    python3 << 'PYTHON_SCRIPT'
import subprocess
from Bio import SeqIO

contig_lengths = {}
for record in SeqIO.parse("${contigs}", "fasta"):
    contig_lengths[record.id] = len(record.seq)

idxstats_output = subprocess.check_output(
    ["samtools", "idxstats", "${bam}"],
    universal_newlines=True,
    stderr=subprocess.PIPE
)

contig_reads = {}
total_mapped_reads = 0

for line in idxstats_output.strip().split('\\n'):
    parts = line.split('\\t')
    if len(parts) >= 3:
        contig_name = parts[0]
        mapped_reads = int(parts[2])
        if contig_name != "*":
            contig_reads[contig_name] = mapped_reads
            total_mapped_reads += mapped_reads

with open("${sample}_megahit_abundance.txt", 'w') as out_f:
    out_f.write("Contig_ID\\tLength(bp)\\tMapped_Reads\\tRPM\\tRPKM\\n")
    for contig_name in sorted(contig_lengths.keys()):
        length = contig_lengths[contig_name]
        reads = contig_reads.get(contig_name, 0)
        rpm = (reads / total_mapped_reads * 1e6) if total_mapped_reads > 0 else 0
        rpkm = (reads / (length / 1000) / (total_mapped_reads / 1e6)) if total_mapped_reads > 0 and length > 0 else 0
        out_f.write(f"{contig_name}\\t{length}\\t{reads}\\t{rpm:.4f}\\t{rpkm:.4f}\\n")

with open("${sample}_megahit_abundance_summary.txt", 'w') as sum_f:
    sum_f.write("="*80 + "\\n")
    sum_f.write("MEGAHIT Contigs Abundance Summary\\n")
    sum_f.write("="*80 + "\\n\\n")
    sum_f.write(f"Sample: ${sample}\\n")
    sum_f.write(f"Total contigs: {len(contig_lengths)}\\n")
    sum_f.write(f"Total mapped reads: {total_mapped_reads:,}\\n")
    if len(contig_lengths) > 0:
        sum_f.write(f"Average contig length: {sum(contig_lengths.values()) / len(contig_lengths):.2f} bp\\n")
        sum_f.write(f"Longest contig: {max(contig_lengths.values()):,} bp\\n")
        sum_f.write(f"Shortest contig: {min(contig_lengths.values()):,} bp\\n")
    sum_f.write("\\n" + "="*80 + "\\n")
PYTHON_SCRIPT
    """
}

// Process: Calculate RPM and RPKM for SPAdes contigs
process CALCULATE_ABUNDANCE_SPADES {
    tag "${sample}_SPAdes"
    label 'process_low'
    conda 'conda-forge::python=3.10 bioconda::samtools=1.17 conda-forge::biopython=1.81'
    publishDir "${params.outdir_short}/abundance_spades", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(bam), path(bai), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_spades_abundance.txt"), emit: abundance
    path("${sample}_spades_abundance_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env bash
    
    # Setup symbolic link to resolve libbz2.so.1.0 dependency
    mkdir -p \$HOME/.local/lib_tmp
    for libdir in /usr/lib /usr/lib64 /lib /lib64 /usr/lib/x86_64-linux-gnu; do
        if [ -f "\$libdir/libbz2.so.1" ]; then
            ln -sf "\$libdir/libbz2.so.1" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        elif [ -f "\$libdir/libbz2.so" ]; then
            ln -sf "\$libdir/libbz2.so" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        fi
    done
    export LD_LIBRARY_PATH=\$HOME/.local/lib_tmp:\$LD_LIBRARY_PATH
    
    python3 << 'PYTHON_SCRIPT'
import subprocess
from Bio import SeqIO

contig_lengths = {}
for record in SeqIO.parse("${contigs}", "fasta"):
    contig_lengths[record.id] = len(record.seq)

idxstats_output = subprocess.check_output(
    ["samtools", "idxstats", "${bam}"],
    universal_newlines=True,
    stderr=subprocess.PIPE
)

contig_reads = {}
total_mapped_reads = 0

for line in idxstats_output.strip().split('\\n'):
    parts = line.split('\\t')
    if len(parts) >= 3:
        contig_name = parts[0]
        mapped_reads = int(parts[2])
        if contig_name != "*":
            contig_reads[contig_name] = mapped_reads
            total_mapped_reads += mapped_reads

with open("${sample}_spades_abundance.txt", 'w') as out_f:
    out_f.write("Contig_ID\\tLength(bp)\\tMapped_Reads\\tRPM\\tRPKM\\n")
    for contig_name in sorted(contig_lengths.keys()):
        length = contig_lengths[contig_name]
        reads = contig_reads.get(contig_name, 0)
        rpm = (reads / total_mapped_reads * 1e6) if total_mapped_reads > 0 else 0
        rpkm = (reads / (length / 1000) / (total_mapped_reads / 1e6)) if total_mapped_reads > 0 and length > 0 else 0
        out_f.write(f"{contig_name}\\t{length}\\t{reads}\\t{rpm:.4f}\\t{rpkm:.4f}\\n")

with open("${sample}_spades_abundance_summary.txt", 'w') as sum_f:
    sum_f.write("="*80 + "\\n")
    sum_f.write("SPAdes Contigs Abundance Summary\\n")
    sum_f.write("="*80 + "\\n\\n")
    sum_f.write(f"Sample: ${sample}\\n")
    sum_f.write(f"Total contigs: {len(contig_lengths)}\\n")
    sum_f.write(f"Total mapped reads: {total_mapped_reads:,}\\n")
    if len(contig_lengths) > 0:
        sum_f.write(f"Average contig length: {sum(contig_lengths.values()) / len(contig_lengths):.2f} bp\\n")
        sum_f.write(f"Longest contig: {max(contig_lengths.values()):,} bp\\n")
        sum_f.write(f"Shortest contig: {min(contig_lengths.values()):,} bp\\n")
    sum_f.write("\\n" + "="*80 + "\\n")
PYTHON_SCRIPT
    """
}

// Process: Kraken2 Classification for MEGAHIT
process KRAKEN2_CLASSIFICATION_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_medium'
    conda 'bioconda::kraken2=2.1.3'
    publishDir "${params.outdir_short}/kraken2_megahit", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs)
    val(kraken2_db)
    
    output:
    tuple val(sample), path("${sample}_megahit_*.txt"), emit: kraken2_megahit
    
    script:
    """
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${sample}_megahit_classification.txt \\
        --report ${sample}_megahit_report.txt \\
        ${contigs}
    """
}

// Process: Kraken2 Classification for SPAdes
process KRAKEN2_CLASSIFICATION_SPADES {
    tag "${sample}_SPAdes"
    label 'process_medium'
    conda 'bioconda::kraken2=2.1.3'
    publishDir "${params.outdir_short}/kraken2_spades", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs)
    val(kraken2_db)
    
    output:
    tuple val(sample), path("${sample}_spades_*.txt"), emit: kraken2_spades
    
    script:
    """
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${sample}_spades_classification.txt \\
        --report ${sample}_spades_report.txt \\
        ${contigs}
    """
}

// Process: Merge Kraken2 Reports for short reads
process MERGE_KRAKEN2_REPORTS_SHORT {
    tag "${sample}"
    label 'process_low'
    publishDir "${params.outdir_short}/merged_reports", mode: 'copy', pattern: "*"
    conda 'pandas=1.5.3 numpy=1.23.5'
    
    input:
    tuple val(sample), path(megahit_reports), path(spades_reports)
    
    output:
    tuple val(sample), path("${sample}_merged_report.txt"), emit: merged_report
    path("${sample}_merged_report.csv"), emit: merged_csv
    path("${sample}_virus_consensus.txt"), emit: virus_consensus_txt
    path("${sample}_virus_consensus.csv"), emit: virus_consensus_csv
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    def parse_kraken2_report(file_path):
        data = []
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\\t')
                if len(parts) >= 6:
                    data.append({
                        'percent': float(parts[0]),
                        'reads': int(parts[1]),
                        'direct_reads': int(parts[2]),
                        'rank': parts[3],
                        'tax_id': parts[4],
                        'name': parts[5].strip()
                    })
        return pd.DataFrame(data)
    
    megahit_report = [f for f in "${megahit_reports}".split() if f.endswith('_report.txt')][0]
    spades_report = [f for f in "${spades_reports}".split() if f.endswith('_report.txt')][0]
    
    megahit_df = parse_kraken2_report(megahit_report)
    spades_df = parse_kraken2_report(spades_report)
    
    spades_df.columns = [f'spades_{col}' if col not in ['tax_id', 'name', 'rank'] else col 
                         for col in spades_df.columns]
    megahit_df.columns = [f'megahit_{col}' if col not in ['tax_id', 'name', 'rank'] else col 
                          for col in megahit_df.columns]
    
    merged = pd.merge(spades_df, megahit_df, on=['tax_id', 'rank'], how='outer', suffixes=('_spades', '_megahit'))
    merged = merged.fillna(0)
    
    if 'name_spades' in merged.columns and 'name_megahit' in merged.columns:
        merged['name'] = merged['name_spades'].where(merged['name_spades'] != 0, merged['name_megahit'])
        merged = merged.drop(['name_spades', 'name_megahit'], axis=1)
    
    with open("${sample}_merged_report.txt", 'w') as f:
        f.write("="*80 + "\\n")
        f.write("Kraken2 Analysis Report - MEGAHIT vs SPAdes\\n")
        f.write("="*80 + "\\n\\n")
        f.write(f"Sample: ${sample}\\n")
        f.write(f"Total SPAdes contigs: {merged['spades_reads'].sum():.1f}\\n")
        f.write(f"Total MEGAHIT contigs: {merged['megahit_reads'].sum():.1f}\\n")
        f.write("="*80 + "\\n")
    
    output_df = merged[['tax_id', 'rank', 'name', 'spades_reads', 'spades_percent', 
                        'megahit_reads', 'megahit_percent']].copy()
    output_df['total_reads'] = output_df['spades_reads'] + output_df['megahit_reads']
    output_df = output_df.sort_values('total_reads', ascending=False)
    output_df.to_csv("${sample}_merged_report.csv", index=False)
    
    # === Consensus Virus Analysis ===
    # Filter only viral classifications (excluding unclassified and root)
    virus_df = output_df[
        (output_df['name'].str.contains('virus|phage|viral|Virus|Phage|Viral', case=False, na=False)) &
        (~output_df['name'].str.contains('unclassified', case=False, na=False)) &
        (output_df['tax_id'] != '0')
    ].copy()
    
    if len(virus_df) > 0:
        # Categorize viruses based on detection in assemblers
        virus_df['detection'] = 'Unknown'
        virus_df.loc[(virus_df['spades_reads'] > 0) & (virus_df['megahit_reads'] > 0), 'detection'] = 'Consensus (Both)'
        virus_df.loc[(virus_df['spades_reads'] > 0) & (virus_df['megahit_reads'] == 0), 'detection'] = 'SPAdes only'
        virus_df.loc[(virus_df['spades_reads'] == 0) & (virus_df['megahit_reads'] > 0), 'detection'] = 'MEGAHIT only'
        
        # Calculate agreement for consensus viruses
        consensus = virus_df[virus_df['detection'] == 'Consensus (Both)'].copy()
        if len(consensus) > 0:
            consensus['agreement_ratio'] = consensus.apply(
                lambda row: min(row['spades_reads'], row['megahit_reads']) / max(row['spades_reads'], row['megahit_reads'])
                if max(row['spades_reads'], row['megahit_reads']) > 0 else 0,
                axis=1
            )
        
        # Save virus consensus report
        with open("${sample}_virus_consensus.txt", 'w') as f:
            f.write("="*80 + "\\n")
            f.write("Viral Consensus Analysis - MEGAHIT vs SPAdes\\n")
            f.write("="*80 + "\\n\\n")
            f.write(f"Sample: ${sample}\\n\\n")
            
            consensus_count = len(virus_df[virus_df['detection'] == 'Consensus (Both)'])
            spades_only_count = len(virus_df[virus_df['detection'] == 'SPAdes only'])
            megahit_only_count = len(virus_df[virus_df['detection'] == 'MEGAHIT only'])
            
            f.write(f"Total viral classifications: {len(virus_df)}\\n")
            f.write(f"  ✅ Consensus viruses (detected by BOTH): {consensus_count}\\n")
            f.write(f"  ⚠️  SPAdes only: {spades_only_count}\\n")
            f.write(f"  ⚠️  MEGAHIT only: {megahit_only_count}\\n")
            f.write(f"\\nConsensus rate: {consensus_count/len(virus_df)*100:.1f}%\\n")
            f.write("="*80 + "\\n\\n")
            
            if len(consensus) > 0:
                f.write("HIGH CONFIDENCE VIRUSES (Consensus - Detected by Both Assemblers):\\n")
                f.write("-"*80 + "\\n")
                for idx, row in consensus.sort_values('total_reads', ascending=False).iterrows():
                    f.write(f"\\n{row['name']}\\n")
                    f.write(f"  Tax ID: {row['tax_id']}\\n")
                    f.write(f"  Rank: {row['rank']}\\n")
                    f.write(f"  SPAdes: {int(row['spades_reads'])} contigs ({row['spades_percent']:.2f}%)\\n")
                    f.write(f"  MEGAHIT: {int(row['megahit_reads'])} contigs ({row['megahit_percent']:.2f}%)\\n")
                    f.write(f"  Total: {int(row['total_reads'])} contigs\\n")
                    f.write(f"  Agreement: {row['agreement_ratio']:.2f}\\n")
                f.write("\\n" + "="*80 + "\\n")
        
        # Save detailed CSV
        virus_df_out = virus_df[['tax_id', 'rank', 'name', 'spades_reads', 'megahit_reads', 
                                  'total_reads', 'detection']].copy()
        if 'agreement_ratio' in consensus.columns and len(consensus) > 0:
            virus_df_out = virus_df_out.merge(
                consensus[['tax_id', 'agreement_ratio']], 
                on='tax_id', 
                how='left'
            )
        virus_df_out = virus_df_out.sort_values(['detection', 'total_reads'], ascending=[True, False])
        virus_df_out.to_csv("${sample}_virus_consensus.csv", index=False)
    else:
        # No viruses found
        with open("${sample}_virus_consensus.txt", 'w') as f:
            f.write("No viral classifications found in either assembler.\\n")
        pd.DataFrame().to_csv("${sample}_virus_consensus.csv", index=False)
    """
}

// ================================================================================
// Process Definitions - LONG-READ PROCESSES
// ================================================================================

// Process: metaFlye Assembly (general metagenome)
process FLYE_ASSEMBLY {
    tag "${sample}_metaFlye"
    label 'process_high'
    container 'docker://quay.io/biocontainers/flye:2.9.2--py39h6935b12_1'
    publishDir "${params.outdir_long}/flye_assembly", mode: 'copy', pattern: "flye_output", saveAs: { "${sample}_flye_assembly" }
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("flye_contigs.fa"), emit: contigs
    tuple val(sample), path("flye_output"), emit: assembly_dir
    
    script:
    def read_type_arg = params.long_read_type == 'nanopore' ? '--nano-raw' : 
                        params.long_read_type == 'pacbio' ? '--pacbio-raw' :
                        params.long_read_type == 'pacbio-hifi' ? '--pacbio-hifi' : '--nano-raw'
    """
    echo "=== metaFlye Assembly: ${sample} ==="
    echo "Read type: ${params.long_read_type}"
    
    flye \\
        ${read_type_arg} ${reads} \\
        --out-dir flye_output \\
        --threads ${task.cpus} \\
        --genome-size ${params.flye_genome_size} \\
        --meta \\
        --min-overlap ${params.flye_min_overlap}
    
    cp flye_output/assembly.fasta flye_contigs.fa
    
    echo "metaFlye: Generated \$(grep -c ">" flye_contigs.fa) contigs"
    echo "Assembly directory saved for viralFlye analysis"
    """
}

// Process: viralFlye - Identify viral contigs from metaFlye results
process VIRALFLYE_IDENTIFY {
    tag "${sample}_viralFlye"
    label 'process_high'
    publishDir "${params.outdir_long}/viralflye", mode: 'copy', 
               saveAs: { filename -> filename.replaceAll('viralflye_output/', '') }
    
    input:
    tuple val(sample), path(reads), path(assembly_dir)
    
    output:
    tuple val(sample), path("viralflye_output/linears_viralFlye.fasta"), emit: linear_contigs, optional: true
    tuple val(sample), path("viralflye_output/circulars_viralFlye.fasta"), emit: circular_contigs, optional: true
    tuple val(sample), path("viralflye_output/components_viralFlye.fasta"), emit: component_contigs, optional: true
    path("viralflye_output/*"), emit: all_outputs
    
    script:
    if (!params.viralflye_hmm) {
        error "viralFlye requires --viralflye_hmm parameter (path to Pfam-A.hmm.gz). Please provide the HMM database path or set run_viralflye=false in config."
    }
    """
    #!/usr/bin/env bash
    set +u  # Allow unbound variables to avoid conda deactivate issues
    
    echo "=== viralFlye: Identifying viral contigs from metaFlye results: ${sample} ==="
    
    # Load conda base and activate viralFlye_env
    CONDA_BASE=\$(conda info --base)
    source "\$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate viralFlye_env
    
    # Verify viralFlye is available
    echo "Checking viralFlye installation..."
    VIRALFLYE_PATH=\$(which viralFlye.py)
    echo "viralFlye.py found at: \$VIRALFLYE_PATH"
    
    # Find viralFlye module location and add to PYTHONPATH
    # viralFlye.py is typically in bin/, the module is in lib/pythonX.X/site-packages/
    CONDA_ENV_PATH=\$(dirname \$(dirname \$VIRALFLYE_PATH))
    echo "Conda env path: \$CONDA_ENV_PATH"
    
    # Add potential viralFlye locations to PYTHONPATH
    for PYTHON_VERSION in \$CONDA_ENV_PATH/lib/python*/site-packages; do
        if [ -d "\$PYTHON_VERSION" ]; then
            export PYTHONPATH="\$PYTHON_VERSION:\$PYTHONPATH"
            echo "Added to PYTHONPATH: \$PYTHON_VERSION"
        fi
    done
    
    # Also check for viralFlye in common locations
    for VIRALFLYE_DIR in \$HOME/viralFlye /scratch/*/viralFlye /tmp/viralFlye; do
        if [ -d "\$VIRALFLYE_DIR" ]; then
            export PYTHONPATH="\$VIRALFLYE_DIR:\$PYTHONPATH"
            echo "Added viralFlye source to PYTHONPATH: \$VIRALFLYE_DIR"
        fi
    done
    
    echo "Current PYTHONPATH: \$PYTHONPATH"
    
    # Verify Python version from conda environment
    echo "Python version check:"
    which python
    python --version
    
    # Use conda environment's Python explicitly
    CONDA_PYTHON="\$CONDA_ENV_PATH/bin/python"
    echo "Using conda Python: \$CONDA_PYTHON"
    \$CONDA_PYTHON --version
    
    # Test viralFlye module import using conda's Python
    \$CONDA_PYTHON -c "from viralflye.main import main; print('✅ viralFlye module loaded successfully')" || {
        echo "❌ ERROR: viralFlye module not found"
        echo "Attempting to locate viralFlye module..."
        \$CONDA_PYTHON -c "import sys; print('Python path:'); [print(p) for p in sys.path]"
        echo ""
        echo "Please ensure viralFlye is installed in viralFlye_env"
        exit 1
    }
    
    mkdir -p viralflye_output
    
    echo "Running viralFlye analysis..."
    echo "HMM database: ${params.viralflye_hmm}"
    
    # Use conda environment's viralFlye.py with explicit Python
    \$CONDA_PYTHON \$VIRALFLYE_PATH \\
        --dir ${assembly_dir} \\
        --reads ${reads} \\
        --hmm ${params.viralflye_hmm} \\
        --outdir viralflye_output \\
        --min_viral_length ${params.viralflye_min_length} \\
        --completeness ${params.viralflye_completeness} \\
        --threads ${params.viralflye_threads}
    
    VIRALFLYE_EXIT=\$?
    echo "viralFlye exit status: \$VIRALFLYE_EXIT"
    
    # Count identified viral contigs
    if [ -f "viralflye_output/linears_viralFlye.fasta" ]; then
        LINEAR_COUNT=\$(grep -c ">" viralflye_output/linears_viralFlye.fasta || echo 0)
        echo "Linear viral contigs: \$LINEAR_COUNT"
    else
        echo "No linear viral contigs identified"
        touch viralflye_output/linears_viralFlye.fasta
    fi
    
    if [ -f "viralflye_output/circulars_viralFlye.fasta" ]; then
        CIRCULAR_COUNT=\$(grep -c ">" viralflye_output/circulars_viralFlye.fasta || echo 0)
        echo "Circular viral contigs: \$CIRCULAR_COUNT"
    else
        echo "No circular viral contigs identified"
        touch viralflye_output/circulars_viralFlye.fasta
    fi
    
    if [ -f "viralflye_output/components_viralFlye.fasta" ]; then
        COMPONENT_COUNT=\$(grep -c ">" viralflye_output/components_viralFlye.fasta || echo 0)
        echo "Component viral contigs: \$COMPONENT_COUNT"
    else
        echo "No component viral contigs identified"
        touch viralflye_output/components_viralFlye.fasta
    fi
    
    echo "viralFlye: Viral identification completed"
    
    # Exit with viralFlye's exit code (but allow 0 even if no viruses found)
    exit 0
    """
}

// Process: Align long reads to metaFlye contigs
process MINIMAP2_ALIGN_FLYE {
    tag "${sample}_metaFlye"
    label 'process_high'
    conda 'bioconda::minimap2=2.26 bioconda::samtools=1.17'
    
    input:
    tuple val(sample), path(reads), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_flye.sorted.bam"), path("${sample}_flye.sorted.bam.bai"), emit: bam
    
    script:
    def minimap_preset = params.long_read_type == 'nanopore' ? 'map-ont' : 
                         params.long_read_type == 'pacbio' ? 'map-pb' :
                         params.long_read_type == 'pacbio-hifi' ? 'asm20' : 'map-ont'
    """
    #!/usr/bin/env bash
    
    echo "=== Aligning long reads to metaFlye contigs: ${sample} ==="
    echo "Minimap2 preset: ${minimap_preset}"
    
    # Setup symbolic link to resolve libbz2.so.1.0 dependency
    mkdir -p \$HOME/.local/lib_tmp
    for libdir in /usr/lib /usr/lib64 /lib /lib64 /usr/lib/x86_64-linux-gnu; do
        if [ -f "\$libdir/libbz2.so.1" ]; then
            ln -sf "\$libdir/libbz2.so.1" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        elif [ -f "\$libdir/libbz2.so" ]; then
            ln -sf "\$libdir/libbz2.so" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        fi
    done
    export LD_LIBRARY_PATH=\$HOME/.local/lib_tmp:\$LD_LIBRARY_PATH
    
    minimap2 \\
        -ax ${minimap_preset} \\
        -t ${task.cpus} \\
        ${contigs} \\
        ${reads} \\
        | samtools view -bS - \\
        | samtools sort -@ ${task.cpus} -o ${sample}_flye.sorted.bam
    
    samtools index ${sample}_flye.sorted.bam
    
    echo "Alignment to metaFlye contigs completed"
    """
}

// Process: Align long reads to viralFlye linear viral contigs
process MINIMAP2_ALIGN_VIRALFLYE_LINEAR {
    tag "${sample}_viralFlye_linear"
    label 'process_medium'
    conda 'bioconda::minimap2=2.26 bioconda::samtools=1.17'
    
    input:
    tuple val(sample), path(reads), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_viralflye_linear.sorted.bam"), path("${sample}_viralflye_linear.sorted.bam.bai"), emit: bam
    
    script:
    def minimap_preset = params.long_read_type == 'nanopore' ? 'map-ont' : 
                         params.long_read_type == 'pacbio' ? 'map-pb' :
                         params.long_read_type == 'pacbio-hifi' ? 'asm20' : 'map-ont'
    """
    #!/usr/bin/env bash
    
    echo "=== Aligning reads to linear viral contigs: ${sample} ==="
    
    # Setup symbolic link
    mkdir -p \$HOME/.local/lib_tmp
    for libdir in /usr/lib /usr/lib64 /lib /lib64 /usr/lib/x86_64-linux-gnu; do
        if [ -f "\$libdir/libbz2.so.1" ]; then
            ln -sf "\$libdir/libbz2.so.1" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        elif [ -f "\$libdir/libbz2.so" ]; then
            ln -sf "\$libdir/libbz2.so" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        fi
    done
    export LD_LIBRARY_PATH=\$HOME/.local/lib_tmp:\$LD_LIBRARY_PATH
    
    minimap2 -ax ${minimap_preset} -t ${task.cpus} ${contigs} ${reads} \\
        | samtools view -bS - \\
        | samtools sort -@ ${task.cpus} -o ${sample}_viralflye_linear.sorted.bam
    
    samtools index ${sample}_viralflye_linear.sorted.bam
    """
}

// Process: Align long reads to viralFlye circular viral contigs
process MINIMAP2_ALIGN_VIRALFLYE_CIRCULAR {
    tag "${sample}_viralFlye_circular"
    label 'process_medium'
    conda 'bioconda::minimap2=2.26 bioconda::samtools=1.17'
    
    input:
    tuple val(sample), path(reads), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_viralflye_circular.sorted.bam"), path("${sample}_viralflye_circular.sorted.bam.bai"), emit: bam
    
    script:
    def minimap_preset = params.long_read_type == 'nanopore' ? 'map-ont' : 
                         params.long_read_type == 'pacbio' ? 'map-pb' :
                         params.long_read_type == 'pacbio-hifi' ? 'asm20' : 'map-ont'
    """
    #!/usr/bin/env bash
    
    echo "=== Aligning reads to circular viral contigs: ${sample} ==="
    
    # Setup symbolic link
    mkdir -p \$HOME/.local/lib_tmp
    for libdir in /usr/lib /usr/lib64 /lib /lib64 /usr/lib/x86_64-linux-gnu; do
        if [ -f "\$libdir/libbz2.so.1" ]; then
            ln -sf "\$libdir/libbz2.so.1" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        elif [ -f "\$libdir/libbz2.so" ]; then
            ln -sf "\$libdir/libbz2.so" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        fi
    done
    export LD_LIBRARY_PATH=\$HOME/.local/lib_tmp:\$LD_LIBRARY_PATH
    
    minimap2 -ax ${minimap_preset} -t ${task.cpus} ${contigs} ${reads} \\
        | samtools view -bS - \\
        | samtools sort -@ ${task.cpus} -o ${sample}_viralflye_circular.sorted.bam
    
    samtools index ${sample}_viralflye_circular.sorted.bam
    """
}

// Process: Calculate RPM and RPKM for metaFlye contigs
process CALCULATE_ABUNDANCE_FLYE {
    tag "${sample}_metaFlye"
    label 'process_low'
    conda 'conda-forge::python=3.10 bioconda::samtools=1.17 conda-forge::biopython=1.81'
    publishDir "${params.outdir_long}/abundance_flye", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(bam), path(bai), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_flye_abundance.txt"), emit: abundance
    path("${sample}_flye_abundance_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env bash
    
    # Setup symbolic link to resolve libbz2.so.1.0 dependency
    mkdir -p \$HOME/.local/lib_tmp
    for libdir in /usr/lib /usr/lib64 /lib /lib64 /usr/lib/x86_64-linux-gnu; do
        if [ -f "\$libdir/libbz2.so.1" ]; then
            ln -sf "\$libdir/libbz2.so.1" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        elif [ -f "\$libdir/libbz2.so" ]; then
            ln -sf "\$libdir/libbz2.so" \$HOME/.local/lib_tmp/libbz2.so.1.0
            break
        fi
    done
    export LD_LIBRARY_PATH=\$HOME/.local/lib_tmp:\$LD_LIBRARY_PATH
    
    python3 << 'PYTHON_SCRIPT'
import subprocess
from Bio import SeqIO

contig_lengths = {}
for record in SeqIO.parse("${contigs}", "fasta"):
    contig_lengths[record.id] = len(record.seq)

idxstats_output = subprocess.check_output(
    ["samtools", "idxstats", "${bam}"],
    universal_newlines=True,
    stderr=subprocess.PIPE
)

contig_reads = {}
total_mapped_reads = 0

for line in idxstats_output.strip().split('\\n'):
    parts = line.split('\\t')
    if len(parts) >= 3:
        contig_name = parts[0]
        mapped_reads = int(parts[2])
        if contig_name != "*":
            contig_reads[contig_name] = mapped_reads
            total_mapped_reads += mapped_reads

with open("${sample}_flye_abundance.txt", 'w') as out_f:
    out_f.write("Contig_ID\\tLength(bp)\\tMapped_Reads\\tRPM\\tRPKM\\n")
    for contig_name in sorted(contig_lengths.keys()):
        length = contig_lengths[contig_name]
        reads = contig_reads.get(contig_name, 0)
        rpm = (reads / total_mapped_reads * 1e6) if total_mapped_reads > 0 else 0
        rpkm = (reads / (length / 1000) / (total_mapped_reads / 1e6)) if total_mapped_reads > 0 and length > 0 else 0
        out_f.write(f"{contig_name}\\t{length}\\t{reads}\\t{rpm:.4f}\\t{rpkm:.4f}\\n")

with open("${sample}_flye_abundance_summary.txt", 'w') as sum_f:
    sum_f.write("="*80 + "\\n")
    sum_f.write("metaFlye Contigs Abundance Summary\\n")
    sum_f.write("="*80 + "\\n\\n")
    sum_f.write(f"Sample: ${sample}\\n")
    sum_f.write(f"Total contigs: {len(contig_lengths)}\\n")
    sum_f.write(f"Total mapped reads: {total_mapped_reads:,}\\n")
    if len(contig_lengths) > 0:
        sum_f.write(f"Average contig length: {sum(contig_lengths.values()) / len(contig_lengths):.2f} bp\\n")
        sum_f.write(f"Longest contig: {max(contig_lengths.values()):,} bp\\n")
        sum_f.write(f"Shortest contig: {min(contig_lengths.values()):,} bp\\n")
    sum_f.write("\\n" + "="*80 + "\\n")
PYTHON_SCRIPT
    """
}

// Process: Calculate RPM and RPKM for viralFlye linear viral contigs
process CALCULATE_ABUNDANCE_VIRALFLYE_LINEAR {
    tag "${sample}_viralFlye_linear"
    label 'process_low'
    conda 'conda-forge::python=3.10 bioconda::samtools=1.17 conda-forge::biopython=1.81'
    publishDir "${params.outdir_long}/abundance_viralflye_linear", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(bam), path(bai), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_viralflye_linear_abundance.txt"), emit: abundance
    path("${sample}_viralflye_linear_abundance_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env bash
    mkdir -p \$HOME/.local/lib_tmp
    for libdir in /usr/lib /usr/lib64 /lib /lib64 /usr/lib/x86_64-linux-gnu; do
        [ -f "\$libdir/libbz2.so.1" ] && ln -sf "\$libdir/libbz2.so.1" \$HOME/.local/lib_tmp/libbz2.so.1.0 && break
        [ -f "\$libdir/libbz2.so" ] && ln -sf "\$libdir/libbz2.so" \$HOME/.local/lib_tmp/libbz2.so.1.0 && break
    done
    export LD_LIBRARY_PATH=\$HOME/.local/lib_tmp:\$LD_LIBRARY_PATH
    
    python3 << 'PYTHON_SCRIPT'
import subprocess
from Bio import SeqIO

contig_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse("${contigs}", "fasta")}
idxstats = subprocess.check_output(["samtools", "idxstats", "${bam}"], universal_newlines=True, stderr=subprocess.PIPE)
contig_reads = {}
total_mapped_reads = 0
for line in idxstats.strip().split('\\n'):
    parts = line.split('\\t')
    if len(parts) >= 3 and parts[0] != "*":
        contig_reads[parts[0]] = int(parts[2])
        total_mapped_reads += int(parts[2])

with open("${sample}_viralflye_linear_abundance.txt", 'w') as f:
    f.write("Contig_ID\\tLength(bp)\\tMapped_Reads\\tRPM\\tRPKM\\n")
    for cname in sorted(contig_lengths.keys()):
        l, r = contig_lengths[cname], contig_reads.get(cname, 0)
        rpm = (r / total_mapped_reads * 1e6) if total_mapped_reads > 0 else 0
        rpkm = (r / (l / 1000) / (total_mapped_reads / 1e6)) if total_mapped_reads > 0 and l > 0 else 0
        f.write(f"{cname}\\t{l}\\t{r}\\t{rpm:.4f}\\t{rpkm:.4f}\\n")

with open("${sample}_viralflye_linear_abundance_summary.txt", 'w') as f:
    f.write("="*80 + "\\nviralFlye Linear Viral Contigs Summary\\n" + "="*80 + "\\n\\n")
    f.write(f"Sample: ${sample}\\nTotal viral contigs: {len(contig_lengths)}\\nTotal mapped reads: {total_mapped_reads:,}\\n")
    if contig_lengths:
        f.write(f"Avg length: {sum(contig_lengths.values())/len(contig_lengths):.2f} bp\\n")
    f.write("\\n" + "="*80 + "\\n")
PYTHON_SCRIPT
    """
}

// Process: Calculate RPM and RPKM for viralFlye circular viral contigs
process CALCULATE_ABUNDANCE_VIRALFLYE_CIRCULAR {
    tag "${sample}_viralFlye_circular"
    label 'process_low'
    conda 'conda-forge::python=3.10 bioconda::samtools=1.17 conda-forge::biopython=1.81'
    publishDir "${params.outdir_long}/abundance_viralflye_circular", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(bam), path(bai), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_viralflye_circular_abundance.txt"), emit: abundance
    path("${sample}_viralflye_circular_abundance_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env bash
    mkdir -p \$HOME/.local/lib_tmp
    for libdir in /usr/lib /usr/lib64 /lib /lib64 /usr/lib/x86_64-linux-gnu; do
        [ -f "\$libdir/libbz2.so.1" ] && ln -sf "\$libdir/libbz2.so.1" \$HOME/.local/lib_tmp/libbz2.so.1.0 && break
        [ -f "\$libdir/libbz2.so" ] && ln -sf "\$libdir/libbz2.so" \$HOME/.local/lib_tmp/libbz2.so.1.0 && break
    done
    export LD_LIBRARY_PATH=\$HOME/.local/lib_tmp:\$LD_LIBRARY_PATH
    
    python3 << 'PYTHON_SCRIPT'
import subprocess
from Bio import SeqIO

contig_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse("${contigs}", "fasta")}
idxstats = subprocess.check_output(["samtools", "idxstats", "${bam}"], universal_newlines=True, stderr=subprocess.PIPE)
contig_reads = {}
total_mapped_reads = 0
for line in idxstats.strip().split('\\n'):
    parts = line.split('\\t')
    if len(parts) >= 3 and parts[0] != "*":
        contig_reads[parts[0]] = int(parts[2])
        total_mapped_reads += int(parts[2])

with open("${sample}_viralflye_circular_abundance.txt", 'w') as f:
    f.write("Contig_ID\\tLength(bp)\\tMapped_Reads\\tRPM\\tRPKM\\n")
    for cname in sorted(contig_lengths.keys()):
        l, r = contig_lengths[cname], contig_reads.get(cname, 0)
        rpm = (r / total_mapped_reads * 1e6) if total_mapped_reads > 0 else 0
        rpkm = (r / (l / 1000) / (total_mapped_reads / 1e6)) if total_mapped_reads > 0 and l > 0 else 0
        f.write(f"{cname}\\t{l}\\t{r}\\t{rpm:.4f}\\t{rpkm:.4f}\\n")

with open("${sample}_viralflye_circular_abundance_summary.txt", 'w') as f:
    f.write("="*80 + "\\nviralFlye Circular Viral Contigs Summary\\n" + "="*80 + "\\n\\n")
    f.write(f"Sample: ${sample}\\nTotal viral contigs: {len(contig_lengths)}\\nTotal mapped reads: {total_mapped_reads:,}\\n")
    if contig_lengths:
        f.write(f"Avg length: {sum(contig_lengths.values())/len(contig_lengths):.2f} bp\\n")
    f.write("\\n" + "="*80 + "\\n")
PYTHON_SCRIPT
    """
}

// Process: Kraken2 Classification for metaFlye contigs
process KRAKEN2_CLASSIFICATION_FLYE {
    tag "${sample}_metaFlye"
    label 'process_medium'
    conda 'bioconda::kraken2=2.1.3'
    publishDir "${params.outdir_long}/kraken2_flye", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs)
    val(kraken2_db)
    
    output:
    tuple val(sample), path("${sample}_flye_*.txt"), emit: kraken2_flye
    
    script:
    """
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${sample}_flye_classification.txt \\
        --report ${sample}_flye_report.txt \\
        ${contigs}
    """
}

// Process: Kraken2 Classification for viralFlye linear viral contigs
process KRAKEN2_CLASSIFICATION_VIRALFLYE_LINEAR {
    tag "${sample}_viralFlye_linear"
    label 'process_medium'
    conda 'bioconda::kraken2=2.1.3'
    publishDir "${params.outdir_long}/kraken2_viralflye_linear", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs)
    val(kraken2_db)
    
    output:
    tuple val(sample), path("${sample}_viralflye_linear_*.txt"), emit: kraken2_linear
    
    script:
    """
    echo "=== Kraken2 classification for linear viral contigs: ${sample} ==="
    
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${sample}_viralflye_linear_classification.txt \\
        --report ${sample}_viralflye_linear_report.txt \\
        ${contigs}
    
    echo "Kraken2 classification completed for linear viral contigs"
    """
}

// Process: Kraken2 Classification for viralFlye circular viral contigs
process KRAKEN2_CLASSIFICATION_VIRALFLYE_CIRCULAR {
    tag "${sample}_viralFlye_circular"
    label 'process_medium'
    conda 'bioconda::kraken2=2.1.3'
    publishDir "${params.outdir_long}/kraken2_viralflye_circular", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs)
    val(kraken2_db)
    
    output:
    tuple val(sample), path("${sample}_viralflye_circular_*.txt"), emit: kraken2_circular
    
    script:
    """
    echo "=== Kraken2 classification for circular viral contigs: ${sample} ==="
    
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${sample}_viralflye_circular_classification.txt \\
        --report ${sample}_viralflye_circular_report.txt \\
        ${contigs}
    
    echo "Kraken2 classification completed for circular viral contigs"
    """
}

// Workflow completion message
workflow.onComplete {
    log.info """
    ==========================================
    🎯 Hybrid Metagenome Assembly Results
    ==========================================
    Pipeline completed successfully!
    
    ${params.input_short ? "Short-read results: ${params.outdir_short}" : ""}
    ${params.input_long ? "Long-read results: ${params.outdir_long}" : ""}
    ==========================================
    """
}

workflow.onError {
    log.error """
    ==========================================
    ❌ Hybrid Metagenome Workflow Failed
    ==========================================
    Error: ${workflow.errorMessage}
    ==========================================
    """
}
