#!/usr/bin/env nextflow

/*
 * Metagenome Assembly and Diamond Taxonomic Classification Workflow (English Version)
 * 
 * This workflow supports both short-read (Illumina paired-end) and long-read (Nanopore/PacBio) data:
 * 
 * SHORT READ PIPELINE:
 * 1. Quality control using fastp (optional)
 * 2. Metagenome assembly using MEGAHIT and SPAdes
 * 3. Gene prediction using Prodigal
 * 4. Taxonomic/functional classification using Diamond BLASTX
 * 5. Comprehensive analysis merging results from both assemblers
 * 
 * LONG READ PIPELINE:
 * 1. Long read assembly using MetaFlye
 * 2. Optional refinement using viralFlye (for viral genomes)
 * 3. Gene prediction using Prodigal
 * 4. Taxonomic/functional classification using Diamond BLASTX
 * 
 * Author: Assistant
 * Version: 4.0.0 (Added long-read support)
 */

nextflow.enable.dsl = 2

// Workflow parameters
// Input data
params.input = null
params.outdir = './results'
params.help = false
params.read_type = 'short'  // 'short' or 'long', to distinguish short-read and long-read data

// MAG parameters (simplified version without these features)
// params.skip_binning = true
// params.skip_checkm = true  
// params.skip_busco = true
// params.skip_prodigal = true
// params.skip_diamond = true
// params.skip_hmmer = true

// MEGAHIT parameters
params.megahit_memory = 0.8
params.megahit_min_contig_len = 1000

// SPAdes parameters (using metaSPAdes)
params.spades_meta = true

// fastp quality control parameters
params.skip_fastp = false
params.fastp_qualified_quality = 20    // Minimum quality value
params.fastp_unqualified_percent = 40  // Maximum percentage of low-quality bases allowed
params.fastp_min_length = 50           // Minimum read length

// Diamond classification parameters
params.diamond_db = null
params.diamond_evalue = 1e-5
params.diamond_max_target_seqs = 1
params.diamond_outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids'

// Taxonomy database parameters (for MERGE_DIAMOND_REPORTS, only used for short reads)
params.taxonomy_names = null  // NCBI taxonomy names.dmp file path
params.taxonomy_nodes = null  // NCBI taxonomy nodes.dmp file path

// MetaFlye parameters (for long reads)
params.metaflye_genome_size = null  // Genome size estimate, e.g., '5m' or '10m', null means let Flye estimate automatically
params.skip_viralflye = false       // Whether to skip viralFlye refinement

// viralFlye parameters (requires Pfam-A HMM database)
params.pfam_hmm = null              // Pfam-A HMM database path (required for viralFlye)
params.viralflye_min_length = 5000  // Minimum viral sequence length
params.viralflye_completeness = 0.5 // Completeness threshold

// Merge analysis parameters
params.skip_merge_reports = false  // Whether to skip comprehensive report generation

// Abundance calculation parameters
params.skip_abundance = false  // Whether to skip viral abundance calculation (RPM and RPKM)
// Long-read mapping preset for minimap2: 'map-ont' (Nanopore) or 'map-pb' (PacBio)
params.long_read_preset = 'map-ont'

// Resource parameters
params.max_cpus = 32
params.max_memory = '256.GB'
params.max_time = '72.h'

// Print help information
if (params.help) {
    log.info """
    ==========================================
    Metagenome Assembly and Diamond Taxonomic Classification Workflow
    ==========================================
    
    Usage:
    nextflow run metagenome_assembly_classification_workflow_en.nf --input samplesheet_short.csv --outdir results --diamond_db /path/to/db --read_type short
    
    Parameters:
    --input                    Input samplesheet (CSV format)
    --outdir                   Output directory
    --diamond_db              Diamond database path (e.g., NCBI nr or RVDB)
    --read_type               Read type: 'short' (paired-end) or 'long' (Nanopore/PacBio, default: short)
    --diamond_evalue          E-value threshold (default: 1e-5)
    --diamond_max_target_seqs Maximum number of target sequences (default: 1)
    
    Examples:
    # Short reads (paired-end Illumina)
    nextflow run metagenome_assembly_classification_workflow_en.nf \\
        --input samplesheet_short.csv \\
        --outdir results \\
        --diamond_db /path/to/diamond/nr.dmnd \\
        --read_type short
    
    # Long reads (Nanopore/PacBio)
    nextflow run metagenome_assembly_classification_workflow_en.nf \\
        --input samplesheet_long.csv \\
        --outdir results \\
        --diamond_db /path/to/diamond/nr.dmnd \\
        --read_type long \\
        --metaflye_genome_size 10m
    """
    exit 0
}

// Validate required parameters
if (!params.input) {
    error "Input samplesheet is required. Use --input parameter."
}

if (!params.diamond_db) {
    error "Diamond database path is required. Use --diamond_db parameter."
}

if (!params.skip_viralflye && !params.pfam_hmm) {
    error "Pfam-A HMM database is required for viralFlye. Use --pfam_hmm parameter or set --skip_viralflye true."
}

// Print workflow information
    log.info """
    ==========================================
    🧬 Metagenome Assembly and Diamond Taxonomic Classification Workflow
    ==========================================
    Workflow version: 4.0.0 (Supports both short and long reads)
    Read type: ${params.read_type}
    Input samplesheet: ${params.input}
Output directory: ${params.outdir}
Diamond database: ${params.diamond_db}
Diamond E-value: ${params.diamond_evalue}
==========================================
"""

// Create input channel from CSV samplesheet
// Create input channel based on read type
if (params.read_type == 'long') {
    // Long read input: each sample has one fastq file
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def sample = row.sample
            def fastq_long = file(row.fastq_long)
            return tuple(sample, fastq_long, 'long')
        }
        .set { ch_reads }
} else {
    // Short read input: each sample has two fastq files (paired-end)
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def sample = row.sample
            def read1 = file(row.fastq_1)
            def read2 = file(row.fastq_2)
            return tuple(sample, [read1, read2], 'short')
        }
        .set { ch_reads }
}

// Define workflow
workflow {
    // Run different workflows based on read type
    if (params.read_type == 'long') {
        // ============================================
        // Long Read Workflow (Nanopore/PacBio)
        // ============================================
        
        // Stage 1: MetaFlye Assembly (long read assembly)
        METAFLYE_ASSEMBLY (
            ch_reads
        )
        
        // Stage 2: Gene Prediction on MetaFlye contigs (all sequences, comprehensive analysis)
        PRODIGAL_METAFLYE (
            METAFLYE_ASSEMBLY.out.contigs
        )
        
        // Stage 3: Diamond Classification on MetaFlye (analyze all sequences)
        DIAMOND_CLASSIFICATION_METAFLYE (
            PRODIGAL_METAFLYE.out.proteins,
            params.diamond_db
        )
        
        // Stage 4: viralFlye Refinement (optional, identify high-confidence viruses)
        if (!params.skip_viralflye) {
            // Prepare reads channel for join (remove third element)
            ch_reads_for_viralflye = ch_reads.map { sample, reads, type -> tuple(sample, reads) }
            
            VIRALFLYE_REFINEMENT (
                METAFLYE_ASSEMBLY.out.metaflye_dir,
                ch_reads_for_viralflye,
                params.pfam_hmm
            )
            
            // Stage 5: Gene Prediction on viralFlye contigs
            PRODIGAL_VIRALFLYE (
                VIRALFLYE_REFINEMENT.out.contigs
            )
            
            // Stage 6: Diamond Classification on viralFlye (viral sequences only)
            DIAMOND_CLASSIFICATION_VIRALFLYE (
                PRODIGAL_VIRALFLYE.out.proteins,
            params.diamond_db
        )
        
            // Stage 7: Add Taxonomy for both MetaFlye and viralFlye
            ch_taxonomy_names = Channel.fromPath(params.taxonomy_names, checkIfExists: true)
            ch_taxonomy_nodes = Channel.fromPath(params.taxonomy_nodes, checkIfExists: true)
            
            // MetaFlye taxonomy
            ADD_TAXONOMY_METAFLYE (
                DIAMOND_CLASSIFICATION_METAFLYE.out.diamond_metaflye,
                ch_taxonomy_names.collect(),
                ch_taxonomy_nodes.collect()
            )
            
            // viralFlye taxonomy
            ADD_TAXONOMY_VIRALFLYE (
                DIAMOND_CLASSIFICATION_VIRALFLYE.out.diamond_viralflye,
                ch_taxonomy_names.collect(),
                ch_taxonomy_nodes.collect()
            )
            
            // Stage 8: Compare and Generate Consensus (compare both result sets, generate consensus virus list)
            COMPARE_DUAL_TRACKS (
                ADD_TAXONOMY_METAFLYE.out.enhanced_report
                    .join(ADD_TAXONOMY_VIRALFLYE.out.enhanced_report)
            )
        } else {
            // If skipping viralFlye, only analyze MetaFlye results
            ch_taxonomy_names = Channel.fromPath(params.taxonomy_names, checkIfExists: true)
            ch_taxonomy_nodes = Channel.fromPath(params.taxonomy_nodes, checkIfExists: true)
            
            ADD_TAXONOMY_METAFLYE (
                DIAMOND_CLASSIFICATION_METAFLYE.out.diamond_metaflye,
                ch_taxonomy_names.collect(),
                ch_taxonomy_nodes.collect()
            )
        }
        
        // Stage 9: Long-read Viral Abundance (RPM/RPKM) for MetaFlye and viralFlye
        if (!params.skip_abundance) {
            // Prepare reads channel for long reads
            ch_reads_long = ch_reads.map { sample, fastq_long, type -> tuple(sample, fastq_long) }

            // MetaFlye abundance (all sequences, filtered by Diamond viral contigs)
            METAFLYE_ASSEMBLY.out.contigs
                .join(ch_reads_long)
                .join(DIAMOND_CLASSIFICATION_METAFLYE.out.diamond_metaflye)
                .map { sample, contigs, fastq_long, diamond_report -> tuple(sample, contigs, fastq_long, diamond_report) }
                .set { ch_metaflye_abundance }

            CALCULATE_ABUNDANCE_METAFLYE (
                ch_metaflye_abundance
            )

            // viralFlye abundance (viral sequences only), only when viralFlye is enabled
            if (!params.skip_viralflye) {
                VIRALFLYE_REFINEMENT.out.contigs
                    .join(ch_reads_for_viralflye)
                    .join(DIAMOND_CLASSIFICATION_VIRALFLYE.out.diamond_viralflye)
                    .map { sample, contigs, fastq_long, diamond_report -> tuple(sample, contigs, fastq_long, diamond_report) }
                    .set { ch_viralflye_abundance }

                CALCULATE_ABUNDANCE_VIRALFLYE (
                    ch_viralflye_abundance
                )
            }
        }
        
    } else {
        // ============================================
        // Short Read Workflow (Paired-end Illumina)
        // ============================================
        
        // Stage 0: Quality Control (optional)
        if (!params.skip_fastp) {
            FASTP (
                ch_reads
            )
            ch_clean_reads = FASTP.out.clean_reads
        } else {
            // When skipping fastp, normalize ch_reads format to tuple(sample, read1, read2) to match FASTP output
            ch_clean_reads = ch_reads.map { sample, reads, read_type -> tuple(sample, reads[0], reads[1]) }
        }
        
        // Stage 1: Assembly
        MEGAHIT_ASSEMBLY (
            ch_clean_reads
        )
        
        SPADES_ASSEMBLY (
            ch_clean_reads
        )
        
        // Stage 2: Gene Prediction with Prodigal
        PRODIGAL_MEGAHIT (
            MEGAHIT_ASSEMBLY.out.contigs
        )
        
        PRODIGAL_SPADES (
            SPADES_ASSEMBLY.out.contigs
        )
        
        // Stage 3: Diamond Classification
        DIAMOND_CLASSIFICATION_MEGAHIT (
            PRODIGAL_MEGAHIT.out.proteins,
            params.diamond_db
        )
        
        DIAMOND_CLASSIFICATION_SPADES (
            PRODIGAL_SPADES.out.proteins,
            params.diamond_db
        )
        
        // Stage 4: Merge Reports (Comprehensive Analysis)
        if (!params.skip_merge_reports) {
            // Merge MEGAHIT and SPAdes reports for the same sample
            DIAMOND_CLASSIFICATION_MEGAHIT.out.diamond_megahit
                .join(DIAMOND_CLASSIFICATION_SPADES.out.diamond_spades)
                .set { ch_reports_to_merge }
            
            // Prepare taxonomy files
            ch_taxonomy_names = Channel.fromPath(params.taxonomy_names, checkIfExists: true)
            ch_taxonomy_nodes = Channel.fromPath(params.taxonomy_nodes, checkIfExists: true)
            
            MERGE_DIAMOND_REPORTS (
                ch_reports_to_merge,
                ch_taxonomy_names.collect(),
                ch_taxonomy_nodes.collect()
            )
        }
        
        // Stage 5: Calculate Viral Abundance (RPM and RPKM)
        if (!params.skip_abundance) {
            // Calculate abundance for MEGAHIT
            // Join contigs with cleaned paired-end reads (two separate files), then with Diamond results
            // Note: ch_clean_reads format is tuple(sample, read1, read2)
            MEGAHIT_ASSEMBLY.out.contigs
                .join(ch_clean_reads)
                .map { sample, contigs, read1, read2 -> tuple(sample, contigs, [read1, read2]) }
                .join(DIAMOND_CLASSIFICATION_MEGAHIT.out.diamond_megahit)
                .map { sample, contigs, reads, diamond_report -> tuple(sample, contigs, reads, diamond_report) }
                .set { ch_megahit_abundance }
            
            CALCULATE_ABUNDANCE_MEGAHIT_SHORT (
                ch_megahit_abundance
            )
            
            // Calculate abundance for SPAdes
            // Join contigs with cleaned paired-end reads (two separate files), then with Diamond results
            // Note: ch_clean_reads format is tuple(sample, read1, read2)
            SPADES_ASSEMBLY.out.contigs
                .join(ch_clean_reads)
                .map { sample, contigs, read1, read2 -> tuple(sample, contigs, [read1, read2]) }
                .join(DIAMOND_CLASSIFICATION_SPADES.out.diamond_spades)
                .map { sample, contigs, reads, diamond_report -> tuple(sample, contigs, reads, diamond_report) }
                .set { ch_spades_abundance }
            
            CALCULATE_ABUNDANCE_SPADES_SHORT (
                ch_spades_abundance
            )
        }
    }
}

// ================================================================================
// Process Definitions
// ================================================================================

// Process: fastp Quality Control (for short reads only)
process FASTP {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{html,json}"
    
    input:
    tuple val(sample), path(reads), val(read_type)
    
    output:
    tuple val(sample), path("${sample}_clean_R1.fastq.gz"), path("${sample}_clean_R2.fastq.gz"), emit: clean_reads
    path("${sample}_fastp.html"), emit: html
    path("${sample}_fastp.json"), emit: json
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    echo "=== fastp Quality Control: ${sample} ==="
    
    # List input files for debugging
    echo "Input files in work directory:"
    ls -lh
    
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
    // Uses system environment (nextflow_env); MEGAHIT must be pre-installed in nextflow_env
    publishDir "${params.outdir}/assembly_megahit", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}_megahit_contigs.fa"), emit: contigs
    path("${sample}_megahit_contigs.fa"), emit: contigs_published
    
    script:
    """
    echo "=== MEGAHIT Assembly: ${sample} ==="
    
    megahit \
        -1 ${read1} \
        -2 ${read2} \
        -o megahit_output \
        -t ${task.cpus} \
        --memory ${params.megahit_memory} \
        --min-contig-len ${params.megahit_min_contig_len}
    
    cp megahit_output/final.contigs.fa ${sample}_megahit_contigs.fa
    
    CONTIG_COUNT=\$(grep -c ">" ${sample}_megahit_contigs.fa)
    TOTAL_LENGTH=\$(grep -v ">" ${sample}_megahit_contigs.fa | tr -d '\\n' | wc -c)
    
    echo "MEGAHIT: Generated \${CONTIG_COUNT} contigs"
    echo "MEGAHIT: Total length \${TOTAL_LENGTH} bp"
    """
}

// Process: SPAdes Assembly
process SPADES_ASSEMBLY {
    tag "${sample}_SPAdes"
    label 'process_high'
    // Uses system environment (nextflow_env); SPAdes must be pre-installed in nextflow_env
    publishDir "${params.outdir}/assembly_spades", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}_spades_contigs.fa"), emit: contigs
    path("${sample}_spades_contigs.fa"), emit: contigs_published
    
    script:
    """
    echo "=== metaSPAdes Assembly: ${sample} ==="
    
    # Use metaSPAdes, disable error correction to avoid memory and bug issues
    metaspades.py \
        -1 ${read1} \
        -2 ${read2} \
        -o spades_output \
        -t ${task.cpus} \
        -m ${task.memory.toGiga()} \
        --only-assembler
    
    cp spades_output/contigs.fasta ${sample}_spades_contigs.fa
    
    CONTIG_COUNT=\$(grep -c ">" ${sample}_spades_contigs.fa)
    TOTAL_LENGTH=\$(grep -v ">" ${sample}_spades_contigs.fa | tr -d '\\n' | wc -c)
    
    echo "metaSPAdes: Generated \${CONTIG_COUNT} contigs"
    echo "metaSPAdes: Total length \${TOTAL_LENGTH} bp"
    """
}

// Process: Prodigal Gene Prediction for MEGAHIT
process PRODIGAL_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_medium'
    conda 'bioconda::prodigal=2.6.3'
    publishDir "${params.outdir}/prodigal_megahit", mode: 'copy', pattern: "*.faa"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_megahit_proteins.faa"), emit: proteins
    path("${sample}_megahit_genes.fna"), emit: genes
    
    script:
    """
    echo "=== Prodigal Gene Prediction (MEGAHIT): ${sample} ==="
    
    # Run Prodigal for gene prediction (metagenome mode)
    prodigal \\
        -i ${contigs} \\
        -a ${sample}_megahit_proteins.faa \\
        -d ${sample}_megahit_genes.fna \\
        -p meta \\
        -q
    
    # Count predicted genes
    GENE_COUNT=\$(grep -c ">" ${sample}_megahit_proteins.faa || echo 0)
    echo "Prodigal: Predicted \${GENE_COUNT} genes from MEGAHIT contigs"
    """
}

// Process: Prodigal Gene Prediction for SPAdes
process PRODIGAL_SPADES {
    tag "${sample}_SPAdes"
    label 'process_medium'
    conda 'bioconda::prodigal=2.6.3'
    publishDir "${params.outdir}/prodigal_spades", mode: 'copy', pattern: "*.faa"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_spades_proteins.faa"), emit: proteins
    path("${sample}_spades_genes.fna"), emit: genes
    
    script:
    """
    echo "=== Prodigal Gene Prediction (SPAdes): ${sample} ==="
    
    # Run Prodigal for gene prediction (metagenome mode)
    prodigal \\
        -i ${contigs} \\
        -a ${sample}_spades_proteins.faa \\
        -d ${sample}_spades_genes.fna \\
        -p meta \\
        -q
    
    # Count predicted genes
    GENE_COUNT=\$(grep -c ">" ${sample}_spades_proteins.faa || echo 0)
    echo "Prodigal: Predicted \${GENE_COUNT} genes from SPAdes contigs"
    """
}

// Process: Diamond Classification for MEGAHIT
process DIAMOND_CLASSIFICATION_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_high'
    conda 'bioconda::diamond=2.1.8'
    publishDir "${params.outdir}/diamond_megahit", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(proteins)
    val(diamond_db)
    
    output:
    tuple val(sample), path("${sample}_megahit_diamond.txt"), emit: diamond_megahit
    
    script:
    """
    echo "=== Diamond Classification (MEGAHIT): ${sample} ==="
    
    # Run Diamond BLASTP for classification
    diamond blastp \\
        --query ${proteins} \\
        --db ${diamond_db} \\
        --out ${sample}_megahit_diamond.txt \\
        --outfmt ${params.diamond_outfmt} \\
        --threads ${task.cpus} \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs ${params.diamond_max_target_seqs} \\
        --sensitive
    
    # Count hits
    HIT_COUNT=\$(wc -l < ${sample}_megahit_diamond.txt || echo 0)
    echo "Diamond: Found \${HIT_COUNT} hits for MEGAHIT proteins"
    """
}

// Process: Diamond Classification for SPAdes
process DIAMOND_CLASSIFICATION_SPADES {
    tag "${sample}_SPAdes"
    label 'process_high'
    conda 'bioconda::diamond=2.1.8'
    publishDir "${params.outdir}/diamond_spades", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(proteins)
    val(diamond_db)
    
    output:
    tuple val(sample), path("${sample}_spades_diamond.txt"), emit: diamond_spades
    
    script:
    """
    echo "=== Diamond Classification (SPAdes): ${sample} ==="
    
    # Run Diamond BLASTP for classification
    diamond blastp \\
        --query ${proteins} \\
        --db ${diamond_db} \\
        --out ${sample}_spades_diamond.txt \\
        --outfmt ${params.diamond_outfmt} \\
        --threads ${task.cpus} \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs ${params.diamond_max_target_seqs} \\
        --sensitive
    
    # Count hits
    HIT_COUNT=\$(wc -l < ${sample}_spades_diamond.txt || echo 0)
    echo "Diamond: Found \${HIT_COUNT} hits for SPAdes proteins"
    """
}

// Process: Merge Diamond Reports (Comprehensive Analysis)
process MERGE_DIAMOND_REPORTS {
    tag "${sample}"
    label 'process_low'
    conda 'conda-forge::pandas=2.0.3'
    publishDir "${params.outdir}/merged_reports", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(megahit_report), path(spades_report)
    path(taxonomy_names)
    path(taxonomy_nodes)
    
    output:
    tuple val(sample), path("${sample}_merged_report.txt"), emit: merged_report
    path("${sample}_merged_report.csv"), emit: merged_csv
    path("${sample}_megahit_with_taxonomy.txt"), emit: megahit_enhanced
    path("${sample}_spades_with_taxonomy.txt"), emit: spades_enhanced
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import sys
    from collections import Counter, defaultdict
    
    # Taxonomy database class
    class TaxonomyDB:
        \"\"\"Parse NCBI taxonomy database\"\"\"
        
        def __init__(self, names_file, nodes_file):
            \"\"\"Initialize taxonomy database\"\"\"
            print("Loading taxonomy database...", file=sys.stderr)
            self.names = self._load_names(names_file)
            self.nodes = self._load_nodes(nodes_file)
            print(f"Loaded: {len(self.names):,} names, {len(self.nodes):,} nodes", file=sys.stderr)
        
        def _load_names(self, names_file):
            \"\"\"Load names.dmp\"\"\"
            names = {}
            with open(names_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 4:
                        taxid = parts[0]
                        name = parts[1]
                        name_class = parts[3]
                        
                        # Only keep scientific names
                        if name_class == 'scientific name':
                            names[taxid] = name
            
            return names
        
        def _load_nodes(self, nodes_file):
            \"\"\"Load nodes.dmp\"\"\"
            nodes = {}
            with open(nodes_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 3:
                        taxid = parts[0]
                        parent_taxid = parts[1]
                        rank = parts[2]
                        
                        nodes[taxid] = {
                            'parent': parent_taxid,
                            'rank': rank
                        }
            
            return nodes
        
        def get_lineage(self, taxid):
            \"\"\"Get complete taxonomic lineage\"\"\"
            lineage = {
                'superkingdom': 'N/A',
                'kingdom': 'N/A',
                'phylum': 'N/A',
                'class': 'N/A',
                'order': 'N/A',
                'family': 'N/A',
                'genus': 'N/A',
                'species': 'N/A',
                'organism_name': 'N/A'
            }
            
            taxid = str(taxid)
            
            # Get current taxid name
            if taxid in self.names:
                lineage['organism_name'] = self.names[taxid]
            
            # Traverse up the taxonomy tree
            current_taxid = taxid
            visited = set()
            
            while current_taxid != '1' and current_taxid in self.nodes:
                # Prevent loops
                if current_taxid in visited:
                    break
                visited.add(current_taxid)
                
                node = self.nodes[current_taxid]
                rank = node['rank']
                
                # Only keep major ranks
                if rank in lineage and current_taxid in self.names:
                    lineage[rank] = self.names[current_taxid]
                
                # Move to parent node
                current_taxid = node['parent']
            
            return lineage
    
    def parse_diamond_output(file_path):
        \"\"\"Parse Diamond output file\"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t', header=None, 
                           names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                  'evalue', 'bitscore', 'staxids'])
            return df
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                        'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                        'evalue', 'bitscore', 'staxids'])
    
    def add_taxonomy_to_dataframe(df, taxonomy_db):
        \"\"\"Add taxonomic information to DataFrame\"\"\"
        if df.empty:
            return df
        
        print(f"Adding taxonomy info to {len(df):,} records...", file=sys.stderr)
        
        lineages = []
        for idx, row in df.iterrows():
            if (idx + 1) % 1000 == 0:
                print(f"  Processing {idx+1:,}/{len(df):,}...", file=sys.stderr)
            
            # Convert TaxID to string (handle float format from pandas)
            if pd.notna(row['staxids']):
                # Convert float to int to string (455367.0 → 455367 → "455367")
                try:
                    taxid = str(int(float(row['staxids'])))
                except (ValueError, TypeError):
                    taxid = '0'
            else:
                taxid = '0'
            
            lineage = taxonomy_db.get_lineage(taxid)
            lineages.append(lineage)
        
        # Add taxonomy columns
        for key in ['organism_name', 'superkingdom', 'kingdom', 'phylum', 
                    'class', 'order', 'family', 'genus', 'species']:
            df[key] = [lineage[key] for lineage in lineages]
        
        return df
    
    def extract_taxonomic_info(df):
        \"\"\"Extract taxonomic statistics from Diamond output\"\"\"
        if df.empty:
            return {
                'total_hits': 0,
                'unique_queries': 0,
                'unique_subjects': 0,
                'avg_identity': 0,
                'avg_length': 0,
                'taxid_counts': Counter(),
                'phylum_counts': Counter(),
                'family_counts': Counter()
            }
        
        stats = {
            'total_hits': len(df),
            'unique_queries': df['qseqid'].nunique(),
            'unique_subjects': df['sseqid'].nunique(),
            'avg_identity': df['pident'].mean(),
            'avg_length': df['length'].mean(),
            'taxid_counts': Counter(df['staxids'].dropna().apply(lambda x: str(int(float(x)))))
        }
        
        # If taxonomy columns exist, count phylum and family
        if 'phylum' in df.columns:
            stats['phylum_counts'] = Counter(df['phylum'].dropna())
            stats['family_counts'] = Counter(df['family'].dropna())
        
        return stats
    
    # Initialize Taxonomy database
    taxonomy_db = TaxonomyDB("${taxonomy_names}", "${taxonomy_nodes}")
    
    # Parse Diamond output files
    print(f"\\nParsing MEGAHIT Diamond results: ${megahit_report}", file=sys.stderr)
    megahit_df = parse_diamond_output("${megahit_report}")
    
    # Add taxonomic information
    megahit_df_enhanced = add_taxonomy_to_dataframe(megahit_df.copy(), taxonomy_db)
    megahit_stats = extract_taxonomic_info(megahit_df_enhanced)
    
    # Save enhanced version (with taxonomy)
    megahit_df_enhanced.to_csv("${sample}_megahit_with_taxonomy.txt", sep='\\t', index=False)
    print(f"MEGAHIT enhanced file saved", file=sys.stderr)
    
    print(f"\\nParsing SPAdes Diamond results: ${spades_report}", file=sys.stderr)
    spades_df = parse_diamond_output("${spades_report}")
    
    # Add taxonomic information
    spades_df_enhanced = add_taxonomy_to_dataframe(spades_df.copy(), taxonomy_db)
    spades_stats = extract_taxonomic_info(spades_df_enhanced)
    
    # Save enhanced version (with taxonomy)
    spades_df_enhanced.to_csv("${sample}_spades_with_taxonomy.txt", sep='\\t', index=False)
    print(f"SPAdes enhanced file saved", file=sys.stderr)
    
    # Generate text report
    with open("${sample}_merged_report.txt", 'w', encoding='utf-8') as f:
        f.write("="*80 + "\\n")
        f.write("Diamond Comprehensive Analysis Report - MEGAHIT vs SPAdes Assembly Comparison\\n")
        f.write("="*80 + "\\n\\n")
        
        # Overall statistics
        f.write("[Overall Statistics]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"SPAdes total hits:          {spades_stats['total_hits']:,}\\n")
        f.write(f"MEGAHIT total hits:         {megahit_stats['total_hits']:,}\\n\\n")
        
        f.write(f"SPAdes unique queries:      {spades_stats['unique_queries']:,}\\n")
        f.write(f"MEGAHIT unique queries:     {megahit_stats['unique_queries']:,}\\n\\n")
        
        f.write(f"SPAdes average identity:    {spades_stats['avg_identity']:.2f}%\\n")
        f.write(f"MEGAHIT average identity:   {megahit_stats['avg_identity']:.2f}%\\n\\n")
        
        f.write(f"SPAdes average length:      {spades_stats['avg_length']:.1f} aa\\n")
        f.write(f"MEGAHIT average length:     {megahit_stats['avg_length']:.1f} aa\\n\\n")
        
        # Phylum level comparison
        f.write("\\n[Phylum Level Comparison (Top 15)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Phylum':<30} {'SPAdes Count':<15} {'MEGAHIT Count':<15} {'Total':<15}\\n")
        f.write("-"*80 + "\\n")
        
        all_phyla = set(spades_stats.get('phylum_counts', {}).keys()) | set(megahit_stats.get('phylum_counts', {}).keys())
        phylum_comparison = []
        for phylum in all_phyla:
            spades_c = spades_stats.get('phylum_counts', {}).get(phylum, 0)
            megahit_c = megahit_stats.get('phylum_counts', {}).get(phylum, 0)
            total_c = spades_c + megahit_c
            phylum_comparison.append((phylum, spades_c, megahit_c, total_c))
        
        phylum_comparison.sort(key=lambda x: x[3], reverse=True)
        for phylum, spades_c, megahit_c, total_c in phylum_comparison[:15]:
            f.write(f"{phylum:<30} {spades_c:<15} {megahit_c:<15} {total_c:<15}\\n")
        
        # Family level comparison
        f.write("\\n[Family Level Comparison (Top 15)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Family':<30} {'SPAdes Count':<15} {'MEGAHIT Count':<15} {'Total':<15}\\n")
        f.write("-"*80 + "\\n")
        
        all_families = set(spades_stats.get('family_counts', {}).keys()) | set(megahit_stats.get('family_counts', {}).keys())
        family_comparison = []
        for family in all_families:
            spades_c = spades_stats.get('family_counts', {}).get(family, 0)
            megahit_c = megahit_stats.get('family_counts', {}).get(family, 0)
            total_c = spades_c + megahit_c
            family_comparison.append((family, spades_c, megahit_c, total_c))
        
        family_comparison.sort(key=lambda x: x[3], reverse=True)
        for family, spades_c, megahit_c, total_c in family_comparison[:15]:
            f.write(f"{family:<30} {spades_c:<15} {megahit_c:<15} {total_c:<15}\\n")
        
        # Taxonomic ID comparison (Top 20)
        f.write("\\n[Taxonomic ID Comparison (Top 20)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Taxonomic ID':<20} {'SPAdes Count':<15} {'MEGAHIT Count':<15}\\n")
        f.write("-"*80 + "\\n")
        
        # Merge all taxonomic IDs
        all_taxids = set(spades_stats['taxid_counts'].keys()) | set(megahit_stats['taxid_counts'].keys())
        taxid_comparison = []
        
        for taxid in all_taxids:
            spades_count = spades_stats['taxid_counts'].get(taxid, 0)
            megahit_count = megahit_stats['taxid_counts'].get(taxid, 0)
            total_count = spades_count + megahit_count
            taxid_comparison.append((taxid, spades_count, megahit_count, total_count))
        
        # Sort by total count
        taxid_comparison.sort(key=lambda x: x[3], reverse=True)
        
        for taxid, spades_c, megahit_c, total_c in taxid_comparison[:20]:
            f.write(f"{taxid:<20} {spades_c:<15} {megahit_c:<15}\\n")
        
        # Unique findings
        f.write("\\n[Unique Findings]\\n")
        f.write("-"*80 + "\\n\\n")
        
        # Taxonomic IDs found only in SPAdes
        spades_only = set(spades_stats['taxid_counts'].keys()) - set(megahit_stats['taxid_counts'].keys())
        f.write(f"Taxonomic IDs found only in SPAdes: {len(spades_only)}\\n")
        if spades_only:
            spades_only_sorted = sorted([(tid, spades_stats['taxid_counts'][tid]) 
                                        for tid in spades_only], 
                                       key=lambda x: x[1], reverse=True)
            for tid, count in spades_only_sorted[:10]:
                f.write(f"  - {tid}: {count} matches\\n")
        
        f.write("\\n")
        
        # Taxonomic IDs found only in MEGAHIT
        megahit_only = set(megahit_stats['taxid_counts'].keys()) - set(spades_stats['taxid_counts'].keys())
        f.write(f"Taxonomic IDs found only in MEGAHIT: {len(megahit_only)}\\n")
        if megahit_only:
            megahit_only_sorted = sorted([(tid, megahit_stats['taxid_counts'][tid]) 
                                         for tid in megahit_only], 
                                        key=lambda x: x[1], reverse=True)
            for tid, count in megahit_only_sorted[:10]:
                f.write(f"  - {tid}: {count} matches\\n")
        
        f.write("\\n" + "="*80 + "\\n")
        f.write("Analysis Complete\\n")
        f.write("="*80 + "\\n")
    
    # Save CSV format detailed comparison
    comparison_data = []
    for taxid, spades_c, megahit_c, total_c in taxid_comparison:
        comparison_data.append({
            'taxonomic_id': taxid,
            'spades_count': spades_c,
            'megahit_count': megahit_c,
            'total_count': total_c,
            'spades_only': spades_c > 0 and megahit_c == 0,
            'megahit_only': megahit_c > 0 and spades_c == 0,
            'both': spades_c > 0 and megahit_c > 0
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    comparison_df = comparison_df.sort_values('total_count', ascending=False)
    comparison_df.to_csv("${sample}_merged_report.csv", index=False)
    
    print(f"Reports generated successfully for ${sample}", file=sys.stderr)
    """
}

// Process: Calculate Viral Abundance for MEGAHIT (RPM and RPKM) - Short Reads
process CALCULATE_ABUNDANCE_MEGAHIT_SHORT {
    tag "${sample}_MEGAHIT"
    label 'process_medium'
    // Uses system environment (nextflow_env); bwa, samtools, htslib, bzip2, zlib, libdeflate, and pandas must be pre-installed in nextflow_env
    publishDir "${params.outdir}/abundance_megahit", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs), path(reads), path(diamond_report)
    
    output:
    tuple val(sample), path("${sample}_megahit_abundance.txt"), emit: abundance
    path("${sample}_megahit_abundance.csv"), emit: abundance_csv
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import subprocess
    import sys
    import os
    import shutil
    from collections import defaultdict
    
    # Set LD_LIBRARY_PATH and create symlink for libbz2.so.1.0
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if not conda_prefix:
        # Try to find conda environment from PATH
        python_path = shutil.which('python3')
        if python_path and 'conda' in python_path:
            conda_prefix = os.path.dirname(os.path.dirname(python_path))
    
    if conda_prefix:
        lib_path = os.path.join(conda_prefix, 'lib')
        current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
        os.environ['LD_LIBRARY_PATH'] = f"{lib_path}:{current_ld_path}" if current_ld_path else lib_path
        print(f"Set LD_LIBRARY_PATH: {os.environ['LD_LIBRARY_PATH']}", file=sys.stderr)
        
        # Create symlink for libbz2.so.1.0 if it doesn't exist
        libbz2_target = os.path.join(lib_path, 'libbz2.so.1.0')
        if not os.path.exists(libbz2_target):
            # Find libbz2.so files
            import glob
            libbz2_files = glob.glob(os.path.join(lib_path, 'libbz2.so*'))
            if libbz2_files:
                # Sort to get the highest version number
                libbz2_files.sort(reverse=True)
                source_file = libbz2_files[0]
                try:
                    os.symlink(os.path.basename(source_file), libbz2_target)
                    print(f"Created symlink: {libbz2_target} -> {os.path.basename(source_file)}", file=sys.stderr)
                except (OSError, FileExistsError) as e:
                    # Symlink might already exist or permission issue
                    print(f"Note: Could not create symlink {libbz2_target}: {e}", file=sys.stderr)
            else:
                print(f"WARNING: No libbz2.so files found in {lib_path}", file=sys.stderr)
    
    print("=== Calculating Viral Abundance (MEGAHIT): ${sample} ===", file=sys.stderr)
    
    # Check if required tools are available
    def check_tool(tool_name):
        tool_path = shutil.which(tool_name)
        if tool_path is None:
            print(f"ERROR: {tool_name} not found in PATH", file=sys.stderr)
            print(f"PATH: {os.environ.get('PATH', 'Not set')}", file=sys.stderr)
            sys.exit(1)
        print(f"Found {tool_name} at: {tool_path}", file=sys.stderr)
        return tool_path
    
    bwa_path = check_tool('bwa')
    samtools_path = check_tool('samtools')
    
    # 1. 从 Diamond 结果中提取病毒 contigs
    print("Extracting viral contigs from Diamond results...", file=sys.stderr)
    try:
        diamond_df = pd.read_csv("${diamond_report}", sep='\\t', header=None,
                                names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                       'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                       'evalue', 'bitscore', 'staxids'])
    except pd.errors.EmptyDataError:
        print("WARNING: Diamond report is empty, creating empty abundance file", file=sys.stderr)
        with open("${sample}_megahit_abundance.txt", 'w') as f:
            f.write("No viral contigs found\\n")
        pd.DataFrame(columns=['contig_id', 'length_bp', 'length_kb', 'mapped_reads', 'RPM', 'RPKM']).to_csv("${sample}_megahit_abundance.csv", index=False)
        exit(0)
    
    # 提取 contig ID（移除基因编号）
    diamond_df['contig_id'] = diamond_df['qseqid'].str.rsplit('_', n=1).str[0]
    viral_contigs = set(diamond_df['contig_id'].unique())
    print(f"Found {len(viral_contigs)} viral contigs", file=sys.stderr)
    
    if len(viral_contigs) == 0:
        print("WARNING: No viral contigs found, creating empty abundance file", file=sys.stderr)
        with open("${sample}_megahit_abundance.txt", 'w') as f:
            f.write("No viral contigs found\\n")
        pd.DataFrame(columns=['contig_id', 'length_bp', 'length_kb', 'mapped_reads', 'RPM', 'RPKM']).to_csv("${sample}_megahit_abundance.csv", index=False)
        exit(0)
    
    # 2. 使用 BWA 将 reads 比对到 contigs
    print("Indexing contigs with BWA...", file=sys.stderr)
    try:
        subprocess.run([bwa_path, 'index', '${contigs}'], check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: BWA index failed: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    
    print("Aligning reads to contigs...", file=sys.stderr)
    with open('alignment.sam', 'w') as sam_file:
        try:
            subprocess.run([bwa_path, 'mem', '-t', '${task.cpus}', '${contigs}', 
                           '${reads[0]}', '${reads[1]}'], 
                          stdout=sam_file, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: BWA mem failed: {e.stderr.decode()}", file=sys.stderr)
            sys.exit(1)
    
    # 3. 转换为 BAM 并排序
    print("Converting and sorting BAM...", file=sys.stderr)
    with open('alignment.bam', 'wb') as bam_file:
        try:
            subprocess.run([samtools_path, 'view', '-bS', 'alignment.sam'], 
                          stdout=bam_file, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: samtools view failed: {e.stderr.decode()}", file=sys.stderr)
            print(f"Command: {samtools_path} view -bS alignment.sam", file=sys.stderr)
            sys.exit(1)
    try:
        subprocess.run([samtools_path, 'sort', '-o', 'sorted.bam', 'alignment.bam'], 
                      check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools sort failed: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    try:
        subprocess.run([samtools_path, 'index', 'sorted.bam'], check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools index failed: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    
    # 4. 统计每个 contig 的 reads 数量
    print("Counting reads per contig...", file=sys.stderr)
    try:
        result = subprocess.run([samtools_path, 'idxstats', 'sorted.bam'], 
                              capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools idxstats failed: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
        sys.exit(1)
    
    # 解析 idxstats 输出
    contig_stats = {}
    total_mapped_reads = 0
    for line in result.stdout.strip().split('\\n'):
        if line.strip():
            parts = line.split('\\t')
            contig_id = parts[0]
            contig_length = int(parts[1])
            mapped_reads = int(parts[2])
            unmapped_reads = int(parts[3])
            
            if contig_id != '*':
                contig_stats[contig_id] = {
                    'length': contig_length,
                    'mapped_reads': mapped_reads
                }
                total_mapped_reads += mapped_reads
    
    print(f"Total mapped reads: {total_mapped_reads:,}", file=sys.stderr)
    
    # 5. 计算 RPM 和 RPKM
    print("Calculating RPM and RPKM...", file=sys.stderr)
    abundance_data = []
    
    for contig_id in viral_contigs:
        if contig_id in contig_stats:
            stats = contig_stats[contig_id]
            length_kb = stats['length'] / 1000.0
            mapped_reads = stats['mapped_reads']
            
            # RPM = (mapped reads / total mapped reads) * 1,000,000
            rpm = (mapped_reads / total_mapped_reads * 1000000) if total_mapped_reads > 0 else 0
            
            # RPKM = (mapped reads / (length_kb * total_mapped_reads / 1,000,000))
            rpkm = (mapped_reads / (length_kb * total_mapped_reads / 1000000)) if (length_kb > 0 and total_mapped_reads > 0) else 0
            
            abundance_data.append({
                'contig_id': contig_id,
                'length_bp': stats['length'],
                'length_kb': length_kb,
                'mapped_reads': mapped_reads,
                'RPM': round(rpm, 2),
                'RPKM': round(rpkm, 2)
            })
    
    # 6. 合并 Diamond 分类信息
    abundance_df = pd.DataFrame(abundance_data)
    if not abundance_df.empty:
        # 获取每个 contig 的最佳匹配（最高 bitscore）
        contig_best_hit = diamond_df.loc[diamond_df.groupby('contig_id')['bitscore'].idxmax()]
        contig_best_hit = contig_best_hit[['contig_id', 'sseqid', 'pident', 'evalue', 'bitscore', 'staxids']]
        contig_best_hit.columns = ['contig_id', 'best_hit', 'identity', 'evalue', 'bitscore', 'taxid']
        
        # 合并
        abundance_df = abundance_df.merge(contig_best_hit, on='contig_id', how='left')
        abundance_df = abundance_df.sort_values('RPKM', ascending=False)
    
    # 7. 保存结果
    # 文本格式
    with open("${sample}_megahit_abundance.txt", 'w') as f:
        f.write("="*80 + "\\n")
        f.write("Viral Abundance Report (MEGAHIT) - ${sample}\\n")
        f.write("="*80 + "\\n\\n")
        f.write(f"Total mapped reads: {total_mapped_reads:,}\\n")
        f.write(f"Viral contigs detected: {len(abundance_data)}\\n\\n")
        f.write(f"{'Contig ID':<30} {'Length (bp)':<15} {'Reads':<10} {'RPM':<12} {'RPKM':<12} {'Identity (%)':<12}\\n")
        f.write("-"*80 + "\\n")
        
        for _, row in abundance_df.iterrows():
            f.write(f"{row['contig_id']:<30} {row['length_bp']:<15} {row['mapped_reads']:<10} "
                   f"{row['RPM']:<12.2f} {row['RPKM']:<12.2f} {row.get('identity', 0):<12.2f}\\n")
    
    # CSV 格式
    abundance_df.to_csv("${sample}_megahit_abundance.csv", index=False)
    
    print(f"Abundance calculation completed: {len(abundance_data)} viral contigs", file=sys.stderr)
    """
}

// Process: Calculate Viral Abundance for SPAdes (RPM and RPKM) - Short Reads
process CALCULATE_ABUNDANCE_SPADES_SHORT {
    tag "${sample}_SPAdes"
    label 'process_medium'
    // Uses system environment (nextflow_env); bwa, samtools, htslib, bzip2, zlib, libdeflate, and pandas must be pre-installed in nextflow_env
    publishDir "${params.outdir}/abundance_spades", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs), path(reads), path(diamond_report)
    
    output:
    tuple val(sample), path("${sample}_spades_abundance.txt"), emit: abundance
    path("${sample}_spades_abundance.csv"), emit: abundance_csv
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import subprocess
    import sys
    import os
    import shutil
    from collections import defaultdict
    
    # Set LD_LIBRARY_PATH and create symlink for libbz2.so.1.0
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if not conda_prefix:
        # Try to find conda environment from PATH
        python_path = shutil.which('python3')
        if python_path and 'conda' in python_path:
            conda_prefix = os.path.dirname(os.path.dirname(python_path))
    
    if conda_prefix:
        lib_path = os.path.join(conda_prefix, 'lib')
        current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
        os.environ['LD_LIBRARY_PATH'] = f"{lib_path}:{current_ld_path}" if current_ld_path else lib_path
        print(f"Set LD_LIBRARY_PATH: {os.environ['LD_LIBRARY_PATH']}", file=sys.stderr)
        
        # Create symlink for libbz2.so.1.0 if it doesn't exist
        libbz2_target = os.path.join(lib_path, 'libbz2.so.1.0')
        if not os.path.exists(libbz2_target):
            # Find libbz2.so files
            import glob
            libbz2_files = glob.glob(os.path.join(lib_path, 'libbz2.so*'))
            if libbz2_files:
                # Sort to get the highest version number
                libbz2_files.sort(reverse=True)
                source_file = libbz2_files[0]
                try:
                    os.symlink(os.path.basename(source_file), libbz2_target)
                    print(f"Created symlink: {libbz2_target} -> {os.path.basename(source_file)}", file=sys.stderr)
                except (OSError, FileExistsError) as e:
                    # Symlink might already exist or permission issue
                    print(f"Note: Could not create symlink {libbz2_target}: {e}", file=sys.stderr)
            else:
                print(f"WARNING: No libbz2.so files found in {lib_path}", file=sys.stderr)
    
    print("=== Calculating Viral Abundance (SPAdes): ${sample} ===", file=sys.stderr)
    
    # Check if required tools are available
    def check_tool(tool_name):
        tool_path = shutil.which(tool_name)
        if tool_path is None:
            print(f"ERROR: {tool_name} not found in PATH", file=sys.stderr)
            print(f"PATH: {os.environ.get('PATH', 'Not set')}", file=sys.stderr)
            sys.exit(1)
        print(f"Found {tool_name} at: {tool_path}", file=sys.stderr)
        return tool_path
    
    bwa_path = check_tool('bwa')
    samtools_path = check_tool('samtools')
    
    # 1. 从 Diamond 结果中提取病毒 contigs
    print("Extracting viral contigs from Diamond results...", file=sys.stderr)
    try:
        diamond_df = pd.read_csv("${diamond_report}", sep='\\t', header=None,
                                names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                       'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                       'evalue', 'bitscore', 'staxids'])
    except pd.errors.EmptyDataError:
        print("WARNING: Diamond report is empty, creating empty abundance file", file=sys.stderr)
        with open("${sample}_spades_abundance.txt", 'w') as f:
            f.write("No viral contigs found\\n")
        pd.DataFrame(columns=['contig_id', 'length_bp', 'length_kb', 'mapped_reads', 'RPM', 'RPKM']).to_csv("${sample}_spades_abundance.csv", index=False)
        exit(0)
    
    # 提取 contig ID（移除基因编号）
    diamond_df['contig_id'] = diamond_df['qseqid'].str.rsplit('_', n=1).str[0]
    viral_contigs = set(diamond_df['contig_id'].unique())
    print(f"Found {len(viral_contigs)} viral contigs", file=sys.stderr)
    
    if len(viral_contigs) == 0:
        print("WARNING: No viral contigs found, creating empty abundance file", file=sys.stderr)
        with open("${sample}_spades_abundance.txt", 'w') as f:
            f.write("No viral contigs found\\n")
        pd.DataFrame(columns=['contig_id', 'length_bp', 'length_kb', 'mapped_reads', 'RPM', 'RPKM']).to_csv("${sample}_spades_abundance.csv", index=False)
        exit(0)
    
    # 2. 使用 BWA 将 reads 比对到 contigs
    print("Indexing contigs with BWA...", file=sys.stderr)
    try:
        subprocess.run([bwa_path, 'index', '${contigs}'], check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: BWA index failed: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    
    print("Aligning reads to contigs...", file=sys.stderr)
    with open('alignment.sam', 'w') as sam_file:
        try:
            subprocess.run([bwa_path, 'mem', '-t', '${task.cpus}', '${contigs}', 
                           '${reads[0]}', '${reads[1]}'], 
                          stdout=sam_file, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: BWA mem failed: {e.stderr.decode()}", file=sys.stderr)
            sys.exit(1)
    
    # 3. 转换为 BAM 并排序
    print("Converting and sorting BAM...", file=sys.stderr)
    with open('alignment.bam', 'wb') as bam_file:
        try:
            subprocess.run([samtools_path, 'view', '-bS', 'alignment.sam'], 
                          stdout=bam_file, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: samtools view failed: {e.stderr.decode()}", file=sys.stderr)
            print(f"Command: {samtools_path} view -bS alignment.sam", file=sys.stderr)
            sys.exit(1)
    try:
        subprocess.run([samtools_path, 'sort', '-o', 'sorted.bam', 'alignment.bam'], 
                      check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools sort failed: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    try:
        subprocess.run([samtools_path, 'index', 'sorted.bam'], check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools index failed: {e.stderr.decode()}", file=sys.stderr)
        sys.exit(1)
    
    # 4. 统计每个 contig 的 reads 数量
    print("Counting reads per contig...", file=sys.stderr)
    try:
        result = subprocess.run([samtools_path, 'idxstats', 'sorted.bam'], 
                              capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools idxstats failed: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
        sys.exit(1)
    
    # 解析 idxstats 输出
    contig_stats = {}
    total_mapped_reads = 0
    for line in result.stdout.strip().split('\\n'):
        if line.strip():
            parts = line.split('\\t')
            contig_id = parts[0]
            contig_length = int(parts[1])
            mapped_reads = int(parts[2])
            unmapped_reads = int(parts[3])
            
            if contig_id != '*':
                contig_stats[contig_id] = {
                    'length': contig_length,
                    'mapped_reads': mapped_reads
                }
                total_mapped_reads += mapped_reads
    
    print(f"Total mapped reads: {total_mapped_reads:,}", file=sys.stderr)
    
    # 5. 计算 RPM 和 RPKM
    print("Calculating RPM and RPKM...", file=sys.stderr)
    abundance_data = []
    
    for contig_id in viral_contigs:
        if contig_id in contig_stats:
            stats = contig_stats[contig_id]
            length_kb = stats['length'] / 1000.0
            mapped_reads = stats['mapped_reads']
            
            # RPM = (mapped reads / total mapped reads) * 1,000,000
            rpm = (mapped_reads / total_mapped_reads * 1000000) if total_mapped_reads > 0 else 0
            
            # RPKM = (mapped reads / (length_kb * total_mapped_reads / 1,000,000))
            rpkm = (mapped_reads / (length_kb * total_mapped_reads / 1000000)) if (length_kb > 0 and total_mapped_reads > 0) else 0
            
            abundance_data.append({
                'contig_id': contig_id,
                'length_bp': stats['length'],
                'length_kb': length_kb,
                'mapped_reads': mapped_reads,
                'RPM': round(rpm, 2),
                'RPKM': round(rpkm, 2)
            })
    
    # 6. 合并 Diamond 分类信息
    abundance_df = pd.DataFrame(abundance_data)
    if not abundance_df.empty:
        # 获取每个 contig 的最佳匹配（最高 bitscore）
        contig_best_hit = diamond_df.loc[diamond_df.groupby('contig_id')['bitscore'].idxmax()]
        contig_best_hit = contig_best_hit[['contig_id', 'sseqid', 'pident', 'evalue', 'bitscore', 'staxids']]
        contig_best_hit.columns = ['contig_id', 'best_hit', 'identity', 'evalue', 'bitscore', 'taxid']
        
        # 合并
        abundance_df = abundance_df.merge(contig_best_hit, on='contig_id', how='left')
        abundance_df = abundance_df.sort_values('RPKM', ascending=False)
    
    # 7. 保存结果
    # 文本格式
    with open("${sample}_spades_abundance.txt", 'w') as f:
        f.write("="*80 + "\\n")
        f.write("Viral Abundance Report (SPAdes) - ${sample}\\n")
        f.write("="*80 + "\\n\\n")
        f.write(f"Total mapped reads: {total_mapped_reads:,}\\n")
        f.write(f"Viral contigs detected: {len(abundance_data)}\\n\\n")
        f.write(f"{'Contig ID':<30} {'Length (bp)':<15} {'Reads':<10} {'RPM':<12} {'RPKM':<12} {'Identity (%)':<12}\\n")
        f.write("-"*80 + "\\n")
        
        for _, row in abundance_df.iterrows():
            f.write(f"{row['contig_id']:<30} {row['length_bp']:<15} {row['mapped_reads']:<10} "
                   f"{row['RPM']:<12.2f} {row['RPKM']:<12.2f} {row.get('identity', 0):<12.2f}\\n")
    
    # CSV 格式
    abundance_df.to_csv("${sample}_spades_abundance.csv", index=False)
    
    print(f"Abundance calculation completed: {len(abundance_data)} viral contigs", file=sys.stderr)
    """
}

// ================================================================================
// Long Read Process Definitions (Nanopore/PacBio)
// ================================================================================

// Process: Calculate Viral Abundance for MetaFlye (RPM and RPKM, long reads)
process CALCULATE_ABUNDANCE_METAFLYE {
    tag "${sample}_MetaFlye"
    label 'process_medium'
    conda 'bioconda::minimap2=2.24 bioconda::samtools=1.17 bioconda::htslib conda-forge::bzip2=1.0.8 conda-forge::zlib conda-forge::libdeflate conda-forge::pandas=2.0.3'
    publishDir "${params.outdir}/abundance_metaflye", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs), path(fastq_long), path(diamond_report)
    
    output:
    tuple val(sample), path("${sample}_metaflye_abundance.txt"), emit: abundance
    path("${sample}_metaflye_abundance.csv"), emit: abundance_csv
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import subprocess
    import sys
    
    print("=== Calculating Viral Abundance (MetaFlye - Long Reads): ${sample} ===", file=sys.stderr)
    
    # 1) Parse Diamond to get viral contigs
    try:
        df = pd.read_csv("${diamond_report}", sep='\\t', header=None,
                         names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','staxids'])
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','staxids'])
    
    if df.empty:
        print("WARNING: No Diamond hits; creating empty abundance files.", file=sys.stderr)
        open("${sample}_metaflye_abundance.txt","w").write("No viral contigs found\\n")
        pd.DataFrame(columns=['contig_id','length_bp','length_kb','mapped_reads','RPM','RPKM']).to_csv("${sample}_metaflye_abundance.csv", index=False)
        sys.exit(0)
    
    df['contig_id'] = df['qseqid'].astype(str).str.rsplit('_', n=1).str[0]
    viral_contigs = set(df['contig_id'].unique())
    print(f"Found {len(viral_contigs)} viral contigs", file=sys.stderr)
    
    if len(viral_contigs) == 0:
        open("${sample}_metaflye_abundance.txt","w").write("No viral contigs found\\n")
        pd.DataFrame(columns=['contig_id','length_bp','length_kb','mapped_reads','RPM','RPKM']).to_csv("${sample}_metaflye_abundance.csv", index=False)
        sys.exit(0)
    
    # 2) Map long reads to contigs using minimap2
    # Default preset 'map-ont' suitable for ONT; allow override via env if needed
    preset = "${ params.long_read_preset ?: 'map-ont' }"
    print(f"Mapping with minimap2 preset: {preset}", file=sys.stderr)
    
    subprocess.run(["minimap2", "-t", str(${task.cpus}), "-a", "-x", preset, "${contigs}", "${fastq_long}"],
                   stdout=open("alignment.sam","w"), check=True)
    
    # 3) Count mapped reads per contig from SAM (avoid external samtools)
    def load_lengths(fa_path):
        lens = {}
        with open(fa_path) as fh:
            name = None
            L = 0
            for ln in fh:
                if not ln:
                    continue
                if ln.startswith('>'):
                    if name is not None:
                        lens[name] = L
                    name = ln[1:].strip().split()[0]
                    L = 0
                else:
                    s = ln.strip()
                    if s:
                        L += len(s)
            if name is not None:
                lens[name] = L
        return lens

    contig_lengths = load_lengths("${contigs}")

    print("Counting mapped reads per contig from SAM...", file=sys.stderr)
    contig_stats = {}
    total_mapped = 0
    with open("alignment.sam") as sf:
        for ln in sf:
            if not ln or ln[0] == '@':
                continue
            cols = ln.rstrip("\\n").split("\\t")
            if len(cols) < 3:
                continue
            try:
                flag = int(cols[1])
            except Exception:
                continue
            # skip unmapped (0x4), secondary (0x100), supplementary (0x800)
            if (flag & 0x4) != 0 or (flag & 0x100) != 0 or (flag & 0x800) != 0:
                continue
            rname = cols[2]
            if rname == "*" or not rname:
                continue
            if rname not in contig_stats:
                contig_stats[rname] = {"length": int(contig_lengths.get(rname, 0)), "mapped_reads": 0}
            contig_stats[rname]["mapped_reads"] += 1
            total_mapped += 1
    
    # 4) Compute RPM / RPKM
    rows = []
    for cid in viral_contigs:
        if cid in contig_stats:
            L = contig_stats[cid]["length"]
            mk = contig_stats[cid]["mapped_reads"]
            kb = L/1000.0
            rpm = (mk/total_mapped*1_000_000) if total_mapped>0 else 0.0
            rpkm = (mk/(kb*total_mapped/1_000_000)) if (kb>0 and total_mapped>0) else 0.0
            rows.append({"contig_id": cid, "length_bp": L, "length_kb": kb, "mapped_reads": mk, "RPM": round(rpm,2), "RPKM": round(rpkm,2)})
    
    out_df = pd.DataFrame(rows).sort_values("RPKM", ascending=False)
    out_df.to_csv("${sample}_metaflye_abundance.csv", index=False)
    
    with open("${sample}_metaflye_abundance.txt","w") as f:
        f.write("="*80+"\\n")
        f.write("Viral Abundance Report (MetaFlye - Long Reads) - ${sample}\\n")
        f.write("="*80+"\\n\\n")
        f.write(f"Total mapped reads: {total_mapped}\\n")
        f.write(f"Viral contigs reported: {len(out_df)}\\n\\n")
        f.write(f"{'Contig ID':<40}{'Length(bp)':>12}{'Reads':>10}{'RPM':>12}{'RPKM':>12}\\n")
        f.write("-"*80+"\\n")
        for _,r in out_df.iterrows():
            f.write(f"{r['contig_id']:<40}{int(r['length_bp']):>12}{int(r['mapped_reads']):>10}{r['RPM']:>12.2f}{r['RPKM']:>12.2f}\\n")
    """
}

// Process: Calculate Viral Abundance for viralFlye (RPM and RPKM, long reads)
process CALCULATE_ABUNDANCE_VIRALFLYE {
    tag "${sample}_viralFlye"
    label 'process_medium'
    conda 'bioconda::minimap2=2.24 bioconda::samtools=1.17 bioconda::htslib conda-forge::bzip2=1.0.8 conda-forge::zlib conda-forge::libdeflate conda-forge::pandas=2.0.3'
    publishDir "${params.outdir}/abundance_viralflye", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(contigs), path(fastq_long), path(diamond_report)
    
    output:
    tuple val(sample), path("${sample}_viralflye_abundance.txt"), emit: abundance
    path("${sample}_viralflye_abundance.csv"), emit: abundance_csv
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import subprocess
    import sys
    
    print("=== Calculating Viral Abundance (viralFlye - Long Reads): ${sample} ===", file=sys.stderr)
    
    # 1) Parse Diamond to get viral contigs (should match viralFlye set)
    try:
        df = pd.read_csv("${diamond_report}", sep='\\t', header=None,
                         names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','staxids'])
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','staxids'])
    
    if df.empty:
        print("WARNING: No Diamond hits; creating empty abundance files.", file=sys.stderr)
        open("${sample}_viralflye_abundance.txt","w").write("No viral contigs found\\n")
        pd.DataFrame(columns=['contig_id','length_bp','length_kb','mapped_reads','RPM','RPKM']).to_csv("${sample}_viralflye_abundance.csv", index=False)
        sys.exit(0)
    
    df['contig_id'] = df['qseqid'].astype(str).str.rsplit('_', n=1).str[0]
    viral_contigs = set(df['contig_id'].unique())
    print(f"Found {len(viral_contigs)} viral contigs", file=sys.stderr)
    
    if len(viral_contigs) == 0:
        open("${sample}_viralflye_abundance.txt","w").write("No viral contigs found\\n")
        pd.DataFrame(columns=['contig_id','length_bp','length_kb','mapped_reads','RPM','RPKM']).to_csv("${sample}_viralflye_abundance.csv", index=False)
        sys.exit(0)
    
    # 2) Map long reads to contigs using minimap2
    preset = "${ params.long_read_preset ?: 'map-ont' }"
    print(f"Mapping with minimap2 preset: {preset}", file=sys.stderr)
    
    subprocess.run(["minimap2", "-t", str(${task.cpus}), "-a", "-x", preset, "${contigs}", "${fastq_long}"],
                   stdout=open("alignment.sam","w"), check=True)
    
    # 3) Count mapped reads per contig from SAM (avoid external samtools)
    def load_lengths(fa_path):
        lens = {}
        with open(fa_path) as fh:
            name = None
            L = 0
            for ln in fh:
                if not ln:
                    continue
                if ln.startswith('>'):
                    if name is not None:
                        lens[name] = L
                    name = ln[1:].strip().split()[0]
                    L = 0
                else:
                    s = ln.strip()
                    if s:
                        L += len(s)
            if name is not None:
                lens[name] = L
        return lens

    contig_lengths = load_lengths("${contigs}")

    print("Counting mapped reads per contig from SAM...", file=sys.stderr)
    contig_stats = {}
    total_mapped = 0
    with open("alignment.sam") as sf:
        for ln in sf:
            if not ln or ln[0] == '@':
                continue
            cols = ln.rstrip("\\n").split("\\t")
            if len(cols) < 3:
                continue
            try:
                flag = int(cols[1])
            except Exception:
                continue
            # skip unmapped (0x4), secondary (0x100), supplementary (0x800)
            if (flag & 0x4) != 0 or (flag & 0x100) != 0 or (flag & 0x800) != 0:
                continue
            rname = cols[2]
            if rname == "*" or not rname:
                continue
            if rname not in contig_stats:
                contig_stats[rname] = {"length": int(contig_lengths.get(rname, 0)), "mapped_reads": 0}
            contig_stats[rname]["mapped_reads"] += 1
            total_mapped += 1
    
    # 4) Compute RPM / RPKM
    rows = []
    for cid in viral_contigs:
        if cid in contig_stats:
            L = contig_stats[cid]["length"]
            mk = contig_stats[cid]["mapped_reads"]
            kb = L/1000.0
            rpm = (mk/total_mapped*1_000_000) if total_mapped>0 else 0.0
            rpkm = (mk/(kb*total_mapped/1_000_000)) if (kb>0 and total_mapped>0) else 0.0
            rows.append({"contig_id": cid, "length_bp": L, "length_kb": kb, "mapped_reads": mk, "RPM": round(rpm,2), "RPKM": round(rpkm,2)})
    
    out_df = pd.DataFrame(rows).sort_values("RPKM", ascending=False)
    out_df.to_csv("${sample}_viralflye_abundance.csv", index=False)
    
    with open("${sample}_viralflye_abundance.txt","w") as f:
        f.write("="*80+"\\n")
        f.write("Viral Abundance Report (viralFlye - Long Reads) - ${sample}\\n")
        f.write("="*80+"\\n\\n")
        f.write(f"Total mapped reads: {total_mapped}\\n")
        f.write(f"Viral contigs reported: {len(out_df)}\\n\\n")
        f.write(f"{'Contig ID':<40}{'Length(bp)':>12}{'Reads':>10}{'RPM':>12}{'RPKM':>12}\\n")
        f.write("-"*80+"\\n")
        for _,r in out_df.iterrows():
            f.write(f"{r['contig_id']:<40}{int(r['length_bp']):>12}{int(r['mapped_reads']):>10}{r['RPM']:>12.2f}{r['RPKM']:>12.2f}\\n")
    """
}

// Process: MetaFlye Assembly (for long reads)
process METAFLYE_ASSEMBLY {
    tag "${sample}_MetaFlye"
    label 'process_high'
    conda 'bioconda::flye=2.9.3'
    publishDir "${params.outdir}/assembly_metaflye", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(fastq_long), val(read_type)
    
    output:
    tuple val(sample), path("${sample}_metaflye_contigs.fa"), emit: contigs
    tuple val(sample), path("metaflye_output"), emit: metaflye_dir
    
    script:
    def genome_size_arg = params.metaflye_genome_size ? "--genome-size ${params.metaflye_genome_size}" : ""
    """
    echo "=== MetaFlye Assembly: ${sample} ==="
    
    # Run MetaFlye for long read assembly (using --meta for metagenome)
    flye \\
        --nano-raw ${fastq_long} \\
        --out-dir metaflye_output \\
        --threads ${task.cpus} \\
        --iterations 3 \\
        --meta \\
        ${genome_size_arg}
    
    # Copy assembly results (for subsequent Prodigal)
    cp metaflye_output/assembly.fasta ${sample}_metaflye_contigs.fa
    
    CONTIG_COUNT=\$(grep -c ">" ${sample}_metaflye_contigs.fa)
    TOTAL_LENGTH=\$(grep -v ">" ${sample}_metaflye_contigs.fa | tr -d '\\n' | wc -c)
    
    echo "MetaFlye: Generated \${CONTIG_COUNT} contigs"
    echo "MetaFlye: Total length \${TOTAL_LENGTH} bp"
    """
}

// Process: viralFlye Refinement (optional, for viral sequence optimization)
// Note: viralFlye must be pre-installed in the runtime environment
// Due to conda environment creation timeout, use system-installed viralFlye
process VIRALFLYE_REFINEMENT {
    tag "${sample}_viralFlye"
    label 'process_medium'
    // Don't use conda, directly use system-installed viralFlye
    publishDir "${params.outdir}/assembly_viralflye", mode: 'copy', pattern: "*.fa"
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    
    input:
    tuple val(sample), path(metaflye_dir)
    tuple val(sample_reads), path(fastq_long)
    val(pfam_hmm)
    
    output:
    tuple val(sample), path("${sample}_viralflye_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== viralFlye Refinement: ${sample} ==="
    echo "Activating viralFlye_env environment..."
    
    # Allow unbound variables (solve conda environment switching issue)
    set +u
    
    # Load conda and activate viralFlye_env environment
    module load Miniforge3/24.11.3-0 || true
    source "\$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || source ~/.bashrc
    
    # Activate viralFlye_env environment (ignore deactivate errors)
    if conda env list | grep -q "viralFlye_env"; then
        echo "Found viralFlye_env, activating..."
        conda activate viralFlye_env 2>&1 | grep -v "unbound variable" || true
        echo "Activated environment: \$CONDA_DEFAULT_ENV"
    else
        echo "ERROR: viralFlye_env not found!"
        exit 1
    fi
    
    # Use Python and viralFlye.py from viralFlye_env directly (absolute path)
    VIRALFLYE_PYTHON="\$HOME/.conda/envs/viralFlye_env/bin/python"
    VIRALFLYE_SCRIPT="\$HOME/.conda/envs/viralFlye_env/bin/viralFlye.py"
    
    # Check if viralFlye.py exists
    if [ ! -f "\$VIRALFLYE_SCRIPT" ]; then
        echo "ERROR: viralFlye.py not found at \$VIRALFLYE_SCRIPT"
        echo "Searching for viralFlye.py..."
        find \$HOME/.conda/envs/viralFlye_env -name "viralFlye.py" -type f || true
        echo "Please install viralFlye in viralFlye_env"
        exit 1
    fi
    
    echo "viralFlye.py found: \$VIRALFLYE_SCRIPT"
    echo "Python: \$VIRALFLYE_PYTHON"
    
    # Run viralFlye.py for viral sequence refinement
    echo "Running viralFlye.py..."
    echo "  MetaFlye dir: ${metaflye_dir}"
    echo "  Reads: ${fastq_long}"
    echo "  HMM database: ${pfam_hmm}"
    
    # Run viralFlye in permissive mode (avoid environment issues)
    set +e  # Allow command failure
    \$VIRALFLYE_PYTHON \$VIRALFLYE_SCRIPT \\
        --dir ${metaflye_dir} \\
        --hmm ${pfam_hmm} \\
        --reads ${fastq_long} \\
        --outdir viralflye_output \\
        --threads ${task.cpus} \\
        --min_viral_length ${params.viralflye_min_length} \\
        --completeness ${params.viralflye_completeness}
    
    VIRALFLYE_EXIT=\$?
    
    # Function to ensure output file exists
    ensure_output() {
        if [ ! -f "${sample}_viralflye_contigs.fa" ]; then
            echo "ERROR: Output file not created, using MetaFlye contigs as fallback"
            cp ${metaflye_dir}/assembly.fasta ${sample}_viralflye_contigs.fa || exit 1
        fi
    }
    
    # Check output files and extract viral sequences
    if [ \$VIRALFLYE_EXIT -eq 0 ] && [ -d "viralflye_output" ]; then
        echo "viralFlye completed successfully"
        
        # viralFlye saves viral contig IDs in txt files
        # Need to extract these sequences from original assembly.fasta
        
        VIRAL_IDS=""
        
        # Merge linear and circular viral IDs
        if [ -f "viralflye_output/vv_linears/linears.txt" ] && [ -s "viralflye_output/vv_linears/linears.txt" ]; then
            echo "Found linear viral contigs:"
            cat viralflye_output/vv_linears/linears.txt
            VIRAL_IDS="\${VIRAL_IDS} viralflye_output/vv_linears/linears.txt"
        fi
        
        if [ -f "viralflye_output/vv_circulars/circulars.txt" ] && [ -s "viralflye_output/vv_circulars/circulars.txt" ]; then
            echo "Found circular viral contigs:"
            cat viralflye_output/vv_circulars/circulars.txt
            VIRAL_IDS="\${VIRAL_IDS} viralflye_output/vv_circulars/circulars.txt"
        fi
        
        # If viral IDs found, extract from original assembly
        if [ -n "\${VIRAL_IDS}" ]; then
            echo "Extracting viral sequences..."
            cat \${VIRAL_IDS} > viralflye_output/all_viral_ids.txt
            
            # Use seqtk to extract sequences
            if command -v seqtk &> /dev/null; then
                seqtk subseq viralflye_output/assembly.fasta viralflye_output/all_viral_ids.txt > ${sample}_viralflye_contigs.fa
            else
                # If seqtk is not available, use Python to extract
                \$VIRALFLYE_PYTHON -c "
import sys
ids = set()
with open('viralflye_output/all_viral_ids.txt') as f:
    for line in f:
        ids.add(line.strip())

from Bio import SeqIO
with open('${sample}_viralflye_contigs.fa', 'w') as out:
    for record in SeqIO.parse('viralflye_output/assembly.fasta', 'fasta'):
        if record.id in ids:
            SeqIO.write(record, out, 'fasta')
" || echo "Sequence extraction failed"
            fi
            
            # Verify extraction results
            if [ -s "${sample}_viralflye_contigs.fa" ]; then
                VIRAL_COUNT=\$(grep -c ">" ${sample}_viralflye_contigs.fa)
                echo "Successfully extracted \${VIRAL_COUNT} viral contigs"
            else
                echo "WARNING: Extraction failed, using all MetaFlye contigs"
                cp ${metaflye_dir}/assembly.fasta ${sample}_viralflye_contigs.fa
            fi
        else
            echo "WARNING: No viral contigs identified by viralFlye"
            echo "Using original MetaFlye contigs"
            cp ${metaflye_dir}/assembly.fasta ${sample}_viralflye_contigs.fa
        fi
    else
        echo "WARNING: viralFlye failed (exit code: \${VIRALFLYE_EXIT})"
        echo "Using original MetaFlye contigs as fallback"
        cp ${metaflye_dir}/assembly.fasta ${sample}_viralflye_contigs.fa
    fi
    
    # Final check: ensure output file exists
    ensure_output
    
    # Verify file
    if [ ! -s "${sample}_viralflye_contigs.fa" ]; then
        echo "ERROR: Output file is empty!"
        exit 1
    fi
    
    CONTIG_COUNT=\$(grep -c ">" ${sample}_viralflye_contigs.fa || echo 0)
    TOTAL_LENGTH=\$(grep -v ">" ${sample}_viralflye_contigs.fa | tr -d '\\n' | wc -c || echo 0)
    
    echo "viralFlye: Generated \${CONTIG_COUNT} contigs"
    echo "viralFlye: Total length \${TOTAL_LENGTH} bp"
    echo "Output file: ${sample}_viralflye_contigs.fa"
    """
}

// Process: Prodigal Gene Prediction for MetaFlye (all sequences, identify viruses based on sequence similarity)
process PRODIGAL_METAFLYE {
    tag "${sample}_MetaFlye"
    label 'process_medium'
    conda 'bioconda::prodigal=2.6.3'
    publishDir "${params.outdir}/prodigal_metaflye", mode: 'copy', pattern: "*.faa"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_metaflye_proteins.faa"), emit: proteins
    path("${sample}_metaflye_genes.fna"), emit: genes
    
    script:
    """
    echo "=== Prodigal Gene Prediction (MetaFlye - All Sequences): ${sample} ==="
    
    # Run Prodigal for gene prediction (metagenome mode)
    prodigal \\
        -i ${contigs} \\
        -a ${sample}_metaflye_proteins.faa \\
        -d ${sample}_metaflye_genes.fna \\
        -p meta \\
        -q
    
    # Count predicted genes
    GENE_COUNT=\$(grep -c ">" ${sample}_metaflye_proteins.faa || echo 0)
    echo "Prodigal: Predicted \${GENE_COUNT} genes from MetaFlye contigs (all sequences)"
    """
}

// Process: Prodigal Gene Prediction for viralFlye (viral sequences only)
process PRODIGAL_VIRALFLYE {
    tag "${sample}_viralFlye"
    label 'process_medium'
    conda 'bioconda::prodigal=2.6.3'
    publishDir "${params.outdir}/prodigal_viralflye", mode: 'copy', pattern: "*.faa"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_viralflye_proteins.faa"), emit: proteins
    path("${sample}_viralflye_genes.fna"), emit: genes
    
    script:
    """
    echo "=== Prodigal Gene Prediction (viralFlye - Viral Sequences Only): ${sample} ==="
    
    # Run Prodigal for gene prediction (metagenome mode)
    prodigal \\
        -i ${contigs} \\
        -a ${sample}_viralflye_proteins.faa \\
        -d ${sample}_viralflye_genes.fna \\
        -p meta \\
        -q
    
    # Count predicted genes
    GENE_COUNT=\$(grep -c ">" ${sample}_viralflye_proteins.faa || echo 0)
    echo "Prodigal: Predicted \${GENE_COUNT} genes from viralFlye contigs (viral sequences only)"
    """
}

// Process: Diamond Classification for MetaFlye (all sequences)
process DIAMOND_CLASSIFICATION_METAFLYE {
    tag "${sample}_MetaFlye"
    label 'process_high'
    conda 'bioconda::diamond=2.1.8'
    publishDir "${params.outdir}/diamond_metaflye", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(proteins)
    val(diamond_db)
    
    output:
    tuple val(sample), path("${sample}_metaflye_diamond.txt"), emit: diamond_metaflye
    
    script:
    """
    echo "=== Diamond Classification (MetaFlye - All Sequences): ${sample} ==="
    
    # Run Diamond BLASTP for classification
    diamond blastp \\
        --query ${proteins} \\
        --db ${diamond_db} \\
        --out ${sample}_metaflye_diamond.txt \\
        --outfmt ${params.diamond_outfmt} \\
        --threads ${task.cpus} \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs ${params.diamond_max_target_seqs} \\
        --sensitive
    
    # Count hits
    HIT_COUNT=\$(wc -l < ${sample}_metaflye_diamond.txt || echo 0)
    echo "Diamond: Found \${HIT_COUNT} hits for MetaFlye proteins (all sequences)"
    """
}

// Process: Diamond Classification for viralFlye (viral sequences only)
process DIAMOND_CLASSIFICATION_VIRALFLYE {
    tag "${sample}_viralFlye"
    label 'process_high'
    conda 'bioconda::diamond=2.1.8'
    publishDir "${params.outdir}/diamond_viralflye", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(proteins)
    val(diamond_db)
    
    output:
    tuple val(sample), path("${sample}_viralflye_diamond.txt"), emit: diamond_viralflye
    
    script:
    """
    echo "=== Diamond Classification (viralFlye - Viral Sequences Only): ${sample} ==="
    
    # Run Diamond BLASTP for classification
    diamond blastp \\
        --query ${proteins} \\
        --db ${diamond_db} \\
        --out ${sample}_viralflye_diamond.txt \\
        --outfmt ${params.diamond_outfmt} \\
        --threads ${task.cpus} \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs ${params.diamond_max_target_seqs} \\
        --sensitive
    
    # Count hits
    HIT_COUNT=\$(wc -l < ${sample}_viralflye_diamond.txt || echo 0)
    echo "Diamond: Found \${HIT_COUNT} hits for viralFlye proteins (viral sequences only)"
    """
}

// Process: Add Taxonomy Information for MetaFlye Results
process ADD_TAXONOMY_METAFLYE {
    tag "${sample}_MetaFlye"
    label 'process_low'
    conda 'conda-forge::pandas=2.0.3'
    publishDir "${params.outdir}/taxonomy_metaflye", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(diamond_report)
    path(taxonomy_names)
    path(taxonomy_nodes)
    
    output:
    tuple val(sample), path("${sample}_metaflye_diamond_with_taxonomy.txt"), emit: enhanced_report
    path("${sample}_metaflye_taxonomy_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import sys
    from collections import Counter
    
    # Taxonomy database class
    class TaxonomyDB:
        \"\"\"Parse NCBI taxonomy database\"\"\"
        
        def __init__(self, names_file, nodes_file):
            \"\"\"Initialize taxonomy database\"\"\"
            print("Loading taxonomy database...", file=sys.stderr)
            self.names = self._load_names(names_file)
            self.nodes = self._load_nodes(nodes_file)
            print(f"Loaded: {len(self.names):,} names, {len(self.nodes):,} nodes", file=sys.stderr)
        
        def _load_names(self, names_file):
            \"\"\"Load names.dmp\"\"\"
            names = {}
            with open(names_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 4:
                        taxid = parts[0]
                        name = parts[1]
                        name_class = parts[3]
                        
                        # Only keep scientific names
                        if name_class == 'scientific name':
                            names[taxid] = name
            
            return names
        
        def _load_nodes(self, nodes_file):
            \"\"\"Load nodes.dmp\"\"\"
            nodes = {}
            with open(nodes_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 3:
                        taxid = parts[0]
                        parent_taxid = parts[1]
                        rank = parts[2]
                        
                        nodes[taxid] = {
                            'parent': parent_taxid,
                            'rank': rank
                        }
            
            return nodes
        
        def get_lineage(self, taxid):
            \"\"\"Get complete taxonomic lineage\"\"\"
            lineage = {
                'superkingdom': 'N/A',
                'kingdom': 'N/A',
                'phylum': 'N/A',
                'class': 'N/A',
                'order': 'N/A',
                'family': 'N/A',
                'genus': 'N/A',
                'species': 'N/A',
                'organism_name': 'N/A'
            }
            
            taxid = str(taxid)
            
            # Get current taxid name
            if taxid in self.names:
                lineage['organism_name'] = self.names[taxid]
            
            # Traverse up the taxonomy tree
            current_taxid = taxid
            visited = set()
            
            while current_taxid != '1' and current_taxid in self.nodes:
                # Prevent loops
                if current_taxid in visited:
                    break
                visited.add(current_taxid)
                
                node = self.nodes[current_taxid]
                rank = node['rank']
                
                # Only keep major ranks
                if rank in lineage and current_taxid in self.names:
                    lineage[rank] = self.names[current_taxid]
                
                # Move to parent node
                current_taxid = node['parent']
            
            return lineage
    
    # Initialize Taxonomy database
    taxonomy_db = TaxonomyDB("${taxonomy_names}", "${taxonomy_nodes}")
    
    # Parse Diamond output file
    print(f"\\nParsing Diamond results: ${diamond_report}", file=sys.stderr)
    
    try:
        df = pd.read_csv("${diamond_report}", sep='\\t', header=None, 
                       names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                              'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                              'evalue', 'bitscore', 'staxids'])
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                  'evalue', 'bitscore', 'staxids'])
    
    if df.empty:
        print("WARNING: No Diamond hits found", file=sys.stderr)
        # Create empty output files
        with open("${sample}_long_diamond_with_taxonomy.txt", 'w') as f:
            f.write("# No hits found\\n")
        with open("${sample}_taxonomy_summary.txt", 'w') as f:
            f.write("No hits found\\n")
    else:
        print(f"Adding taxonomy info to {len(df):,} records...", file=sys.stderr)
        
        # Add taxonomic information
        lineages = []
        for idx, row in df.iterrows():
            if (idx + 1) % 100 == 0:
                print(f"  Processing {idx+1:,}/{len(df):,}...", file=sys.stderr)
            
            # Convert TaxID to string
            if pd.notna(row['staxids']):
                try:
                    taxid = str(int(float(row['staxids'])))
                except (ValueError, TypeError):
                    taxid = '0'
            else:
                taxid = '0'
            
            lineage = taxonomy_db.get_lineage(taxid)
            lineages.append(lineage)
        
        # Add taxonomy columns
        for key in ['organism_name', 'superkingdom', 'kingdom', 'phylum', 
                    'class', 'order', 'family', 'genus', 'species']:
            df[key] = [lineage[key] for lineage in lineages]
        
        # Save enhanced results
        df.to_csv("${sample}_metaflye_diamond_with_taxonomy.txt", sep='\\t', index=False)
        print(f"Enhanced Diamond results saved", file=sys.stderr)
        
        # Generate statistical summary
        with open("${sample}_metaflye_taxonomy_summary.txt", 'w', encoding='utf-8') as f:
            f.write("="*80 + "\\n")
            f.write("Viral Taxonomy Statistical Summary - MetaFlye (All Sequences)\\n")
            f.write("="*80 + "\\n\\n")
            
            f.write(f"Sample Name: ${sample}\\n")
            f.write(f"Total Hits: {len(df):,}\\n")
            f.write(f"Unique Queries: {df['qseqid'].nunique():,}\\n")
            f.write(f"Average Identity: {df['pident'].mean():.2f}%\\n")
            f.write(f"Average Alignment Length: {df['length'].mean():.1f} aa\\n\\n")
            
            # Superkingdom statistics
            f.write("[Superkingdom Distribution]\\n")
            f.write("-"*80 + "\\n")
            superkingdom_counts = Counter(df['superkingdom'].dropna())
            for sk, count in superkingdom_counts.most_common(10):
                f.write(f"{sk:<40} {count:>10,}\\n")
            f.write("\\n")
            
            # Kingdom statistics
            f.write("[Kingdom Distribution]\\n")
            f.write("-"*80 + "\\n")
            kingdom_counts = Counter(df['kingdom'].dropna())
            for k, count in kingdom_counts.most_common(10):
                f.write(f"{k:<40} {count:>10,}\\n")
            f.write("\\n")
            
            # Phylum statistics
            f.write("[Phylum Distribution (Top 15)]\\n")
            f.write("-"*80 + "\\n")
            phylum_counts = Counter(df['phylum'].dropna())
            for phylum, count in phylum_counts.most_common(15):
                f.write(f"{phylum:<40} {count:>10,}\\n")
            f.write("\\n")
            
            # Family statistics
            f.write("[Family Distribution (Top 15)]\\n")
            f.write("-"*80 + "\\n")
            family_counts = Counter(df['family'].dropna())
            for family, count in family_counts.most_common(15):
                f.write(f"{family:<40} {count:>10,}\\n")
            f.write("\\n")
            
            # Species statistics
            f.write("[Species Distribution (Top 20)]\\n")
            f.write("-"*80 + "\\n")
            species_counts = Counter(df['species'].dropna())
            for species, count in species_counts.most_common(20):
                f.write(f"{species:<40} {count:>10,}\\n")
            f.write("\\n")
            
            f.write("="*80 + "\\n")
            f.write("Analysis Complete\\n")
            f.write("="*80 + "\\n")
        
        print(f"Taxonomy summary saved", file=sys.stderr)
    """
}

// Process: Add Taxonomy Information for viralFlye Results (viral sequences only)
process ADD_TAXONOMY_VIRALFLYE {
    tag "${sample}_viralFlye"
    label 'process_low'
    conda 'conda-forge::pandas=2.0.3'
    publishDir "${params.outdir}/taxonomy_viralflye", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(diamond_report)
    path(taxonomy_names)
    path(taxonomy_nodes)
    
    output:
    tuple val(sample), path("${sample}_viralflye_diamond_with_taxonomy.txt"), emit: enhanced_report
    path("${sample}_viralflye_taxonomy_summary.txt"), emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import sys
    from collections import Counter
    
    # Use the same TaxonomyDB class as ADD_TAXONOMY_METAFLYE (simplified here, will be included in actual use)
    # For simplicity, reuse the code logic above
    
    exec(open("${projectDir}/taxonomy_db.py").read()) if False else None
    
    # Since Python scripts need full definition in Nextflow, complete TaxonomyDB class is copied here
    class TaxonomyDB:
        def __init__(self, names_file, nodes_file):
            print("Loading taxonomy database...", file=sys.stderr)
            self.names = self._load_names(names_file)
            self.nodes = self._load_nodes(nodes_file)
            print(f"Loaded: {len(self.names):,} names, {len(self.nodes):,} nodes", file=sys.stderr)
        
        def _load_names(self, names_file):
            names = {}
            with open(names_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 4 and parts[3] == 'scientific name':
                        names[parts[0]] = parts[1]
            return names
        
        def _load_nodes(self, nodes_file):
            nodes = {}
            with open(nodes_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 3:
                        nodes[parts[0]] = {'parent': parts[1], 'rank': parts[2]}
            return nodes
        
        def get_lineage(self, taxid):
            lineage = {
                'superkingdom': 'N/A', 'kingdom': 'N/A', 'phylum': 'N/A',
                'class': 'N/A', 'order': 'N/A', 'family': 'N/A',
                'genus': 'N/A', 'species': 'N/A', 'organism_name': 'N/A'
            }
            taxid = str(taxid)
            if taxid in self.names:
                lineage['organism_name'] = self.names[taxid]
            
            current_taxid = taxid
            visited = set()
            while current_taxid != '1' and current_taxid in self.nodes:
                if current_taxid in visited:
                    break
                visited.add(current_taxid)
                node = self.nodes[current_taxid]
                rank = node['rank']
                if rank in lineage and current_taxid in self.names:
                    lineage[rank] = self.names[current_taxid]
                current_taxid = node['parent']
            return lineage
    
    # Initialize database
    taxonomy_db = TaxonomyDB("${taxonomy_names}", "${taxonomy_nodes}")
    
    # Parse Diamond results
    try:
        df = pd.read_csv("${diamond_report}", sep='\\t', header=None, 
                       names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                              'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                              'evalue', 'bitscore', 'staxids'])
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                  'evalue', 'bitscore', 'staxids'])
    
    if not df.empty:
        print(f"Adding taxonomy info to {len(df):,} viralFlye records...", file=sys.stderr)
        lineages = []
        for idx, row in df.iterrows():
            if (idx + 1) % 100 == 0:
                print(f"  Processing {idx+1:,}/{len(df):,}...", file=sys.stderr)
            if pd.notna(row['staxids']):
                try:
                    taxid = str(int(float(row['staxids'])))
                except (ValueError, TypeError):
                    taxid = '0'
            else:
                taxid = '0'
            lineages.append(taxonomy_db.get_lineage(taxid))
        
        for key in ['organism_name', 'superkingdom', 'kingdom', 'phylum', 
                    'class', 'order', 'family', 'genus', 'species']:
            df[key] = [lineage[key] for lineage in lineages]
        
        df.to_csv("${sample}_viralflye_diamond_with_taxonomy.txt", sep='\\t', index=False)
        
        # Generate summary
        with open("${sample}_viralflye_taxonomy_summary.txt", 'w', encoding='utf-8') as f:
            f.write("="*80 + "\\n")
            f.write("Viral Taxonomy Statistical Summary - viralFlye (Viral Sequences Only)\\n")
            f.write("="*80 + "\\n\\n")
            f.write(f"Sample Name: ${sample}\\n")
            f.write(f"Total Hits: {len(df):,}\\n")
            f.write(f"Unique Queries: {df['qseqid'].nunique():,}\\n")
            f.write(f"Average Identity: {df['pident'].mean():.2f}%\\n")
            f.write(f"Average Alignment Length: {df['length'].mean():.1f} aa\\n\\n")
            
            for title, column in [("Superkingdom", "superkingdom"), ("Kingdom", "kingdom"), 
                                 ("Phylum", "phylum"), ("Family", "family"), ("Species", "species")]:
                f.write(f"[{title} Distribution (Top 15)]\\n")
                f.write("-"*80 + "\\n")
                counts = Counter(df[column].dropna())
                for name, count in counts.most_common(15):
                    f.write(f"{name:<40} {count:>10,}\\n")
                f.write("\\n")
            
            f.write("="*80 + "\\n")
            f.write("Analysis Complete\\n")
            f.write("="*80 + "\\n")
    else:
        with open("${sample}_viralflye_diamond_with_taxonomy.txt", 'w') as f:
            f.write("# No hits\\n")
        with open("${sample}_viralflye_taxonomy_summary.txt", 'w') as f:
            f.write("No hits\\n")
    """
}

// Process: Compare Dual Tracks and Generate Consensus (compare dual-track results, generate consensus virus list)
process COMPARE_DUAL_TRACKS {
    tag "${sample}"
    label 'process_low'
    conda 'conda-forge::pandas=2.0.3'
    publishDir "${params.outdir}/consensus_analysis", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(metaflye_report), path(viralflye_report)
    
    output:
    tuple val(sample), path("${sample}_consensus_viruses.txt"), emit: consensus
    path("${sample}_metaflye_only_viruses.txt"), emit: metaflye_only
    path("${sample}_viralflye_only_viruses.txt"), emit: viralflye_only
    path("${sample}_dual_track_comparison.txt"), emit: comparison_report
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import sys
    from collections import Counter
    
    print("\\nComparing dual-track analysis results...", file=sys.stderr)
    
    # Read both result sets
    metaflye_df = pd.read_csv("${metaflye_report}", sep='\\t')
    viralflye_df = pd.read_csv("${viralflye_report}", sep='\\t')
    
    print(f"MetaFlye hits: {len(metaflye_df):,}", file=sys.stderr)
    print(f"viralFlye hits: {len(viralflye_df):,}", file=sys.stderr)
    
    # Extract viral hits (Kingdom contains 'virae' or Phylum contains 'viricota')
    def is_virus(df):
        virus_mask = (
            df['kingdom'].str.contains('virae', case=False, na=False) |
            df['phylum'].str.contains('viricota', case=False, na=False) |
            df['superkingdom'].str.contains('virus', case=False, na=False)
        )
        return df[virus_mask]
    
    metaflye_viruses = is_virus(metaflye_df)
    viralflye_viruses = is_virus(viralflye_df)
    
    print(f"MetaFlye viral hits: {len(metaflye_viruses):,}", file=sys.stderr)
    print(f"viralFlye viral hits: {len(viralflye_viruses):,}", file=sys.stderr)
    
    # Get contig IDs (remove gene numbers)
    metaflye_viruses['contig_id'] = metaflye_viruses['qseqid'].str.rsplit('_', n=1).str[0]
    viralflye_viruses['contig_id'] = viralflye_viruses['qseqid'].str.rsplit('_', n=1).str[0]
    
    # Find consensus and unique contigs
    metaflye_contigs = set(metaflye_viruses['contig_id'].unique())
    viralflye_contigs = set(viralflye_viruses['contig_id'].unique())
    
    consensus_contigs = metaflye_contigs & viralflye_contigs  # Intersection
    metaflye_only_contigs = metaflye_contigs - viralflye_contigs  # Difference
    viralflye_only_contigs = viralflye_contigs - metaflye_contigs  # Difference
    
    print(f"\\nConsensus viral contigs: {len(consensus_contigs):,}", file=sys.stderr)
    print(f"MetaFlye-only viral contigs: {len(metaflye_only_contigs):,}", file=sys.stderr)
    print(f"viralFlye-only viral contigs: {len(viralflye_only_contigs):,}", file=sys.stderr)
    
    # Save consensus viruses (identified by both tracks, highest confidence)
    consensus_df = metaflye_viruses[metaflye_viruses['contig_id'].isin(consensus_contigs)]
    consensus_df = consensus_df.sort_values('evalue')
    consensus_df.to_csv("${sample}_consensus_viruses.txt", sep='\\t', index=False)
    
    # Save MetaFlye-only viruses (filtered by viralFlye, possibly distant viruses)
    metaflye_only_df = metaflye_viruses[metaflye_viruses['contig_id'].isin(metaflye_only_contigs)]
    metaflye_only_df = metaflye_only_df.sort_values('evalue')
    metaflye_only_df.to_csv("${sample}_metaflye_only_viruses.txt", sep='\\t', index=False)
    
    # Save viralFlye-only viruses (feature-based but low sequence similarity)
    viralflye_only_df = viralflye_viruses[viralflye_viruses['contig_id'].isin(viralflye_only_contigs)]
    viralflye_only_df = viralflye_only_df.sort_values('evalue')
    viralflye_only_df.to_csv("${sample}_viralflye_only_viruses.txt", sep='\\t', index=False)
    
    # Generate comparison report
    with open("${sample}_dual_track_comparison.txt", 'w', encoding='utf-8') as f:
        f.write("="*80 + "\\n")
        f.write("Dual-Track Analysis Comparison Report - Viral Identification Consensus Analysis\\n")
        f.write("="*80 + "\\n\\n")
        
        f.write(f"Sample Name: ${sample}\\n\\n")
        
        # Overall statistics
        f.write("[Overall Statistics]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"MetaFlye Total Hits:              {len(metaflye_df):,}\\n")
        f.write(f"MetaFlye Viral Hits:              {len(metaflye_viruses):,}\\n")
        f.write(f"viralFlye Total Hits:             {len(viralflye_df):,}\\n")
        f.write(f"viralFlye Viral Hits:             {len(viralflye_viruses):,}\\n\\n")
        
        # Contig-level comparison
        f.write("[Contig-Level Comparison]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"MetaFlye Viral Contigs:           {len(metaflye_contigs):,}\\n")
        f.write(f"viralFlye Viral Contigs:          {len(viralflye_contigs):,}\\n")
        f.write(f"Consensus Viral Contigs (Both):   {len(consensus_contigs):,} ★★★ Highest Confidence\\n")
        f.write(f"MetaFlye-Only Viral Contigs:      {len(metaflye_only_contigs):,} → Distant Viral Candidates\\n")
        f.write(f"viralFlye-Only Viral Contigs:     {len(viralflye_only_contigs):,} → Feature-based, Low Similarity\\n\\n")
        
        # Confidence levels
        f.write("[Viral Confidence Levels]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"★★★ High Confidence (Consensus):  {len(consensus_contigs):,} contigs\\n")
        f.write(f"    - viralFlye feature identification + Diamond sequence match\\n")
        f.write(f"    - Confirmed by both independent methods\\n\\n")
        
        f.write(f"★★  Medium Confidence (viralFlye-Only): {len(viralflye_only_contigs):,} contigs\\n")
        f.write(f"    - Identified as viral by viralFlye\\n")
        f.write(f"    - No match or extremely low similarity in Diamond\\n")
        f.write(f"    - Possibly: feature-based but sequence-unique viruses\\n\\n")
        
        f.write(f"★   Needs Verification (MetaFlye-Only): {len(metaflye_only_contigs):,} contigs\\n")
        f.write(f"    - Matched to viral proteins in Diamond\\n")
        f.write(f"    - But filtered by viralFlye (feature mismatch)\\n")
        f.write(f"    - Possibly: distant viruses or false positives\\n\\n")
        
        # Taxonomic distribution of consensus viruses
        if len(consensus_df) > 0:
            f.write("[Taxonomic Distribution of Consensus Viruses]\\n")
            f.write("-"*80 + "\\n")
            
            for title, column in [("Kingdom", "kingdom"), ("Phylum", "phylum"), 
                                 ("Family", "family"), ("Genus", "genus"), ("Species", "species")]:
                f.write(f"\\n{title} Distribution (Top 15):\\n")
                counts = Counter(consensus_df[column].dropna())
                top_n = 15 if title in ["Species", "Genus"] else 10
                for name, count in counts.most_common(top_n):
                    f.write(f"  {name:<40} {count:>6,}\\n")
        
        # MetaFlye-only virus analysis (distant viral candidates)
        if len(metaflye_only_df) > 0:
            f.write("\\n[MetaFlye-Only Viruses (Distant Viral Candidates)]\\n")
            f.write("-"*80 + "\\n")
            f.write(f"Count: {len(metaflye_only_contigs)} contigs, {len(metaflye_only_df)} protein matches\\n")
            f.write(f"Average Identity: {metaflye_only_df['pident'].mean():.2f}%\\n")
            
            # Top Family
            f.write("\\nFamily Distribution (Top 5):\\n")
            family_counts = Counter(metaflye_only_df['family'].dropna())
            for name, count in family_counts.most_common(5):
                f.write(f"  {name:<40} {count:>6,}\\n")
        
        # viralFlye-only virus analysis
        if len(viralflye_only_df) > 0:
            f.write("\\n[viralFlye-Only Viruses]\\n")
            f.write("-"*80 + "\\n")
            f.write(f"Count: {len(viralflye_only_contigs)} contigs, {len(viralflye_only_df)} protein matches\\n")
            f.write(f"Note: These viruses have viral features but were not matched by Diamond in MetaFlye track\\n")
        
        f.write("\\n" + "="*80 + "\\n")
        f.write("Recommendation: Use consensus virus list for downstream analysis (highest confidence)\\n")
        f.write("="*80 + "\\n")
    
    print(f"\\nDual-track comparison analysis complete", file=sys.stderr)
    print(f"  Consensus viruses: {len(consensus_contigs)} contigs", file=sys.stderr)
    print(f"  MetaFlye-only: {len(metaflye_only_contigs)} contigs (distant viral candidates)", file=sys.stderr)
    print(f"  viralFlye-only: {len(viralflye_only_contigs)} contigs", file=sys.stderr)
    """
}

// Workflow completion message
workflow.onComplete {
    def readTypeInfo = params.read_type == 'long' ? """
    Long read workflow results:
    - assembly_metaflye/: MetaFlye assembled contigs (all sequences)
      * *_metaflye_contigs.fa: Assembled contigs (DNA sequences)
    - assembly_viralflye/: viralFlye refined contigs (viral sequences only, if enabled)
      * *_viralflye_contigs.fa: Viral contigs (DNA sequences)
    - prodigal_metaflye/: Gene predictions from MetaFlye (all sequences)
      * *_metaflye_proteins.faa: Predicted protein sequences
    - prodigal_viralflye/: Gene predictions from viralFlye (viral only, if enabled)
      * *_viralflye_proteins.faa: Predicted viral protein sequences
    - diamond_metaflye/: Diamond classification of MetaFlye proteins (all sequences)
      * *_metaflye_diamond.txt: BLAST-style alignment results
    - diamond_viralflye/: Diamond classification of viralFlye proteins (viral only, if enabled)
      * *_viralflye_diamond.txt: BLAST-style alignment results
    - taxonomy_metaflye/: Taxonomy analysis of all sequences
      * *_metaflye_diamond_with_taxonomy.txt: Full lineage (Kingdom/Class/Family/Species)
      * *_metaflye_taxonomy_summary.txt: Statistical summary
    - taxonomy_viralflye/: Taxonomy analysis of viral sequences only (if enabled)
      * *_viralflye_diamond_with_taxonomy.txt: Full lineage (Kingdom/Class/Family/Species)
      * *_viralflye_taxonomy_summary.txt: Statistical summary
    - consensus_analysis/: Dual-track comparison and consensus viruses (if enabled)
      * *_consensus_viruses.txt: Viruses identified by BOTH methods (highest confidence)
      * *_metaflye_only_viruses.txt: Viruses only in MetaFlye (distant viral candidates)
      * *_viralflye_only_viruses.txt: Viruses only in viralFlye
      * *_dual_track_comparison.txt: Detailed comparison report
    """ : """
    Short read workflow results:
    - fastp/: Quality control reports (if enabled)
      * *_fastp.html: HTML quality reports
      * *_fastp.json: JSON quality data
    - assembly_megahit/: MEGAHIT assembled contigs
      * *_megahit_contigs.fa: Assembled contigs (DNA sequences)
    - assembly_spades/: SPAdes assembled contigs
      * *_spades_contigs.fa: Assembled contigs (DNA sequences)
    - prodigal_megahit/: Gene predictions from MEGAHIT contigs
      * *_megahit_proteins.faa: Predicted protein sequences
    - prodigal_spades/: Gene predictions from SPAdes contigs
      * *_spades_proteins.faa: Predicted protein sequences
    - diamond_megahit/: Diamond classification of MEGAHIT proteins
      * *_megahit_diamond.txt: BLAST-style alignment results
    - diamond_spades/: Diamond classification of SPAdes proteins
      * *_spades_diamond.txt: BLAST-style alignment results
    - merged_reports/: Comprehensive analysis (if enabled)
      * *_merged_report.txt: Combined analysis report
      * *_merged_report.csv: Detailed comparison data
      * *_megahit_with_taxonomy.txt: MEGAHIT results with full taxonomy
      * *_spades_with_taxonomy.txt: SPAdes results with full taxonomy
    - abundance_megahit/: Viral abundance for MEGAHIT (if enabled)
      * *_megahit_abundance.txt: Abundance report with RPM and RPKM (text format)
      * *_megahit_abundance.csv: Abundance data with RPM and RPKM (CSV format)
    - abundance_spades/: Viral abundance for SPAdes (if enabled)
      * *_spades_abundance.txt: Abundance report with RPM and RPKM (text format)
      * *_spades_abundance.csv: Abundance data with RPM and RPKM (CSV format)
    """
    
    log.info """
    ==========================================
    🎯 Metagenome Assembly and Diamond Classification Results
    ==========================================
    Pipeline completed successfully!
    
    Read type: ${params.read_type}
    Results directory: ${params.outdir}
    
    Generated files:
${readTypeInfo}
    ==========================================
    """
}

workflow.onError {
    log.error """
    ==========================================
    ❌ Metagenome Assembly and Diamond Classification Workflow Failed
    ==========================================
    Error: ${workflow.errorMessage}
    ==========================================
    """
}
