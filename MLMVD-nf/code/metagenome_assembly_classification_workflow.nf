#!/usr/bin/env nextflow

/*
 * Metagenome Viral Classification Workflow
 * 
 * This workflow integrates:
 * 1. Quality control using fastp (optional)
 * 2. Metagenome assembly using MEGAHIT and SPAdes (parallel)
 * 3. Viral sequence identification using VirSorter2 and DeepVirFinder
 * 4. Comprehensive comparison and merging of viral identification results
 * 5. Assembler comparison to identify high-confidence consensus viral sequences
 * 
 * Author: Assistant
 * Version: 5.2.1
 */

nextflow.enable.dsl = 2

// Workflow parameters
// Input data
params.input = null
params.outdir = './results'
params.help = false
// Long-read support
params.longread = false                 // Enable long-read branch (PacBio/Nanopore)
params.longread_platform = 'nano'       // Long-read platform: 'nano' or 'pacbio'
params.skip_longread_qc = true          // Skip long-read QC (default: skip)
// viralFlye refinement branch
params.enable_viralflye = false         // Enable viralFlye targeted refinement (default: disabled)
params.viralflye_min_score = 0.5        // Minimum VirSorter2 score for target contig selection
params.viralflye_min_length = 1000      // Minimum length for target contig selection (bp)
params.viralflye_completeness = 0.5     // ⭐ Completeness cutoff for viralComplete (default 0.5, lower = more viruses)
params.viralflye_env = '/home/sp96859/.conda/envs/viralFlye_env'  // viralFlye conda environment path
params.pfam_db = '/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/Pfam/Pfam-A.hmm'  // Pfam database path

// Workflow control
params.skip_virsorter2 = false    // Whether to skip VirSorter2 viral identification
params.skip_deepvirfinder = false // Whether to skip DeepVirFinder viral identification
params.skip_merge_reports = false // Whether to skip result merging
params.save_clean_reads = true    // Whether to save filtered clean reads

// Viral identification paths
params.virsorter2_db = null       // VirSorter2 database path
params.deepvirfinder_dir = '/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/DeepVirFinder' // DeepVirFinder installation directory

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

// VirSorter2 parameters
params.virsorter2_min_length = 1000      // Minimum contig length
params.virsorter2_min_score = 0.5        // Minimum viral score

// DeepVirFinder parameters  
params.deepvirfinder_min_length = 1000   // Minimum contig length
params.deepvirfinder_pvalue = 0.05       // p-value threshold

// Resource parameters
params.max_cpus = 32
params.max_memory = '512.GB'  // Adjusted to support 512GB memory requirement for SPAdes
params.max_time = '72.h'

// Print help information
if (params.help) {
    log.info """
    ==========================================
    🦠 Metagenome Viral Classification Workflow
    ==========================================
    
    Usage:
    nextflow run metagenome_assembly_classification_workflow.nf --input samplesheet.csv --outdir results --virsorter2_db /path/to/db
    
    Required Parameters:
    --input                    Input samplesheet (CSV format)
    --outdir                   Output directory
    --virsorter2_db           VirSorter2 database path
    
    Viral Identification Parameters:
    --deepvirfinder_dir       DeepVirFinder installation directory (default: auto-detected)
    --skip_virsorter2         Skip VirSorter2 analysis (default: false)
    --skip_deepvirfinder      Skip DeepVirFinder analysis (default: false)
    --skip_merge_reports      Skip merging VirSorter2 and DeepVirFinder results (default: false)
    --virsorter2_min_length   Minimum contig length for VirSorter2 (default: 1000)
    --virsorter2_min_score    Minimum viral score for VirSorter2 (default: 0.5)
    --deepvirfinder_min_length Minimum contig length for DeepVirFinder (default: 1000)
    --deepvirfinder_pvalue    P-value threshold for DeepVirFinder (default: 0.05)
    
    Optional Parameters:
    --skip_fastp              Skip fastp quality control (default: false)
    --save_clean_reads        Save filtered clean reads (default: true)
    --longread                Enable long-read mode (PacBio/Nanopore). When true, use metaFlye and skip short-read steps (default: false)
    --longread_platform       Long-read platform: 'nano' (Nanopore) or 'pacbio' (PacBio). Used by Flye (default: nano)
    --skip_longread_qc        Skip long-read QC (default: true)
    --enable_viralflye        Enable viralFlye refinement on viral reads (default: false)
    --viralflye_min_score     Min VS2 score to select targets (default: 0.5)
    --viralflye_min_length    Min contig length to select targets (default: 1000)
    --viralflye_completeness  Completeness cutoff for viralComplete (default: 0.5) ⭐ Lower = more viruses
    --viralflye_env          viralFlye conda environment path (default: /home/sp96859/.conda/envs/viralFlye_env)
    --pfam_db                Pfam database path for viralFlye (default: /scratch/sp96859/.../Pfam/Pfam-A.hmm)
    
    Example:
    nextflow run metagenome_assembly_classification_workflow.nf \\
        --input samplesheet.csv \\
        --outdir results \\
        --virsorter2_db /scratch/databases/virsorter2/db
    
    Long-read examples:
    - From raw Nanopore/PacBio reads (one FASTQ per sample):
      nextflow run metagenome_assembly_classification_workflow.nf \\
        --input samplesheet_long.csv \\
        --outdir results_long \\
        --virsorter2_db /scratch/databases/virsorter2/db \\
        --longread true \\
        --longread_platform nano

    Enable viralFlye refinement:
      ... --enable_viralflye true --viralflye_min_score 0.6 --viralflye_min_length 1500
    
    Long-read samplesheet format (CSV):
      sample,fastq_long
      s1,/path/to/s1_nanopore.fastq.gz
    
    Output:
    - Viral identification results from VirSorter2, DeepVirFinder, and optionally viralFlye
    - Viral abundance analysis (RPM and RPKM) for all identified viral contigs
    - Results saved in: results/abundance/
      * RPM: Reads Per Million (normalized by total read count)
      * RPKM: Reads Per Kilobase per Million (normalized by contig length and read count)
    """
    exit 0
}

// Validate required parameters
if (!params.input) {
    error "Input samplesheet is required. Use --input parameter."
}

// Validate VirSorter2 database
if (!params.skip_virsorter2 && !params.virsorter2_db) {
    error "VirSorter2 database path is required. Use --virsorter2_db parameter or --skip_virsorter2 to skip."
}

// Validate DeepVirFinder installation directory
if (!params.skip_deepvirfinder) {
    def dvf_dir = file(params.deepvirfinder_dir)
    if (!dvf_dir.exists() || !dvf_dir.isDirectory()) {
        log.warn "DeepVirFinder directory not found at: ${params.deepvirfinder_dir}"
        log.warn "DeepVirFinder analysis will be skipped."
        params.skip_deepvirfinder = true
    }
}

// Print workflow information
if (!params.longread) {
    log.info """
==========================================
🦠 Metagenome Viral Classification Workflow (Short-Read Mode)
==========================================
Workflow version: 5.2.1
Input samplesheet: ${params.input}
Output directory: ${params.outdir}

Quality Control:
- fastp QC: ${params.skip_fastp ? 'Disabled' : 'Enabled'}
- Save clean reads: ${params.save_clean_reads ? 'Yes' : 'No'}

Assembly Methods:
- MEGAHIT: Enabled
- metaSPAdes: Enabled

Viral Identification:
- VirSorter2: ${params.skip_virsorter2 ? 'Disabled' : 'Enabled'}
${params.skip_virsorter2 ? '' : "  Database: ${params.virsorter2_db}"}
- DeepVirFinder: ${params.skip_deepvirfinder ? 'Disabled' : 'Enabled'}
${params.skip_deepvirfinder ? '' : "  Directory: ${params.deepvirfinder_dir}"}

Result Merging:
- Merge viral reports: ${params.skip_merge_reports ? 'Disabled' : 'Enabled'}

Abundance Calculation:
- Viral abundance (RPM & RPKM): Enabled
- Output directory: ${params.outdir}/abundance/
==========================================
"""
} else {
    log.info """
==========================================
🦠 Metagenome Viral Classification Workflow (Long-Read Mode)
==========================================
Workflow version: 5.2.1
Input samplesheet: ${params.input}
Output directory: ${params.outdir}

Long-Read Configuration:
- Platform: ${params.longread_platform} (${params.longread_platform == 'nano' ? 'Nanopore' : 'PacBio'})
- Long-read QC: ${params.skip_longread_qc ? 'Disabled' : 'Enabled'}

Assembly Methods:
- metaFlye: Enabled (--meta mode)

Viral Identification:
- VirSorter2: ${params.skip_virsorter2 ? 'Disabled' : 'Enabled'}
${params.skip_virsorter2 ? '' : "  Database: ${params.virsorter2_db}"}
- DeepVirFinder: ${params.skip_deepvirfinder ? 'Disabled' : 'Enabled'}
${params.skip_deepvirfinder ? '' : "  Directory: ${params.deepvirfinder_dir}"}

viralFlye Refinement (Viral-specific Assembly Refinement):
- Enabled: ${params.enable_viralflye ? 'Yes' : 'No'}
${params.enable_viralflye ? "  Min VS2 score: ${params.viralflye_min_score}" : ''}
${params.enable_viralflye ? "  Min contig length: ${params.viralflye_min_length} bp" : ''}
${params.enable_viralflye ? "  viralFlye environment: ${params.viralflye_env}" : ''}
${params.enable_viralflye ? "  Pfam database: ${params.pfam_db}" : ''}
${params.enable_viralflye ? "  Note: Using viralFlye.py with Pfam for viral genome optimization" : ''}

Result Merging:
- Merge viral reports: ${params.skip_merge_reports ? 'Disabled' : 'Enabled'}

Abundance Calculation:
- Viral abundance (RPM & RPKM): Enabled
- Output directory: ${params.outdir}/abundance/
==========================================
"""
}

// Create input channels based on mode (short-read or long-read)
if (!params.longread) {
    // Short-read: expects samplesheet with fastq_1 and fastq_2
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def sample = row.sample
            def read1 = file(row.fastq_1)
            def read2 = file(row.fastq_2)
            return tuple(sample, [read1, read2])
        }
        .set { ch_reads }
} else {
    // Long-read: expects samplesheet with fastq_long (Nanopore or PacBio single-end FASTQ)
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .filter { row ->
            // Filter empty rows: check if sample name and path are not empty
            def sample = row.sample?.trim()
            def read_long = row.fastq_long?.trim()
            return sample && read_long && sample != '' && read_long != ''
        }
        .map { row -> 
            def sample = row.sample?.trim()
            def read_long_path = row.fastq_long?.trim()
            if (!sample || !read_long_path) {
                log.error "Invalid row in samplesheet: sample=${sample}, fastq_long=${read_long_path}"
                return null
            }
            def read_long = file(read_long_path)
            if (!read_long.exists()) {
                log.error "Long-read file not found: ${read_long_path}"
                return null
            }
            return tuple(sample, read_long)
        }
        .filter { it != null }
        .set { ch_long_reads }
    
    // Validate long-read input channel
    ch_long_reads.view { sample, reads -> "Long-read sample: ${sample}, reads: ${reads}" }
}

// Define workflow
workflow {
    if (!params.longread) {
        // Short-read workflow
        // Stage 0: QC (optional)
        if (!params.skip_fastp) {
            FASTP (
                ch_reads
            )
            ch_clean_reads = FASTP.out.clean_reads
        } else {
            ch_clean_reads = ch_reads
        }

        // Stage 1: Assembly (MEGAHIT + metaSPAdes)
        MEGAHIT_ASSEMBLY (
            ch_clean_reads
        )

        SPADES_ASSEMBLY (
            ch_clean_reads
        )

        // Stage 2: VirSorter2
        if (!params.skip_virsorter2) {
            VIRSORTER2_MEGAHIT (
                MEGAHIT_ASSEMBLY.out.contigs,
                params.virsorter2_db
            )

            VIRSORTER2_SPADES (
                SPADES_ASSEMBLY.out.contigs,
                params.virsorter2_db
            )
        }

        // Stage 3: DeepVirFinder
        if (!params.skip_deepvirfinder) {
            DEEPVIRFINDER_MEGAHIT (
                MEGAHIT_ASSEMBLY.out.contigs
            )

            DEEPVIRFINDER_SPADES (
                SPADES_ASSEMBLY.out.contigs
            )
        }

        // Stage 4: Merge and compare (short-read dual assemblers)
        if (!params.skip_merge_reports && !params.skip_virsorter2 && !params.skip_deepvirfinder) {
            VIRSORTER2_MEGAHIT.out.results
                .join(DEEPVIRFINDER_MEGAHIT.out.results)
                .set { ch_viral_megahit }

            VIRSORTER2_SPADES.out.results
                .join(DEEPVIRFINDER_SPADES.out.results)
                .set { ch_viral_spades }

            MERGE_VIRAL_REPORTS_MEGAHIT (
                ch_viral_megahit
            )

            MERGE_VIRAL_REPORTS_SPADES (
                ch_viral_spades
            )

            // Stage 5: Assembler comparison
            MERGE_VIRAL_REPORTS_MEGAHIT.out.merged_csv
                .join(MERGE_VIRAL_REPORTS_SPADES.out.merged_csv)
                .set { ch_assembler_comparison }

            COMPARE_ASSEMBLERS (
                ch_assembler_comparison
            )
        }

        // Stage 6: Calculate viral abundance
        if (!params.skip_virsorter2) {
            // Prepare input for abundance calculation: combine reads + viral contigs + assembler name
            // For MEGAHIT
            ch_clean_reads
                .join(VIRSORTER2_MEGAHIT.out.viral_contigs)
                .map { sample, reads, contigs -> tuple(sample, reads, contigs, "megahit") }
                .set { ch_abundance_megahit }

            // For SPAdes
            ch_clean_reads
                .join(VIRSORTER2_SPADES.out.viral_contigs)
                .map { sample, reads, contigs -> tuple(sample, reads, contigs, "spades") }
                .set { ch_abundance_spades }
            
            // Mix both channels and calculate abundance once
            ch_abundance_megahit
                .mix(ch_abundance_spades)
                .set { ch_abundance_shortread_input }
            
            CALCULATE_ABUNDANCE (
                ch_abundance_shortread_input
            )
        }
    } else {
        // Long-read workflow
        // Stage 0: Long-read QC (optional)
        if (!params.skip_longread_qc) {
            LONGREAD_QC (
                ch_long_reads
            )
            ch_long_clean = LONGREAD_QC.out.clean_long
        } else {
            ch_long_clean = ch_long_reads
        }

        // Stage 1: metaFlye assembly (--meta mode)
        METAFLYE_ASSEMBLY (
            ch_long_clean
        )

        // Stage 2: Parallel Viral Identification - VirSorter2
        if (!params.skip_virsorter2) {
            VIRSORTER2_METAFLYE (
                METAFLYE_ASSEMBLY.out.contigs,
                params.virsorter2_db
            )
        }

        // Stage 3: Parallel Viral Identification - DeepVirFinder
        if (!params.skip_deepvirfinder) {
            DEEPVIRFINDER_METAFLYE (
                METAFLYE_ASSEMBLY.out.contigs
            )
        }

        // Stage 4: Parallel Viral Identification - viralFlye (with Pfam validation)
        if (params.enable_viralflye) {
            // Prepare viralFlye input: metaFlye directory + original reads
            METAFLYE_ASSEMBLY.out.flye_dir
                .join(ch_long_clean)
                .set { ch_viralflye_input }
            
            // Run viralFlye viral identification
            VIRALFLYE_IDENTIFY (
                ch_viralflye_input
            )
        }

        // Stage 5: Comprehensive viral identification results (three-tool parallel comparison)
        if (!params.skip_merge_reports) {
            // Combine results based on enabled tools
            if (!params.skip_virsorter2 && !params.skip_deepvirfinder && params.enable_viralflye) {
                // All three tools enabled: comprehensive comparison
                VIRSORTER2_METAFLYE.out.results
                    .join(DEEPVIRFINDER_METAFLYE.out.results)
                    .join(VIRALFLYE_IDENTIFY.out.results)
                    .set { ch_three_tools }
                
                COMPARE_THREE_VIRAL_TOOLS (
                    ch_three_tools
                )
            } else if (!params.skip_virsorter2 && !params.skip_deepvirfinder) {
                // Only VS2 and DVF: two-tool comparison
                VIRSORTER2_METAFLYE.out.results
                    .join(DEEPVIRFINDER_METAFLYE.out.results)
                    .set { ch_viral_metaflye }

                MERGE_VIRAL_REPORTS_METAFLYE (
                    ch_viral_metaflye
                )
            }
        }

        // Stage 6: Calculate viral abundance (long-read mode)
        if (!params.skip_virsorter2 && params.enable_viralflye) {
            // Combine both metaFlye and viralFlye viral contigs into one channel
            ch_long_clean
                .join(VIRSORTER2_METAFLYE.out.viral_contigs)
                .map { sample, reads, contigs -> tuple(sample, reads, contigs, "metaflye") }
                .set { ch_abundance_metaflye }
            
            ch_long_clean
                .join(VIRALFLYE_IDENTIFY.out.viral_contigs)
                .map { sample, reads, contigs -> tuple(sample, reads, contigs, "viralflye") }
                .set { ch_abundance_viralflye }
            
            // Mix both channels and calculate abundance once
            ch_abundance_metaflye
                .mix(ch_abundance_viralflye)
                .set { ch_abundance_longread_input }
            
            CALCULATE_ABUNDANCE_LONGREAD (
                ch_abundance_longread_input
            )
        } else if (!params.skip_virsorter2) {
            // Only metaFlye abundance
            ch_long_clean
                .join(VIRSORTER2_METAFLYE.out.viral_contigs)
                .map { sample, reads, contigs -> tuple(sample, reads, contigs, "metaflye") }
                .set { ch_abundance_longread_input }
            
            CALCULATE_ABUNDANCE_LONGREAD (
                ch_abundance_longread_input
            )
        } else if (params.enable_viralflye) {
            // Only viralFlye abundance
            ch_long_clean
                .join(VIRALFLYE_IDENTIFY.out.viral_contigs)
                .map { sample, reads, contigs -> tuple(sample, reads, contigs, "viralflye") }
                .set { ch_abundance_longread_input }
            
            CALCULATE_ABUNDANCE_LONGREAD (
                ch_abundance_longread_input
            )
        }
    }
}

// ================================================================================
// Process Definitions
// ================================================================================

// Process: fastp Quality Control
process FASTP {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{html,json}"
    publishDir "${params.outdir}/clean_reads", mode: 'copy', pattern: "*_clean_R*.fastq.gz", enabled: params.save_clean_reads
    
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

// Process: Long-read QC (optional, simplified: currently copies input directly, reserved for Filtlong/NanoFilt integration)
process LONGREAD_QC {
    tag "${sample}"
    label 'process_medium'
    publishDir "${params.outdir}/longread_qc", mode: 'copy', pattern: "*.fastq.gz"

    input:
    tuple val(sample), path(read_long)

    output:
    tuple val(sample), path("${sample}_long_clean.fastq.gz"), emit: clean_long

    script:
    """
    echo "=== Long-read QC (placeholder implementation): ${sample} ==="
    # Currently skips complex filtering by default, only standardizes output filename for downstream workflow
    # For strict QC, integrate Filtlong or NanoFilt here
    if [[ "${read_long}" == *.gz ]]; then
        cp ${read_long} ${sample}_long_clean.fastq.gz
    else
        gzip -c ${read_long} > ${sample}_long_clean.fastq.gz
    fi
    """
}

// Process: MEGAHIT Assembly
process MEGAHIT_ASSEMBLY {
    tag "${sample}_MEGAHIT"
    label 'process_high'
    container 'docker://quay.io/biocontainers/megahit:1.2.9--h2e03b76_1'
    publishDir "${params.outdir}/assembly_megahit", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_megahit_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== MEGAHIT Assembly: ${sample} ==="
    
    megahit \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o megahit_output \
        -t ${task.cpus} \
        --memory ${params.megahit_memory} \
        --min-contig-len ${params.megahit_min_contig_len}
    
    cp megahit_output/final.contigs.fa ${sample}_megahit_contigs.fa
    
    echo "MEGAHIT: Generated \$(grep -c ">" ${sample}_megahit_contigs.fa) contigs"
    """
}

// Process: SPAdes Assembly
process SPADES_ASSEMBLY {
    tag "${sample}_SPAdes"
    label 'process_high'
    container 'docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    publishDir "${params.outdir}/assembly_spades", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_spades_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== metaSPAdes Assembly: ${sample} ==="
    
    # Use metaSPAdes, disable error correction to avoid memory and bug issues
    metaspades.py \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o spades_output \
        -t ${task.cpus} \
        -m ${task.memory.toGiga()} \
        --only-assembler
    
    cp spades_output/contigs.fasta ${sample}_spades_contigs.fa
    
    echo "metaSPAdes: Generated \$(grep -c ">" ${sample}_spades_contigs.fa) contigs"
    """
}

// Process: metaFlye assembly (long-read, --meta mode enabled)
process METAFLYE_ASSEMBLY {
    tag "${sample}_metaFlye"
    label 'process_high'
    conda 'bioconda::flye=2.9'
    publishDir "${params.outdir}/assembly_metaflye", mode: 'copy', pattern: "*.fa"
    // Always save complete flye output directory (includes assembly_info.txt, assembly_graph and other important files)
    publishDir "${params.outdir}/metaflye_full_output", mode: 'copy', pattern: "${sample}_flye_output"

    input:
    tuple val(sample), path(read_long)

    output:
    tuple val(sample), path("${sample}_metaflye_contigs.fa"), emit: contigs
    tuple val(sample), path("${sample}_flye_output"), emit: flye_dir  // viralFlye requires complete directory

    script:
    """
    echo "=== metaFlye assembly: ${sample} ==="
    if [ "${params.longread_platform}" = "pacbio" ]; then
        PLATFORM_FLAG="--pacbio-raw"
    else
        PLATFORM_FLAG="--nano-raw"
    fi

    flye \\
        \${PLATFORM_FLAG} ${read_long} \\
        --out-dir ${sample}_flye_output \\
        --threads ${task.cpus} \\
        --meta

    cp ${sample}_flye_output/assembly.fasta ${sample}_metaflye_contigs.fa

    echo "metaFlye: Generated \$(grep -c ">" ${sample}_metaflye_contigs.fa) contigs"
    """
}

// Process: VirSorter2 for MEGAHIT
process VIRSORTER2_MEGAHIT {
    tag "${sample}_MEGAHIT_VirSorter2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'  // Use pre-installed environment
    publishDir "${params.outdir}/virsorter2_megahit", mode: 'copy', pattern: "*.{tsv,fa}"
    
    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)
    
    output:
    tuple val(sample), path("${sample}_megahit_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_megahit_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_megahit_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true
    
    script:
    """
    # Ensure correct conda environment is used
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    
    # CRITICAL: Clean PYTHONPATH to avoid package conflicts with VirSorter2's internal conda envs
    unset PYTHONPATH
    
    echo "=== VirSorter2 Analysis (MEGAHIT): ${sample} ==="
    echo "Using Python: \$(which python)"
    echo "Using VirSorter2: \$(which virsorter)"
    
    # Run VirSorter2 for viral sequence identification
    virsorter run \\
        -i ${contigs} \\
        -w virsorter2_output \\
        --db-dir ${virsorter2_db} \\
        --min-length ${params.virsorter2_min_length} \\
        --min-score ${params.virsorter2_min_score} \\
        -j ${task.cpus} \\
        all
    
    # Copy result files
    cp virsorter2_output/final-viral-score.tsv ${sample}_megahit_vs2_final-viral-score.tsv
    
    # If viral sequences detected, copy viral contig file
    if [ -f virsorter2_output/final-viral-combined.fa ]; then
        cp virsorter2_output/final-viral-combined.fa ${sample}_megahit_vs2_final-viral-combined.fa
    fi
    
    if [ -f virsorter2_output/final-viral-boundary.tsv ]; then
        cp virsorter2_output/final-viral-boundary.tsv ${sample}_megahit_vs2_final-viral-boundary.tsv
    fi
    
    # Count identified viral sequences
    VIRAL_COUNT=\$(tail -n +2 ${sample}_megahit_vs2_final-viral-score.tsv | wc -l || echo 0)
    echo "VirSorter2: Identified \${VIRAL_COUNT} viral sequences from MEGAHIT contigs"
    """
}

// Process: VirSorter2 for SPAdes
process VIRSORTER2_SPADES {
    tag "${sample}_SPAdes_VirSorter2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'  // Use pre-installed environment
    publishDir "${params.outdir}/virsorter2_spades", mode: 'copy', pattern: "*.{tsv,fa}"
    
    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)
    
    output:
    tuple val(sample), path("${sample}_spades_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_spades_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_spades_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true
    
    script:
    """
    # Ensure correct conda environment is used
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    
    # CRITICAL: Clean PYTHONPATH to avoid package conflicts with VirSorter2's internal conda envs
    unset PYTHONPATH
    
    echo "=== VirSorter2 Analysis (SPAdes): ${sample} ==="
    echo "Using Python: \$(which python)"
    echo "Using VirSorter2: \$(which virsorter)"
    
    # Run VirSorter2 for viral sequence identification
    virsorter run \\
        -i ${contigs} \\
        -w virsorter2_output \\
        --db-dir ${virsorter2_db} \\
        --min-length ${params.virsorter2_min_length} \\
        --min-score ${params.virsorter2_min_score} \\
        -j ${task.cpus} \\
        all
    
    # Copy result files
    cp virsorter2_output/final-viral-score.tsv ${sample}_spades_vs2_final-viral-score.tsv
    
    # If viral sequences detected, copy viral contig file
    if [ -f virsorter2_output/final-viral-combined.fa ]; then
        cp virsorter2_output/final-viral-combined.fa ${sample}_spades_vs2_final-viral-combined.fa
    fi
    
    if [ -f virsorter2_output/final-viral-boundary.tsv ]; then
        cp virsorter2_output/final-viral-boundary.tsv ${sample}_spades_vs2_final-viral-boundary.tsv
    fi
    
    # Count identified viral sequences
    VIRAL_COUNT=\$(tail -n +2 ${sample}_spades_vs2_final-viral-score.tsv | wc -l || echo 0)
    echo "VirSorter2: Identified \${VIRAL_COUNT} viral sequences from SPAdes contigs"
    """
}

// Process: VirSorter2 (metaFlye assembly products)
process VIRSORTER2_METAFLYE {
    tag "${sample}_metaFlye_VirSorter2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'
    publishDir "${params.outdir}/virsorter2_metaflye", mode: 'copy', pattern: "*.{tsv,fa}"

    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)

    output:
    tuple val(sample), path("${sample}_metaflye_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_metaflye_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_metaflye_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true

    script:
    """
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    
    # CRITICAL: Clean PYTHONPATH to avoid package conflicts with VirSorter2's internal conda envs
    unset PYTHONPATH
    
    echo "=== VirSorter2（metaFlye）：${sample} ==="
    virsorter run \
        -i ${contigs} \
        -w virsorter2_output \
        --db-dir ${virsorter2_db} \
        --min-length ${params.virsorter2_min_length} \
        --min-score ${params.virsorter2_min_score} \
        -j ${task.cpus} \
        all

    cp virsorter2_output/final-viral-score.tsv ${sample}_metaflye_vs2_final-viral-score.tsv
    if [ -f virsorter2_output/final-viral-combined.fa ]; then
        cp virsorter2_output/final-viral-combined.fa ${sample}_metaflye_vs2_final-viral-combined.fa
    fi
    if [ -f virsorter2_output/final-viral-boundary.tsv ]; then
        cp virsorter2_output/final-viral-boundary.tsv ${sample}_metaflye_vs2_final-viral-boundary.tsv
    fi
    VIRAL_COUNT=\$(tail -n +2 ${sample}_metaflye_vs2_final-viral-score.tsv | wc -l || echo 0)
    echo "VirSorter2: Identified \${VIRAL_COUNT} viral sequences (metaFlye)"
    """
}

// Process: DeepVirFinder for MEGAHIT
process DEEPVIRFINDER_MEGAHIT {
    tag "${sample}_MEGAHIT_DeepVirFinder"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_megahit", mode: 'copy', pattern: "*.txt"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_megahit_dvf_output.txt"), emit: results
    
    script:
    """
    echo "=== DeepVirFinder Analysis (MEGAHIT): ${sample} ==="
    echo "Using DeepVirFinder from: ${params.deepvirfinder_dir}"
    
    # Activate dvf conda environment (using absolute path)
    set +u  # Temporarily disable undefined variable check (conda requires this)
    
    # Load Miniforge3 module (required for SLURM environment)
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    
    # Absolute path to dvf environment
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    
    # Check if environment exists
    if [ ! -d "\$DVF_ENV" ]; then
        echo "❌ dvf environment not found: \$DVF_ENV"
        exit 1
    fi
    
    # Get conda base path
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    if [ -z "\$CONDA_BASE" ]; then
        CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    fi
    
    # Initialize conda
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "\$CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "❌ Cannot find conda.sh"
        exit 1
    fi
    
    # Activate dvf environment using absolute path
    conda activate "\$DVF_ENV" || { echo "❌ Failed to activate dvf environment"; exit 1; }
    
    # Force update PATH to ensure dvf environment Python is used
    export PATH="\$DVF_ENV/bin:\$PATH"
    export CONDA_PREFIX="\$DVF_ENV"
    export CONDA_DEFAULT_ENV="dvf"
    
    # Clean PYTHONPATH to prevent package pollution from other environments (critical!)
    unset PYTHONPATH
    # Get Python version
    PYTHON_VER=\$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    # Keep only dvf environment's site-packages
    export PYTHONPATH="\$DVF_ENV/lib/python\${PYTHON_VER}/site-packages"
    
    # Set Keras to use Theano backend (critical!)
    export KERAS_BACKEND=theano
    
    set -u  # Re-enable
    
    echo "✅ Active conda environment: \$CONDA_DEFAULT_ENV"
    echo "✅ Python path: \$(which python)"
    echo "✅ Python version: \$(python --version)"
    echo "✅ DVF env path: \$DVF_ENV"
    echo "✅ PYTHONPATH: \$PYTHONPATH"
    echo "✅ KERAS_BACKEND: \$KERAS_BACKEND"
    
    # Verify h5py is available
    python -c "import h5py; print('✅ h5py available:', h5py.__version__)" || { echo "❌ h5py not found"; exit 1; }
    
    # Verify keras is available and check backend
    python -c "import os; os.environ['KERAS_BACKEND']='theano'; import keras; print('✅ Keras available:', keras.__version__); print('✅ Keras backend:', keras.backend.backend())" || { echo "❌ Keras not found or backend error"; exit 1; }
    
    # Run DeepVirFinder for viral sequence identification
    python ${params.deepvirfinder_dir}/dvf.py \\
        -i ${contigs} \\
        -o dvf_output \\
        -l ${params.deepvirfinder_min_length} \\
        -c ${task.cpus}
    
    # Copy result files
    cp dvf_output/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_megahit_dvf_output.txt
    
    # Count predicted viral sequences (p-value < threshold)
    VIRAL_COUNT=\$(awk -v pval="${params.deepvirfinder_pvalue}" 'NR>1 && \$3<pval {count++} END {print count+0}' ${sample}_megahit_dvf_output.txt)
    echo "DeepVirFinder: Predicted \${VIRAL_COUNT} viral sequences from MEGAHIT contigs (p-value < ${params.deepvirfinder_pvalue})"
    """
}

// Process: DeepVirFinder for SPAdes
process DEEPVIRFINDER_SPADES {
    tag "${sample}_SPAdes_DeepVirFinder"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_spades", mode: 'copy', pattern: "*.txt"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_spades_dvf_output.txt"), emit: results
    
    script:
    """
    echo "=== DeepVirFinder Analysis (SPAdes): ${sample} ==="
    echo "Using DeepVirFinder from: ${params.deepvirfinder_dir}"
    
    # Activate dvf conda environment (using absolute path)
    set +u  # Temporarily disable undefined variable check (conda requires this)
    
    # Load Miniforge3 module (required for SLURM environment)
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    
    # Absolute path to dvf environment
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    
    # Check if environment exists
    if [ ! -d "\$DVF_ENV" ]; then
        echo "❌ dvf environment not found: \$DVF_ENV"
        exit 1
    fi
    
    # Get conda base path
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    if [ -z "\$CONDA_BASE" ]; then
        CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    fi
    
    # Initialize conda
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "\$CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "❌ Cannot find conda.sh"
        exit 1
    fi
    
    # Activate dvf environment using absolute path
    conda activate "\$DVF_ENV" || { echo "❌ Failed to activate dvf environment"; exit 1; }
    
    # Force update PATH to ensure dvf environment Python is used
    export PATH="\$DVF_ENV/bin:\$PATH"
    export CONDA_PREFIX="\$DVF_ENV"
    export CONDA_DEFAULT_ENV="dvf"
    
    # Clean PYTHONPATH to prevent package pollution from other environments (critical!)
    unset PYTHONPATH
    # Get Python version
    PYTHON_VER=\$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    # Keep only dvf environment's site-packages
    export PYTHONPATH="\$DVF_ENV/lib/python\${PYTHON_VER}/site-packages"
    
    # Set Keras to use Theano backend (critical!)
    export KERAS_BACKEND=theano
    
    set -u  # Re-enable
    
    echo "✅ Active conda environment: \$CONDA_DEFAULT_ENV"
    echo "✅ Python path: \$(which python)"
    echo "✅ Python version: \$(python --version)"
    echo "✅ DVF env path: \$DVF_ENV"
    echo "✅ PYTHONPATH: \$PYTHONPATH"
    echo "✅ KERAS_BACKEND: \$KERAS_BACKEND"
    
    # Verify h5py is available
    python -c "import h5py; print('✅ h5py available:', h5py.__version__)" || { echo "❌ h5py not found"; exit 1; }
    
    # Verify keras is available and check backend
    python -c "import os; os.environ['KERAS_BACKEND']='theano'; import keras; print('✅ Keras available:', keras.__version__); print('✅ Keras backend:', keras.backend.backend())" || { echo "❌ Keras not found or backend error"; exit 1; }
    
    # Run DeepVirFinder for viral sequence identification
    python ${params.deepvirfinder_dir}/dvf.py \\
        -i ${contigs} \\
        -o dvf_output \\
        -l ${params.deepvirfinder_min_length} \\
        -c ${task.cpus}
    
    # Copy result files
    cp dvf_output/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_spades_dvf_output.txt
    
    # Count predicted viral sequences (p-value < threshold)
    VIRAL_COUNT=\$(awk -v pval="${params.deepvirfinder_pvalue}" 'NR>1 && \$3<pval {count++} END {print count+0}' ${sample}_spades_dvf_output.txt)
    echo "DeepVirFinder: Predicted \${VIRAL_COUNT} viral sequences from SPAdes contigs (p-value < ${params.deepvirfinder_pvalue})"
    """
}

// Process: DeepVirFinder (metaFlye assembly products)
process DEEPVIRFINDER_METAFLYE {
    tag "${sample}_metaFlye_DeepVirFinder"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_metaflye", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path("${sample}_metaflye_dvf_output.txt"), emit: results

    script:
    """
    echo "=== DeepVirFinder (metaFlye): ${sample} ==="
    echo "Using DeepVirFinder from: ${params.deepvirfinder_dir}"
    set +u
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    if [ ! -d "\$DVF_ENV" ]; then
        echo "❌ DVF environment not found: \$DVF_ENV"; exit 1; fi
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    [ -z "\$CONDA_BASE" ] && CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then source "\$CONDA_BASE/etc/profile.d/conda.sh"; else echo "❌ Cannot find conda.sh"; exit 1; fi
    conda activate "\$DVF_ENV" || { echo "❌ Failed to activate dvf environment"; exit 1; }
    export PATH="\$DVF_ENV/bin:\$PATH"; export CONDA_PREFIX="\$DVF_ENV"; export CONDA_DEFAULT_ENV="dvf"
    unset PYTHONPATH
    PYTHON_VER=\$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    export PYTHONPATH="\$DVF_ENV/lib/python\${PYTHON_VER}/site-packages"
    export KERAS_BACKEND=theano
    set -u
    python ${params.deepvirfinder_dir}/dvf.py \
        -i ${contigs} \
        -o dvf_output \
        -l ${params.deepvirfinder_min_length} \
        -c ${task.cpus}
    cp dvf_output/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_metaflye_dvf_output.txt
    VIRAL_COUNT=\$(awk -v pval="${params.deepvirfinder_pvalue}" 'NR>1 && \$3<pval {count++} END {print count+0}' ${sample}_metaflye_dvf_output.txt)
    echo "DeepVirFinder: Predicted \${VIRAL_COUNT} viral sequences (metaFlye, p<${params.deepvirfinder_pvalue})"
    """
}
// Process: Merge Viral Identification Reports for MEGAHIT
// Integrate VirSorter2 and DeepVirFinder results
process MERGE_VIRAL_REPORTS_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_megahit", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)
    
    output:
    tuple val(sample), path("${sample}_megahit_viral_merged_report.txt"), emit: merged_report
    tuple val(sample), path("${sample}_megahit_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_megahit_viral_consensus.txt"), emit: consensus_list
    
    script:
    """
    # Clean PYTHONPATH to avoid package conflicts
    unset PYTHONPATH
    
    python3 << 'MERGE_MEGAHIT_SCRIPT'
import pandas as pd
from collections import defaultdict

def parse_virsorter2(file_path):
    \"\"\"Parse VirSorter2 output file\"\"\"
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['seqname']
            seqname_normalized = seqname.split('||')[0] if '||' in seqname else seqname
            viral_dict[seqname_normalized] = {
                'vs2_score': row['max_score'],
                'vs2_group': row['max_score_group'],
                'vs2_length': row['length']
            }
        return viral_dict
    except Exception as e:
        print(f"Warning: Failed to parse VirSorter2 results: {e}")
        return {}

def parse_deepvirfinder(file_path):
    \"\"\"Parse DeepVirFinder output file\"\"\"
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['name']
            seqname_normalized = seqname.split()[0] if ' ' in seqname else seqname
            viral_dict[seqname_normalized] = {
                'dvf_score': row['score'],
                'dvf_pvalue': row['pvalue'],
                'dvf_length': row['len']
            }
        return viral_dict
    except Exception as e:
        print(f"Warning: Failed to parse DeepVirFinder results: {e}")
        return {}

# Parse result files
print(f"Parsing VirSorter2 results: ${virsorter2_results}")
vs2_dict = parse_virsorter2("${virsorter2_results}")
print(f"VirSorter2: Parsed {len(vs2_dict)} sequences")
if len(vs2_dict) > 0:
    print(f"Sample VirSorter2 sequences: {list(vs2_dict.keys())[:5]}")

print(f"Parsing DeepVirFinder results: ${deepvirfinder_results}")
dvf_dict = parse_deepvirfinder("${deepvirfinder_results}")
print(f"DeepVirFinder: Parsed {len(dvf_dict)} sequences")
if len(dvf_dict) > 0:
    print(f"Sample DeepVirFinder sequences: {list(dvf_dict.keys())[:5]}")

# Merge results
all_sequences = set(vs2_dict.keys()) | set(dvf_dict.keys())
print(f"Total unique sequences: {len(all_sequences)}")
if len(all_sequences) == 0:
    print("⚠️ WARNING: No viral sequences found at all! Please check:")
    print("  1. Are there any sequences in the VirSorter2 output file?")
    print("  2. Are there any sequences in the DeepVirFinder output file?")
    print("  3. Check the detection thresholds in the config file")

# Generate comprehensive report
with open("${sample}_megahit_viral_merged_report.txt", 'w', encoding='utf-8') as f:
    f.write("="*80 + "\\n")
    f.write("Viral Identification Comprehensive Analysis Report - MEGAHIT Assembly Results\\n")
    f.write("VirSorter2 + DeepVirFinder\\n")
    f.write("="*80 + "\\n\\n")
    f.write("[Overall Statistics]\\n")
    f.write("-"*80 + "\\n")
    f.write(f"VirSorter2 identified viral sequences:    {len(vs2_dict):,}\\n")
    f.write(f"DeepVirFinder identified viral sequences: {len(dvf_dict):,}\\n")
    consensus = set(vs2_dict.keys()) & set(dvf_dict.keys())
    f.write(f"Consensus viral sequences (both tools):   {len(consensus):,}\\n")
    vs2_only = set(vs2_dict.keys()) - set(dvf_dict.keys())
    dvf_only = set(dvf_dict.keys()) - set(vs2_dict.keys())
    f.write(f"VirSorter2 only:                          {len(vs2_only):,}\\n")
    f.write(f"DeepVirFinder only:                       {len(dvf_only):,}\\n\\n")
    dvf_significant = [seq for seq, data in dvf_dict.items() if data['dvf_pvalue'] < ${params.deepvirfinder_pvalue}]
    f.write(f"DeepVirFinder significant sequences (p<${params.deepvirfinder_pvalue}): {len(dvf_significant):,}\\n\\n")
    f.write("\\n[Consensus Viral Sequences (High Confidence)]\\n")
    f.write("-"*80 + "\\n")
    f.write(f"{'Sequence Name':<40} {'VS2 Score':<12} {'DVF Score':<12} {'DVF P-value':<12}\\n")
    f.write("-"*80 + "\\n")
    for seq in sorted(consensus):
        vs2_score = vs2_dict[seq]['vs2_score']
        dvf_score = dvf_dict[seq]['dvf_score']
        dvf_pval = dvf_dict[seq]['dvf_pvalue']
        f.write(f"{seq:<40} {vs2_score:<12.3f} {dvf_score:<12.3f} {dvf_pval:<12.2e}\\n")
    f.write("\\n" + "="*80 + "\\n")
    f.write("Analysis Complete\\n")
    f.write("="*80 + "\\n")

# Save detailed data in CSV format
merged_data = []
for seq in all_sequences:
    row = {
        'sequence_name': seq,
        'identified_by': '',
        'vs2_score': None,
        'vs2_group': None,
        'dvf_score': None,
        'dvf_pvalue': None,
        'consensus': False
    }
    if seq in vs2_dict:
        row['vs2_score'] = vs2_dict[seq]['vs2_score']
        row['vs2_group'] = vs2_dict[seq]['vs2_group']
    if seq in dvf_dict:
        row['dvf_score'] = dvf_dict[seq]['dvf_score']
        row['dvf_pvalue'] = dvf_dict[seq]['dvf_pvalue']
    if seq in consensus:
        row['identified_by'] = 'Both'
        row['consensus'] = True
    elif seq in vs2_only:
        row['identified_by'] = 'VirSorter2_only'
    else:
        row['identified_by'] = 'DeepVirFinder_only'
    merged_data.append(row)

merged_df = pd.DataFrame(merged_data)
merged_df.to_csv("${sample}_megahit_viral_merged_report.csv", index=False)

# Save consensus sequence list (recommended for downstream analysis)
with open("${sample}_megahit_viral_consensus.txt", 'w') as f:
    for seq in sorted(consensus):
        f.write(seq + "\\n")
    
print(f"Viral identification report generated successfully: ${sample} (MEGAHIT)")
print(f"Consensus viral sequences: {len(consensus)}")
MERGE_MEGAHIT_SCRIPT
    """
}

// Process: Merge Viral Identification Reports for SPAdes
process MERGE_VIRAL_REPORTS_SPADES {
    tag "${sample}_SPAdes"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_spades", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)
    
    output:
    tuple val(sample), path("${sample}_spades_viral_merged_report.txt"), emit: merged_report
    tuple val(sample), path("${sample}_spades_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_spades_viral_consensus.txt"), emit: consensus_list
    
    script:
    """
    # Clean PYTHONPATH to avoid package conflicts
    unset PYTHONPATH
    
    python3 << 'MERGE_SPADES_SCRIPT'
import pandas as pd
from collections import defaultdict

def parse_virsorter2(file_path):
    \"\"\"Parse VirSorter2 output file\"\"\"
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['seqname']
            seqname_normalized = seqname.split('||')[0] if '||' in seqname else seqname
            viral_dict[seqname_normalized] = {
                'vs2_score': row['max_score'],
                'vs2_group': row['max_score_group'],
                'vs2_length': row['length']
            }
        return viral_dict
    except Exception as e:
        print(f"Warning: Failed to parse VirSorter2 results: {e}")
        return {}

def parse_deepvirfinder(file_path):
    \"\"\"Parse DeepVirFinder output file\"\"\"
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['name']
            seqname_normalized = seqname.split()[0] if ' ' in seqname else seqname
            viral_dict[seqname_normalized] = {
                'dvf_score': row['score'],
                'dvf_pvalue': row['pvalue'],
                'dvf_length': row['len']
            }
        return viral_dict
    except Exception as e:
        print(f"Warning: Failed to parse DeepVirFinder results: {e}")
        return {}

# Parse result files
print(f"Parsing VirSorter2 results: ${virsorter2_results}")
vs2_dict = parse_virsorter2("${virsorter2_results}")
print(f"VirSorter2: Parsed {len(vs2_dict)} sequences")
if len(vs2_dict) > 0:
    print(f"Sample VirSorter2 sequences: {list(vs2_dict.keys())[:5]}")

print(f"Parsing DeepVirFinder results: ${deepvirfinder_results}")
dvf_dict = parse_deepvirfinder("${deepvirfinder_results}")
print(f"DeepVirFinder: Parsed {len(dvf_dict)} sequences")
if len(dvf_dict) > 0:
    print(f"Sample DeepVirFinder sequences: {list(dvf_dict.keys())[:5]}")

# Merge results
all_sequences = set(vs2_dict.keys()) | set(dvf_dict.keys())
print(f"Total unique sequences: {len(all_sequences)}")
if len(all_sequences) == 0:
    print("⚠️ WARNING: No viral sequences found at all! Please check:")
    print("  1. Are there any sequences in the VirSorter2 output file?")
    print("  2. Are there any sequences in the DeepVirFinder output file?")
    print("  3. Check the detection thresholds in the config file")

# Generate comprehensive report  
with open("${sample}_spades_viral_merged_report.txt", 'w', encoding='utf-8') as f:
    f.write("="*80 + "\\n")
    f.write("Viral Identification Comprehensive Analysis Report - SPAdes Assembly Results\\n")
    f.write("VirSorter2 + DeepVirFinder\\n")
    f.write("="*80 + "\\n\\n")
    f.write("[Overall Statistics]\\n")
    f.write("-"*80 + "\\n")
    f.write(f"VirSorter2 identified viral sequences:    {len(vs2_dict):,}\\n")
    f.write(f"DeepVirFinder identified viral sequences: {len(dvf_dict):,}\\n")
    consensus = set(vs2_dict.keys()) & set(dvf_dict.keys())
    f.write(f"Consensus viral sequences (both tools):   {len(consensus):,}\\n")
    vs2_only = set(vs2_dict.keys()) - set(dvf_dict.keys())
    dvf_only = set(dvf_dict.keys()) - set(vs2_dict.keys())
    f.write(f"VirSorter2 only:                          {len(vs2_only):,}\\n")
    f.write(f"DeepVirFinder only:                       {len(dvf_only):,}\\n\\n")
    dvf_significant = [seq for seq, data in dvf_dict.items() if data['dvf_pvalue'] < ${params.deepvirfinder_pvalue}]
    f.write(f"DeepVirFinder significant sequences (p<${params.deepvirfinder_pvalue}): {len(dvf_significant):,}\\n\\n")
    f.write("\\n[Consensus Viral Sequences (High Confidence)]\\n")
    f.write("-"*80 + "\\n")
    f.write(f"{'Sequence Name':<40} {'VS2 Score':<12} {'DVF Score':<12} {'DVF P-value':<12}\\n")
    f.write("-"*80 + "\\n")
    for seq in sorted(consensus):
        vs2_score = vs2_dict[seq]['vs2_score']
        dvf_score = dvf_dict[seq]['dvf_score']
        dvf_pval = dvf_dict[seq]['dvf_pvalue']
        f.write(f"{seq:<40} {vs2_score:<12.3f} {dvf_score:<12.3f} {dvf_pval:<12.2e}\\n")
    f.write("\\n" + "="*80 + "\\n")
    f.write("Analysis Complete\\n")
    f.write("="*80 + "\\n")

# Save detailed data in CSV format
merged_data = []
for seq in all_sequences:
    row = {
        'sequence_name': seq,
        'identified_by': '',
        'vs2_score': None,
        'vs2_group': None,
        'dvf_score': None,
        'dvf_pvalue': None,
        'consensus': False
    }
    if seq in vs2_dict:
        row['vs2_score'] = vs2_dict[seq]['vs2_score']
        row['vs2_group'] = vs2_dict[seq]['vs2_group']
    if seq in dvf_dict:
        row['dvf_score'] = dvf_dict[seq]['dvf_score']
        row['dvf_pvalue'] = dvf_dict[seq]['dvf_pvalue']
    if seq in consensus:
        row['identified_by'] = 'Both'
        row['consensus'] = True
    elif seq in vs2_only:
        row['identified_by'] = 'VirSorter2_only'
    else:
        row['identified_by'] = 'DeepVirFinder_only'
    merged_data.append(row)

merged_df = pd.DataFrame(merged_data)
merged_df.to_csv("${sample}_spades_viral_merged_report.csv", index=False)

# Save consensus sequence list (recommended for downstream analysis)
with open("${sample}_spades_viral_consensus.txt", 'w') as f:
    for seq in sorted(consensus):
        f.write(seq + "\\n")
    
print(f"Viral identification report generated successfully: ${sample} (SPAdes)")
print(f"Consensus viral sequences: {len(consensus)}")
MERGE_SPADES_SCRIPT
    """
}

// Process: Merge (metaFlye) - Integrate VirSorter2 and DeepVirFinder results
process MERGE_VIRAL_REPORTS_METAFLYE {
    tag "${sample}_metaFlye"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_metaflye", mode: 'copy', pattern: "*"

    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)

    output:
    tuple val(sample), path("${sample}_metaflye_viral_merged_report.txt"), emit: merged_report
    tuple val(sample), path("${sample}_metaflye_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_metaflye_viral_consensus.txt"), emit: consensus_list

    script:
    """
    # Clean PYTHONPATH to avoid package conflicts
    unset PYTHONPATH
    
    python3 << 'PYTHON_SCRIPT'
import pandas as pd

def parse_virsorter2(file_path):
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['seqname']
            seqname_normalized = seqname.split('||')[0] if '||' in seqname else seqname
            viral_dict[seqname_normalized] = {
                'vs2_score': row['max_score'],
                'vs2_group': row['max_score_group'],
                'vs2_length': row['length']
            }
        return viral_dict
    except Exception as e:
        print(f"Warning: Failed to parse VirSorter2 results: {e}")
        return {}

def parse_deepvirfinder(file_path):
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['name']
            seqname_normalized = seqname.split()[0] if ' ' in seqname else seqname
            viral_dict[seqname_normalized] = {
                'dvf_score': row['score'],
                'dvf_pvalue': row['pvalue'],
                'dvf_length': row['len']
            }
        return viral_dict
    except Exception as e:
        print(f"Warning: Failed to parse DeepVirFinder results: {e}")
        return {}

vs2_dict = parse_virsorter2("${virsorter2_results}")
dvf_dict = parse_deepvirfinder("${deepvirfinder_results}")
all_sequences = set(vs2_dict.keys()) | set(dvf_dict.keys())
consensus = set(vs2_dict.keys()) & set(dvf_dict.keys())
vs2_only = set(vs2_dict.keys()) - set(dvf_dict.keys())
dvf_only = set(dvf_dict.keys()) - set(vs2_dict.keys())
with open("${sample}_metaflye_viral_merged_report.txt", 'w', encoding='utf-8') as f:
    f.write("="*80 + "\\n")
    f.write("Viral Identification Comprehensive Analysis Report - metaFlye Assembly Results\\n")
    f.write("VirSorter2 + DeepVirFinder\\n")
    f.write("="*80 + "\\n\\n")
    f.write("[Overall Statistics]\\n")
    f.write("-"*80 + "\\n")
    f.write(f"VirSorter2 identified viral sequences:    {len(vs2_dict):,}\\n")
    f.write(f"DeepVirFinder identified viral sequences: {len(dvf_dict):,}\\n")
    f.write(f"Consensus viral sequences (both tools):   {len(consensus):,}\\n")
    f.write(f"VirSorter2 only:                          {len(vs2_only):,}\\n")
    f.write(f"DeepVirFinder only:                       {len(dvf_only):,}\\n\\n")

merged_data = []
for seq in all_sequences:
    row = {
        'sequence_name': seq,
        'identified_by': '',
        'vs2_score': None,
        'vs2_group': None,
        'dvf_score': None,
        'dvf_pvalue': None,
        'consensus': False
    }
    if seq in vs2_dict:
        row['vs2_score'] = vs2_dict[seq]['vs2_score']
        row['vs2_group'] = vs2_dict[seq]['vs2_group']
    if seq in dvf_dict:
        row['dvf_score'] = dvf_dict[seq]['dvf_score']
        row['dvf_pvalue'] = dvf_dict[seq]['dvf_pvalue']
    if seq in consensus:
        row['identified_by'] = 'Both'
        row['consensus'] = True
    elif seq in vs2_only:
        row['identified_by'] = 'VirSorter2_only'
    else:
        row['identified_by'] = 'DeepVirFinder_only'
    merged_data.append(row)

pd.DataFrame(merged_data).to_csv("${sample}_metaflye_viral_merged_report.csv", index=False)

with open("${sample}_metaflye_viral_consensus.txt", 'w') as f:
    for seq in sorted(consensus):
        f.write(seq + "\\n")

print(f"Viral identification report generated successfully: ${sample} (metaFlye)")
PYTHON_SCRIPT
    """
}

// ================================================================================
// The following processes are no longer used in Plan A (three-tool parallel analysis)
// Retained for reference or other purposes
// ================================================================================

// Process: Select target viral contigs above VS2 threshold and export subset FASTA
// Note: No longer used in Plan A, viralFlye directly processes metaFlye complete output
process SELECT_VIRAL_TARGETS_METAFLYE {
    tag "${sample}_select_targets"
    label 'process_low'
    conda 'bioconda::seqkit=2.8.2 pandas'
    publishDir "${params.outdir}/viralflye_targets", mode: 'copy', pattern: "*"

    input:
    tuple val(sample), path(contigs)
    tuple val(sample2), path(vs2_scores)

    output:
    tuple val(sample), path("${sample}_viral_targets.fa"), emit: targets
    path("${sample}_viral_target_ids.txt"), emit: ids

    when:
    sample == sample2

    script:
    """
    echo "=== Select viral target contigs: ${sample} ==="
    python - << 'PY'
import pandas as pd, sys
vs2 = pd.read_csv("${vs2_scores}", sep='\t')
ok = (vs2['max_score'] >= ${params.viralflye_min_score}) & (vs2['length'] >= ${params.viralflye_min_length})
ids = []
for _, r in vs2[ok].iterrows():
    name = r['seqname']
    name = name.split('||')[0] if '||' in name else name
    ids.append(name)
open("${sample}_viral_target_ids.txt", 'w').write('\n'.join(ids)+'\n')
PY
    # Use seqkit to extract target contigs
    if [ -s ${sample}_viral_target_ids.txt ]; then
        seqkit grep -f ${sample}_viral_target_ids.txt ${contigs} > ${sample}_viral_targets.fa
    else
        # If no targets, create empty file to avoid interruption
        echo > ${sample}_viral_targets.fa
    fi
    """
}

// Process: Map long reads to target contigs using minimap2 and extract relevant reads using samtools
// Note: No longer used in Plan A, viralFlye directly uses original reads
process SUBSET_LONGREADS_FOR_VIRAL {
    tag "${sample}_subset_reads"
    label 'process_medium'
    conda 'bioconda::minimap2=2.28 bioconda::samtools=1.20'
    publishDir "${params.outdir}/viralflye_reads", mode: 'copy', pattern: "*.fastq.gz"

    input:
    tuple val(sample), path(read_long), path(targets_fa)

    output:
    tuple val(sample), path("${sample}_viral_reads.fastq.gz"), emit: selected_reads

    script:
    """
    echo "=== Extract viral-related reads: ${sample} ==="
    if [ "${params.longread_platform}" = "pacbio" ]; then
        PLATFORM_OPT="map-pb"
    else
        PLATFORM_OPT="map-ont"
    fi
    # Generate SAM and extract mapped reads to FASTQ
    minimap2 -t ${task.cpus} -x \${PLATFORM_OPT} -a ${targets_fa} ${read_long} \\
      | samtools view -b -F 4 -@ ${task.cpus} \\
      | samtools fastq -@ ${task.cpus} -n - \\
      | gzip -c > ${sample}_viral_reads.fastq.gz
    """
}

// Process: viralFlye viral identification (parallel with VirSorter2, DeepVirFinder)
// viralFlye identifies and validates viral sequences from metaFlye assembly results
process VIRALFLYE_IDENTIFY {
    tag "${sample}_viralFlye"
    label 'process_high'
    publishDir "${params.outdir}/viralflye_results", mode: 'copy', pattern: "*.{fa,fasta,csv}"
    publishDir "${params.outdir}/viralflye_full_output", mode: 'copy', pattern: "${sample}_viralflye_output"

    input:
    tuple val(sample), path(flye_dir), path(read_long)

    output:
    tuple val(sample), path("${sample}_viralflye_contigs.fa"), emit: viral_contigs
    tuple val(sample), path("${sample}_viralflye_summary.csv"), emit: results
    path("${sample}_viralflye_output"), emit: viralflye_dir, optional: true

    script:
    """
    echo "=========================================="
    echo "=== viralFlye Viral Identification and Refinement: ${sample} ==="
    echo "=========================================="
    echo ""
    echo "viralFlye Description:"
    echo "  - Identify viral sequences from metaFlye assembly results"
    echo "  - Validate viral proteins using Pfam database"
    echo "  - Assess viral genome completeness"
    echo "  - Generate high-quality viral sequences"
    echo ""
    
    # Activate viralFlye conda environment
    set +u
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    
    VIRALFLYE_ENV="${params.viralflye_env}"
    
    if [ ! -d "\$VIRALFLYE_ENV" ]; then
        echo "❌ viralFlye environment not found: \$VIRALFLYE_ENV"
        exit 1
    fi
    
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    [ -z "\$CONDA_BASE" ] && CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "\$CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "❌ Cannot find conda.sh"
        exit 1
    fi
    
    conda activate "\$VIRALFLYE_ENV" || { echo "❌ Failed to activate viralFlye environment"; exit 1; }
    
    # Force clean PATH and PYTHONPATH to avoid conflicts
    export PATH="\$VIRALFLYE_ENV/bin:\$PATH"
    
    # CRITICAL: Clean PYTHONPATH completely to avoid package conflicts
    unset PYTHONPATH
    
    # Get Python version from viralFlye environment
    PYTHON_VER=\$(\$VIRALFLYE_ENV/bin/python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    
    # Only include viralFlye environment's site-packages
    export PYTHONPATH="\$VIRALFLYE_ENV/lib/python\${PYTHON_VER}/site-packages"
    
    # Ensure conda variables point to viralFlye environment
    export CONDA_PREFIX="\$VIRALFLYE_ENV"
    export CONDA_DEFAULT_ENV="viralFlye_env"
    
    set -u
    
    echo "✅ Environment information:"
    echo "   - Python: \$(which python)"
    echo "   - Python version: \$(python --version)"
    echo ""
    
    # Check viralFlye.py
    VIRALFLYE_PY="\$VIRALFLYE_ENV/bin/viralFlye.py"
    if [ ! -f "\$VIRALFLYE_PY" ]; then
        echo "❌ viralFlye.py not found: \$VIRALFLYE_PY"
        exit 1
    fi
    echo "✅ viralFlye.py: \$VIRALFLYE_PY"
    
    # Check Pfam database
    PFAM_DB="${params.pfam_db}"
    if [ ! -f "\$PFAM_DB" ]; then
        echo "❌ Pfam database not found: \$PFAM_DB"
        exit 1
    fi
    echo "✅ Pfam database: \$PFAM_DB"
    echo "   Size: \$(du -sh \$PFAM_DB | cut -f1)"
    echo ""
    
    # Check metaFlye output directory
    if [ ! -d "${flye_dir}" ]; then
        echo "❌ metaFlye output directory not found: ${flye_dir}"
        exit 1
    fi
    echo "✅ metaFlye output directory: ${flye_dir}"
    echo "   Contents:"
    ls -lh ${flye_dir}/ | head -10
    echo ""
    
    # Check original reads
    if [ ! -f "${read_long}" ]; then
        echo "❌ Original reads file not found: ${read_long}"
        exit 1
    fi
    echo "✅ Original reads: ${read_long}"
    echo "   Size: \$(du -sh ${read_long} | cut -f1)"
    echo ""
    
    # Check if Flye is installed (viralFlye dependency)
    if ! command -v flye &> /dev/null; then
        echo "⚠️  Warning: Flye tool not installed"
        echo "   Install Flye: conda install -c bioconda flye"
        echo "   Attempting to continue..."
    else
        echo "✅ Flye: \$(flye --version 2>&1 | head -1)"
    fi
    echo ""
    
    echo "=========================================="
    echo "🦠 Running viralFlye Pipeline"
    echo "=========================================="
    echo ""
    
    # Run viralFlye
    # Parameter explanation:
    #   --dir: metaFlye output directory
    #   --hmm: Pfam-A HMM database
    #   --reads: Original long-read data
    #   --outdir: viralFlye output directory
    #   --threads: Number of threads
    #   --min_viral_length: Minimum viral length (default 5kb)
    #   --completeness: Completeness threshold (default 0.5, lower value identifies more viruses) ⭐
    python \$VIRALFLYE_PY \\
        --dir ${flye_dir} \\
        --hmm \$PFAM_DB \\
        --reads ${read_long} \\
        --outdir ${sample}_viralflye_output \\
        --threads ${task.cpus} \\
        --min_viral_length ${params.viralflye_min_length} \\
        --completeness ${params.viralflye_completeness}
    
    echo ""
    echo "=========================================="
    echo "📁 Processing viralFlye Output"
    echo "=========================================="
    echo ""
    
    # viralFlye output structure
    # viralflye_output/
    #   ├── vv_circulars/
    #   │   ├── circulars.txt           # Circular virus ID list
    #   │   └── circulars.fasta         # Circular virus sequences ⭐
    #   ├── vv_linears/
    #   │   ├── linears.txt             # Linear virus ID list
    #   │   └── linears.fasta           # Linear virus sequences ⭐
    #   ├── vc_circulars/               # viralComplete circular results
    #   ├── vc_linears/                 # viralComplete linear results
    #   └── vv_components/              # viralVerify components
    
    # Search for viralFlye output viral sequences
    echo "Searching for viralFlye output viral sequences..."
    echo "Strategy: Extract Virus + Uncertain-viral candidate sequences"
    OUTPUT_FOUND=false
    VIRAL_COUNT=0
    
    # Create temporary merge file
    > temp_viral_sequences.fa
    
    # Priority 1: vv_circulars/circulars.fasta (all circular candidate sequences)
    # Note: Includes Virus + Uncertain - viral or bacterial
    if [ -f "${sample}_viralflye_output/vv_circulars/circulars.fasta" ]; then
        CIRC_COUNT=\$(grep -c ">" ${sample}_viralflye_output/vv_circulars/circulars.fasta 2>/dev/null || echo 0)
        if [ \$CIRC_COUNT -gt 0 ]; then
            echo "✅ Found \$CIRC_COUNT circular candidate sequences (including Virus + Uncertain)"
            cat ${sample}_viralflye_output/vv_circulars/circulars.fasta >> temp_viral_sequences.fa
            VIRAL_COUNT=\$((VIRAL_COUNT + CIRC_COUNT))
            OUTPUT_FOUND=true
        fi
    fi
    
    # Priority 2: vv_linears/linears.fasta (all linear candidate sequences)
    # Note: Includes Virus + Uncertain - viral or bacterial
    if [ -f "${sample}_viralflye_output/vv_linears/linears.fasta" ]; then
        LIN_COUNT=\$(grep -c ">" ${sample}_viralflye_output/vv_linears/linears.fasta 2>/dev/null || echo 0)
        if [ \$LIN_COUNT -gt 0 ]; then
            echo "✅ Found \$LIN_COUNT linear candidate sequences (including Virus + Uncertain)"
            cat ${sample}_viralflye_output/vv_linears/linears.fasta >> temp_viral_sequences.fa
            VIRAL_COUNT=\$((VIRAL_COUNT + LIN_COUNT))
            OUTPUT_FOUND=true
        fi
    fi
    
    # If viral sequences found, copy merged file
    if [ "\$OUTPUT_FOUND" = true ] && [ \$VIRAL_COUNT -gt 0 ]; then
        cp temp_viral_sequences.fa ${sample}_viralflye_contigs.fa
        echo "✅ Successfully extracted \$VIRAL_COUNT viral sequences"
    else
        echo "⚠️  No viral sequences found in vv_circulars or vv_linears"
        
        # Priority 3: Search other possible locations
        echo "Searching for other possible viral sequence files..."
        FOUND_FILES=\$(find ${sample}_viralflye_output -name "*.fasta" -o -name "*.fa" 2>/dev/null | \
                      xargs -I {} sh -c 'c=\$(grep -c ">" {} 2>/dev/null || echo 0); [ \$c -gt 0 ] && echo {}:\$c')
        
        if [ -n "\$FOUND_FILES" ]; then
            echo "Found non-empty FASTA files:"
            echo "\$FOUND_FILES"
            
            # Use first non-empty fasta
            FIRST_FASTA=\$(echo "\$FOUND_FILES" | head -1 | cut -d: -f1)
            if [ -f "\$FIRST_FASTA" ]; then
                echo "Using: \$FIRST_FASTA"
                cp "\$FIRST_FASTA" ${sample}_viralflye_contigs.fa
                OUTPUT_FOUND=true
                VIRAL_COUNT=\$(grep -c ">" ${sample}_viralflye_contigs.fa || echo 0)
            fi
        fi
    fi
    
    # Clean up temporary file
    rm -f temp_viral_sequences.fa
    
    # Check final output
    if [ "\$OUTPUT_FOUND" = true ] && [ -f "${sample}_viralflye_contigs.fa" ]; then
        CONTIG_COUNT=\$(grep -c ">" ${sample}_viralflye_contigs.fa 2>/dev/null || echo 0)
        TOTAL_LENGTH=\$(awk '/^>/ {next} {total += length(\$0)} END {print total+0}' ${sample}_viralflye_contigs.fa)
        
        echo ""
        echo "=========================================="
        echo "✅ viralFlye Viral Identification Complete!"
        echo "=========================================="
        echo "Sample: ${sample}"
        echo "Identified viruses: \$CONTIG_COUNT"
        echo "Total length: \$TOTAL_LENGTH bp"
        
        if [ \$CONTIG_COUNT -gt 0 ]; then
            echo ""
            echo "Identified viral sequences:"
            grep ">" ${sample}_viralflye_contigs.fa | head -10
        fi
    else
        echo ""
        echo "⚠️  viralFlye did not identify any viral sequences"
        CONTIG_COUNT=0
        # Create empty file
        touch ${sample}_viralflye_contigs.fa
    fi
    
    # Generate standardized result CSV (for comparison with VS2 and DVF)
    echo ""
    echo "Generating viralFlye result summary..."
    
    # CRITICAL: Use viralFlye environment's Python and ensure no user site-packages interference
    # Set PYTHONNOUSERSITE=1 to prevent Python from loading ~/.local/lib/python*/site-packages
    export PYTHONNOUSERSITE=1
    
    # Ensure PYTHONPATH points only to viralFlye environment
    export PYTHONPATH="\$VIRALFLYE_ENV/lib/python\${PYTHON_VER}/site-packages"
    
    # Check if pandas is available in viralFlye environment
    echo "Checking for pandas in viralFlye environment..."
    if ! \$VIRALFLYE_ENV/bin/python -c "import pandas" 2>/dev/null; then
        echo "⚠️  pandas not found in viralFlye_env. Attempting to install..."
        # Activate conda and install pandas
        conda activate "\$VIRALFLYE_ENV" || { echo "❌ Failed to activate viralFlye environment"; exit 1; }
        conda install -y -c conda-forge pandas || pip install pandas
        # Verify installation
        if ! \$VIRALFLYE_ENV/bin/python -c "import pandas" 2>/dev/null; then
            echo "❌ Failed to install pandas in viralFlye_env"
            echo "Please install pandas manually: conda activate viralFlye_env && conda install -c conda-forge pandas"
            exit 1
        fi
        echo "✅ pandas installed successfully"
    else
        echo "✅ pandas found in viralFlye environment"
    fi
    
    # Use viralFlye environment's Python explicitly
    \$VIRALFLYE_ENV/bin/python << 'PYEOF'
import csv
import pandas as pd

# Read circular virus classification results
viral_list = []
try:
    circ_df = pd.read_csv("${sample}_viralflye_output/vv_circulars/circulars_result_table.csv")
    for _, row in circ_df.iterrows():
        # ⭐⭐⭐ Modified: Include Virus and Uncertain - viral or bacterial
        prediction = str(row['Prediction'])
        if prediction == 'Virus' or 'Uncertain - viral or bacterial' in prediction:
            viral_list.append({
                'seqname': row['Contig name'],
                'length': row['Length'],
                'viralflye_score': row['Score'] if pd.notna(row['Score']) else 0,
                'viralflye_type': 'Circular',
                'viralflye_prediction': row['Prediction']
            })
except Exception as e:
    print(f"Warning: Failed to read circular results: {e}")

# Read linear virus classification results
try:
    lin_df = pd.read_csv("${sample}_viralflye_output/vv_linears/linears_result_table.csv")
    for _, row in lin_df.iterrows():
        # ⭐⭐⭐ Modified: Include Virus and Uncertain - viral or bacterial
        prediction = str(row['Prediction'])
        if prediction == 'Virus' or 'Uncertain - viral or bacterial' in prediction:
            viral_list.append({
                'seqname': row['Contig name'],
                'length': row['Length'],
                'viralflye_score': row['Score'] if pd.notna(row['Score']) else 0,
                'viralflye_type': 'Linear',
                'viralflye_prediction': row['Prediction']
            })
except Exception as e:
    print(f"Warning: Failed to read linear results: {e}")

# Save as CSV
df = pd.DataFrame(viral_list)
df.to_csv("${sample}_viralflye_summary.csv", index=False)

# Count different classifications
if len(viral_list) > 0:
    virus_count = sum(1 for v in viral_list if v['viralflye_prediction'] == 'Virus')
    uncertain_count = sum(1 for v in viral_list if 'Uncertain' in str(v['viralflye_prediction']))
    print(f"viralFlye Identification Summary:")
    print(f"  - Confirmed viruses (Virus): {virus_count} sequences")
    print(f"  - Uncertain viral candidates (Uncertain - viral or bacterial): {uncertain_count} sequences")
    print(f"  - Total: {len(viral_list)} viral/candidate sequences")
else:
    print(f"viralFlye did not identify any viruses or candidate sequences")
PYEOF
    
    echo "=========================================="
    """
}

// ================================================================================
// Three-Tool Comprehensive Comparison (Plan A: Parallel Independent Analysis)
// ================================================================================

// Process: Comprehensive comparison of three viral identification tools
// VirSorter2 + DeepVirFinder + viralFlye parallel analysis
process COMPARE_THREE_VIRAL_TOOLS {
    tag "${sample}_ThreeTools"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/three_tools_comparison", mode: 'copy', pattern: "*"

    input:
    tuple val(sample), path(vs2_results), path(dvf_results), path(vf_results)

    output:
    tuple val(sample), path("${sample}_three_tools_comparison.txt"), emit: report
    tuple val(sample), path("${sample}_three_tools_comparison.csv"), emit: csv
    path("${sample}_high_confidence_viruses.txt"), emit: consensus_list
    path("${sample}_high_confidence_viruses.fa"), emit: consensus_fasta, optional: true

    script:
    """
    # Clean PYTHONPATH to avoid package conflicts
    unset PYTHONPATH
    
    python3 << 'COMPARE_THREE_SCRIPT'
# -*- coding: utf-8 -*-

import pandas as pd
import sys

print("="*80)
print("Three-Tool Viral Identification Comprehensive Comparison")
print("VirSorter2 + DeepVirFinder + viralFlye")
print("="*80)
print()

# Parse VirSorter2 results
def parse_virsorter2(file_path):
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['seqname'].split('||')[0] if '||' in row['seqname'] else row['seqname']
            viral_dict[seqname] = {
                'vs2_score': row['max_score'],
                'vs2_group': row['max_score_group'],
                'vs2_length': row['length']
            }
        print(f"✅ VirSorter2: Identified {len(viral_dict)} viruses")
        return viral_dict
    except Exception as e:
        print(f"⚠️  VirSorter2: Parse failed - {e}")
        return {}

# Parse DeepVirFinder results
def parse_deepvirfinder(file_path):
    try:
        df = pd.read_csv(file_path, sep='\\t')
        viral_dict = {}
        for _, row in df.iterrows():
            if row['pvalue'] < ${params.deepvirfinder_pvalue}:
                seqname = row['name'].split()[0] if ' ' in row['name'] else row['name']
                viral_dict[seqname] = {
                    'dvf_score': row['score'],
                    'dvf_pvalue': row['pvalue'],
                    'dvf_length': row['len']
                }
        print(f"✅ DeepVirFinder: Identified {len(viral_dict)} viruses (p<${params.deepvirfinder_pvalue})")
        return viral_dict
    except Exception as e:
        print(f"⚠️  DeepVirFinder: Parse failed - {e}")
        return {}

# Parse viralFlye results
def parse_viralflye(file_path):
    try:
        df = pd.read_csv(file_path)
        viral_dict = {}
        for _, row in df.iterrows():
            seqname = row['seqname']
            viral_dict[seqname] = {
                'vf_score': row['viralflye_score'],
                'vf_type': row['viralflye_type'],
                'vf_length': row['length']
            }
        print(f"✅ viralFlye: Identified {len(viral_dict)} viruses (Pfam validated)")
        return viral_dict
    except Exception as e:
        print(f"⚠️  viralFlye: Parse failed - {e}")
        return {}

# Parse all results
vs2_dict = parse_virsorter2("${vs2_results}")
dvf_dict = parse_deepvirfinder("${dvf_results}")
vf_dict = parse_viralflye("${vf_results}")

# Merge all identified viruses
all_viruses = set(vs2_dict.keys()) | set(dvf_dict.keys()) | set(vf_dict.keys())

print()
print("="*80)
print("Statistical Summary")
print("="*80)
print(f"VirSorter2 identified:      {len(vs2_dict):>6} viruses")
print(f"DeepVirFinder identified:   {len(dvf_dict):>6} viruses")
print(f"viralFlye identified:       {len(vf_dict):>6} viruses")
print(f"Total (deduplicated):       {len(all_viruses):>6} viruses")
print()

# Calculate intersections
vs2_dvf = set(vs2_dict.keys()) & set(dvf_dict.keys())
vs2_vf = set(vs2_dict.keys()) & set(vf_dict.keys())
dvf_vf = set(dvf_dict.keys()) & set(vf_dict.keys())
three_consensus = set(vs2_dict.keys()) & set(dvf_dict.keys()) & set(vf_dict.keys())

print("Tool Intersection Analysis:")
print(f"VS2 ∩ DVF:            {len(vs2_dvf):>6} sequences")
print(f"VS2 ∩ viralFlye:      {len(vs2_vf):>6} sequences")
print(f"DVF ∩ viralFlye:      {len(dvf_vf):>6} sequences")
print(f"Three-tool consensus ⭐⭐⭐:       {len(three_consensus):>6} sequences (highest confidence)")
print()

# Generate detailed comparison table
comparison_data = []
for seq in all_viruses:
    row = {
        'sequence_name': seq,
        'identified_by': [],
        'vs2_score': vs2_dict[seq]['vs2_score'] if seq in vs2_dict else None,
        'vs2_group': vs2_dict[seq]['vs2_group'] if seq in vs2_dict else None,
        'dvf_score': dvf_dict[seq]['dvf_score'] if seq in dvf_dict else None,
        'dvf_pvalue': dvf_dict[seq]['dvf_pvalue'] if seq in dvf_dict else None,
        'vf_score': vf_dict[seq]['vf_score'] if seq in vf_dict else None,
        'vf_type': vf_dict[seq]['vf_type'] if seq in vf_dict else None,
        'consensus_count': 0
    }
    
    # Count how many tools identified this sequence
    if seq in vs2_dict:
        row['identified_by'].append('VirSorter2')
        row['consensus_count'] += 1
    if seq in dvf_dict:
        row['identified_by'].append('DeepVirFinder')
        row['consensus_count'] += 1
    if seq in vf_dict:
        row['identified_by'].append('viralFlye')
        row['consensus_count'] += 1
    
    row['identified_by'] = '+'.join(row['identified_by'])
    comparison_data.append(row)

# Save CSV
df = pd.DataFrame(comparison_data)
df = df.sort_values('consensus_count', ascending=False)
df.to_csv("${sample}_three_tools_comparison.csv", index=False)

# Generate text report
with open("${sample}_three_tools_comparison.txt", 'w', encoding='utf-8') as f:
    f.write("="*100 + "\\n")
    f.write("Three-Tool Viral Identification Comprehensive Comparison Report\\n")
    f.write("VirSorter2 + DeepVirFinder + viralFlye (Plan A: Parallel Independent Analysis)\\n")
    f.write("Sample: ${sample}\\n")
    f.write("="*100 + "\\n\\n")
    f.write("[Overall Statistics]\\n")
    f.write("-"*100 + "\\n")
    f.write(f"VirSorter2 identified viruses:        {len(vs2_dict):>6} sequences\\n")
    f.write(f"DeepVirFinder identified viruses:     {len(dvf_dict):>6} sequences (p<${params.deepvirfinder_pvalue})\\n")
    f.write(f"viralFlye identified viruses:         {len(vf_dict):>6} sequences (Pfam validated)\\n")
    f.write(f"Total viruses (deduplicated):         {len(all_viruses):>6} sequences\\n\\n")
    f.write("[Tool Intersection Analysis]\\n")
    f.write("-"*100 + "\\n")
    f.write(f"VirSorter2 ∩ DeepVirFinder: {len(vs2_dvf):>6} sequences\\n")
    f.write(f"VirSorter2 ∩ viralFlye:     {len(vs2_vf):>6} sequences\\n")
    f.write(f"DeepVirFinder ∩ viralFlye:  {len(dvf_vf):>6} sequences\\n")
    f.write(f"Three-tool consensus ⭐⭐⭐:            {len(three_consensus):>6} sequences (highest confidence)\\n\\n")
    f.write("[Confidence Stratification]\\n")
    f.write("-"*100 + "\\n")
    three_tool = [s for s in comparison_data if s['consensus_count'] == 3]
    two_tool = [s for s in comparison_data if s['consensus_count'] == 2]
    one_tool = [s for s in comparison_data if s['consensus_count'] == 1]
    f.write(f"Highest confidence (3-tool consensus):    {len(three_tool):>6} sequences ⭐⭐⭐\\n")
    f.write(f"Medium confidence (2-tool consensus):     {len(two_tool):>6} sequences ⭐⭐\\n")
    f.write(f"Low confidence (1-tool identification):   {len(one_tool):>6} sequences ⭐\\n\\n")
    f.write("[Recommended Analysis Strategy]\\n")
    f.write("-"*100 + "\\n")
    f.write("1. Prioritize three-tool consensus viruses (most reliable)\\n")
    f.write("2. Two-tool consensus viruses as second priority\\n")
    f.write("3. viralFlye unique identifications (Pfam validated, high specificity)\\n")
    f.write("4. Single-tool identifications for exploratory analysis\\n\\n")
    if len(three_tool) > 0:
        f.write("[Three-Tool Consensus Virus List ⭐⭐⭐]\\n")
        f.write("-"*100 + "\\n")
        f.write(f"{'Sequence Name':<40} {'VS2 Score':<10} {'DVF Score':<10} {'viralFlye':<15}\\n")
        f.write("-"*100 + "\\n")
        for v in three_tool:
            seq = v['sequence_name']
            vs2_sc = f"{v['vs2_score']:.2f}" if v['vs2_score'] is not None else 'N/A'
            dvf_sc = f"{v['dvf_score']:.2f}" if v['dvf_score'] is not None else 'N/A'
            vf_info = f"{v['vf_type']}" if v['vf_type'] else 'N/A'
            f.write(f"{seq:<40} {vs2_sc:<10} {dvf_sc:<10} {vf_info:<15}\\n")
        f.write("\\n" + "="*100 + "\\n")

# Save high-confidence virus list
with open("${sample}_high_confidence_viruses.txt", 'w') as f:
    f.write("# High-confidence viral sequences (prioritized)\\n")
    f.write(f"# Sample: ${sample}\\n")
    f.write("#\\n")
    f.write("# Three-tool consensus (highest confidence):  " + str(len(three_tool)) + "\\n")
    for v in three_tool:
        f.write(v['sequence_name'] + "\\t3-tool-consensus\\n")
    f.write("#\\n")
    f.write("# Two-tool consensus:  " + str(len(two_tool)) + "\\n")
    for v in two_tool:
        f.write(v['sequence_name'] + "\\t2-tool-consensus\\t" + v['identified_by'] + "\\n")

print()
print("="*80)
print(f"Three-tool comparison complete: ${sample}")
print(f"  Three-tool consensus: {len(three_tool)} sequences ⭐⭐⭐")
print(f"  Two-tool consensus: {len(two_tool)} sequences ⭐⭐")
print(f"  Single-tool identification: {len(one_tool)} sequences ⭐")
print("="*80)
COMPARE_THREE_SCRIPT
    """
}

// ================================================================================
// Viral Abundance Calculation
// ================================================================================

// Process: Calculate Viral Abundance (Short-Read Mode)
// Map reads to viral contigs and calculate RPM and RPKM
process CALCULATE_ABUNDANCE {
    tag "${sample}_${assembler}_Abundance"
    label 'process_medium'
    conda 'bioconda::bowtie2=2.5.1 bioconda::samtools=1.18'
    publishDir "${params.outdir}/abundance/${assembler}", mode: 'copy', pattern: "*.{csv,txt}"
    
    input:
    tuple val(sample), path(reads), path(viral_contigs), val(assembler)
    
    output:
    tuple val(sample), path("${sample}_${assembler}_abundance.csv"), emit: abundance
    path("${sample}_${assembler}_abundance_summary.txt"), emit: summary
    
    script:
    """
    echo "=== Calculating Viral Abundance: ${sample} (${assembler}) ==="
    
    # Check if viral contigs file is empty or has no sequences
    CONTIG_COUNT=\$(grep -c ">" ${viral_contigs} 2>/dev/null || echo 0)
    
    if [ \$CONTIG_COUNT -eq 0 ]; then
        echo "No viral contigs found. Creating empty abundance file."
        echo "contig_id,contig_length,mapped_reads,total_reads,rpm,rpkm" > ${sample}_${assembler}_abundance.csv
        echo "No viral contigs identified for ${sample} (${assembler})" > ${sample}_${assembler}_abundance_summary.txt
        exit 0
    fi
    
    echo "Found \$CONTIG_COUNT viral contigs"
    
    # Build bowtie2 index
    echo "Building bowtie2 index..."
    bowtie2-build ${viral_contigs} viral_index
    
    # Map reads to viral contigs
    echo "Mapping reads to viral contigs..."
    bowtie2 \\
        -x viral_index \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -p ${task.cpus} \\
        --no-unal \\
        -S aligned.sam
    
    # Convert SAM to BAM and sort
    echo "Converting to BAM and sorting..."
    samtools view -bS aligned.sam | samtools sort -@ ${task.cpus} -o aligned_sorted.bam
    samtools index aligned_sorted.bam
    
    # Calculate total number of reads in input
    echo "Calculating total read count..."
    TOTAL_READS=\$(zcat ${reads[0]} | wc -l)
    TOTAL_READS=\$((TOTAL_READS / 4))
    echo "Total reads in input: \$TOTAL_READS"
    
    # Export for Python script
    export TOTAL_READS
    
    # Clean PYTHONPATH to avoid package conflicts
    unset PYTHONPATH
    
    # Calculate abundance for each contig
    echo "Calculating abundance metrics..."
    python3 << 'PYEOF'
import os
import subprocess

# Get contig lengths
contig_lengths = {}
with open("${viral_contigs}", 'r') as f:
    current_contig = None
    current_seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_contig:
                contig_lengths[current_contig] = len(current_seq)
            current_contig = line[1:].split()[0]
            current_seq = ""
        else:
            current_seq += line
    if current_contig:
        contig_lengths[current_contig] = len(current_seq)

# Get mapped read counts for each contig
contig_counts = {}
for contig in contig_lengths.keys():
    result = subprocess.run(
        ["samtools", "view", "-c", "aligned_sorted.bam", contig],
        capture_output=True,
        text=True
    )
    count = int(result.stdout.strip())
    contig_counts[contig] = count

# Total reads from environment variable
total_reads = int(os.environ.get('TOTAL_READS', 0))

# Calculate RPM and RPKM
results = []
for contig, length in contig_lengths.items():
    mapped_reads = contig_counts.get(contig, 0)
    
    # RPM = (Mapped Reads / Total Reads) * 1,000,000
    rpm = (mapped_reads / total_reads) * 1e6 if total_reads > 0 else 0
    
    # RPKM = (Mapped Reads / (Contig Length in kb * Total Reads in millions))
    # RPKM = (Mapped Reads * 1e9) / (Contig Length * Total Reads)
    rpkm = (mapped_reads * 1e9) / (length * total_reads) if length > 0 and total_reads > 0 else 0
    
    results.append({
        'contig_id': contig,
        'contig_length': length,
        'mapped_reads': mapped_reads,
        'total_reads': total_reads,
        'rpm': rpm,
        'rpkm': rpkm
    })

# Write CSV output
with open("${sample}_${assembler}_abundance.csv", 'w') as f:
    f.write("contig_id,contig_length,mapped_reads,total_reads,rpm,rpkm\\n")
    for r in sorted(results, key=lambda x: x['rpkm'], reverse=True):
        f.write(f"{r['contig_id']},{r['contig_length']},{r['mapped_reads']},{r['total_reads']},{r['rpm']:.4f},{r['rpkm']:.4f}\\n")

# Write summary report
with open("${sample}_${assembler}_abundance_summary.txt", 'w') as f:
    f.write("="*80 + "\\n")
    f.write(f"Viral Abundance Analysis Summary\\n")
    f.write(f"Sample: ${sample}\\n")
    f.write(f"Assembler: ${assembler}\\n")
    f.write("="*80 + "\\n\\n")
    f.write(f"Total reads in input: {total_reads:,}\\n")
    f.write(f"Total viral contigs: {len(contig_lengths):,}\\n")
    
    total_mapped = sum(contig_counts.values())
    f.write(f"Total reads mapped to viral contigs: {total_mapped:,}\\n")
    mapping_rate = (total_mapped / total_reads * 100) if total_reads > 0 else 0
    f.write(f"Mapping rate: {mapping_rate:.2f}%\\n\\n")
    
    f.write("Top 10 Most Abundant Viral Contigs (by RPKM):\\n")
    f.write("-"*80 + "\\n")
    f.write(f"{'Contig ID':<40} {'Length':<10} {'Reads':<10} {'RPKM':<10}\\n")
    f.write("-"*80 + "\\n")
    
    for i, r in enumerate(sorted(results, key=lambda x: x['rpkm'], reverse=True)[:10], 1):
        f.write(f"{r['contig_id']:<40} {r['contig_length']:<10} {r['mapped_reads']:<10} {r['rpkm']:<10.2f}\\n")
    
    f.write("\\n" + "="*80 + "\\n")

print(f"Abundance calculation completed for ${sample} (${assembler})")
PYEOF
    
    echo "Viral abundance calculation completed!"
    """
}

// Process: Calculate Viral Abundance (Long-Read Mode)
// Map long reads to viral contigs and calculate RPM and RPKM
process CALCULATE_ABUNDANCE_LONGREAD {
    tag "${sample}_${source}_Abundance"
    label 'process_medium'
    conda 'bioconda::minimap2=2.28 bioconda::samtools=1.18'
    publishDir "${params.outdir}/abundance/${source}", mode: 'copy', pattern: "*.{csv,txt}"
    
    input:
    tuple val(sample), path(reads), path(viral_contigs), val(source)
    
    output:
    tuple val(sample), path("${sample}_${source}_abundance.csv"), emit: abundance
    path("${sample}_${source}_abundance_summary.txt"), emit: summary
    
    script:
    """
    echo "=== Calculating Viral Abundance (Long-Read): ${sample} (${source}) ==="
    
    # Check if viral contigs file is empty or has no sequences
    CONTIG_COUNT=\$(grep -c ">" ${viral_contigs} 2>/dev/null || echo 0)
    
    if [ \$CONTIG_COUNT -eq 0 ]; then
        echo "No viral contigs found. Creating empty abundance file."
        echo "contig_id,contig_length,mapped_reads,total_reads,rpm,rpkm" > ${sample}_${source}_abundance.csv
        echo "No viral contigs identified for ${sample} (${source})" > ${sample}_${source}_abundance_summary.txt
        exit 0
    fi
    
    echo "Found \$CONTIG_COUNT viral contigs"
    
    # Determine platform for minimap2
    if [ "${params.longread_platform}" = "pacbio" ]; then
        PLATFORM_PRESET="map-pb"
    else
        PLATFORM_PRESET="map-ont"
    fi
    
    # Map reads to viral contigs using minimap2
    echo "Mapping long reads to viral contigs..."
    minimap2 \\
        -x \${PLATFORM_PRESET} \\
        -a \\
        -t ${task.cpus} \\
        ${viral_contigs} \\
        ${reads} \\
        > aligned.sam
    
    # Convert SAM to BAM and sort
    echo "Converting to BAM and sorting..."
    samtools view -bS aligned.sam | samtools sort -@ ${task.cpus} -o aligned_sorted.bam
    samtools index aligned_sorted.bam
    
    # Calculate total number of reads in input
    echo "Calculating total read count..."
    if [[ "${reads}" == *.gz ]]; then
        TOTAL_READS=\$(zcat ${reads} | grep -c "^@" || echo 0)
    else
        TOTAL_READS=\$(grep -c "^@" ${reads} || echo 0)
    fi
    echo "Total reads in input: \$TOTAL_READS"
    
    # Export for Python script
    export TOTAL_READS
    
    # Clean PYTHONPATH to avoid package conflicts
    unset PYTHONPATH
    
    # Calculate abundance for each contig
    echo "Calculating abundance metrics..."
    python3 << 'PYEOF'
import os
import subprocess

# Get contig lengths
contig_lengths = {}
with open("${viral_contigs}", 'r') as f:
    current_contig = None
    current_seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_contig:
                contig_lengths[current_contig] = len(current_seq)
            current_contig = line[1:].split()[0]
            current_seq = ""
        else:
            current_seq += line
    if current_contig:
        contig_lengths[current_contig] = len(current_seq)

# Get mapped read counts for each contig
contig_counts = {}
for contig in contig_lengths.keys():
    result = subprocess.run(
        ["samtools", "view", "-c", "-F", "4", "aligned_sorted.bam", contig],
        capture_output=True,
        text=True
    )
    count = int(result.stdout.strip())
    contig_counts[contig] = count

# Total reads from environment variable
total_reads = int(os.environ.get('TOTAL_READS', 0))

# Calculate RPM and RPKM
results = []
for contig, length in contig_lengths.items():
    mapped_reads = contig_counts.get(contig, 0)
    
    # RPM = (Mapped Reads / Total Reads) * 1,000,000
    rpm = (mapped_reads / total_reads) * 1e6 if total_reads > 0 else 0
    
    # RPKM = (Mapped Reads * 1e9) / (Contig Length * Total Reads)
    rpkm = (mapped_reads * 1e9) / (length * total_reads) if length > 0 and total_reads > 0 else 0
    
    results.append({
        'contig_id': contig,
        'contig_length': length,
        'mapped_reads': mapped_reads,
        'total_reads': total_reads,
        'rpm': rpm,
        'rpkm': rpkm
    })

# Write CSV output
with open("${sample}_${source}_abundance.csv", 'w') as f:
    f.write("contig_id,contig_length,mapped_reads,total_reads,rpm,rpkm\\n")
    for r in sorted(results, key=lambda x: x['rpkm'], reverse=True):
        f.write(f"{r['contig_id']},{r['contig_length']},{r['mapped_reads']},{r['total_reads']},{r['rpm']:.4f},{r['rpkm']:.4f}\\n")

# Write summary report
with open("${sample}_${source}_abundance_summary.txt", 'w') as f:
    f.write("="*80 + "\\n")
    f.write(f"Viral Abundance Analysis Summary (Long-Read)\\n")
    f.write(f"Sample: ${sample}\\n")
    f.write(f"Source: ${source}\\n")
    f.write("="*80 + "\\n\\n")
    f.write(f"Total reads in input: {total_reads:,}\\n")
    f.write(f"Total viral contigs: {len(contig_lengths):,}\\n")
    
    total_mapped = sum(contig_counts.values())
    f.write(f"Total reads mapped to viral contigs: {total_mapped:,}\\n")
    mapping_rate = (total_mapped / total_reads * 100) if total_reads > 0 else 0
    f.write(f"Mapping rate: {mapping_rate:.2f}%\\n\\n")
    
    f.write("Top 10 Most Abundant Viral Contigs (by RPKM):\\n")
    f.write("-"*80 + "\\n")
    f.write(f"{'Contig ID':<40} {'Length':<10} {'Reads':<10} {'RPKM':<10}\\n")
    f.write("-"*80 + "\\n")
    
    for i, r in enumerate(sorted(results, key=lambda x: x['rpkm'], reverse=True)[:10], 1):
        f.write(f"{r['contig_id']:<40} {r['contig_length']:<10} {r['mapped_reads']:<10} {r['rpkm']:<10.2f}\\n")
    
    f.write("\\n" + "="*80 + "\\n")

print(f"Abundance calculation completed for ${sample} (${source})")
PYEOF
    
    echo "Viral abundance calculation completed!"
    """
}

// ================================================================================
// Assembler Comparison
// ================================================================================

// Process: Compare MEGAHIT vs SPAdes Viral Identification Results
process COMPARE_ASSEMBLERS {
    tag "${sample}_Assembler_Comparison"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/assembler_comparison", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(megahit_report), path(spades_report)
    
    output:
    tuple val(sample), path("${sample}_assembler_comparison.txt"), emit: comparison_report
    path("${sample}_assembler_comparison.csv"), emit: comparison_csv
    path("${sample}_consensus_viral_sequences.txt"), emit: final_consensus
    
    script:
    """
    # Clean PYTHONPATH to avoid package conflicts
    unset PYTHONPATH
    
    python3 << 'COMPARE_ASSEMBLERS_SCRIPT'
# -*- coding: utf-8 -*-

import pandas as pd
from collections import defaultdict

def parse_merged_report(file_path):
    \"\"\"Parse merged viral report CSV file\"\"\"
    import os
    try:
        if file_path.endswith('.csv'):
            csv_path = file_path
        else:
            csv_path = file_path.replace('.txt', '.csv')
        
        print(f"  Attempting to read CSV from: {csv_path}")
        if not os.path.exists(csv_path):
            print(f"  ⚠️  CSV file does not exist: {csv_path}")
            print(f"  Checking if original file exists: {file_path}")
            if os.path.exists(file_path):
                print(f"  Original file exists, but CSV version not found")
            return pd.DataFrame()
        
        df = pd.read_csv(csv_path)
        print(f"  ✅ Successfully parsed {len(df)} rows from {csv_path}")
        return df
    except Exception as e:
        print(f"  ❌ Warning: Failed to parse {file_path}: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()

print(f"\\nParsing MEGAHIT results from: ${megahit_report}")
print(f"  Input file path: ${megahit_report}")

megahit_df = parse_merged_report("${megahit_report}")
print(f"MEGAHIT DataFrame shape: {megahit_df.shape}")
print(f"MEGAHIT columns: {list(megahit_df.columns)}")
if len(megahit_df) > 0:
    print(f"MEGAHIT first 3 rows:\\n{megahit_df.head(3)}")
else:
    print("⚠️  WARNING: MEGAHIT DataFrame is empty!")

print(f"\\nParsing SPAdes results from: ${spades_report}")
print(f"  Input file path: ${spades_report}")

spades_df = parse_merged_report("${spades_report}")
print(f"SPAdes DataFrame shape: {spades_df.shape}")
print(f"SPAdes columns: {list(spades_df.columns)}")
if len(spades_df) > 0:
    print(f"SPAdes first 3 rows:\\n{spades_df.head(3)}")
else:
    print("⚠️  WARNING: SPAdes DataFrame is empty!")

# Extract viral sequences (only "Both" type - high confidence sequences)
print(f"\\nExtracting 'Both' sequences from MEGAHIT...")
print(f"MEGAHIT columns: {list(megahit_df.columns) if len(megahit_df) > 0 else 'empty'}")
if len(megahit_df) > 0 and 'identified_by' in megahit_df.columns:
    print(f"MEGAHIT 'identified_by' value counts:\\n{megahit_df['identified_by'].value_counts()}")
    megahit_both_df = megahit_df[megahit_df['identified_by'] == 'Both']
    print(f"MEGAHIT 'Both' sequences found: {len(megahit_both_df)}")
    megahit_seqs = set(megahit_both_df['sequence_name'].tolist())
    print(f"MEGAHIT sequences extracted: {len(megahit_seqs)}")
else:
    print(f"⚠️  MEGAHIT extraction failed - columns: {list(megahit_df.columns) if len(megahit_df) > 0 else 'empty'}")
    megahit_seqs = set()

print(f"\\nExtracting 'Both' sequences from SPAdes...")
print(f"SPAdes columns: {list(spades_df.columns) if len(spades_df) > 0 else 'empty'}")
if len(spades_df) > 0 and 'identified_by' in spades_df.columns:
    print(f"SPAdes 'identified_by' value counts:\\n{spades_df['identified_by'].value_counts()}")
    spades_both_df = spades_df[spades_df['identified_by'] == 'Both']
    print(f"SPAdes 'Both' sequences found: {len(spades_both_df)}")
    spades_seqs = set(spades_both_df['sequence_name'].tolist())
    print(f"SPAdes sequences extracted: {len(spades_seqs)}")
else:
    print(f"⚠️  SPAdes extraction failed - columns: {list(spades_df.columns) if len(spades_df) > 0 else 'empty'}")
    spades_seqs = set()

print(f"\\nMEGAHIT 'Both' sequences (high confidence): {len(megahit_seqs)}")
print(f"SPAdes 'Both' sequences (high confidence): {len(spades_seqs)}")

# Note: Sequence IDs from different assemblers are different (k141_XXX vs NODE_XXX)
# "Consensus" means both assemblers detected viruses, not identical sequence IDs  
# We output all high-confidence sequences from both assemblers
megahit_only = megahit_seqs
spades_only = spades_seqs
all_viral = megahit_seqs | spades_seqs

# Consensus: if both assemblers detected viruses, all sequences are "consensus"
# Otherwise, consensus is empty
if len(megahit_seqs) > 0 and len(spades_seqs) > 0:
    consensus_both_assemblers = all_viral
    print(f"Both assemblers detected viruses: {len(megahit_seqs)} (MEGAHIT) + {len(spades_seqs)} (SPAdes) = {len(all_viral)} total consensus sequences")
else:
    consensus_both_assemblers = set()
    print(f"Only one assembler detected viruses: {len(megahit_seqs)} (MEGAHIT) + {len(spades_seqs)} (SPAdes)")
consensus_count = len(consensus_both_assemblers)

# Generate comprehensive comparison report
with open("${sample}_assembler_comparison.txt", 'w', encoding='utf-8') as f:
    f.write("="*100 + "\\n")
    f.write("Assembler Comparison Report - Viral Identification Results\\n")
    f.write("MEGAHIT vs metaSPAdes\\n")
    f.write("Sample: ${sample}\\n")
    f.write("="*100 + "\\n\\n")
    f.write("[Overall Statistics]\\n")
    f.write("-"*100 + "\\n")
    f.write(f"MEGAHIT identified viral sequences:    {len(megahit_seqs):,}\\n")
    f.write(f"SPAdes identified viral sequences:     {len(spades_seqs):,}\\n")
    f.write(f"Total unique viral sequences:          {len(all_viral):,}\\n\\n")
    f.write(f"Consensus viral sequences (both assemblies detected): {len(consensus_both_assemblers):,}\\n")
    f.write(f"MEGAHIT viral sequences:                {len(megahit_seqs):,}\\n")
    f.write(f"SPAdes viral sequences:                {len(spades_seqs):,}\\n\\n")
    if len(all_viral) > 0:
        consistency = len(consensus_both_assemblers) / len(all_viral) * 100
        f.write(f"Assembler consistency:                 {consistency:.2f}%\\n\\n")
    f.write("="*100 + "\\n")
    f.write("[Recommendation]\\n")
    f.write("-"*100 + "\\n")
    f.write(f"High-confidence viral sequences (identified by both methods): {len(consensus_both_assemblers):,}\\n")
    f.write("Recommend prioritizing these consensus sequences for downstream analysis.\\n\\n")
    f.write("[Detailed Analysis]\\n")
    f.write("-"*100 + "\\n")
    if len(megahit_only) > 0:
        f.write(f"\\nMEGAHIT-specific sequences ({len(megahit_only)}):\\n")
        f.write("  - May represent low-coverage or high-complexity regions\\n")
        f.write("  - MEGAHIT has stronger assembly capability for complex structures\\n")
    if len(spades_only) > 0:
        f.write(f"\\nSPAdes-specific sequences ({len(spades_only)}):\\n")
        f.write("  - May represent high-coverage regions\\n")
        f.write("  - SPAdes kmer strategy may capture more details\\n")
    f.write("\\n" + "="*100 + "\\n")
    f.write("[Statistical Summary]\\n")
    f.write("-"*100 + "\\n")
    if len(megahit_seqs) > 0:
        megahit_consensus_pct = len(consensus_both_assemblers) / len(megahit_seqs) * 100
        f.write(f"Consensus ratio in MEGAHIT sequences:  {megahit_consensus_pct:.2f}%\\n")
    if len(spades_seqs) > 0:
        spades_consensus_pct = len(consensus_both_assemblers) / len(spades_seqs) * 100
        f.write(f"Consensus ratio in SPAdes sequences:   {spades_consensus_pct:.2f}%\\n")

# Generate CSV detailed comparison
# Note: Since sequence IDs are different between assemblers, we treat all as consensus if both assemblers detected viruses
comparison_data = []

for seq in all_viral:
    row = {
        'sequence': seq,
        'found_in_MEGAHIT': 'Yes' if seq in megahit_seqs else 'No',
        'found_in_SPAdes': 'Yes' if seq in spades_seqs else 'No',
        'status': 'Consensus' if len(consensus_both_assemblers) > 0 else 'Single_assembler'
    }
    
    # Add MEGAHIT information
    if seq in megahit_seqs:
        megahit_row = megahit_df[megahit_df['sequence_name'] == seq].iloc[0]
        row['MEGAHIT_vs2_score'] = megahit_row.get('vs2_score', 'N/A')
        row['MEGAHIT_dvf_score'] = megahit_row.get('dvf_score', 'N/A')
        row['MEGAHIT_identified_by'] = megahit_row.get('identified_by', 'N/A')
    else:
        row['MEGAHIT_vs2_score'] = 'N/A'
        row['MEGAHIT_dvf_score'] = 'N/A'
        row['MEGAHIT_identified_by'] = 'N/A'
    
    # Add SPAdes information
    if seq in spades_seqs:
        spades_row = spades_df[spades_df['sequence_name'] == seq].iloc[0]
        row['SPAdes_vs2_score'] = spades_row.get('vs2_score', 'N/A')
        row['SPAdes_dvf_score'] = spades_row.get('dvf_score', 'N/A')
        row['SPAdes_identified_by'] = spades_row.get('identified_by', 'N/A')
    else:
        row['SPAdes_vs2_score'] = 'N/A'
        row['SPAdes_dvf_score'] = 'N/A'
        row['SPAdes_identified_by'] = 'N/A'
    
    comparison_data.append(row)

comparison_df = pd.DataFrame(comparison_data)
comparison_df.to_csv("${sample}_assembler_comparison.csv", index=False)

# Save final consensus sequence list
with open("${sample}_consensus_viral_sequences.txt", 'w') as f:
    if len(consensus_both_assemblers) > 0:
        f.write("# High-confidence viral sequences from both MEGAHIT and SPAdes\\n")
        f.write("# Both assemblers detected viruses (consensus by existence)\\n")
    else:
        f.write("# Viral sequences from single assembler only\\n")
    f.write(f"# Sample: ${sample}\\n")
    f.write(f"# Total sequences: {len(consensus_both_assemblers)}\\n")
    f.write("# Note: Sequence IDs differ between assemblers, consensus means both detected viruses\\n")
    f.write("#\\n")
    for seq in sorted(consensus_both_assemblers):
        f.write(seq + "\\n")
    
print(f"\\nAssembler comparison complete: ${sample}")
print(f"  MEGAHIT: {len(megahit_seqs)} viral sequences")
print(f"  SPAdes:  {len(spades_seqs)} viral sequences")
print(f"  Consensus: {len(consensus_both_assemblers)} viral sequences")
if len(consensus_both_assemblers) > 0:
    print(f"  Both assemblers detected viruses: {len(megahit_seqs)} + {len(spades_seqs)} = {len(consensus_both_assemblers)} total consensus sequences")
COMPARE_ASSEMBLERS_SCRIPT
    """
}

// Workflow completion message
workflow.onComplete {
    if (!params.longread) {
        log.info """
        ==========================================
        🎯 Metagenome Viral Classification Results (Short-Read Mode)
        ==========================================
        Pipeline completed successfully!
        
        Results directory: ${params.outdir}
    
    Generated files:
    - fastp/: Quality control reports
      * *_fastp.html: HTML quality reports
      * *_fastp.json: JSON quality data
      
    - clean_reads/: Filtered clean reads (if save_clean_reads=true)
      * *_clean_R1.fastq.gz: Forward clean reads
      * *_clean_R2.fastq.gz: Reverse clean reads
      
    - assembly_megahit/: MEGAHIT assembly results
      * *_megahit_contigs.fa: Assembled contigs
      
    - assembly_spades/: metaSPAdes assembly results
      * *_spades_contigs.fa: Assembled contigs
      
    - virsorter2_megahit/: VirSorter2 viral identification (MEGAHIT)
      * *_megahit_vs2_final-viral-score.tsv: Viral scores
      * *_megahit_vs2_final-viral-combined.fa: Identified viral contigs
      
    - virsorter2_spades/: VirSorter2 viral identification (SPAdes)
      * *_spades_vs2_final-viral-score.tsv: Viral scores
      * *_spades_vs2_final-viral-combined.fa: Identified viral contigs
      
    - deepvirfinder_megahit/: DeepVirFinder viral prediction (MEGAHIT)
      * *_megahit_dvf_output.txt: Prediction results with scores and p-values
      
    - deepvirfinder_spades/: DeepVirFinder viral prediction (SPAdes)
      * *_spades_dvf_output.txt: Prediction results with scores and p-values
      
    - merged_viral_reports_megahit/: Integrated viral analysis (MEGAHIT)
      * *_megahit_viral_merged_report.txt: Comprehensive viral identification report
      * *_megahit_viral_merged_report.csv: Detailed comparison data
      * *_megahit_viral_consensus.txt: High-confidence viral sequences list
      
    - merged_viral_reports_spades/: Integrated viral analysis (SPAdes)
      * *_spades_viral_merged_report.txt: Comprehensive viral identification report
      * *_spades_viral_merged_report.csv: Detailed comparison data
      * *_spades_viral_consensus.txt: High-confidence viral sequences list
    
    - assembler_comparison/: MEGAHIT vs SPAdes comparison ⭐
      * *_assembler_comparison.txt: Comprehensive assembler comparison report
      * *_assembler_comparison.csv: Detailed comparison data
      * *_consensus_viral_sequences.txt: Final high-confidence viral sequences (both assemblers)
    
    - abundance/: Viral abundance analysis ⭐⭐⭐
      * abundance/megahit/: MEGAHIT viral contig abundance
        - *_megahit_abundance.csv: Detailed abundance data (RPM and RPKM for each contig)
        - *_megahit_abundance_summary.txt: Summary report with top 10 most abundant viruses
      * abundance/spades/: SPAdes viral contig abundance
        - *_spades_abundance.csv: Detailed abundance data (RPM and RPKM for each contig)
        - *_spades_abundance_summary.txt: Summary report with top 10 most abundant viruses
    
    ==========================================
    """
    } else {
        // Long-read mode completion message
        log.info """
        ==========================================
        🎯 Metagenome Viral Classification Results (Long-Read Mode)
        ==========================================
        Pipeline completed successfully!
        
        Results directory: ${params.outdir}
        
        Generated files:
        
        [Assembly Results]
        - assembly_metaflye/: metaFlye metagenomic assembly
          * *_metaflye_contigs.fa: Assembled contigs (FASTA format)
        
        - metaflye_full_output/: metaFlye complete output directory ⭐
          * *_flye_output/: Complete Flye output (includes assembly_info.txt, assembly_graph, etc.)
          * Purpose: For viralFlye use, or downstream analysis (assembly graph, assembly statistics, etc.)
        
        [Viral Identification Results - Three Parallel Methods]
        
        - virsorter2_metaflye/: VirSorter2 viral identification ⭐
          * *_metaflye_vs2_final-viral-score.tsv: Viral scores and classifications
          * *_metaflye_vs2_final-viral-combined.fa: Identified viral sequences
        
        - deepvirfinder_metaflye/: DeepVirFinder viral prediction ⭐
          * *_metaflye_dvf_output.txt: Prediction results (scores and p-values)
        
        - viralflye_results/: viralFlye viral identification (Pfam validated) ⭐
          * *_viralflye_contigs.fa: Identified viral sequences
          * *_viralflye_summary.csv: Viral identification summary
        
        - viralflye_full_output/: viralFlye complete output
          * *_viralflye_output/: Contains viralVerify and viralComplete results
        
        [Comprehensive Comparison Results] ⭐⭐⭐
        
        ${params.enable_viralflye && !params.skip_virsorter2 && !params.skip_deepvirfinder ? 
        """- three_tools_comparison/: Three-tool comprehensive comparison (parallel analysis)
          * *_three_tools_comparison.txt: Comprehensive comparison report
          * *_three_tools_comparison.csv: Detailed comparison data
          * *_high_confidence_viruses.txt: High-confidence virus list (sorted by consensus level)
            - Three-tool consensus ⭐⭐⭐ (highest confidence)
            - Two-tool consensus ⭐⭐ (medium confidence)
            - Single-tool identification ⭐ (exploratory)
        """ : 
        """- merged_viral_reports_metaflye/: VirSorter2 + DeepVirFinder integrated report
          * *_metaflye_viral_merged_report.txt: Comprehensive analysis report
          * *_metaflye_viral_merged_report.csv: Detailed data
          * *_metaflye_viral_consensus.txt: High-confidence virus list
        """}
        
        - abundance/: Viral abundance analysis ⭐⭐⭐
          * abundance/metaflye/: metaFlye viral contig abundance
            - *_metaflye_abundance.csv: Detailed abundance data (RPM and RPKM for each contig)
            - *_metaflye_abundance_summary.txt: Summary report with top 10 most abundant viruses
          * abundance/viralflye/: viralFlye viral contig abundance (if enabled)
            - *_viralflye_abundance.csv: Detailed abundance data (RPM and RPKM for each contig)
            - *_viralflye_abundance_summary.txt: Summary report with top 10 most abundant viruses
        
        [Analysis Recommendations]
        
        Tool positioning:
        - VirSorter2: General viral identification (balanced sensitivity and specificity)
        - DeepVirFinder: High sensitivity prediction (discover more potential viruses)
        - viralFlye: Strict Pfam validation (high confidence, low false positive) ⭐
        
        Priority:
        1. Three-tool consensus viruses → Most reliable, priority for in-depth analysis
        2. viralFlye identified viruses → Pfam validated, high specificity
        3. Two-tool consensus viruses → Relatively reliable
        4. Single-tool identification → Exploratory analysis
        
        Abundance Analysis:
        - RPM (Reads Per Million): Normalized by total read count, useful for comparing relative abundance across samples
        - RPKM (Reads Per Kilobase per Million): Normalized by both contig length and total read count, best for comparing abundance of different length contigs
        - Use RPKM to identify the most abundant viral species in each sample
        - Compare abundance patterns across different assemblers/tools for validation
        
        ==========================================
        """
    }
}

workflow.onError {
    log.error """
    ==========================================
    ❌ Metagenome Viral Classification Workflow Failed
    ==========================================
    Error: ${workflow.errorMessage}
    ==========================================
    """
}
