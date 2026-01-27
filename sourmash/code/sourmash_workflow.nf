/*
 * Sourmash metagenomic analysis workflow
 *
 * Features:
 * - Optional assembly using metaSPAdes or MEGAHIT
 * - Optional direct read sketching (no assembly)
 * - Viral detection/classification using sourmash
 *
 * Usage:
 * nextflow run sourmash_workflow.nf --reads "data/*_{R1,R2}.fastq.gz" --skip_assembly false --assembler megahit
 */

// Workflow parameters
params.reads = null
params.long_reads = null  // Long-read data (single-end)
params.samplesheet = null  // samplesheet.csv path (optional, legacy combined format)
params.samplesheet_short = null  // samplesheet_short.csv path (short reads)
params.samplesheet_long = null  // samplesheet_long.csv path (long reads)
params.skip_assembly = false
params.skip_long_assembly = false  // Whether to skip long-read assembly
params.assembler = 'megahit'  // 'metaspades' or 'megahit' or 'both' (run both)
params.run_both_assemblers = false  // Whether to run both MEGAHIT and metaSPAdes
params.long_assembler = 'metaflye'  // Long-read assembler: 'metaflye'
params.use_local_flye = true  // Use locally installed Flye (in sourmash_env), avoid container image issues
params.use_local_sourmash = true  // Use locally installed sourmash (in sourmash_env), avoid container image issues
params.use_local_megahit = true  // Use locally installed MEGAHIT (in sourmash_env), avoid container image issues
params.use_local_spades = true  // Use locally installed metaSPAdes (in sourmash_env), avoid container image issues
// MEGAHIT parameters (tunable)
// Note: MEGAHIT --memory is a fraction of *physical node RAM*. Under Slurm/cgroup limits, using e.g. 0.9 can OOM-kill the job.
// We use megahit_memory_safety as a safety factor against the job's effective RAM limit and convert it to a physical-RAM fraction at runtime.
params.megahit_memory_safety = 0.85  // Safety factor (0-1). Smaller => less RAM. Recommended 0.7~0.9
params.megahit_min_count = 2
params.megahit_k_min = 21
params.megahit_k_max = 121
params.megahit_k_step = 20
params.megahit_extra_args = ''  // Extra args (optional), e.g. "--prune-level 2"
params.viral_db = '/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/sourmash_db/ncbi-viruses-2025.01.k31.zip'
params.bacterial_db = null
params.k_value = 31
params.scaled = 1000
params.threshold_bp = 500  // Default threshold: 500 bp (good tradeoff for wastewater metagenomes)
params.gather_num_results = null  // Limit number of gather results (null = unlimited)
params.gather_linear = true  // Use linear mode (more sensitive, slower)
params.gather_ignore_abundance = true  // Ignore abundance (may detect more low-abundance matches)
params.outdir = 'results'
params.work_dir = 'work'

// Validate required parameters
if (!params.samplesheet && !params.samplesheet_short && !params.samplesheet_long && !params.reads && !params.long_reads) {
    exit 1, "Error: you must provide one of the following:\n" +
            "  --samplesheet: samplesheet.csv path (legacy combined format)\n" +
            "  --samplesheet_short: samplesheet_short.csv path (short reads)\n" +
            "  --samplesheet_long: samplesheet_long.csv path (long reads)\n" +
            "  --reads: short-read glob (e.g. 'data/*_{R1,R2}.fastq.gz')\n" +
            "  --long_reads: long-read glob (e.g. 'data/*.fastq.gz')"
}

// Print configuration
log.info """
========================================
  Sourmash metagenomic analysis workflow
========================================
Samplesheet (legacy): ${params.samplesheet ?: 'Not set'}
Samplesheet (short reads): ${params.samplesheet_short ?: 'Not set'}
Samplesheet (long reads): ${params.samplesheet_long ?: 'Not set'}
Reads: ${params.reads ?: 'Not set'}
Long reads: ${params.long_reads ?: 'Not set'}
Skip assembly: ${params.skip_assembly}
Skip long-read assembly: ${params.skip_long_assembly}
Assembler: ${params.assembler}
Long-read assembler: ${params.long_assembler}
Use local Flye: ${params.use_local_flye}
Use local sourmash: ${params.use_local_sourmash}
Use local MEGAHIT: ${params.use_local_megahit}
Use local metaSPAdes: ${params.use_local_spades}
Viral DB: ${params.viral_db}
Bacterial DB: ${params.bacterial_db ?: 'Not set'}
k-mer size: ${params.k_value}
Scaled: ${params.scaled}
Output directory: ${params.outdir}
========================================
"""

// Channel definitions
def reads_ch
def long_reads_ch

// Prefer separate samplesheets (recommended)
if (params.samplesheet_short || params.samplesheet_long) {
    // Parse short-read samplesheet
    if (params.samplesheet_short) {
        log.info "Using short-read samplesheet: ${params.samplesheet_short}"
        reads_ch = Channel.fromPath(params.samplesheet_short)
            .splitCsv(header: true, sep: ',')
            // Note: skip blank lines; if required columns (R1/R2) are missing, fail fast to avoid file('') crashes
            .filter { row ->
                def r1 = (row.fastq_1 ?: row.R1 ?: row.read1 ?: row.fastq_r1 ?: '').toString().trim()
                def r2 = (row.fastq_2 ?: row.R2 ?: row.read2 ?: row.fastq_r2 ?: '').toString().trim()
                // Skip blank line (e.g. a trailing empty row in CSV)
                if (!r1 && !r2) return false
                // Non-blank line but missing required columns: provide a clear error
                if (!r1 || !r2) {
                    throw new IllegalArgumentException("Error: samplesheet_short contains a row missing R1/R2. Please check fastq_1/fastq_2 (or R1/R2) columns. Row: ${row}")
                }
                return true
            }
            .map { row ->
                def sample_id = row.sample_id ?: row.sample ?: row.id
                def r1 = (row.fastq_1 ?: row.R1 ?: row.read1 ?: row.fastq_r1 ?: '').toString().trim()
                def r2 = (row.fastq_2 ?: row.R2 ?: row.read2 ?: row.fastq_r2 ?: '').toString().trim()
                
                // If sample_id is empty, infer it from the filename
                if (!sample_id || sample_id.trim().isEmpty()) {
                    if (r1) {
                        def r1_file = file(r1)
                        sample_id = r1_file.baseName.replaceAll('_R1$', '').replaceAll('\\.fastq.*', '')
                    } else {
                        sample_id = "sample_${UUID.randomUUID().toString().substring(0, 8)}"
                    }
                }
                
                def r1_file = file(r1)
                def r2_file = file(r2)
                [sample_id.trim(), [r1_file, r2_file]]
            }
    } else {
        reads_ch = Channel.empty()
    }
    
    // Parse long-read samplesheet
    if (params.samplesheet_long) {
        log.info "Using long-read samplesheet: ${params.samplesheet_long}"
        long_reads_ch = Channel.fromPath(params.samplesheet_long)
            .splitCsv(header: true, sep: ',')
            // Note: skip blank lines; if required column (long_read) is missing, fail fast to avoid file('') crashes
            .filter { row ->
                def long_read = (row.fastq_long ?: row.long_read ?: row.long_reads ?: '').toString().trim()
                // Skip blank line
                if (!long_read) return false
                return true
            }
            .map { row ->
                def sample_id = row.sample_id ?: row.sample ?: row.id
                def long_read = (row.fastq_long ?: row.long_read ?: row.long_reads ?: '').toString().trim()
                
                // If sample_id is empty, infer it from the filename
                if (!sample_id || sample_id.trim().isEmpty()) {
                    if (long_read) {
                        def long_file = file(long_read)
                        sample_id = long_file.baseName.replaceAll('_R[12]$', '').replaceAll('\\.fastq.*', '')
                    } else {
                        sample_id = "sample_${UUID.randomUUID().toString().substring(0, 8)}"
                    }
                }
                
                def long_file = file(long_read)
                [sample_id.trim(), long_file]
            }
    } else {
        long_reads_ch = Channel.empty()
    }
    
} else if (params.samplesheet) {
    // Legacy mode: single combined samplesheet
    log.info "Using samplesheet: ${params.samplesheet}"
    
    // Read samplesheet
    def samples = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            def sample_id = row.sample_id ?: row.sample ?: row.id
            def r1 = row.R1 ?: row.read1 ?: row.fastq_r1 ?: ''
            def r2 = row.R2 ?: row.read2 ?: row.fastq_r2 ?: ''
            def long_read = row.long_read ?: row.long_reads ?: row.fastq_long ?: ''
            
            // If sample_id is empty, infer it from the filename
            if (!sample_id || sample_id.trim().isEmpty()) {
                if (long_read) {
                    def long_file = file(long_read)
                    sample_id = long_file.baseName.replaceAll('_R[12]$', '').replaceAll('\\.fastq.*', '')
                } else if (r1) {
                    def r1_file = file(r1)
                    sample_id = r1_file.baseName.replaceAll('_R1$', '').replaceAll('\\.fastq.*', '')
                } else {
                    sample_id = "sample_${UUID.randomUUID().toString().substring(0, 8)}"
                }
            }
            
            // Return sample record
            [
                sample_id: sample_id.trim(),
                r1: r1.trim(),
                r2: r2.trim(),
                long_read: long_read.trim()
            ]
        }
    
    // Split short-read and long-read inputs
    // Short reads require both R1 and R2
    reads_ch = samples.filter { it.r1 && !it.r1.isEmpty() && it.r2 && !it.r2.isEmpty() }
        .map { 
            def r1_file = file(it.r1)
            def r2_file = file(it.r2)
            [it.sample_id, [r1_file, r2_file]] 
        }
    
    // Long reads require long_read
    long_reads_ch = samples.filter { it.long_read && !it.long_read.isEmpty() }
        .map { 
            def long_file = file(it.long_read)
            [it.sample_id, long_file] 
        }
    
} else {
    // Glob mode
    // Short reads (paired-end)
    if (params.reads) {
        reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    } else {
        reads_ch = Channel.empty()
    }
    
    // Long reads (single-end)
    if (params.long_reads) {
        long_reads_ch = Channel.fromPath(params.long_reads, checkIfExists: true)
            .map { file -> 
                def sample_id = file.baseName.replaceAll('_R[12]$', '').replaceAll('\\.fastq.*', '')
                [sample_id, file] 
            }
    } else {
        long_reads_ch = Channel.empty()
    }
}

// Process 1: Assembly - metaSPAdes
process META_SPADES {
    label 'high_memory'
    
    // Conditional container: if using local metaSPAdes, do not use a container
    if (!params.use_local_spades) {
        container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    }
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    // Note: emit contigs together with assembler to avoid downstream signature/CSV naming collisions
    tuple val(sample_id), val('metaspades'), path("${sample_id}/contigs.fasta"), emit: contigs
    path("${sample_id}"), emit: assembly_dir
    
    publishDir "${params.outdir}/assembly/metaspades", mode: 'copy'
    
    script:
    """
    mkdir -p ${sample_id}
    
    # Decompress reads if needed
    if [[ \$(echo ${reads[0]} | grep -c '\\.gz\$') -eq 1 ]]; then
        zcat ${reads[0]} > ${sample_id}/R1.fastq
        zcat ${reads[1]} > ${sample_id}/R2.fastq
        R1=${sample_id}/R1.fastq
        R2=${sample_id}/R2.fastq
    else
        R1=${reads[0]}
        R2=${reads[1]}
    fi
    
    # Run metaSPAdes
    if [ "${params.use_local_spades}" == "true" ]; then
        # Locate conda
        CONDA_CMD=""
        # Prefer cluster Miniforge absolute path (Slurm child jobs may not inherit module env)
        if [ -f "/apps/eb/Miniforge3/24.11.3-0/bin/conda" ]; then
            CONDA_CMD="/apps/eb/Miniforge3/24.11.3-0/bin/conda"
        elif [ -f "/opt/conda/bin/conda" ]; then
            CONDA_CMD="/opt/conda/bin/conda"
        elif [ -f "\$HOME/miniconda3/bin/conda" ]; then
            CONDA_CMD="\$HOME/miniconda3/bin/conda"
        elif [ -f "\$HOME/anaconda3/bin/conda" ]; then
            CONDA_CMD="\$HOME/anaconda3/bin/conda"
        elif [ -f "/scratch/sp96859/conda/bin/conda" ]; then
            CONDA_CMD="/scratch/sp96859/conda/bin/conda"
        else
            CONDA_CMD=\$(which conda 2>/dev/null || echo "")
        fi
        
        if [ -z "\$CONDA_CMD" ] || [ ! -f "\$CONDA_CMD" ]; then
            echo "Error: could not find conda" >&2
            exit 1
        fi
        
        echo "Using local metaSPAdes (sourmash_env)..."
        \$CONDA_CMD run -n sourmash_env metaspades.py \\
            -1 \$R1 \\
            -2 \$R2 \\
            -o ${sample_id}/assembly \\
            -t ${task.cpus} \\
            -m ${task.memory.toGiga()} \\
            --only-assembler
    else
        # Use containerized metaSPAdes
        echo "Using containerized metaSPAdes..."
        metaspades.py \\
            -1 \$R1 \\
            -2 \$R2 \\
            -o ${sample_id}/assembly \\
            -t ${task.cpus} \\
            -m ${task.memory.toGiga()} \\
            --only-assembler
    fi
    
    # Copy contigs to output directory
    cp ${sample_id}/assembly/contigs.fasta ${sample_id}/contigs.fasta
    
    # Cleanup temporary files
    if [[ -f ${sample_id}/R1.fastq ]]; then
        rm ${sample_id}/R1.fastq ${sample_id}/R2.fastq
    fi
    """
}

// Process 2: Assembly - MEGAHIT
process MEGAHIT_ASSEMBLY {
    label 'high_memory'
    
    // Conditional container: if using local MEGAHIT, do not use a container
    if (!params.use_local_megahit) {
        container 'quay.io/biocontainers/megahit:1.2.9--h5bf99c6_3'
    }
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    // Note: emit contigs together with assembler to avoid downstream signature/CSV naming collisions
    tuple val(sample_id), val('megahit'), path("${sample_id}/contigs.fasta"), emit: contigs
    path("${sample_id}"), emit: assembly_dir
    
    publishDir "${params.outdir}/assembly/megahit", mode: 'copy'
    
    script:
    """
    mkdir -p ${sample_id}
    
    # Note: Nextflow stages inputs into the Slurm work directory, but samplesheets often contain
    # ".fastq" while the actual file is ".fastq.gz". Auto-correct and validate existence here.
    R1_IN="${reads[0]}"
    R2_IN="${reads[1]}"
    
    # Auto-add/remove .gz
    if [ ! -f "\$R1_IN" ] && [ -f "\${R1_IN}.gz" ]; then R1_IN="\${R1_IN}.gz"; fi
    if [ ! -f "\$R2_IN" ] && [ -f "\${R2_IN}.gz" ]; then R2_IN="\${R2_IN}.gz"; fi
    if [ ! -f "\$R1_IN" ] && [[ "\$R1_IN" == *.gz ]] && [ -f "\${R1_IN%.gz}" ]; then R1_IN="\${R1_IN%.gz}"; fi
    if [ ! -f "\$R2_IN" ] && [[ "\$R2_IN" == *.gz ]] && [ -f "\${R2_IN%.gz}" ]; then R2_IN="\${R2_IN%.gz}"; fi
    
    # Final existence check
    if [ ! -f "\$R1_IN" ] || [ ! -f "\$R2_IN" ]; then
        echo "Error: input read files not found (samplesheet path or compression suffix mismatch?)" >&2
        echo "R1_IN=\$R1_IN" >&2
        echo "R2_IN=\$R2_IN" >&2
        echo "Working directory: \$(pwd)" >&2
        echo "Directory listing:" >&2
        ls -la >&2
        exit 1
    fi
    
    # Decompress reads if needed
    if [[ \$(echo "\$R1_IN" | grep -c '\\.gz\$') -eq 1 ]]; then
        zcat "\$R1_IN" > ${sample_id}/R1.fastq
        zcat "\$R2_IN" > ${sample_id}/R2.fastq
        R1=${sample_id}/R1.fastq
        R2=${sample_id}/R2.fastq
    else
        R1="\$R1_IN"
        R2="\$R2_IN"
    fi
    
    # Run MEGAHIT
    if [ "${params.use_local_megahit}" == "true" ]; then
        # Locate conda
        CONDA_CMD=""
        # Prefer cluster Miniforge absolute path (Slurm child jobs may not inherit module env)
        if [ -f "/apps/eb/Miniforge3/24.11.3-0/bin/conda" ]; then
            CONDA_CMD="/apps/eb/Miniforge3/24.11.3-0/bin/conda"
        elif [ -f "/opt/conda/bin/conda" ]; then
            CONDA_CMD="/opt/conda/bin/conda"
        elif [ -f "\$HOME/miniconda3/bin/conda" ]; then
            CONDA_CMD="\$HOME/miniconda3/bin/conda"
        elif [ -f "\$HOME/anaconda3/bin/conda" ]; then
            CONDA_CMD="\$HOME/anaconda3/bin/conda"
        elif [ -f "/scratch/sp96859/conda/bin/conda" ]; then
            CONDA_CMD="/scratch/sp96859/conda/bin/conda"
        else
            CONDA_CMD=\$(which conda 2>/dev/null || echo "")
        fi
        
        if [ -z "\$CONDA_CMD" ] || [ ! -f "\$CONDA_CMD" ]; then
            echo "Error: could not find conda" >&2
            exit 1
        fi
        
        echo "Using local MEGAHIT (sourmash_env)..."
        
        # ---------------------------
        # Key: compute MEGAHIT --memory dynamically
        # - /proc/meminfo: physical node RAM (MEGAHIT uses this by default)
        # - cgroup limit: Slurm job memory limit (more realistic, but may be unavailable on some clusters)
        # - Slurm env: if cgroup is unavailable, infer from SLURM_MEM_PER_NODE/SLURM_MEM_PER_CPU
        # Goal: convert (job RAM limit * safety factor) into a physical-RAM fraction for --memory, to avoid OOM kill
        # ---------------------------
        PHYS_KB=\$(awk '/MemTotal/ {print \$2}' /proc/meminfo 2>/dev/null || echo 0)
        PHYS_BYTES=\$((PHYS_KB * 1024))
        
        # cgroup v2: /sys/fs/cgroup/memory.max；v1: /sys/fs/cgroup/memory/memory.limit_in_bytes
        CGROUP_BYTES=""
        if [ -f /sys/fs/cgroup/memory.max ]; then
            CGROUP_BYTES=\$(cat /sys/fs/cgroup/memory.max 2>/dev/null || echo "")
        elif [ -f /sys/fs/cgroup/memory/memory.limit_in_bytes ]; then
            CGROUP_BYTES=\$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo "")
        fi
        
        # Treat "max" or empty as unlimited
        if [ -z "\$CGROUP_BYTES" ] || [ "\$CGROUP_BYTES" = "max" ]; then
            CGROUP_BYTES=\$PHYS_BYTES
        fi
        
        # Take the smaller one as the effective RAM limit
        ALLOWED_BYTES=\$PHYS_BYTES
        if [ "\$CGROUP_BYTES" -gt 0 ] && [ "\$CGROUP_BYTES" -lt "\$PHYS_BYTES" ]; then
            ALLOWED_BYTES=\$CGROUP_BYTES
        fi
        
        # If cgroup looks like physical RAM (likely no limit detected), infer from Slurm env vars
        # SLURM_MEM_PER_NODE is typically in MB
        if [ "\$ALLOWED_BYTES" -ge "\$PHYS_BYTES" ] || [ "\$ALLOWED_BYTES" -eq 0 ]; then
            if [ -n "\${SLURM_MEM_PER_NODE:-}" ] && [ "\${SLURM_MEM_PER_NODE}" -gt 0 ]; then
                ALLOWED_BYTES=\$((SLURM_MEM_PER_NODE * 1024 * 1024))
            elif [ -n "\${SLURM_MEM_PER_CPU:-}" ] && [ "\${SLURM_MEM_PER_CPU}" -gt 0 ] && [ -n "\${SLURM_CPUS_PER_TASK:-}" ] && [ "\${SLURM_CPUS_PER_TASK}" -gt 0 ]; then
                ALLOWED_BYTES=\$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK * 1024 * 1024))
            fi
        fi
        
        # Compute memory fraction: (effective_limit / physical_ram) * safety; clamp to [0.05, 0.95]
        MEM_FRAC=\$(awk -v a="\$ALLOWED_BYTES" -v p="\$PHYS_BYTES" -v s="${params.megahit_memory_safety}" 'BEGIN{
            if (p<=0) { print "0.25"; exit }
            v=(a/p)*s
            if (v<0.05) v=0.05
            if (v>0.95) v=0.95
            printf "%.3f", v
        }')
        
        echo "MEGAHIT memory calc: phys_bytes=\$PHYS_BYTES cgroup_bytes=\$CGROUP_BYTES allowed_bytes=\$ALLOWED_BYTES safety=${params.megahit_memory_safety} -> --memory \$MEM_FRAC"
        \$CONDA_CMD run -n sourmash_env megahit \\
            -1 \$R1 \\
            -2 \$R2 \\
            -o ${sample_id}/assembly \\
            -t ${task.cpus} \\
            --memory \$MEM_FRAC \\
            --min-count ${params.megahit_min_count} \\
            --k-min ${params.megahit_k_min} \\
            --k-max ${params.megahit_k_max} \\
            --k-step ${params.megahit_k_step} \\
            ${params.megahit_extra_args}
    else
        # Use containerized MEGAHIT
        echo "Using containerized MEGAHIT..."
        
        # Same: compute --memory dynamically to avoid OOM kill under Slurm/cgroup limits
        PHYS_KB=\$(awk '/MemTotal/ {print \$2}' /proc/meminfo 2>/dev/null || echo 0)
        PHYS_BYTES=\$((PHYS_KB * 1024))
        
        CGROUP_BYTES=""
        if [ -f /sys/fs/cgroup/memory.max ]; then
            CGROUP_BYTES=\$(cat /sys/fs/cgroup/memory.max 2>/dev/null || echo "")
        elif [ -f /sys/fs/cgroup/memory/memory.limit_in_bytes ]; then
            CGROUP_BYTES=\$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo "")
        fi
        
        if [ -z "\$CGROUP_BYTES" ] || [ "\$CGROUP_BYTES" = "max" ]; then
            CGROUP_BYTES=\$PHYS_BYTES
        fi
        
        ALLOWED_BYTES=\$PHYS_BYTES
        if [ "\$CGROUP_BYTES" -gt 0 ] && [ "\$CGROUP_BYTES" -lt "\$PHYS_BYTES" ]; then
            ALLOWED_BYTES=\$CGROUP_BYTES
        fi
        
        # If cgroup looks like physical RAM, infer from Slurm env vars
        if [ "\$ALLOWED_BYTES" -ge "\$PHYS_BYTES" ] || [ "\$ALLOWED_BYTES" -eq 0 ]; then
            if [ -n "\${SLURM_MEM_PER_NODE:-}" ] && [ "\${SLURM_MEM_PER_NODE}" -gt 0 ]; then
                ALLOWED_BYTES=\$((SLURM_MEM_PER_NODE * 1024 * 1024))
            elif [ -n "\${SLURM_MEM_PER_CPU:-}" ] && [ "\${SLURM_MEM_PER_CPU}" -gt 0 ] && [ -n "\${SLURM_CPUS_PER_TASK:-}" ] && [ "\${SLURM_CPUS_PER_TASK}" -gt 0 ]; then
                ALLOWED_BYTES=\$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK * 1024 * 1024))
            fi
        fi
        
        MEM_FRAC=\$(awk -v a="\$ALLOWED_BYTES" -v p="\$PHYS_BYTES" -v s="${params.megahit_memory_safety}" 'BEGIN{
            if (p<=0) { print "0.25"; exit }
            v=(a/p)*s
            if (v<0.05) v=0.05
            if (v>0.95) v=0.95
            printf "%.3f", v
        }')
        
        echo "MEGAHIT memory calc: phys_bytes=\$PHYS_BYTES cgroup_bytes=\$CGROUP_BYTES allowed_bytes=\$ALLOWED_BYTES safety=${params.megahit_memory_safety} -> --memory \$MEM_FRAC"
        megahit \\
            -1 \$R1 \\
            -2 \$R2 \\
            -o ${sample_id}/assembly \\
            -t ${task.cpus} \\
            --memory \$MEM_FRAC \\
            --min-count ${params.megahit_min_count} \\
            --k-min ${params.megahit_k_min} \\
            --k-max ${params.megahit_k_max} \\
            --k-step ${params.megahit_k_step} \\
            ${params.megahit_extra_args}
    fi
    
    # Check MEGAHIT exit status
    MEGAHIT_EXIT=\$?
    if [ \$MEGAHIT_EXIT -ne 0 ]; then
        echo "Error: MEGAHIT assembly failed, exit code: \$MEGAHIT_EXIT" >&2
        exit 1
    fi
    
    # Copy contigs to output directory
    if [ -f "${sample_id}/assembly/final.contigs.fa" ]; then
        cp ${sample_id}/assembly/final.contigs.fa ${sample_id}/contigs.fasta
    else
        echo "Error: assembly failed, final.contigs.fa not found" >&2
        echo "Assembly directory listing:" >&2
        ls -la ${sample_id}/assembly/ >&2 || true
        exit 1
    fi
    
    # Cleanup temporary files
    if [[ -f ${sample_id}/R1.fastq ]]; then
        rm ${sample_id}/R1.fastq ${sample_id}/R2.fastq
    fi
    """
}

// Process 3: Long-read assembly - metaFlye
process METAFLYE_ASSEMBLY {
    label 'high_memory'
    
    input:
    tuple val(sample_id), path(read_file)
    
    output:
    // Note: emit contigs together with assembler to avoid downstream signature/CSV naming collisions
    tuple val(sample_id), val('metaflye'), path("${sample_id}/contigs.fasta"), emit: contigs
    path("${sample_id}"), emit: assembly_dir
    
    publishDir "${params.outdir}/assembly/metaflye", mode: 'copy'
    
    script:
    """
    mkdir -p ${sample_id}
    
    # If using locally installed Flye, run via conda run (no activation needed)
    if [ "${params.use_local_flye}" == "true" ]; then
        # Run flye in sourmash_env via conda run
        # Locate conda absolute path (try common locations)
        CONDA_CMD=""
        # Prefer cluster Miniforge absolute path (Slurm child jobs may not inherit module env)
        if [ -f "/apps/eb/Miniforge3/24.11.3-0/bin/conda" ]; then
            CONDA_CMD="/apps/eb/Miniforge3/24.11.3-0/bin/conda"
        elif [ -f "/opt/conda/bin/conda" ]; then
            CONDA_CMD="/opt/conda/bin/conda"
        elif [ -f "\$HOME/miniconda3/bin/conda" ]; then
            CONDA_CMD="\$HOME/miniconda3/bin/conda"
        elif [ -f "\$HOME/anaconda3/bin/conda" ]; then
            CONDA_CMD="\$HOME/anaconda3/bin/conda"
        elif [ -f "/scratch/sp96859/conda/bin/conda" ]; then
            CONDA_CMD="/scratch/sp96859/conda/bin/conda"
        else
            # Try to find conda in PATH
            CONDA_CMD=\$(which conda 2>/dev/null || echo "")
        fi
        
        if [ -z "\$CONDA_CMD" ] || [ ! -f "\$CONDA_CMD" ]; then
            echo "Error: could not find conda" >&2
            exit 1
        fi
        
        \$CONDA_CMD run -n sourmash_env flye \\
            --nano-raw ${read_file} \\
            --out-dir ${sample_id}/assembly \\
            --threads ${task.cpus} \\
            --genome-size 5m \\
            --meta \\
            --iterations 3
    else
        # Use containerized flye
        flye \\
            --nano-raw ${read_file} \\
            --out-dir ${sample_id}/assembly \\
            --threads ${task.cpus} \\
            --genome-size 5m \\
            --meta \\
            --iterations 3
    fi
    
    # Copy contigs to output directory
    if [ -f "${sample_id}/assembly/assembly.fasta" ]; then
        cp ${sample_id}/assembly/assembly.fasta ${sample_id}/contigs.fasta
    else
        echo "Error: assembly failed, assembly.fasta not found" >&2
        exit 1
    fi
    """
}

// Process 4: Sketch signatures from contigs
process SKETCH_CONTIGS {
    tag "${sample_id}_${assembler}"
    label 'medium_memory'
    // Conditional container: if using local sourmash, do not use a container
    if (!params.use_local_sourmash) {
        container 'quay.io/biocontainers/sourmash:latest'
    }
    
    input:
    // Note: must carry sample_id + assembler to produce unique signature filenames and avoid collisions
    tuple val(sample_id), val(assembler), path(contigs)
    
    output:
    path("${sample_id}_${assembler}_contigs.sig"), emit: sig
    
    publishDir "${params.outdir}/signatures", mode: 'copy'
    
    script:
    """
    # Debug info
    echo "=========================================="
    echo "  SKETCH_CONTIGS started"
    echo "=========================================="
    echo "Input: ${contigs}"
    echo "Working directory: \$(pwd)"
    echo "Files:"
    ls -lh ${contigs}* 2>/dev/null || echo "No files found"
    echo ""
    
    # Validate input file (supports path or filename)
    INPUT_FILE="${contigs}"
    if [ ! -f "\$INPUT_FILE" ]; then
        # Try to locate by basename
        CONTIGS_BASENAME=\$(basename "${contigs}")
        INPUT_FILE=\$(find . -name "\$CONTIGS_BASENAME" -type f | head -1)
        if [ -z "\$INPUT_FILE" ] || [ ! -f "\$INPUT_FILE" ]; then
            echo "Error: input file not found: ${contigs}" >&2
            echo "Directory listing:" >&2
            ls -la >&2
            exit 1
        fi
    fi
    
    echo "Using input file: \$INPUT_FILE"
    echo "File size: \$(ls -lh \$INPUT_FILE | awk '{print \$5}')"
    echo ""
    
    # Output filename
    # Note: must be unique (sample_id + assembler) or outputs may overwrite each other
    OUTPUT_SIG="${sample_id}_${assembler}_contigs.sig"
    echo "Output: \$OUTPUT_SIG"
    echo ""
    
    # Run sourmash sketch
    echo "Generating signature..."
    
    # If using local sourmash, run via conda run
    if [ "${params.use_local_sourmash}" == "true" ]; then
        # Locate conda
        CONDA_CMD=""
        # Prefer cluster Miniforge absolute path (Slurm child jobs may not inherit module env)
        if [ -f "/apps/eb/Miniforge3/24.11.3-0/bin/conda" ]; then
            CONDA_CMD="/apps/eb/Miniforge3/24.11.3-0/bin/conda"
        elif [ -f "/opt/conda/bin/conda" ]; then
            CONDA_CMD="/opt/conda/bin/conda"
        elif [ -f "\$HOME/miniconda3/bin/conda" ]; then
            CONDA_CMD="\$HOME/miniconda3/bin/conda"
        elif [ -f "\$HOME/anaconda3/bin/conda" ]; then
            CONDA_CMD="\$HOME/anaconda3/bin/conda"
        elif [ -f "/scratch/sp96859/conda/bin/conda" ]; then
            CONDA_CMD="/scratch/sp96859/conda/bin/conda"
        else
            CONDA_CMD=\$(which conda 2>/dev/null || echo "")
        fi
        
        if [ -z "\$CONDA_CMD" ] || [ ! -f "\$CONDA_CMD" ]; then
            echo "Error: could not find conda" >&2
            exit 1
        fi
        
        echo "Using local sourmash (sourmash_env)..."
        \$CONDA_CMD run -n sourmash_env sourmash sketch dna \\
            \$INPUT_FILE \\
            -p k=${params.k_value},scaled=${params.scaled} \\
            -o \$OUTPUT_SIG \\
            --name "${sample_id}_${assembler}_contigs"
    else
        # Use containerized sourmash
        echo "Using containerized sourmash..."
        sourmash sketch dna \\
            \$INPUT_FILE \\
            -p k=${params.k_value},scaled=${params.scaled} \\
            -o \$OUTPUT_SIG \\
            --name "${sample_id}_${assembler}_contigs"
    fi
    
    EXIT_CODE=\$?
    echo "sourmash exit code: \$EXIT_CODE"
    
    # Validate output
    if [ ! -f "\$OUTPUT_SIG" ]; then
        echo "Error: signature file was not created: \$OUTPUT_SIG" >&2
        echo "Directory listing:" >&2
        ls -la >&2
        exit 1
    fi
    
    echo "Signature created: \$OUTPUT_SIG"
    echo "File size: \$(ls -lh \$OUTPUT_SIG | awk '{print \$5}')"
    """
}

// Process 4 (alt): Sketch signatures from reads (when skipping assembly)
process SKETCH_READS {
    tag "${sample_id}"
    
    // Conditional container: if using local sourmash, do not use a container
    if (!params.use_local_sourmash) {
        container 'quay.io/biocontainers/sourmash:latest'
    }
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    // Note: short-read signatures must use a suffix to avoid collisions with long-read/contigs signatures
    path("${sample_id}_short_reads.sig"), emit: sig
    
    publishDir "${params.outdir}/signatures", mode: 'copy'
    
    script:
    """
    # Concatenate paired-end reads
    if [[ \$(echo ${reads[0]} | grep -c '\\.gz\$') -eq 1 ]]; then
        zcat ${reads[0]} ${reads[1]} > ${sample_id}_combined.fastq
    else
        cat ${reads[0]} ${reads[1]} > ${sample_id}_combined.fastq
    fi
    
    # Generate signature
    if [ "${params.use_local_sourmash}" == "true" ]; then
        # Locate conda
        CONDA_CMD=\$(which conda 2>/dev/null || echo "")
        if [ -z "\$CONDA_CMD" ]; then
            for path in /apps/eb/Miniforge3/24.11.3-0/bin/conda /opt/conda/bin/conda \$HOME/miniconda3/bin/conda \$HOME/anaconda3/bin/conda /scratch/sp96859/conda/bin/conda; do
                if [ -f "\$path" ]; then
                    CONDA_CMD="\$path"
                    break
                fi
            done
        fi
        
        if [ -z "\$CONDA_CMD" ] || [ ! -f "\$CONDA_CMD" ]; then
            echo "Error: could not find conda" >&2
            exit 1
        fi
        
        \$CONDA_CMD run -n sourmash_env sourmash sketch dna \\
            ${sample_id}_combined.fastq \\
            -p k=${params.k_value},scaled=${params.scaled} \\
            -o ${sample_id}_short_reads.sig \\
            --name "${sample_id}_short_reads"
    else
        sourmash sketch dna \\
            ${sample_id}_combined.fastq \\
            -p k=${params.k_value},scaled=${params.scaled} \\
            -o ${sample_id}_short_reads.sig \\
            --name "${sample_id}_short_reads"
    fi
    
    # Cleanup temporary file
    rm ${sample_id}_combined.fastq
    """
}

// Process 4b: Sketch signatures from long reads (single-end)
process SKETCH_LONG_READS {
    tag "${sample_id}"
    
    // Conditional container: if using local sourmash, do not use a container
    if (!params.use_local_sourmash) {
        container 'quay.io/biocontainers/sourmash:latest'
    }
    
    input:
    tuple val(sample_id), path(read_file)
    
    output:
    // Note: long-read signatures must use a suffix to avoid collisions with short-read/contigs signatures
    path("${sample_id}_long_reads.sig"), emit: sig
    
    publishDir "${params.outdir}/signatures", mode: 'copy'
    
    script:
    """
    # Generate signature (long reads, single-end)
    if [ "${params.use_local_sourmash}" == "true" ]; then
        # Locate conda
        CONDA_CMD=\$(which conda 2>/dev/null || echo "")
        if [ -z "\$CONDA_CMD" ]; then
            for path in /apps/eb/Miniforge3/24.11.3-0/bin/conda /opt/conda/bin/conda \$HOME/miniconda3/bin/conda \$HOME/anaconda3/bin/conda /scratch/sp96859/conda/bin/conda; do
                if [ -f "\$path" ]; then
                    CONDA_CMD="\$path"
                    break
                fi
            done
        fi
        
        if [ -z "\$CONDA_CMD" ] || [ ! -f "\$CONDA_CMD" ]; then
            echo "Error: could not find conda" >&2
            exit 1
        fi
        
        \$CONDA_CMD run -n sourmash_env sourmash sketch dna \\
            ${read_file} \\
            -p k=${params.k_value},scaled=${params.scaled} \\
            -o ${sample_id}_long_reads.sig \\
            --name "${sample_id}_long_reads"
    else
        sourmash sketch dna \\
            ${read_file} \\
            -p k=${params.k_value},scaled=${params.scaled} \\
            -o ${sample_id}_long_reads.sig \\
            --name "${sample_id}_long_reads"
    fi
    """
}

// Process 5: Viral detection (sourmash gather)
process VIRAL_GATHER {
    tag "${sig_file.baseName.replaceAll('\\.sig$', '')}"
    
    // Conditional container: if using local sourmash, do not use a container
    if (!params.use_local_sourmash) {
        container 'quay.io/biocontainers/sourmash:latest'
    }
    
    input:
    path sig_file
    
    output:
    path("*_viral.csv"), emit: results
    path("*_viral_summary.txt"), emit: summary
    
    publishDir "${params.outdir}/viral_detection", mode: 'copy'
    
    script:
    // Extract sample ID from sig_file path
    def sig_file_obj = file(sig_file)
    def sample_id = sig_file_obj.baseName.replaceAll('\\.sig$', '')
    if (!sample_id || sample_id.isEmpty()) {
        sample_id = "sample"
    }
    """
    if [ "${params.use_local_sourmash}" == "true" ]; then
        # Locate conda
        CONDA_CMD=\$(which conda 2>/dev/null || echo "")
        if [ -z "\$CONDA_CMD" ]; then
            for path in /apps/eb/Miniforge3/24.11.3-0/bin/conda /opt/conda/bin/conda \$HOME/miniconda3/bin/conda \$HOME/anaconda3/bin/conda /scratch/sp96859/conda/bin/conda; do
                if [ -f "\$path" ]; then
                    CONDA_CMD="\$path"
                    break
                fi
            done
        fi
        
        if [ -z "\$CONDA_CMD" ] || [ ! -f "\$CONDA_CMD" ]; then
            echo "Error: could not find conda" >&2
            exit 1
        fi
        
        \$CONDA_CMD run -n sourmash_env sourmash gather \\
            ${sig_file} \\
            ${params.viral_db} \\
            -o ${sample_id}_viral.csv \\
            --threshold-bp ${params.threshold_bp} \\
            --save-matches ${sample_id}_viral_matches.sig || true
    else
        sourmash gather \\
            ${sig_file} \\
            ${params.viral_db} \\
            -o ${sample_id}_viral.csv \\
            --threshold-bp ${params.threshold_bp} \\
            --save-matches ${sample_id}_viral_matches.sig \\
            ${params.gather_num_results ? "--num-results ${params.gather_num_results}" : ""} \\
            ${params.gather_linear ? "--linear" : ""} \\
            ${params.gather_ignore_abundance ? "--ignore-abundance" : ""} || true
    fi
    
    # Write summary
    echo "Viral detection summary" > ${sample_id}_viral_summary.txt
    echo "Sample: ${sample_id}" >> ${sample_id}_viral_summary.txt
    echo "Time: \$(date)" >> ${sample_id}_viral_summary.txt
    echo "" >> ${sample_id}_viral_summary.txt
    
    # Validate CSV exists (if no matches or below threshold, sourmash gather may not create a CSV)
    if [ -f "${sample_id}_viral.csv" ] && [ -s "${sample_id}_viral.csv" ]; then
        echo "Viral matches found (showing first lines):" >> ${sample_id}_viral_summary.txt
        echo "" >> ${sample_id}_viral_summary.txt
        head -20 ${sample_id}_viral.csv >> ${sample_id}_viral_summary.txt || true
    else
        echo "No viral matches above threshold (threshold: ${params.threshold_bp} bp)." >> ${sample_id}_viral_summary.txt
        echo "This can be normal (insufficient viral sequence detected in the sample)." >> ${sample_id}_viral_summary.txt
        echo "Note: even if there are matches, if overlap is below the threshold, sourmash gather may not create a CSV." >> ${sample_id}_viral_summary.txt
        # Create an empty CSV so downstream steps can proceed
        touch ${sample_id}_viral.csv
    fi
    """
}

// Process 6: Bacterial detection (optional)
process BACTERIAL_GATHER {
    tag "${sig_file.baseName.replaceAll('\\.sig$', '')}"
    
    // Conditional container: if using local sourmash, do not use a container
    if (!params.use_local_sourmash) {
        container 'quay.io/biocontainers/sourmash:latest'
    }
    
    input:
    path sig_file
    
    output:
    path("*_bacterial.csv"), emit: results
    path("*_bacterial_summary.txt"), emit: summary
    
    when:
    params.bacterial_db != null
    
    publishDir "${params.outdir}/bacterial_detection", mode: 'copy'
    
    script:
    // Extract sample ID from sig_file path
    def sig_file_obj = file(sig_file)
    def sample_id = sig_file_obj.baseName.replaceAll('\\.sig$', '')
    if (!sample_id || sample_id.isEmpty()) {
        sample_id = "sample"
    }
    """
    if [ "${params.use_local_sourmash}" == "true" ]; then
        # Locate conda
        CONDA_CMD=\$(which conda 2>/dev/null || echo "")
        if [ -z "\$CONDA_CMD" ]; then
            for path in /apps/eb/Miniforge3/24.11.3-0/bin/conda /opt/conda/bin/conda \$HOME/miniconda3/bin/conda \$HOME/anaconda3/bin/conda /scratch/sp96859/conda/bin/conda; do
                if [ -f "\$path" ]; then
                    CONDA_CMD="\$path"
                    break
                fi
            done
        fi
        
        if [ -z "\$CONDA_CMD" ] || [ ! -f "\$CONDA_CMD" ]; then
            echo "Error: could not find conda" >&2
            exit 1
        fi
        
        \$CONDA_CMD run -n sourmash_env sourmash gather \\
            ${sig_file} \\
            ${params.bacterial_db} \\
            -o ${sample_id}_bacterial.csv \\
            --threshold-bp ${params.threshold_bp} \\
            --save-matches ${sample_id}_bacterial_matches.sig \\
            ${params.gather_num_results ? "--num-results ${params.gather_num_results}" : ""} \\
            ${params.gather_linear ? "--linear" : ""} \\
            ${params.gather_ignore_abundance ? "--ignore-abundance" : ""} || true
    else
        sourmash gather \\
            ${sig_file} \\
            ${params.bacterial_db} \\
            -o ${sample_id}_bacterial.csv \\
            --threshold-bp ${params.threshold_bp} \\
            --save-matches ${sample_id}_bacterial_matches.sig \\
            ${params.gather_num_results ? "--num-results ${params.gather_num_results}" : ""} \\
            ${params.gather_linear ? "--linear" : ""} \\
            ${params.gather_ignore_abundance ? "--ignore-abundance" : ""} || true
    fi
    
    # Write summary
    echo "Bacterial detection summary" > ${sample_id}_bacterial_summary.txt
    echo "Sample: ${sample_id}" >> ${sample_id}_bacterial_summary.txt
    echo "Time: \$(date)" >> ${sample_id}_bacterial_summary.txt
    echo "" >> ${sample_id}_bacterial_summary.txt
    
    # Validate CSV exists (if no matches or below threshold, sourmash gather may not create a CSV)
    if [ -f "${sample_id}_bacterial.csv" ] && [ -s "${sample_id}_bacterial.csv" ]; then
        echo "Bacterial matches found (showing first lines):" >> ${sample_id}_bacterial_summary.txt
        echo "" >> ${sample_id}_bacterial_summary.txt
        head -20 ${sample_id}_bacterial.csv >> ${sample_id}_bacterial_summary.txt || true
    else
        echo "No bacterial matches above threshold (threshold: ${params.threshold_bp} bp)." >> ${sample_id}_bacterial_summary.txt
        echo "This can be normal (insufficient bacterial sequence detected in the sample)." >> ${sample_id}_bacterial_summary.txt
        echo "Note: even if there are matches, if overlap is below the threshold, sourmash gather may not create a CSV." >> ${sample_id}_bacterial_summary.txt
        # Create an empty CSV so downstream steps can proceed
        touch ${sample_id}_bacterial.csv
    fi
    """
}

// Process 7: Generate report
process GENERATE_REPORT {
    tag "report"
    
    input:
    // Note: stage all input CSVs into a subdirectory to further reduce same-name collision risk
    path viral_results, stageAs: 'viral_results/*'
    
    output:
    path("analysis_report.html"), emit: report
    
    publishDir "${params.outdir}", mode: 'copy'
    
    script:
    """
    cat > analysis_report.html << 'EOF'
    <!DOCTYPE html>
    <html>
    <head>
        <title>Sourmash analysis report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            h1 { color: #2c3e50; }
            h2 { color: #34495e; }
            table { border-collapse: collapse; width: 100%; margin: 20px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #3498db; color: white; }
            tr:nth-child(even) { background-color: #f2f2f2; }
        </style>
    </head>
    <body>
        <h1>Sourmash metagenomic analysis report</h1>
        <h2>Run parameters</h2>
        <ul>
            <li>Assembler: ${params.assembler}</li>
            <li>Skip assembly: ${params.skip_assembly}</li>
            <li>k-mer size: ${params.k_value}</li>
            <li>Scaled: ${params.scaled}</li>
            <li>Viral DB: ${params.viral_db}</li>
        </ul>
        <h2>Completed</h2>
        <p>See the CSV files in the output directory for detailed results.</p>
        <p>Generated at: \$(date)</p>
    </body>
    </html>
    EOF
    """
}

// Main workflow
workflow {
    // Print workflow start
    log.info "Starting workflow..."
    
    // Merge all signature channels
    def all_sig_ch = Channel.empty()
    
    // Collect all contigs (short + long) for unified processing
    def all_contigs_ch = Channel.empty()
    
    // Short-read inputs
    // Nextflow will skip processes for empty channels automatically
    def short_sig_ch = Channel.empty()
    
    // Only process if short reads were provided
    def has_short_reads = params.samplesheet_short || params.samplesheet || params.reads
    
    if (has_short_reads) {
        if (params.skip_assembly) {
            // Use reads directly (no assembly)
            log.info "Skipping assembly; sketching short reads directly"
            short_sig_ch = SKETCH_READS(reads_ch).sig
        } else {
            // Assemble first
            if (params.run_both_assemblers || params.assembler == 'both') {
                // Run both assemblers
                log.info "Running both MEGAHIT and metaSPAdes for short-read assembly"
                all_contigs_ch = all_contigs_ch.mix(META_SPADES(reads_ch).contigs)
                all_contigs_ch = all_contigs_ch.mix(MEGAHIT_ASSEMBLY(reads_ch).contigs)
            } else {
                // Run a single assembler
                log.info "Running ${params.assembler} for short-read assembly"
                
                if (params.assembler == 'metaspades') {
                    all_contigs_ch = all_contigs_ch.mix(META_SPADES(reads_ch).contigs)
                } else if (params.assembler == 'megahit') {
                    all_contigs_ch = all_contigs_ch.mix(MEGAHIT_ASSEMBLY(reads_ch).contigs)
                } else {
                    exit 1, "Error: unsupported assembler '${params.assembler}'. Choose 'metaspades', 'megahit', or 'both'."
                }
            }
        }
    }
    
    // Long-read inputs
    // Nextflow will skip processes for empty channels automatically
    def long_sig_ch = Channel.empty()
    
    // Check if long reads were provided
    def has_long_reads = params.samplesheet_long || params.samplesheet || params.long_reads
    
    if (has_long_reads) {
        if (params.skip_long_assembly) {
            // Use long reads directly (no assembly)
            log.info "Sketching long reads directly (no assembly)"
            long_sig_ch = SKETCH_LONG_READS(long_reads_ch).sig
        } else {
            // Assemble long reads
            log.info "Running ${params.long_assembler} for long-read assembly"
            
            if (params.long_assembler == 'metaflye') {
                all_contigs_ch = all_contigs_ch.mix(METAFLYE_ASSEMBLY(long_reads_ch).contigs)
            } else {
                exit 1, "Error: unsupported long-read assembler '${params.long_assembler}'. Choose 'metaflye'."
            }
        }
    }
    
    // Sketch all contigs (short + long) via a single SKETCH_CONTIGS call
    // Nextflow will skip the process for an empty channel automatically
    def contigs_sig_ch = SKETCH_CONTIGS(all_contigs_ch).sig
    
    // Merge all signatures
    all_sig_ch = all_sig_ch.mix(short_sig_ch)
    all_sig_ch = all_sig_ch.mix(long_sig_ch)
    all_sig_ch = all_sig_ch.mix(contigs_sig_ch)
    
    // Viral detection
    // If channel is empty, the process won't run
    def viral_results_ch = VIRAL_GATHER(all_sig_ch).results
    
    // Bacterial detection (optional)
    def bacterial_results_ch = Channel.empty()
    if (params.bacterial_db) {
        bacterial_results_ch = BACTERIAL_GATHER(all_sig_ch).results
    }
    
    // Generate report
    // Collect all viral result files (report only needs filenames; it does not parse contents)
    def all_viral_results = viral_results_ch.collect()
    
    // Only generate report if there are viral result files
    GENERATE_REPORT(all_viral_results)
}

// -----------------------------
// Note: the workflow block assembles the dataflow and submits tasks.
// Printing "completed" at the end of the workflow block can be misleading (tasks may only be queued).
// Use onComplete/onError callbacks to report final status after execution finishes.
// -----------------------------
workflow.onComplete {
    // Printed only when the workflow has actually completed
    log.info "Workflow completed. Results saved to: ${params.outdir}"
}

workflow.onError {
    // Printed on failure (Nextflow log contains detailed error/stack trace)
    log.error "Workflow failed. Check .nextflow.log and each task's .command.err/.command.out. Output directory: ${params.outdir}"
}

