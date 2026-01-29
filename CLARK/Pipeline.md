# CLARK Virus Classification Pipeline

This document describes the virus classification pipeline implemented in `slurm_virus_classification_paired.sh`: its workflow, the software used, and the final result files.


[![Apptainer](https://img.shields.io/badge/Apptainer-%3D%201.3-42b983)](https://apptainer.org/)
[![Python](https://img.shields.io/badge/python-%3D%203.10-blue.svg)](https://www.python.org/)

---

## Pipeline Workflow

The pipeline takes paired-end FASTQ reads (R1 and R2), classifies them against a viral reference database with CLARK, and optionally estimates taxonomic abundance. The flow is linear with one conditional step.

### Workflow Diagram

```
                    Paired-end FASTQ (R1, R2)
                    optionally gzip-compressed
                                │
                                ▼
    ┌───────────────────────────────────────────────────────────┐
    │  Step 1: Preprocessing & validation                       │
    │  • Validate CLARK database (db_central_k*.tsk.*)          │
    │  • Validate input files (existence, size, gzip integrity) │
    │  • Decompress .gz to output dir if needed; cleanup later  │
    └───────────────────────────────────────────────────────────┘
                                │
                                ▼
    ┌─────────────────────────────────────────────────────────┐
    │  Step 2: Virus classification (CLARK)                   │
    │  • classify_metagenome.sh -P R1 R2 -R <prefix> ...      │
    │  • k-mer–based assignment to viral taxa                 │
    │  • Retry with fewer threads on segmentation fault       │
    └─────────────────────────────────────────────────────────┘
                                │
                                ▼
                    Classification successful?
                          /           \
                        No             Yes
                         │               │
                         ▼               ▼
                    Job exits    ┌─────────────────────────────────────┐
                                 │  Step 3: Abundance estimation (opt) │
                                 │  • estimate_abundance.sh if present │
                                 │  • else Python script on CSV        │
                                 │  • Adds proportion/abundance columns│
                                 └─────────────────────────────────────┘
                                                │
                                                ▼
                                         Final result files
```

### Workflow in Short

1. **Preprocessing and validation**  
   The script checks that the CLARK database and input FASTQ files are valid. If inputs are gzip-compressed, they are decompressed into the output directory and removed after the run.

2. **Virus classification**  
   CLARK is run in paired-end mode (`-P`) with the given k-mer size and mode. It assigns reads to viral taxa using the reference database. If the run fails with a segmentation fault, the pipeline retries with fewer threads.

3. **Abundance estimation (optional)**  
   If enabled (`ENABLE_ABUNDANCE=true`) and classification succeeded, the pipeline either runs CLARK’s `estimate_abundance.sh` on the classification results or, if that script is missing, uses an embedded Python step to compute proportions and write an abundance table.

---

## Preprocessing and Preprocessing Software

Step 1 of the pipeline performs **preprocessing and validation** only. There is no quality trimming, adapter removal, or read filtering.

| Preprocessing step      | What is done | Software / tool |
|-------------------------|--------------|------------------|
| **Database validation** | Check that `CLARK_DB_DIR` exists and contains database files `db_central_k*.tsk.*` for the chosen k-mer size. | Bash, `find` |
| **Input file validation** | Check that R1 and R2 exist, are readable, non-empty; run `gzip -t` on `.gz` files; basic FASTQ format check (e.g. header starts with `@`). | Bash, `gzip -t`, `stat`, `head` |
| **Decompression**       | If inputs are `.fastq.gz`, decompress to plain FASTQ in `OUTPUT_DIR` so CLARK can read them; temporary decompressed files are removed at the end. | **gunzip** (gzip/gunzip) |

**Included:** Validation (Bash + standard Unix tools) and decompression (**gunzip**).  
**Not included:** No dedicated preprocessing tools such as fastp, Trimmomatic, cutadapt, or BBDuk; the pipeline assumes input FASTQ is already (optionally gzip-compressed) and does not perform quality or adapter trimming.

---

## Software Used

| Software / component      | Role in the pipeline |
|---------------------------|----------------------|
| **CLARK**                 | Metagenomic classifier. Assigns reads to viral taxa using k-mer matching against a custom database. Used via the `classify_metagenome.sh` wrapper, which calls the `exe/CLARK` binary. |
| **classify_metagenome.sh**| CLARK’s classification script. Run with `-P` (paired-end), input R1/R2, output prefix `-R`, mode `-m`, k-mer `-k`, and thread count `-n`. |
| **estimate_abundance.sh** | CLARK’s abundance script (optional). Reads the classification CSV and writes abundance/proportion statistics. Used only if the file exists in the CLARK directory. |
| **Python 3**              | Used when `estimate_abundance.sh` is not available: reads the classification CSV, computes counts and proportions (e.g. Proportion_All(%), Proportion_Classified(%)), and writes `<OUTPUT_PREFIX>_abundance.csv`. |
| **Bash**                  | Scripting language for the pipeline logic, file checks, and job control. |
| **SLURM**                 | Job scheduler. The pipeline is submitted with `sbatch`; resources (CPUs, memory, time) are set via SBATCH directives. |
| **gunzip / gzip**         | Preprocessing: decompression of `.fastq.gz` inputs when needed; decompressed files are written to the output directory and cleaned up at the end. |

The pipeline runs on the host (or current environment). It does **not** use Apptainer or Singularity; any container would need to be started outside this script.

---

## Abundance Estimation (Step 3)

When **abundance estimation** is enabled (`ENABLE_ABUNDANCE=true`) and classification has succeeded, the pipeline computes **relative abundance** (proportions) per taxon.

### What it does

- **Input:** The classification result file `<OUTPUT_PREFIX>.csv` (taxon name, taxon ID, lineage, **count**, etc.).
- **Output:** `<OUTPUT_PREFIX>_abundance.csv` — same taxa with extra columns:
  - **Proportion_All(%)** — proportion of that taxon’s count over **all** reads (or total in file).
  - **Proportion_Classified(%)** — proportion over **classified** reads only.

### Two methods (only one is used per run)

| Method | When used | Software |
|--------|-----------|----------|
| **CLARK script** | If `estimate_abundance.sh` exists under `CLARK_DIR`. | **estimate_abundance.sh** (CLARK), with database `-D`, input CSV `-F`, and optional filters. |
| **Built-in calculation** | If `estimate_abundance.sh` is **not** found. | **Python 3** (inline script): reads CSV, finds Count/Name columns, computes proportions, writes `<OUTPUT_PREFIX>_abundance.csv`. |

### Optional parameters (when using estimate_abundance.sh)

- `CONFIDENCE_THRESHOLD` (-c), `GAMMA_THRESHOLD` (-g), `ABUNDANCE_THRESHOLD` (-a)
- `HIGH_CONFIDENCE` (--highconfidence), `MPA_FORMAT` (--mpa for MetaPhlAn-style output)

These do **not** apply to the Python fallback; the Python step only adds Proportion_All(%) and Proportion_Classified(%).

---

## Final Result Files

All pipeline outputs are written under the directory given by `OUTPUT_DIR`. The prefix for result names is `OUTPUT_PREFIX` (default: `virus_results`).

### Primary result files

| Result file                     | Description |
|---------------------------------|-------------|
| **`<OUTPUT_PREFIX>.csv`**       | **Classification results.** One row per assigned taxon (e.g. virus), with columns such as name, taxon ID, lineage, and **count** of reads assigned. Produced after step 2 succeeds. Example: `virus_results.csv`. |
| **`<OUTPUT_PREFIX>_abundance.csv`** | **Abundance results.** Same taxa as the classification table, with added columns for **Proportion_All(%)** and **Proportion_Classified(%)** (and optionally other filters). Produced in step 3 when abundance estimation is enabled. Example: `virus_results_abundance.csv`. |

### Job logs

| Log file | Description |
|---------|-------------|
| **`CLARK_Virus_Classification_Paired_<jobid>.out`** | SLURM standard output (echoes, tool output). |
| **`CLARK_Virus_Classification_Paired_<jobid>.err`** | SLURM standard error (warnings, errors). |

### Where they are written

- Classification: `OUTPUT_DIR/<OUTPUT_PREFIX>.csv`
- Abundance: `OUTPUT_DIR/<OUTPUT_PREFIX>_abundance.csv`
- Logs: current working directory at submit time (or as set by SBATCH `--output` and `--error`).

### Summary of final outputs

- **Always (if classification succeeds):**  
  `OUTPUT_DIR/<OUTPUT_PREFIX>.csv` — viral classification counts.

- **If abundance estimation is enabled and step 3 runs:**  
  `OUTPUT_DIR/<OUTPUT_PREFIX>_abundance.csv` — viral classification with abundance/proportion columns.

- **Always:**  
  SLURM `.out` and `.err` log files for the job.
