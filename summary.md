# Metagenomic Viral Detection Summary (Pipelines + Tool Assessment)

## 1) Executive Summary (What to use, when)

- **Broad screening + standardized abundance tables (production-friendly)**: use **MetaTaxProfiler (nf-core/taxprofiler)**, optionally followed by confirmatory validation (GOTTCHA2 / sourmash / assembly evidence).
- **Discovery + genome/contig evidence chain (especially large DNA viruses / NCLDV)**: use **rvdb-viral-metagenome-nf** (assembly + protein homology via DIAMOND/RVDB). Add **MLMVD-nf** when you need a smaller, higher-confidence candidate set.
- **Assembly-before-classification (contig-level Kraken2/Bracken)**: use **KrakenMetaReads-nf**; strong for interpretable contig evidence, but may lose ultra-low-abundance taxa during assembly.
- **Fast similarity screening for assembled contigs**: use **sourmash**.
- **Non-Kraken2 confirmatory profiling options**: **GOTTCHA2** (signature-based, conservative) and **CLARK/CLARK-S** (discriminative k-mers; best for in-database targets, needs follow-up validation).

> **Key principle** (from the assessment): do not treat a “detection list” as a conclusion. Viral detection should be interpreted as an **evidence chain** (coverage, contig/assembly support, multi-tool consistency, and thresholds).

---

## 2) Catalog: 8 Pipelines in This Folder (At a Glance)

| Pipeline | Primary goal | Data types | Assembly | Core method | Key databases | Main outputs |
|---|---|---|---|---|---|---|
| **Sourmash** | Fast MinHash similarity/containment search | Assembled contigs | External (MetaFlye/SPAdes/MEGAHIT) | `sketch` + `gather` | NCBI viruses (k=31) | `.sig`, `gather_results*.csv` |
| **rvdb-viral-metagenome-nf** | Discovery via assembly + protein homology; dual comparisons | Short + long reads | ✅ | Prodigal + DIAMOND vs RVDB + taxonomy comparison | RVDB + NCBI taxonomy (+ Pfam-A for viralFlye) | taxonomy reports; consensus/diff lists; optional RPM/RPKM |
| **mag-nf (nf-core/mag approach)** | Assembly + Kraken2 (reads + contigs) | Short + long reads | ✅ | Kraken2 on reads/contigs | Kraken2 viral DB | `results_short/`, `results_long/` |
| **TaxProfiler-nf (MetaTaxProfiler)** | Standardized profiling + auto RPM/RPKM | Short + long reads | ❌ | Kraken2 (+ optional Bracken) + reporting | Kraken2 DB (+ Bracken DB) | `abundance/` tables + reports |
| **MLMVD-nf** | Novel/distant virus mining with multi-tool consensus | Short + long reads | ✅ | VirSorter2 + DeepVirFinder (+ viralFlye for long reads) | VS2 DB + DVF models + Pfam-A | consensus tiers; high-confidence lists; optional RPM/RPKM |
| **KrakenMetaReads-nf** | Assembly-before-classification + RPM/RPKM | Short + long reads | ✅ | Kraken2 (+ Bracken for short reads) on contigs | Kraken2 DB (+ Bracken DB) | abundance folders; merged reports |
| **GOTTCHA2** | Conservative profiling via unique signatures | Short + long reads | ❌ | GOTTCHA2 signature profiling | GOTTCHA2 `.mmi` | `*.tsv`, `*.summary.tsv`, `*.full.tsv` |
| **CLARK** | Fast discriminative k-mer classification + abundance | SE/PE reads | ❌ | CLARK / CLARK-S / CLARK-l | RefSeq (viruses-only targets possible) | classification CSV + abundance scripts |

---

## 3) Assessment Framework (Why results differ across tools)

The assessment explicitly distinguishes two output concepts to prevent contradictory claims:
- **Default output**: tool’s standard/default run results (typically higher sensitivity, more noise).
- **Consensus/threshold output**: results after applying **multi-tool cross-validation**, **host removal**, **coverage/contig thresholds**, etc. (typically higher specificity, lower sensitivity).

It also groups strategies into four practical categories:
1. **Direct read classification**: taxprofiler (and direct Kraken2/CLARK/GOTTCHA2 usage).
2. **Contig classification after assembly**: KrakenMetaReads-nf; contig interpretation from nf-core/mag outputs.
3. **Protein/alignment evidence chains**: rvdb-viral (DIAMOND/RVDB), etc.
4. **Viral feature + consensus screening**: MLMVD-nf (multi-tool voting/consensus).

---

## 4) Deployability Matters (HPC reality check)

From the CDC-list assessment under **UGA Sapelo2** constraints (no sudo; Apptainer security restrictions):
- **Runnable and assessed (5 from CDC list)**: nf-core/taxprofiler, GOTTCHA2, sourmash, CLARK, nf-core/mag
- **Not runnable (deployment failure recorded; not benchmarked)**: CZID (IDseq), SURPI+, DHO Lab, NAO MGS, TaxTriage  

**Conclusion**: *Deployability and reproducibility are first-order selection criteria*, not optional “engineering details”.

---

## 5) Key Findings (Sensitivity vs Specificity, and what to trust)

### 5.1 Sensitivity–Specificity trade-off (high-level)

| Tool/Workflow | Default sensitivity | Default specificity | Specificity after consensus/threshold | Notes |
|---|---|---|---|---|
| taxprofiler | High | Medium | High | Broad screening; requires thresholds + evidence to control false positives |
| GOTTCHA2 | Medium–Low | High | Medium | Conservative; strong as a confirmatory tool |
| sourmash | Medium | Medium–High | Medium | Good for containment/consistency validation |
| CLARK (CLARK-S) | Medium | High | Medium | Strong for in-database targets; long reads prefer CLARK-S + validation |
| nf-core/mag | Not detection-focused | — | High (via MAG QC) | Value is genome/MAG evidence chains (CheckM/GUNC/coverage, etc.) |
| rvdb-viral | High | Medium | High | Strong discovery; consensus sets are most reliable |
| MLMVD-nf | Low | Very high | (built-in) | Designed as strict filter: fewer calls, higher confidence |
| KrakenMetaReads-nf | Medium | Medium | Medium–High | Assembly-first reduces ambiguity; may lose low-abundance taxa |

### 5.2 Resource/engineering implications
- **Resource-heavy (CPU/time)**: nf-core/mag, rvdb-viral, MLMVD-nf (assembly and multi-tool inference dominate).
- **Memory-heavy**: Kraken2/CLARK family (DB index size is the main driver).
- **Lightweight**: sourmash (sketching and search are relatively small/fast).
- **Production maturity**: nf-core + Nextflow workflows generally offer better batch execution, logging, and resume/checkpointing than standalone tools.

### 5.3 “Database is the invisible deciding factor”
- RefSeq-like standard DBs often underrepresent environmental virus/phage diversity.
- RVDB improves viral coverage (especially environmental/distant viruses) but still has annotation bias risks.

---

## 6) False Positive Control (Minimum recommended evidence chain)

Common false-positive sources: host contamination, kitome/reagent contamination, short-fragment homology, DB misannotation, low-complexity/repeats.

### 6.1 Minimum essential strategy (strongly recommended)
1. **Host removal**: Bowtie2 for short reads; Minimap2 for long reads.
2. **Threshold filtering**: minimum read counts / relative abundance thresholds for read-level outputs (scale with sequencing depth).

### 6.2 Recommended multi-evidence consensus
- **Multi-tool support (≥2 tools)**: e.g., taxprofiler Kraken2 plus an additional independent signal (Kaiju/DIAMOND module if enabled, and/or GOTTCHA2/CLARK).
- **Assembly evidence priority for key findings**: require at least one medium-to-long contig (e.g., >5–10 kb; adjust by virus type) with supportive evidence:
  - hallmark genes / gene-feature plausibility
  - DIAMOND/RVDB or VirSorter2/CheckV support
  - coverage distribution consistency (avoid single-spot artifacts)

### 6.3 Interpretation rule
Anything that is **single-tool only**, **very low abundance**, and **lacks assembly/coverage evidence** should be labeled:
**“Low-confidence candidate (requires validation)”**, not a confirmed finding.

---

## 7) Decision Matrix (Scores 1–5; from the assessment)

| Dimension | taxprofiler | rvdb-viral | MLMVD-nf | nf-core/mag | GOTTCHA2 | CLARK | sourmash | KrakenMetaReads-nf |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Known virus sensitivity (default) | 5 | 4 | 2 | 2 | 3 | 3 | 3 | 3 |
| Distant/environmental virus discovery | 3 | 5 | 3 | 4 | 1 | 1 | 2 | 3 |
| Specificity (default) | 3 | 3 | 5 | 4* | 5 | 4 | 4 | 3 |
| Result reliability (interpretable evidence chain) | 4 | 4 | 4 | 5 | 4 | 3 | 4 | 4 |
| Analysis speed | 3 | 1 | 2 | 1 | 4 | 5 | 5 | 3 |
| Memory efficiency | 2 | 3 | 3 | 2 | 4 | 2 | 5 | 2 |
| Genome reconstruction capability | 1 | 5 | 2 | 5 | 1 | 1 | 1 | 4 |
| Engineering maturity (HPC reproducible) | 5 | 5 | 4 | 5 | 3 | 3 | 4 | 5 |

\* nf-core/mag’s high specificity is mainly due to MAG QC evidence chains (CheckM/GUNC/coverage), not read-level classification.

---

## 8) Recommended Deployable Solutions (Decision-oriented)

### 8.1 High-confidence consensus + benchmark-style reporting
- **Recommended stack**: taxprofiler (primary screen) + (key targets) GOTTCHA2 / sourmash / assembly evidence
- **Use when**: wastewater/environmental monitoring; pre-clinical screening; you need clear “screen-positive vs confirm-positive” labeling.

### 8.2 Emerging/outbreak discovery (discovery & evidence chain)
- **Recommended stack**: rvdb-viral (assembly + protein homology) + MLMVD-nf (strict candidate set) + optional nf-core/mag (deep reconstruction)
- **Use when**: environmental unknown virus surveys; outbreak tracing requiring high-quality contig/genome evidence.

### 8.3 Wastewater surveillance (high background/noise)
- **Recommended stack**: taxprofiler (batch) → GOTTCHA2 confirmation → (key samples) KrakenMetaReads-nf / rvdb-viral for assembly validation
- **Reporting rule**: explicitly separate **screen-positive** from **confirm-positive**.

### 8.4 Ultra-diverse, large-scale monitoring (cost control)
- **Recommended stack**: sourmash (fast coarse screening/dedup) + taxprofiler (fine screening) + (key samples) nf-core/mag / rvdb-viral deep mining

---

## 9) Appendix: Terms (short)

- **Sensitivity**: avoid misses (more calls).
- **Specificity**: avoid false alarms (fewer, cleaner calls).
- **Consensus/threshold output**: apply rules (multi-tool, thresholds, coverage/assembly evidence) to raise specificity.
- **NCLDV**: Nucleocytoviricota (large dsDNA viruses; e.g., Mimiviridae-related signals often require contig/gene evidence).
- **Contig**: assembled contiguous sequence.
- **MAG**: metagenome-assembled genome (after binning).

---

**Notes**
- For full details (tool-by-tool dependency stacks, failure reasons for non-runnable items, and empirical observations), see `Metagenomic Viral Detection Tool Assessment Report.md`.
- For the full per-pipeline overview (as originally summarized), see `PIPELINES_SUMMARY_EN.md`.

