# Metagenomic Viral Detection Tool Assessment
## Presentation for CDC Expert Review (NWSS)

**Author**: Sihua Peng  
**Assessment environment**: UGA Sapelo2 (HPC cluster; no sudo; Apptainer restrictions)  
**Last updated**: 2026-01-28  

---

## Slide 1 — Title

# Metagenomic Viral Detection Tool Assessment

- Assessment scope: CDC-listed tools + deployable alternatives/workflows
- Focus: deployability, reproducibility, evidence-chain interpretability, false-positive control

---

## Slide 2 — Goals (What this review answers)

- **Can it run on HPC without sudo?** (deployment + reproducibility)
- **What evidence chain does it provide?** (read hits vs contigs vs protein homology vs multi-tool consensus)
- **How does it behave on short vs long reads?**
- **How to control false positives in wastewater-style noisy samples?**
- **What is the recommended “stack” for different operational scenarios?**

---

## Slide 3 — Assessment Framework (Avoiding contradictory conclusions)

To prevent statements like “most sensitive and most specific” at the same time:

- **Default output**: tool’s standard/default result set  
  - typically higher sensitivity, higher noise
- **Consensus/threshold output**: default output after applying rules  
  - multi-tool cross-validation, host removal, coverage/contig thresholds, etc.  
  - typically higher specificity, lower sensitivity

---

## Slide 4 — Assessment Objects (8 runnable items)

### CDC list (runnable; 5)
- **nf-core/taxprofiler (MetaTaxProfiler)**
- **GOTTCHA2**
- **sourmash**
- **CLARK (CLARK / CLARK-S / CLARK-l)**
- **nf-core/mag**

### Author-developed workflows (3)
- **rvdb-viral-metagenome-nf**
- **MLMVD-nf**
- **KrakenMetaReads-nf**

---

## Slide 5 — Non-runnable CDC items (deployment blockers)

| Software | Primary blocker on HPC (no sudo / Apptainer policy) |
|---|---|
| CZID (IDseq) | cloud-first; local path not well adapted to cluster constraints |
| SURPI+ | requires sudo/container build patterns incompatible with HPC policy |
| DHO Lab | sudo/container + complex dependency chain |
| NAO MGS | sudo/container + complex dependency chain |
| TaxTriage | sudo/container + complex dependency chain |

**Key point**: deployability is a first-order selection criterion for public-health workflows.

---

## Slide 6 — Two runs: datasets and scaling (simulated vs real)

### Run #1 (Simulated CDC dataset)
- Short reads: `llnl_66ce4dde_R1.fastq.gz` + `llnl_66ce4dde_R2.fastq.gz`  
  - 670.52 MB + 693.43 MB = **1.36 GB**
- Long reads: `llnl_66d1047e.fastq.gz` = **263.93 MB**

### Run #2 (Real dataset)
- Short reads: `SRR35987572_1.fastq.gz` + `SRR35987572_1.fastq.gz`*  
  - 145 GB + 144 GB = **289 GB**
- Long reads: **not available**
- Scale: **~212.5×** larger than Run #1 (short-read volume)

\* Filename kept verbatim from notes.

---

## Slide 7 — Empirical compute footprint (8 pipelines, two runs)

| Run | Data type | Relative scale | Wall time (all 8 pipelines) | Peak memory (max across 8) | Main memory driver |
|---|---|---:|---|---:|---|
| #1 | Simulated (CDC) | 1× | ≤ 48 hours (≤ 2 days) | 256 GB | overall cap; assembly steps typically peak |
| #2 | Real | 212.5× | 2–7 days | 768 GB | genome assembly software dominates |

**Operational implication**: as input scales, **assembly becomes the critical memory/time bottleneck**.

---

## Slide 8 — Four strategy archetypes (how tools generate evidence)

1. **Direct read classification**
   - taxprofiler; also direct Kraken2/CLARK/GOTTCHA2 usage
2. **Contig classification after assembly**
   - KrakenMetaReads-nf; contig interpretation from nf-core/mag outputs
3. **Protein/alignment evidence chain**
   - rvdb-viral (DIAMOND vs RVDB), protein-level discovery
4. **Viral feature + multi-tool consensus screening**
   - MLMVD-nf (VirSorter2 + DeepVirFinder + optional viralFlye voting)

---

## Slide 9 — Pipelines catalog (what is in this repository)

| Pipeline | Core purpose | Assembly | Core method (1 line) |
|---|---|---|---|
| Sourmash | fast similarity/containment screening | external | MinHash sketch + gather vs NCBI viruses |
| MetaTaxProfiler | standardized profiling + abundance | no | Kraken2 (+ optional Bracken) + reporting |
| mag-nf (nf-core/mag) | genome reconstruction | yes | assembly + bins/MAG QC; optional taxonomy |
| KrakenMetaReads-nf | assembly-before-classification | yes | Kraken2/Bracken on contigs + RPM/RPKM |
| rvdb-viral-metagenome-nf | discovery via protein homology | yes | Prodigal + DIAMOND vs RVDB + taxonomy |
| MLMVD-nf | strict consensus candidate set | yes | VS2 + DVF (+ viralFlye long-read) voting |
| GOTTCHA2 | conservative confirmatory profiling | no | unique-signature profiling |
| CLARK | fast k-mer classification | no | discriminative k-mers (CLARK-S for long reads) |

---

## Slide 10 — taxprofiler (MetaTaxProfiler): why it is the default “front door”

- **Best role**: broad screening + rapid overview + standardized outputs
- **Strength**: engineering maturity (nf-core/Nextflow; batch + logging + resume)
- **Risk**: long-tail low-read hits need thresholds + validation (not “confirmed”)
- **Best practice**: treat taxprofiler as **screen-positive generator**, not final truth

---

## Slide 11 — GOTTCHA2 / CLARK / sourmash: confirmatory tools (different flavors)

- **GOTTCHA2**: high-specificity unique signatures  
  - strong confirmation; conservative (may miss distant/low-abundance targets)
- **CLARK / CLARK-S**: discriminative k-mers  
  - fast; strong for in-database targets; long reads prefer CLARK-S + validation
- **sourmash**: containment-style MinHash metrics  
  - lightweight; good for contig-level consistency/coverage-style validation

---

## Slide 12 — Assembly-centric workflows: nf-core/mag vs KrakenMetaReads-nf

### nf-core/mag
- **Not a detection-first tool**: primary value is **MAG reconstruction + QC evidence chain**
- Strong when you need genome-level interpretation (coverage, completeness, contamination)

### KrakenMetaReads-nf
- **Assembly-before-classification** (contig-level Kraken2/Bracken)
- Benefits: longer contigs reduce ambiguity; good for large-virus positioning
- Risk: **low-abundance taxa may be lost during assembly**

---

## Slide 13 — Discovery-centric workflow: rvdb-viral-metagenome-nf

- **Core objective**: maximize discovery, especially environmental/distant viruses
- Evidence chain: **assembly → Prodigal → DIAMOND protein search vs RVDB → taxonomy**
- Strong for: long contigs + large DNA viruses (e.g., NCLDV-related signals)
- Best practice: prioritize **consensus sets** (dual-track overlap) over single-path hits

---

## Slide 14 — Strict candidate set workflow: MLMVD-nf

- **Positioning**: high-specificity screening via multi-tool voting
- Long-read mode: **VirSorter2 ∥ DeepVirFinder ∥ viralFlye** (Pfam validation) → consensus tiers
- Benefit: small, high-confidence list for follow-up validation
- Cost: reduced sensitivity (may miss true but weakly supported contigs)

---

## Slide 15 — Key finding: sensitivity vs specificity trade-off (high-level)

| Tool/Workflow | Default sensitivity | Default specificity | What improves specificity |
|---|---|---|---|
| taxprofiler | High | Medium | thresholds + host removal + multi-evidence validation |
| GOTTCHA2 | Medium–Low | High | signature coverage confirmation |
| sourmash | Medium | Medium–High | containment thresholds; contig context |
| CLARK (CLARK-S) | Medium | High | best for in-db targets; follow-up validation |
| rvdb-viral | High | Medium | consensus sets (track overlap) + additional validation |
| MLMVD-nf | Low | Very high | built-in multi-tool consensus |
| KrakenMetaReads-nf | Medium | Medium | assembly evidence + (optional) Bracken + thresholds |
| nf-core/mag | not detection-first | — | MAG QC evidence chain |

---

## Slide 16 — False positive control: minimum required safeguards

**Do not treat a detection list as a conclusion.**

Minimum essentials:
1. **Host removal**
   - short reads: Bowtie2
   - long reads: Minimap2
2. **Threshold filtering**
   - minimum read count / relative abundance thresholds

Interpretation rule:
- single-tool only + very low abundance + no coverage/assembly evidence ⇒ **Low-confidence candidate**

---

## Slide 17 — Multi-evidence chain: what “confirmed” should look like

Recommended for key findings (especially public-health reporting):

- **≥2 independent tools** support the target, *and/or*
- **Assembly evidence** (≥5–10 kb contig; virus-type dependent) with:
  - plausible gene features / hallmark genes
  - DIAMOND/RVDB or VirSorter2/CheckV support
  - reasonable coverage distribution (avoid single-spot artifacts)

Labeling recommendation:
- **screen-positive** / **supportive** / **confirmed**

---

## Slide 18 — Decision matrix (scores 1–5; assessment summary)

| Dimension | taxprofiler | rvdb-viral | MLMVD-nf | nf-core/mag | GOTTCHA2 | CLARK | sourmash | KrakenMetaReads-nf |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Known virus sensitivity (default) | 5 | 4 | 2 | 2 | 3 | 3 | 3 | 3 |
| Distant/environmental discovery | 3 | 5 | 3 | 4 | 1 | 1 | 2 | 3 |
| Specificity (default) | 3 | 3 | 5 | 4* | 5 | 4 | 4 | 3 |
| Evidence-chain interpretability | 4 | 4 | 4 | 5 | 4 | 3 | 4 | 4 |
| Speed | 3 | 1 | 2 | 1 | 4 | 5 | 5 | 3 |
| Memory efficiency | 2 | 3 | 3 | 2 | 4 | 2 | 5 | 2 |
| Genome reconstruction | 1 | 5 | 2 | 5 | 1 | 1 | 1 | 4 |
| Engineering maturity (HPC) | 5 | 5 | 4 | 5 | 3 | 3 | 4 | 5 |

\* nf-core/mag specificity is primarily driven by MAG QC evidence chains (CheckM/GUNC/coverage), not read-level classification.

---

## Slide 19 — Recommended deployable solutions (decision-oriented)

### A) High-confidence consensus reporting (benchmark-style)
- **taxprofiler** (primary screen)
- confirm key targets with **GOTTCHA2 / sourmash / assembly evidence**

### B) Discovery & outbreak tracing (evidence-first)
- **rvdb-viral** (discovery + protein homology evidence chain)
- add **MLMVD-nf** for strict candidate set
- optional **nf-core/mag** for deep reconstruction/ecology/host association

### C) Wastewater surveillance (high noise)
- taxprofiler batch → GOTTCHA2 confirmation → key samples to KrakenMetaReads-nf / rvdb-viral
- report **screen-positive vs confirm-positive**

---

## Slide 20 — Final conclusions

- **No single tool meets all needs**: use layered strategy  
  **rapid screening → confirmation → deep evidence**
- **Deployability is a hard gate** for public-health operationalization
- **Assembly dominates compute at scale** (memory/time bottleneck)
- **False-positive control is essential**: thresholds + host removal + multi-evidence chain

**Thank you. Q&A**

