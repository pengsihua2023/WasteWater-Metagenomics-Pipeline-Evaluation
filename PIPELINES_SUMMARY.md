# 8 个病毒宏基因组 Pipeline 总结（总览版）

本文汇总当前文件夹中 8 个 README（每个对应一个 pipeline/工具链）的核心目标、输入输出、关键步骤、主要工具与数据库，以及适用场景，便于快速选型与对比。

## 总览对比表

| Pipeline | 主要目标 | 支持数据类型 | 是否组装 | 核心分类/鉴定方法 | 关键数据库 | 主要输出 |
|---|---|---|---|---|---|---|
| **Sourmash** | 基于 MinHash 的快速相似性检索（contigs vs 参考库） | 组装后的 contigs（长读/短读均可） | 依赖外部组装结果（MetaFlye/SPAdes/MEGAHIT） | `sourmash sketch` + `sourmash gather` | NCBI viruses（k=31，zip） | `.sig` 签名；`gather_results*.csv` 命中表 |
| **rvdb-viral-metagenome-nf** | 病毒蛋白相似性鉴定 +（短读双组装对比 / 长读双轨道对比） | Illumina（短读）+ Nanopore/PacBio（长读） | ✅（短读：MEGAHIT+SPAdes；长读：MetaFlye） | Prodigal 预测蛋白 + Diamond BLASTP（vs RVDB）+ 完整分类谱系对比 | RVDB 蛋白库（Diamond index）+ NCBI taxonomy；长读还用 Pfam-A（viralFlye） | 7 层级分类对比报告；双轨道共识/差异列表；可选 RPM/RPKM |
| **mag-nf（nf-core/mag 方案）** | 组装 + Kraken2 分类（读段与 contigs） | 短读（nf-core/mag）；长读（直跑容器） | ✅（短读：MEGAHIT+SPAdes；长读：Flye） | Kraken2（读段 + contigs） | Kraken2 viral DB | `results_short/` 与 `results_long/`（含 kraken2 report、contigs 分类结果等） |
| **TaxProfiler-nf（MetaTaxProfiler）** | 标准化分类 + 自动计算病毒丰度（RPM/RPKM） | Illumina + Nanopore/PacBio | ❌（以分类为主） | Kraken2；短读可选 Bracken 校正；自动汇总与报表 | Kraken2 DB（短读建议配 Bracken DB） | `results_viral_short/abundance/` 或 `results_viral_long/abundance/`（含汇总表、Top viruses） |
| **MLMVD-nf（ML 增强新病毒发现）** | 面向“新/远缘病毒”的多工具并行鉴定与交叉验证 | 短读 + 长读 | ✅（短读：MEGAHIT+SPAdes；长读：metaFlye） | VirSorter2（混合 ML）+ DeepVirFinder（深度学习）+ viralFlye（Pfam 功能域验证，长读）+ 三工具交集分层 | VirSorter2 DB；DeepVirFinder 模型；Pfam-A（viralFlye） | 三工具对比报告（1/2/3-tool 共识分层）；高置信列表；可选 RPM/RPKM |
| **KrakenMetaReads-nf** | nf-core/taxprofiler 为核心的分类 +（可选组装/多组装）+ 丰度 | 短读 + 长读 | ✅（短读：MEGAHIT+SPAdes；长读：Flye+ViralFlye） | Kraken2；短读可配 Bracken；自动批处理 + RPM/RPKM；可输出 BIOM | Kraken2 DB +（短读）Bracken DB | `results_viral_short/` 与 `results_viral_long/`（含 abundance_*、merged reports 等） |
| **GOTTCHA2** | 用 GOTTCHA2 做宏基因组分类/病毒分类 | 短读 + 长读 | ❌（read-based 为主） | GOTTCHA2（k-mer signature profiling） | GOTTCHA2 DB（`.mmi`） | `*.tsv` / `*.summary.tsv` / `*.full.tsv` |
| **CLARK** | 以判别性 k-mer 做快速分类（含病毒专用 DB 构建与丰度估计） | 单端/双端均可（更偏 reads） | ❌（工具本身不负责组装） | CLARK / CLARK-l / CLARK-S 分类；自带 abundance 估计脚本或自定义脚本 | NCBI RefSeq（可构建 viruses 目标库） | 分类结果 CSV；可用脚本生成丰度表 |

## 每个 Pipeline 的要点总结

### 1) `sourmash README.md`（Sourmash Metagenome Analysis Pipeline）
- **定位**：把 contigs 做成 MinHash 签名（sketch），再去参考库里做包含性检索（gather），适合“快速粗筛/相似性检索”。
- **输入**：不同组装器产出的 contigs（示例：MetaFlye/SPAdes/MEGAHIT）。
- **关键参数**：
  - sketch：k=21/31/51，scaled=1000，开启 abundance（`abund`）。
  - gather：使用 k=31（与数据库一致），阈值 300 bp。
- **输出**：
  - `.sig`：签名文件。
  - `gather_results*.csv`：含 overlap（`intersect_bp`）、containment ANI、abundance 等。
- **适用场景**：已有组装结果，想快速判断“更像哪些已知病毒参考序列/库条目”，或者做快速比较与筛查。

### 2) `rvdb-viral-metagenome-nf README.md`（RVDB + 双策略病毒宏基因组）
- **定位**：面向“病毒蛋白相似性”鉴定，用 RVDB（蛋白库）+ Diamond；短读用 **双组装器对比**，长读用 **viralFlye 特征过滤 + Diamond 两轨并行对比**，强调覆盖率与互补性。
- **短读流程**：fastp →（MEGAHIT ∥ SPAdes）→ Prodigal → Diamond vs RVDB → 7 层分类对比 →（可选）BWA 计算 RPM/RPKM。
- **长读流程**：MetaFlye → viralFlye（Pfam + viralVerify）→ 两轨并行：
  - Track1：全 contigs → Prodigal → Diamond → taxonomy
  - Track2：viralFlye 过滤后的病毒 contigs → Prodigal → Diamond → taxonomy
  - 输出共识/差异（如“Diamond 找到但 viralFlye 过滤掉”的远缘候选）。
- **输出亮点**：短读按 Kingdom→Species 的并列对比报告；长读给出共识病毒、MetaFlye-only（远缘候选）等分层文件。
- **适用场景**：要做“短读组装器差异”与“长读特征法 vs 相似性法互补”分析，且分类统计要完整可解释。

### 3) `mag-nf README.md`（nf-core/mag + 自定义 contig 分类；长读直跑）
- **定位**：把 nf-core/mag（短读）当作稳定的工程化流水线，并补上它默认缺失的 **contigs 分类**；长读则用 Flye + Kraken2 走更直接的方案。
- **短读**：fastp、PhiX 去除（Bowtie2）、MEGAHIT+SPAdes、读段 Kraken2（nf-core/mag 自带）+ contigs Kraken2（自定义）、MultiQC。
- **长读**：Flye 组装后对 contigs 与 raw reads 分别跑 Kraken2。
- **适用场景**：希望“短读走 nf-core 标准化、同时又要 contigs 级别分类”；长读追求简单快速落地。

### 4) `TaxProfiler-nf README.md`（MetaTaxProfiler：分类到标准化丰度输出）
- **定位**：把“Kraken2 分类 +（短读可选 Bracken 校正）+ 自动算 RPM/RPKM + 报表”打包成一键化工作流，产出标准化的丰度表。
- **输入**：samplesheet（短读/长读格式统一），数据库配置 `databases.csv`。
- **输出**：
  - `abundance/` 下的单样本详细表、全样本汇总表、Top viruses（默认 RPM≥10）。
- **适用场景**：主要目标是“得到可直接用于统计/画图的病毒丰度矩阵”，而不是深度的新病毒发现或复杂交叉验证。

### 5) `MLMVD-nf README.md`（机器学习增强的新病毒发现 + 多工具共识）
- **定位**：专门面向“新/远缘病毒发现”，用 **VirSorter2 + DeepVirFinder**（短读/长读）再加 **viralFlye（Pfam 功能域验证，长读）**，通过交集分层给出不同置信度集合。
- **长读（核心亮点）**：metaFlye →（VS2 ∥ DVF ∥ viralFlye）并行 → 三工具对比：
  - 3-tool 共识：最高置信
  - 2-tool 共识：中等置信
  - 单工具：探索性集合
- **短读**：MEGAHIT+SPAdes 并行组装；对每个组装结果跑 VS2 与 DVF，并做跨组装器比较。
- **适用场景**：研究重点是“尽可能发现更多病毒候选，并用独立方法交叉验证以提高可信度”，尤其适合探索新环境样本的病毒“暗物质”。

### 6) `KrakenMetaReads-nf README.md`（TaxProfiler 核心的 Kraken2/Bracken + 丰度）
- **定位**：以 nf-core/taxprofiler 为核心，把短读/长读的分类、（可选组装与多组装）与 RPM/RPKM 丰度计算标准化，并支持批处理与（描述中提到）BIOM 输出。
- **短读**：MEGAHIT/SPAdes 产物分别 Kraken2（可结合 Bracken）→ `abundance_megahit/` 与 `abundance_spades/`。
- **长读**：Flye + ViralFlye（区分 circular/linear）后分别 Kraken2 → 对应 abundance 目录。
- **适用场景**：想要“Kraken2/Bracken 体系下的标准化分类 + 丰度”，同时希望短读/长读都用一套 Nextflow 工程化方式管理与批处理。

### 7) `GOTTCHA2 README.md`（GOTTCHA2 分类脚本）
- **定位**：用 GOTTCHA2 做 read-based 的分类/病毒分类，脚本化（SLURM）提交，偏“快速 profile”。
- **输入**：短读 paired-end FASTQ 或长读 single FASTQ。
- **输出**：`.tsv`、`.summary.tsv`、`.full.tsv`。
- **适用场景**：需要一条相对轻量的分类路径，且已有 GOTTCHA2 DB 与运行习惯。

### 8) `CLARK README.md`（CLARK 安装与病毒分类指南）
- **定位**：CLARK 系列（含 CLARK-l/CLARK-S）基于判别性 k-mer 的快速分类；文档侧重“安装、构建 viruses 目标库、分类与丰度估计”。
- **关键点**：
  - 可仅构建 viruses 目标库，降低 DB 规模与时间成本。
  - k-mer 长度影响敏感性：k=20/21 更敏感（病毒检测常用），k=31 更均衡。
  - 丰度可用自带 `estimate_abundance.sh` 或自定义脚本计算。
- **适用场景**：希望用 CLARK 做高速分类，或需要基于 RefSeq 的自建病毒库流程参考。

## 选型建议（怎么挑）

- **只想要标准化病毒丰度表（RPM/RPKM）**：优先 `TaxProfiler-nf` 或 `KrakenMetaReads-nf`（Kraken2/Bracken 体系，输出结构更“统计友好”）。
- **想做“短读双组装器对比”或“长读双轨互补（特征法 vs 相似性法）”**：优先 `rvdb-viral-metagenome-nf`。
- **想最大化发现新/远缘病毒，并分层输出置信度集合**：优先 `MLMVD-nf`（VS2 + DVF + viralFlye 三证据交叉）。
- **已有组装 contigs，想快速做相似性检索/粗筛**：优先 `Sourmash`。
- **想走 k-mer profile/快速分类（非 Kraken2 体系）**：看 `GOTTCHA2`（签名 profile）或 `CLARK`（判别性 k-mer + 自建 viruses 目标库）。

