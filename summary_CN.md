# 宏基因组病毒检测总结（Pipelines + 工具评估整合版）

## 1）执行摘要（什么时候用什么）

- **广谱筛查 + 标准化丰度表（适合生产/批量）**：优先用 **MetaTaxProfiler（nf-core/taxprofiler）**，并可叠加二次验证（GOTTCHA2 / sourmash / 组装证据）。
- **偏发现（discovery）+ 基因组/contig 证据链（尤其适合大 DNA 病毒 / NCLDV）**：优先用 **rvdb-viral-metagenome-nf**（组装 + DIAMOND/RVDB 蛋白同源）。当需要“更小但更可靠”的候选集时，再叠加 **MLMVD-nf** 做严格筛选。
- **先组装再分类（contig 级 Kraken2/Bracken）**：优先用 **KrakenMetaReads-nf**；contig 证据更易解释，但组装阶段可能丢失极低丰度类群。
- **对已组装 contigs 做快速相似性筛查**：优先用 **sourmash**。
- **非 Kraken2 体系的确认型 profiling**：**GOTTCHA2**（基于 unique signatures，较保守）与 **CLARK/CLARK-S**（判别性 k-mer；对“数据库内覆盖的目标”效果好，但通常需要后续验证）。

> **核心原则**（来自评估报告）：不要把“检测列表”当作结论。病毒检测应被解释为一条 **证据链**（覆盖度、contig/组装支撑、多工具一致性、阈值策略等）。

---

## 2）本文件夹 8 个 Pipeline 总览（速览表）

| Pipeline | 主要目标 | 数据类型 | 是否组装 | 核心方法 | 关键数据库 | 主要输出 |
|---|---|---|---|---|---|---|
| **Sourmash** | MinHash 相似性/包含性快速检索 | 组装后的 contigs | 外部组装（MetaFlye/SPAdes/MEGAHIT） | `sketch` + `gather` | NCBI viruses（k=31） | `.sig`、`gather_results*.csv` |
| **rvdb-viral-metagenome-nf** | 组装 + 蛋白同源发现；双对比策略 | 短读 + 长读 | ✅ | Prodigal + DIAMOND vs RVDB + taxonomy 对比 | RVDB + NCBI taxonomy（长读还用 Pfam-A/viralFlye） | 分类报告；共识/差异列表；可选 RPM/RPKM |
| **mag-nf（nf-core/mag 方案）** | 组装 + Kraken2（reads + contigs） | 短读 + 长读 | ✅ | Kraken2（reads/contigs） | Kraken2 viral DB | `results_short/`、`results_long/` |
| **TaxProfiler-nf（MetaTaxProfiler）** | 标准化 profiling + 自动 RPM/RPKM | 短读 + 长读 | ❌ | Kraken2（短读可选 Bracken）+ 报表 | Kraken2 DB（+ Bracken DB） | `abundance/` 表格 + 报告 |
| **MLMVD-nf** | 新/远缘病毒挖掘：多工具共识 | 短读 + 长读 | ✅ | VirSorter2 + DeepVirFinder（长读可加 viralFlye） | VS2 DB + DVF 模型 + Pfam-A | 共识分层；高置信列表；可选 RPM/RPKM |
| **KrakenMetaReads-nf** | 先组装后分类 + RPM/RPKM | 短读 + 长读 | ✅ | Kraken2（短读可加 Bracken）跑 contigs | Kraken2 DB（+ Bracken DB） | abundance 目录；merged reports |
| **GOTTCHA2** | 基于 unique signatures 的保守 profiling | 短读 + 长读 | ❌ | GOTTCHA2 signature profiling | GOTTCHA2 `.mmi` | `*.tsv`、`*.summary.tsv`、`*.full.tsv` |
| **CLARK** | 判别性 k-mer 快速分类 + 丰度 | 短读 + 长读 | ❌ | CLARK / CLARK-S / CLARK-l | RefSeq（可做 viruses-only targets） | 分类 CSV + 丰度脚本输出 |

---

## 3）评估框架（为什么不同工具结果会不一样）

为避免“既最敏感又最特异”等自相矛盾结论，评估明确区分两类输出：
- **默认输出（Default output）**：工具按常用/默认参数跑出的标准结果（通常更敏感，但噪声更大）。
- **共识/阈值输出（Consensus/threshold output）**：在默认输出基础上叠加 **多工具交叉验证**、**去宿主**、**覆盖度/contig 阈值** 等规则得到的结果集（通常特异性更高，但敏感性更低）。

同时将策略划分为 4 类（便于选型组合）：
1. **直接读段分类（read-level）**：taxprofiler（以及直接跑 Kraken2/CLARK/GOTTCHA2）。
2. **组装后做 contig 分类**：KrakenMetaReads-nf；以及对 nf-core/mag 产物的 contig 解释。
3. **蛋白/比对证据链**：rvdb-viral（DIAMOND/RVDB）等。
4. **病毒特征 + 多工具共识筛选**：MLMVD-nf（多工具投票/交集共识）。

---

## 4）可部署性很重要（HPC 现实检查）

在 **UGA Sapelo2**（无 sudo；Apptainer 安全限制）环境下对 CDC 清单工具的部署测试结论：
- **可运行并进入评估（CDC 清单中 5 个）**：nf-core/taxprofiler、GOTTCHA2、sourmash、CLARK、nf-core/mag
- **无法运行（记录失败原因，未进入性能对比）**：CZID（IDseq）、SURPI+、DHO Lab、NAO MGS、TaxTriage

**结论**：*可部署性与可复现性是一级筛选条件*，不是可有可无的“工程细节”。

---

## 5）关键发现（敏感性 vs 特异性：哪些更可信）

### 5.1 敏感性–特异性权衡（高层概览）

| 工具/工作流 | 默认敏感性 | 默认特异性 | 共识/阈值后特异性 | 备注 |
|---|---|---|---|---|
| taxprofiler | 高 | 中 | 高 | 广谱筛查；需要阈值 + 证据链来控假阳 |
| GOTTCHA2 | 中–低 | 高 | 中 | 保守；很适合作为确认工具 |
| sourmash | 中 | 中–高 | 中 | 适合 containment/一致性验证 |
| CLARK（CLARK-S） | 中 | 高 | 中 | 对库内目标强；长读建议 CLARK-S + 验证 |
| nf-core/mag | 非检测导向 | — | 高（来自 MAG QC） | 价值在基因组/MAG 证据链（CheckM/GUNC/coverage 等） |
| rvdb-viral | 高 | 中 | 高 | 发现能力强；共识集更可靠 |
| MLMVD-nf | 低 | 很高 | （策略内置） | 严格过滤器：更少但更可信 |
| KrakenMetaReads-nf | 中 | 中 | 中–高 | 先组装降低歧义；可能丢低丰度 |

### 5.2 资源与工程代价
- **重资源（CPU/时间）**：nf-core/mag、rvdb-viral、MLMVD-nf（主要耗在组装/多工具推断）。
- **吃内存**：Kraken2/CLARK 家族（数据库索引大小是主因）。
- **轻量**：sourmash（sketch 与检索相对小且快）。
- **生产成熟度**：nf-core + Nextflow 工作流通常比 standalone 工具更适合批量、日志、断点续跑与可复现生产化。

### 5.3 “数据库是隐形的决定因素”
- RefSeq 类标准库往往低估环境病毒/噬菌体多样性。
- RVDB 提升病毒覆盖（尤其环境/远缘病毒），但依然存在注释偏倚风险。

---

## 6）假阳性控制（最低建议证据链）

常见假阳性来源：宿主污染、试剂污染（kitome）、短片段同源、数据库误注释、低复杂度/重复序列等。

### 6.1 最低必备策略（强烈推荐）
1. **去宿主（Host removal）**：短读用 Bowtie2；长读用 Minimap2。
2. **阈值过滤**：对 read-level 结果设置最小 reads 数 / 相对丰度阈值（阈值需随测序深度调整）。

### 6.2 推荐的“多证据链共识”
- **多工具支持（≥2 工具）**：例如 taxprofiler 的 Kraken2 结果，再叠加一个独立证据（如启用 Kaiju/DIAMOND 模块，或 GOTTCHA2/CLARK 的支持）。
- **关键发现优先要求组装证据**：至少需要一个中长 contig（例如 >5–10 kb；随病毒类型调整）并具备支持证据：
  - hallmark genes / 基因特征合理性
  - DIAMOND/RVDB 或 VirSorter2/CheckV 支撑
  - 覆盖分布一致性（避免“单点覆盖伪迹”）

### 6.3 解读规则
任何 **仅单工具命中**、**极低丰度**、且 **缺乏组装/覆盖证据** 的条目，应标注为：
**“低置信候选（需要验证）”**，而不是“已确认发现”。

---

## 7）决策矩阵（1–5 分；来自评估报告）

| 维度 | taxprofiler | rvdb-viral | MLMVD-nf | nf-core/mag | GOTTCHA2 | CLARK | sourmash | KrakenMetaReads-nf |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 已知病毒敏感性（默认） | 5 | 4 | 2 | 2 | 3 | 3 | 3 | 3 |
| 远缘/环境病毒发现能力 | 3 | 5 | 3 | 4 | 1 | 1 | 2 | 3 |
| 特异性（默认） | 3 | 3 | 5 | 4* | 5 | 4 | 4 | 3 |
| 结果可靠性（可解释证据链） | 4 | 4 | 4 | 5 | 4 | 3 | 4 | 4 |
| 速度 | 3 | 1 | 2 | 1 | 4 | 5 | 5 | 3 |
| 内存效率 | 2 | 3 | 3 | 2 | 4 | 2 | 5 | 2 |
| 基因组重建能力 | 1 | 5 | 2 | 5 | 1 | 1 | 1 | 4 |
| 工程成熟度（HPC 可复现） | 5 | 5 | 4 | 5 | 3 | 3 | 4 | 5 |

\* nf-core/mag 的高特异性主要来自 MAG QC 证据链（CheckM/GUNC/coverage），而不是 read-level 分类本身。

---

## 8）可部署推荐方案（面向决策）

### 8.1 高置信共识 + 类 benchmark 报告输出
- **推荐组合**：taxprofiler（主筛）+（关键目标）GOTTCHA2 / sourmash / 组装证据
- **适用场景**：污水/环境监测；临床前筛查；需要清晰区分“筛查阳性 vs 确认阳性”。

### 8.2 新发/暴发发现（discovery + 证据链）
- **推荐组合**：rvdb-viral（组装 + 蛋白同源）+ MLMVD-nf（严格候选集）+（可选）nf-core/mag（深入重建）
- **适用场景**：环境未知病毒调查；暴发溯源中需要快速产出高质量 contig/基因组证据链。

### 8.3 污水监测（高背景噪声）
- **推荐组合**：taxprofiler（批量筛）→ GOTTCHA2 二次确认 →（关键样本）KrakenMetaReads-nf / rvdb-viral 做组装验证
- **报告规则**：明确区分 **screen-positive** 与 **confirm-positive**。

### 8.4 超高多样性大规模监测（成本控制）
- **推荐组合**：sourmash（快速粗筛/去冗余）+ taxprofiler（精筛）+（关键样本）nf-core/mag / rvdb-viral 深挖

---

## 9）附录：术语（简版）

- **敏感性（Sensitivity）**：尽量不漏检（更“多报”）。
- **特异性（Specificity）**：尽量不误报（更“少但干净”）。
- **共识/阈值输出**：在默认结果上叠加多工具、阈值、覆盖/组装证据等规则以提高特异性。
- **NCLDV**：Nucleocytoviricota（大型 dsDNA 病毒；例如 Mimiviridae 相关信号通常需要 contig/基因证据支撑）。
- **Contig**：组装得到的连续序列片段。
- **MAG**：宏基因组组装基因组（binning 后的基因组草图）。

---

**备注**
- 需要更完整的工具依赖栈、不可运行条目的失败原因、以及经验观察，请参考 `Metagenomic Viral Detection Tool Assessment Report.md`。
- 需要每条 pipeline 的更细致总览（最初的整理），请参考 `PIPELINES_SUMMARY_EN.md`。

