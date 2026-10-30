# 小样本扩增子统计关联网络探索

**中文** | [English](README.en.md)

**当前开发版本：`0.1.0`**。计划首次开源时使用 Git 标签 `v0.1.0` 并标记为预发布版；目前尚未发布或创建标签。版本变更记录见 [CHANGELOG.md](CHANGELOG.md)。

这是一个免费、本地运行的探索性工具：从 ASV、OTU 或其他特征的丰度表计算 Spearman 关联，生成可审查的节点表和边表。项目重点关注扩增子数据的小样本应用，也可分析其他满足输入要求的数据。

> 统计关联用于提出研究假设，不能直接证明微生物互作或因果。小样本、高维特征、组成效应及实验设计均可能影响结果；本工具仍处于开发和验证阶段，不保证特定数据的假阳性率。

## 文件结构

```text
network_tool/
├── cli.py / gui.py / network_core.py      # 命令行、图形界面和共享分析核心
├── app_version.py / CHANGELOG.md          # 软件版本和变更记录
├── requirements.txt / environment.yml    # Python 与 Conda 依赖
├── README.md / README.en.md               # 中英文说明
├── LICENSE / CONTRIBUTING.md              # 许可与贡献指南
├── .gitignore / .gitattributes           # Git 忽略规则及数据文件换行策略
├── .github/workflows/tests.yml            # GitHub Actions 测试配置
├── docs/                                  # 数据来源、样本清单、许可说明和发布清单
│   ├── prjna1137280_t20_samples.tsv
│   ├── DATA_PROVENANCE.md
│   ├── THIRD_PARTY_NOTICES.md
│   └── RELEASE_CHECKLIST.md
├── examples/                              # 人工构造的示例输入
├── tests/                                 # 自动测试
└── test_data/                             # 真实数据案例输入及来源说明
```

`test_data/` 包含三份真实数据案例输入；原始 FASTQ 和本地运行结果不随源码发布。提交前仍应审查实际文件清单。

## 安装与快速开始

需要 Python 3.10 或更高版本、NumPy、pandas；图形界面还需要 Tkinter 和图形显示环境。已有环境满足要求时可直接使用；如需专用 Conda 环境，在项目目录执行：

```bash
conda env create -f environment.yml
conda activate microbial-network
```

也可在其他 Python 环境中执行 `python -m pip install -r requirements.txt`。在 Windows 中若 `python` 不可用，可尝试 `py -3`。项目自带人工构造示例，不代表真实微生物群落：

```bash
python cli.py --counts examples/abundance.tsv --metadata examples/metadata.tsv --taxonomy-path examples/taxonomy.tsv --output demo_output
```

默认按组分析，不计算 p/q，输出目录必须尚不存在；再次运行请换一个名称。运行 `python gui.py` 可启动图形界面。WSL 2 配置 WSLg 后可运行 GUI；Linux 需要可用图形会话，无显示的服务器可以直接使用 CLI，不强制使用 Conda。参见 [Microsoft WSL 图形应用指南](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gui-apps)。

通过 `python cli.py --version` 查看软件版本；新生成的 `report.json` 也会记录 `software_version`。此前已生成的报告不会自动添加版本字段。

## 输入

- **丰度表**：特征为行、样本为列，首列为唯一特征 ID，表头为唯一样本 ID；TSV 或 CSV。支持 QIIME 2 BIOM 导出文件开头的 `# Constructed from biom file` 注释。值应为非负有限数，每个样本总丰度大于零。
- **元数据（可选）**：TSV/CSV，默认样本列 `sample-id`、分组列 `group`，可通过 CLI 修改；样本必须与丰度表一一匹配。不提供时所有样本作为一组。
- **分类注释（可选）**：QIIME 2 导出的 TSV，包含唯一 `Feature ID` 和 `Taxon`；可选 `Confidence`。未匹配的注释字段为空。

工具将每个样本输入值除以样本总和；若输入已做抽平、转换或其他处理，使用者需判断这种归一化是否适用。每个分析组至少需要 4 个样本。

## GUI 和 CLI

GUI 依次选择丰度表、可选元数据和分类表，设置参数及输出目录，然后点击“开始分析”。默认新建输出目录；也可选已有目录，但如报告或分组目录同名则拒绝覆盖。运行失败可能留下不完整目录，请检查后更换目录重试。

```bash
# 分组分析；p/q 留空，仅使用相关系数和稳定性筛选
python cli.py --counts abundance.tsv --metadata metadata.tsv --taxonomy-path taxonomy.tsv --output grouped_results

# 合并全部样本，默认不计算 p/q
python cli.py --counts abundance.tsv --metadata metadata.tsv --output pooled_results --use-all-samples yes

# 在置换假设成立时，以原始 p 值筛边；q 值仍计算并输出
python cli.py --counts abundance.tsv --metadata metadata.tsv --output p_results --independent --p-max 0.05
```

如需按 BH 校正后的值筛选，可加入 `--q-max 0.05`；也可以同时设置 p 和 q。p/q 阈值为严格小于，且须启用 `--independent`。这项确认由用户根据实验设计判断，软件不能自动验证。GUI 的“合并全部样本”默认不勾选；CLI 使用 `--use-all-samples y` 或 `yes` 开启，`n` 或 `no` 关闭。

## 参数与计算

| 参数（CLI） | 默认值 | 作用 |
|---|---:|---|
| `--sample-column` / `--group-column` | `sample-id` / `group` | 元数据中的样本及分组列名 |
| `--use-all-samples` | `no` | `yes` 或 `y` 合并全部样本 |
| `--min-prevalence` | 0.5 | 特征在组内非零样本的最小比例 |
| `--min-abundance` | 0.001 | 特征在组内的最小平均相对丰度，即 0.1% |
| `--top-n` | 50 | 筛选后按组内平均相对丰度保留前 N 个节点；最少 3 个 |
| `--max-features` | 200 | 可设置的节点上限；N 不可超过此值 |
| `--pseudocount` | 0.5 | CLR 前加到 CPM 的伪计数；当前仅 CLI 可修改 |
| `--rho-min` | 0.7 | 候选边的 `|ρ|` 下限，包含等于阈值 |
| `--min-stability` | 0.8 | 重选节点后恢复频率的下限，包含等于阈值 |
| `--p-max` / `--q-max` | 不筛选 | 候选边的原始 p 或 BH q 上限，严格小于；可分别或同时设置 |
| `--permutations` | 9999 | 样本数超过 9 时的随机置换次数 |
| `--repeats` | 36 | 子样本重复次数；0 表示不运行子样本分析 |
| `--seed` | 2026 | 随机种子；置换与重采样使用独立随机序列 |
| `--independent` | 关闭 | 用户确认适用的置换假设后计算 p/q |
| `--existing-output` | 关闭 | 使用已有输出目录，但不覆盖同名结果 |

节点先经过出现率和平均丰度筛选，再按组内平均相对丰度排序取前 N；丰度相同时按特征 ID 排序。未达到 N 时保留合格节点，CLR 后的常数特征会被移除。程序以该组每个样本的总丰度计算 CPM，对选中节点加伪计数并做 CLR，再计算 Spearman 秩相关；CLR 的参照取决于选中的节点。

确认置换条件后，4–9 样本使用所有排列计算双侧精确 p 值；更多样本使用随机置换 `(b+1)/(B+1)`，并对组内全部节点对做 BH 校正。较小的随机 p 值可能受置换次数影响，p/q 的解释还取决于零假设、数据筛选和检验依赖。p/q 阈值留空时不以它们筛边；若确认条件则仍计算并输出，未确认则 p/q 为空。

稳定性包括留一和子样本评估。每个子集取 `max(3, floor(0.8 × 样本数))` 个样本；组合数不超过指定次数时使用全部组合，否则随机抽取（组合可能重复）。重选稳定性会重新进行节点筛选、前 N 排序和 CLR，缺席或无法计算的边按未恢复计入全部子集的分母；频率不是边为真的概率。`--repeats 0` 时没有子样本恢复频率，候选边表为空。

**分组与置换限制**：合并样本得到的是全局关联，组间差异可能驱动相关，不证明在每组内部均稳定。样本独立不自动保证置换所需的可交换性；配对、重复测量、批次、时间和处理组结构需要单独评估。当前不支持分层、配对或受限置换，合并模式的子样本也不按组分层。

## 输出

每个分析组拥有独立子目录；全局模式目录为 `all_samples`。`report.json` 记录软件版本、输入哈希、参数、组别、原分组样本数和 NumPy/pandas 版本。开发阶段版本可能继续变化，运行者应同时保留对应代码或提交记录。

| 文件 | 内容 |
|---|---|
| `edges.csv` | 全部节点对及 `source`、`target`、Spearman `rho`、`p_value`、`q_value`、留一和重采样指标 |
| `candidate_edges.csv` | 满足 `|ρ|`、重选稳定性及用户启用的 p/q 条件的边；不是已证实互作 |
| `nodes.csv` | 最终节点 ID 和分类注释 |
| `selected_taxonomy.csv` | 最终节点的分类、置信度和组内平均相对丰度 |

`rho` 可作为带符号边权重；需要非负权重的可视化软件可使用 `|rho|` 并单独保留符号。CSV 可导入 Gephi 等软件。空候选边表是有效结果，不会自动放宽阈值。

## 真实数据案例：PRJNA1137280

本项目使用 Hoefle 等人的公开数据进行功能与方法探索：Hoefle D, Sperber S, Marzban N, Antonietti M, Hoffmann T, Wicaksono WA, Abdelfattah A, Berg G. 2024. *Artificial Humic Acid Diminishes the Effect of Drought on the Soil Microbiome*. Journal of Sustainable Agriculture and Environment 3(4):e70034. [doi:10.1002/sae2.70034](https://doi.org/10.1002/sae2.70034)。原始数据项目：[PRJNA1137280](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1137280)。

我们选取作者元数据中的 `t20`，共 36 个样本，四组各 9 个：[36 个样本及 SRR accession 对照表](docs/prjna1137280_t20_samples.tsv)。

| 处理 | 样本 | 数量 |
|---|---|---:|
| 腐殖酸 + 干旱 | `S_5A1`–`S_5A9` | 9 |
| 腐殖酸 + 正常供水 | `S_5B1`–`S_5B9` | 9 |
| 无腐殖酸 + 干旱 | `S_5C1`–`S_5C9` | 9 |
| 无腐殖酸 + 正常供水 | `S_5D1`–`S_5D9` | 9 |

原文报告使用 515F/806R 扩增 16S rRNA V4 区，以 Illumina NovaSeq 6000 进行 2 × 150 bp 测序，并在 QIIME 2 中按单端序列处理。**本项目同样采用单端策略，使用 R1 进行 DADA2 去噪，不进行双端拼接。** 原文正文未明确说明所用方向，因此不将 R1 描述为已证实与原文相同。130 bp 截断及下表参数为本项目的实际设置，不代表完全复现原文。

| 步骤 | 本项目实际设置 |
|---|---|
| 去引物 | QIIME 2 2025.7 / Cutadapt 5.1；两端 5′ 锚定，前向 `^GTGYCAGCMGCCGCGGTAA`，反向 `^GGACTACNVGGGTWTCTAAT`；错误率 0.1，丢弃未匹配 reads，4 核 |
| DADA2 | `denoise-single`，R1；`trim-left=0`、`trunc-len=130`、`max-ee=2`、`trunc-q=2`，独立去噪，`consensus` 去嵌合体，16 线程；DADA2 1.30.0 |
| 错误模型 | `learn_min_reads=1,000,000`；日志实际使用 2 个样本的 1,814,611 条 reads、235,899,430 个碱基 |
| 物种注释 | QIIME 2 `classify-sklearn`，SILVA 138.2 SSURef NR99 全长 uniform 分类器；置信度 0.7，方向 `same`，4 进程 |
| 网络输入 | 根据分类注释排除叶绿体和线粒体；45,721 个 ASV 中排除 1,026 个，保留 44,695 个；未抽平，未分类和古菌保留 |

上述去引物序列是本项目实际使用的简并序列，不应在缺少进一步证据时断言其为原研究 PCR 引物的精确序列。原文报告 BLAST 分类及基于阴性对照的污染识别；本项目案例使用 `classify-sklearn`，未在该案例中完成同样的污染识别。130 bp 单端序列的分类分辨率有限。

`test_data/` 包含本案例的派生丰度表、元数据和分类表；原始 FASTQ 不随仓库发布。`examples/` 是人工数据。相关来源、处理及许可参见 [DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md) 与 [test_data/README.md](test_data/README.md)。

## 参考与许可

- Bolyen E, et al. 2019. QIIME 2. *Nature Biotechnology* 37:852–857. [doi:10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9)。
- Callahan BJ, et al. 2016. DADA2. *Nature Methods* 13:581–583. [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)。
- Quast C, et al. 2013. The SILVA rRNA gene database project. *Nucleic Acids Research* 41:D590–D596. [doi:10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)。SILVA 数据使用条件见 [SILVA license](https://www.arb-silva.de/silva-license-information/)。

项目开发与维护：**Xianyu Yin**；开发协作：**Codex（OpenAI）**。代码采用 [MIT License](LICENSE)；第三方依赖及数据来源见 [THIRD_PARTY_NOTICES.md](docs/THIRD_PARTY_NOTICES.md)，贡献方式见 [CONTRIBUTING.md](CONTRIBUTING.md)。在项目目录运行 `python -m unittest discover -s tests -v` 执行自动测试。
