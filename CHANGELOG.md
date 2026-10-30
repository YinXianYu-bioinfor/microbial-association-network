# 更新记录 / Changelog

## 0.1.0 — 未发布 / Unreleased

计划作为首个探索性预发布版本，尚未创建 Git 标签或 GitHub Release。

Planned first exploratory pre-release; no Git tag or GitHub Release has been created.

- GUI 和 CLI 使用同一计算核心，支持分组及全部样本合并分析。
- 输入筛选、组内丰度前 N 节点、CLR、Spearman 关联及边表、节点表和分类表导出。
- 可选置换检验与 p/q 筛选；4–9 个样本使用精确置换，更多样本使用随机置换。
- 留一和子样本稳定性指标、参数和输入哈希记录；新增软件版本记录。
- 提供中英文文档、人工示例和自动测试，以及附来源说明的真实数据案例输入。

The release includes a shared CLI/GUI analysis core, group and pooled modes, feature filtering and top-N selection, CLR and Spearman correlations, optional permutation p/q filtering, stability metrics, input and software-version reporting, bilingual documentation, synthetic examples, and automated tests.

**限制 / Limitations:** 统计关联不证明直接生物互作；合并组间差异可能驱动关联。置换检验要求适当的独立性和可交换性，尚不支持分层或配对置换。高维、小样本及随机置换精度问题仍需进一步验证。

Statistical associations do not prove direct interactions. Pooled group differences can drive associations. Permutation inference requires appropriate assumptions; stratified and paired permutations are not supported. High dimensionality, small-sample behavior, and Monte Carlo precision require further validation.
