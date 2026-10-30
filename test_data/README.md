# 真实数据测试集

本项目使用 PRJNA1137280 第 20 天的 36 个样本（四组各 9 个）进行测试。以下文件由本项目对公开数据重新处理或整理后提供，不是原始测序文件：

- `asv_counts_no_organelle.tsv`：DADA2 单端 R1 截断至 130 bp 后的 ASV 原始计数，已根据分类注释去除标注为叶绿体或线粒体的特征；44,695 个特征。
- `sample_metadata.tsv`：样本与分组信息，TSV 字段外围的双引号是标准格式，程序可正常解析。
- `taxonomy.tsv`：基于 SILVA 138.2 分类的 ASV 注释，130 bp 序列的分类分辨率有限。

原始研究数据参见 [NCBI BioProject PRJNA1137280](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1137280) 和 Hoefle 等（2024），[doi:10.1002/sae2.70034](https://doi.org/10.1002/sae2.70034)。样本 accession、去引物和去噪等处理参数见项目 README。丰度表是 R1 单端去噪后再移除分类标为叶绿体或线粒体的结果；元数据是按 t20 选择和整理的子集；分类表由代表序列经 `classify-sklearn` 分类生成。

分类注释基于 **SILVA 138.2**；SILVA 数据库及分类信息由 SILVA 提供，适用 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)，详见 [SILVA 许可说明](https://www.arb-silva.de/silva-license-information/)。本项目对数据进行的处理包括分类、样本选择及特征过滤；分类输出不是原始 SILVA 数据库文件。本项目软件的 MIT 许可不改变第三方数据的权利和使用条件。

这些文件供功能测试和探索，不能作为真实微生物互作的验证标准。公开序列数据可再分发并不意味着所有相关内容均适用同一许可证；再使用者应保留来源与署名，并检查数据及参考资料的适用条款。
