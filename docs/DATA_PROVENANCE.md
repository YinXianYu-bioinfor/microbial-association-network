# 数据来源与公开范围

`examples/` 中的数据和分类名称是专门为软件使用示例人工构造的，不代表实际物种、样本或研究结论。可用于检查输入格式及程序是否运行。

本地开发时使用过公开项目 [PRJNA1137280](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1137280) 的处理结果。相关论文：Hoefle 等，*Artificial Humic Acid Diminishes the Effect of Drought on the Soil Microbiome*, DOI [10.1002/sae2.70034](https://doi.org/10.1002/sae2.70034)。测试结果还使用 SILVA 138.2 的分类信息；详见 [SILVA 许可说明](https://www.arb-silva.de/silva-license-information/) 和 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)。

`test_data/` 提供基于上述公开数据处理的丰度表、样本元数据和 SILVA 分类结果；本地运行结果和原始 FASTQ 不随源码发布。原始研究数据由 Hoefle 等人提供，请引用论文及 BioProject。`taxonomy.tsv` 是使用 [SILVA 138.2 分类器](https://www.arb-silva.de/archive/release_138_2/QIIME2/2025.7/SSU/full-length/uniform)生成的注释结果；SILVA 数据库与分类信息由 SILVA 提供，适用 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)，须保留来源、许可链接及所做修改的说明。分类器、处理步骤及版本见项目 README 和 `test_data/README.md`。

公开数据可访问不意味着本项目代码的 MIT 许可证涵盖这些数据；再使用者应核实各来源的适用条款。详见 [NCBI 对公开序列数据的说明](https://www.ncbi.nlm.nih.gov/sra/docs/sequence-data-processing/)及 [SILVA 许可说明](https://www.arb-silva.de/silva-license-information/)。
