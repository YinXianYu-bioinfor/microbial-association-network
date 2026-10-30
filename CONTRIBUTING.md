# 贡献指南

欢迎通过 GitHub issue 讨论问题或提出改进，并在提交前说明预期行为与复现步骤。请不要提交个人数据、真实样本数据、测序文件或本地运行结果；数据许可不确定时使用人工生成的最小示例。

在项目根目录执行：

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python cli.py --counts examples/abundance.tsv --metadata examples/metadata.tsv --taxonomy-path examples/taxonomy.tsv --output demo_output
```

`demo_output` 必须不存在，已有结果请先检查并换个目录。统计相关的改动应说明检验假设、参数以及使用了什么数据验证；请保持 GUI 与 CLI 的核心行为一致。提交内容受项目 MIT 许可证约束。
