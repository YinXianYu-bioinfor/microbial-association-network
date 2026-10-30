"""Command-line entry point."""
import argparse
import sys
from app_version import VERSION
from network_core import analyze, EXPLORATION_NOTICE, POOLED_NOTICE


def main():
    p = argparse.ArgumentParser(description='微生物丰度 Spearman 关联网络')
    p.add_argument('--version', action='version', version=f'network-tool {VERSION}')
    p.add_argument('--counts', required=True, help='特征×样本 TSV/CSV')
    p.add_argument('--output', required=True, help='输出目录，默认必须不存在')
    p.add_argument('--existing-output', action='store_true', help='输出到已有目录；拒绝覆盖同名结果')
    p.add_argument('--metadata', help='可选样本元数据 TSV/CSV')
    p.add_argument('--use-all-samples', choices=['y', 'yes', 'n', 'no'], default='no', type=str.lower, metavar='{y,yes,n,no}', help='合并全部样本构建全局网络（默认 no）')
    p.add_argument('--taxonomy-path', help='可选 QIIME 2 taxonomy.tsv')
    p.add_argument('--sample-column', default='sample-id')
    p.add_argument('--group-column', default='group')
    p.add_argument('--min-prevalence', type=float, default=0.5)
    p.add_argument('--min-abundance', type=float, default=0.001)
    p.add_argument('--top-n', type=int, default=50, help='每组按平均相对丰度保留前 N 个节点（默认 50）')
    p.add_argument('--pseudocount', type=float, default=0.5)
    p.add_argument('--rho-min', type=float, default=0.7)
    p.add_argument('--min-stability', type=float, default=0.8)
    p.add_argument('--p-max', type=float, help='可选：候选边原始 p 值必须小于此阈值；需 --independent')
    p.add_argument('--q-max', type=float, help='可选：候选边 BH q 值必须小于此阈值；需 --independent')
    p.add_argument('--permutations', type=int, default=9999, help='超过 9 个样本时的随机置换次数；≤9 使用全部排列')
    p.add_argument('--repeats', type=int, default=36)
    p.add_argument('--seed', type=int, default=2026)
    p.add_argument('--max-features', type=int, default=200)
    p.add_argument('--independent', action='store_true', help='确认样本独立且置换所需可交换性成立')
    args = vars(p.parse_args())
    args['use_all_samples'] = args['use_all_samples'] in ('y', 'yes')
    print(EXPLORATION_NOTICE, file=sys.stderr)
    if args['use_all_samples']:
        print(POOLED_NOTICE, file=sys.stderr)
    try:
        report = analyze(counts_path=args.pop('counts'), output_dir=args.pop('output'), metadata_path=args.pop('metadata'), **args)
    except (ValueError, OSError) as e:
        print(f'分析失败：{e}', file=sys.stderr)
        return 1
    for g in report['groups']:
        print(f"{g['group']}: {g['samples']} 样本，{g['nodes']} 节点，{g['candidate_edges']} 条候选边")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
