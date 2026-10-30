"""Input validation for the local association-network tool."""
from pathlib import Path
import csv
import itertools
import hashlib
import json
import math
import re
from functools import lru_cache
import numpy as np
import pandas as pd
from app_version import VERSION

EXPLORATION_NOTICE = ('本工具面向扩增子等数据的统计关联探索。小样本与大量特征可能使多重校正后的 q 值较大、检出能力降低；'
                      '结果不保证均为真阳性，也不能直接证明生物学互作，请结合实验设计与后续证据谨慎解读。')
POOLED_NOTICE = ('合并全部样本的关联可能受处理组差异影响，不代表各组内均稳定出现；'
                 '合并后的样本可交换性也需单独判断，本工具不支持分层置换。')


def load_inputs(counts_path, metadata_path=None, sample_column='sample-id', group_column='group'):
    path = Path(counts_path)
    delimiter = ',' if path.suffix.lower() == '.csv' else '\t'
    with path.open(encoding='utf-8-sig', newline='') as f:
        lines = f.readlines()
    if lines and lines[0].startswith('# Constructed from biom file'):
        lines = lines[1:]
    rows = list(csv.reader(lines, delimiter=delimiter))
    if len(rows) < 2 or len(rows[0]) < 3:
        raise ValueError('丰度表需要特征列与至少两个样本')
    header = [x.strip() for x in rows[0]]
    samples = header[1:]
    if any(not x for x in samples) or len(set(samples)) != len(samples):
        raise ValueError('样本名称为空或重复')
    if any(len(row) != len(header) for row in rows[1:]):
        raise ValueError('丰度表各行列数不一致')
    features = [r[0].strip() for r in rows[1:]]
    if any(not x for x in features) or len(set(features)) != len(features):
        raise ValueError('特征名称为空或重复')
    try:
        values = np.array([[float(x) for x in r[1:]] for r in rows[1:]], dtype=float)
    except ValueError as e:
        raise ValueError('丰度表包含缺失或非数值数据') from e
    if not np.isfinite(values).all() or (values < 0).any() or (values.sum(axis=0) == 0).any():
        raise ValueError('丰度必须为非负有限数，且样本总丰度不能为零')
    counts = pd.DataFrame(values.T, index=samples, columns=features)
    if metadata_path is None:
        meta = pd.DataFrame({'Group': ['all'] * len(samples)}, index=samples)
    else:
        metadata_path = Path(metadata_path)
        meta = pd.read_csv(metadata_path, sep=',' if metadata_path.suffix.lower() == '.csv' else '\t', dtype=str, keep_default_na=False)
        if sample_column not in meta or group_column not in meta:
            raise ValueError('元数据缺少指定的样本列或分组列')
        meta[sample_column] = meta[sample_column].str.strip()
        meta[group_column] = meta[group_column].str.strip()
        if (meta[sample_column] == '').any() or meta[sample_column].duplicated().any() or (meta[group_column] == '').any():
            raise ValueError('元数据样本或分组为空，或样本重复')
        meta = meta.set_index(sample_column)
        if set(meta.index) != set(samples):
            raise ValueError('元数据样本与丰度表样本不一致')
        meta = meta.loc[samples].rename(columns={group_column: 'Group'})
    return counts, meta


def rank_correlation(features_by_samples):
    values = np.asarray(features_by_samples, dtype=float)
    ranks = np.asarray([pd.Series(row).rank(method='average').to_numpy() for row in values])
    ranks -= ranks.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranks, axis=1)
    unit = np.divide(ranks, norms[:, None], out=np.full_like(ranks, np.nan), where=norms[:, None] > 0)
    return np.clip(unit @ unit.T, -1, 1)


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if not len(p):
        return p.copy()
    if not np.isfinite(p).all() or (p < 0).any() or (p > 1).any():
        raise ValueError('p 值必须在 0 到 1 之间')
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty(len(p)); result[order] = np.minimum(q, 1)
    return result


def exact_permutation_p(x, y):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    if len(x) != len(y) or not 3 <= len(x) <= 9:
        raise ValueError('精确置换要求两向量等长且样本数为 3–9')
    ranks = np.asarray([pd.Series(v).rank(method='average').to_numpy() for v in (x, y)])
    ranks -= ranks.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(ranks, axis=1)
    if (norm == 0).any():
        raise ValueError('常数特征无法计算相关性')
    ranks /= norm[:, None]
    observed = abs(ranks[0] @ ranks[1])
    permutations = np.asarray(list(itertools.permutations(range(len(x)))), dtype=int)
    null = ranks[1][permutations] @ ranks[0]
    return float(np.mean(np.abs(null) >= observed - 1e-12))


@lru_cache(maxsize=1)
def _permutation_indices(n):
    return np.fromiter((i for p in itertools.permutations(range(n)) for i in p), dtype=np.uint8).reshape(-1, n)


@lru_cache(maxsize=512)
def _exact_distribution(a, b):
    # Average ranks multiplied by two are integers, including ties.
    a = np.asarray(a, dtype=np.int32)
    b = np.asarray(b, dtype=np.int32)
    scores = np.abs(b[_permutation_indices(len(a))] @ a)
    values, counts = np.unique(scores, return_counts=True)
    return values, np.cumsum(counts[::-1])[::-1]


def exact_pairwise_p(features_by_samples, i, j):
    x = np.asarray(features_by_samples, dtype=float)
    n = x.shape[1]
    if not 3 <= n <= 9:
        raise ValueError('精确置换仅支持 3–9 个样本')
    ranks = np.asarray([pd.Series(row).rank(method='average').to_numpy() for row in x])
    ranks = np.rint(2 * ranks - (n + 1)).astype(np.int32)
    if (np.sum(ranks * ranks, axis=1) == 0).any():
        raise ValueError('常数特征无法计算相关性')
    patterns = [tuple(sorted(row)) for row in ranks]
    result = np.empty(len(i), dtype=float)
    for k, (left, right) in enumerate(zip(i, j)):
        a, b = sorted((patterns[left], patterns[right]))
        values, tail_counts = _exact_distribution(a, b)
        observed = abs(int(ranks[left] @ ranks[right]))
        result[k] = tail_counts[np.searchsorted(values, observed)] / math.factorial(n)
    return result


def group_folder_name(group):
    name = re.sub(r'[^A-Za-z0-9_-]', '_', str(group)) or 'group'
    if name != group:
        name += '_' + hashlib.sha256(str(group).encode()).hexdigest()[:8]
    return name


def select_top_features(sub, min_prevalence, min_abundance, top_n):
    relative = sub.div(sub.sum(axis=1), axis=0)
    means = relative.mean(axis=0)
    keep = ((sub > 0).sum(axis=0) >= math.ceil(len(sub) * min_prevalence)) & (means >= min_abundance)
    features = sorted(means.index[keep], key=lambda feature: (-means[feature], feature))[:top_n]
    return relative, features


def analyze(counts_path, output_dir, metadata_path=None, sample_column='sample-id', group_column='group',
            min_prevalence=0.5, min_abundance=0.001, pseudocount=0.5, rho_min=0.7,
            permutations=9999, repeats=36, seed=2026, independent=False, max_features=200, taxonomy_path=None, min_stability=0.8, existing_output=False, top_n=50, p_max=None, q_max=None, use_all_samples=False):
    if not 0 <= min_prevalence <= 1 or min_abundance < 0 or pseudocount <= 0:
        raise ValueError('筛选阈值或伪计数无效')
    if not 0 <= rho_min <= 1 or not 0 <= min_stability <= 1 or permutations < 1 or repeats < 0:
        raise ValueError('统计参数无效')
    if not isinstance(top_n, int) or top_n < 3 or max_features < 3 or top_n > max_features:
        raise ValueError('保留节点数 N 至少为 3，且不能超过计算节点上限；如需更多节点，请同时调高上限')
    for label, value in (('p', p_max), ('q', q_max)):
        if value is not None and (not math.isfinite(value) or not 0 < value <= 1):
            raise ValueError(f'{label} 值筛选阈值必须大于 0 且不超过 1')
    if not independent and (p_max is not None or q_max is not None):
        raise ValueError('使用 p/q 筛选前，需确认样本独立且满足置换检验所需条件')
    out = Path(output_dir)
    if existing_output and not out.is_dir():
        raise FileNotFoundError(f'指定的已有输出目录不存在：{out}')
    if not existing_output and out.exists():
        raise FileExistsError(f'输出目录已存在：{out}')
    counts, meta = load_inputs(counts_path, metadata_path, sample_column, group_column)
    original_group_counts = {str(k): int(v) for k, v in meta['Group'].value_counts().sort_index().items()}
    if use_all_samples:
        meta = meta.copy()
        meta['Group'] = 'all_samples'
    taxonomy = {}
    if taxonomy_path:
        tax = pd.read_csv(taxonomy_path, sep='\t', dtype=str, keep_default_na=False)
        if 'Feature ID' not in tax or 'Taxon' not in tax or tax['Feature ID'].duplicated().any():
            raise ValueError('分类注释表须含唯一 Feature ID 和 Taxon 列')
        if 'Confidence' not in tax:
            tax['Confidence'] = ''
        taxonomy = dict(zip(tax['Feature ID'], zip(tax['Taxon'], tax['Confidence'])))
    if existing_output:
        targets = [out / 'report.json'] + [out / group_folder_name(g) for g in meta['Group'].unique()]
        for target in targets:
            if target.exists() or target.is_symlink():
                raise FileExistsError(f'输出路径已存在，无法覆盖：{target}')
    else:
        out.mkdir(parents=True)
    rng = np.random.default_rng(seed)
    permutation_rng = np.random.default_rng(np.random.SeedSequence([seed, 1]))
    reports = []
    for group in sorted(meta['Group'].unique()):
        sub = counts.loc[meta['Group'] == group]
        n = len(sub)
        if n < 4:
            raise ValueError(f'{group} 样本不足 4 个')
        relative, chosen = select_top_features(sub, min_prevalence, min_abundance, top_n)
        selected = sub.loc[:, chosen]
        if selected.shape[1] > max_features:
            raise ValueError(f'{group} 有 {selected.shape[1]} 个节点，超过计算上限 {max_features}；请调整筛选或上限')
        if selected.shape[1] < 3:
            raise ValueError(f'{group} 保留节点少于 3，无法进行可靠的 CLR 分析')
        cpm = selected.div(sub.sum(axis=1), axis=0) * 1e6
        x = np.log(cpm.to_numpy().T + pseudocount)
        x -= x.mean(axis=0, keepdims=True)
        valid = np.std(x, axis=1) > 1e-12
        x = x[valid]; features = selected.columns.to_numpy()[valid]
        if len(features) < 3:
            raise ValueError(f'{group} 可变节点少于 3')
        rho = rank_correlation(x)
        i, j = np.triu_indices(len(features), 1)
        observed = rho[i, j]
        rank = np.asarray([pd.Series(row).rank(method='average').to_numpy() for row in x])
        rank -= rank.mean(axis=1, keepdims=True)
        rank /= np.linalg.norm(rank, axis=1, keepdims=True)
        exceed = np.zeros(len(i), dtype=int)
        if independent and n <= 9:
            pvalues = exact_pairwise_p(x, i, j)
            qvalues = bh_adjust(pvalues)
            permutation_method, permutations_used = 'exact', math.factorial(n)
        elif independent:
            for _ in range(permutations):
                p = permutation_rng.permutation(n)
                null = (rank @ rank[:, p].T)[i, j]
                exceed += np.abs(null) >= np.abs(observed) - 1e-12
            pvalues = (exceed + 1) / (permutations + 1)
            qvalues = bh_adjust(pvalues)
            permutation_method, permutations_used = 'monte_carlo', permutations
        else:
            pvalues = qvalues = np.full(len(i), np.nan)
            permutation_method, permutations_used = 'not_run', 0
        loo = np.array([rank_correlation(np.delete(x, k, axis=1))[i, j] for k in range(n)])
        loo_valid = np.isfinite(loo)
        loo_sign = (loo_valid & (np.sign(loo) == np.sign(observed))).sum(axis=0) / n
        loo_delta = np.max(np.where(loo_valid, np.abs(loo - observed), 0), axis=0)
        loo_delta[~loo_valid.any(axis=0)] = np.nan
        hits = np.zeros(len(i), dtype=int)
        reselected_hits = np.zeros(len(i), dtype=int)
        reselected_eligible = np.zeros(len(i), dtype=int)
        size = max(3, math.floor(n * .8))
        combinations = list(itertools.combinations(range(n), size)) if math.comb(n, size) <= repeats else None
        subsets = combinations if combinations is not None else [rng.choice(n, size, replace=False) for _ in range(repeats)]
        for indices in subsets:
            r = rank_correlation(x[:, indices])[i, j]
            hits += np.isfinite(r) & (np.sign(r) == np.sign(observed)) & (np.abs(r) >= rho_min)
            trial = sub.iloc[list(indices)]
            trial_rel, trial_features = select_top_features(trial, min_prevalence, min_abundance, top_n)
            if len(trial_features) < 3:
                continue
            trial_cpm = trial_rel.loc[:, trial_features].to_numpy().T * 1e6
            log_trial = np.log(trial_cpm + pseudocount)
            center = log_trial.mean(axis=0)
            available = np.isin(features, trial_features)
            common = features[available]
            values = np.log(trial_rel.loc[:, common].to_numpy().T * 1e6 + pseudocount) - center
            if len(common) < 2:
                continue
            rr = rank_correlation(values)
            positions = np.full(len(features), -1, dtype=int)
            positions[np.flatnonzero(available)] = np.arange(len(common))
            e = (positions[i] >= 0) & (positions[j] >= 0)
            target = np.flatnonzero(e)
            vals = rr[positions[i[target]], positions[j[target]]]
            finite = np.isfinite(vals)
            reselected_eligible[target[finite]] += 1
            reselected_hits[target[finite]] += (np.sign(vals[finite]) == np.sign(observed[target[finite]])) & (np.abs(vals[finite]) >= rho_min)
        frequency = hits / len(subsets) if subsets else np.full(len(i), np.nan)
        edges = pd.DataFrame({'source': features[i], 'target': features[j], 'rho': observed,
                              'p_value': pvalues, 'q_value': qvalues, 'stability_hits': hits,
                              'stability_subsets': len(subsets), 'stability_frequency': frequency,
                              'loo_sign_fraction': loo_sign, 'loo_max_delta': loo_delta,
                              'reselected_eligible': reselected_eligible, 'reselected_hits': reselected_hits,
                              'reselected_frequency': reselected_hits / len(subsets) if subsets else np.full(len(i), np.nan)})
        edges['candidate'] = (np.abs(observed) >= rho_min) & (edges['reselected_frequency'] >= min_stability)
        if p_max is not None:
            edges['candidate'] &= edges['p_value'] < p_max
        if q_max is not None:
            edges['candidate'] &= edges['q_value'] < q_max
        folder = out / group_folder_name(group)
        folder.mkdir()
        edges.to_csv(folder / 'edges.csv', index=False)
        edges.loc[edges['candidate']].to_csv(folder / 'candidate_edges.csv', index=False)
        node_table = pd.DataFrame({'feature_id': features, 'taxonomy': [taxonomy.get(f, ('', ''))[0] for f in features]})
        node_table.to_csv(folder / 'nodes.csv', index=False)
        pd.DataFrame({'feature_id': features, 'taxonomy': [taxonomy.get(f, ('', ''))[0] for f in features],
                      'confidence': [taxonomy.get(f, ('', ''))[1] for f in features],
                      'mean_relative_abundance': [relative[f].mean() for f in features]}).to_csv(folder / 'selected_taxonomy.csv', index=False)
        reports.append({'group': group, 'samples': n, 'samples_ids': sub.index.tolist(),
                        'raw_features': sub.shape[1], 'nodes': len(features), 'pairs': len(i),
                        'candidate_edges': int(edges['candidate'].sum()), 'subsets': len(subsets),
                        'permutation_method': permutation_method, 'permutations_used': permutations_used})
    report = {'software_version': VERSION, 'input_sha256': hashlib.sha256(Path(counts_path).read_bytes()).hexdigest(),
              'metadata_sha256': hashlib.sha256(Path(metadata_path).read_bytes()).hexdigest() if metadata_path else None,
              'taxonomy_sha256': hashlib.sha256(Path(taxonomy_path).read_bytes()).hexdigest() if taxonomy_path else None,
              'independence_confirmed': independent, 'analysis_mode': 'pooled' if use_all_samples else 'per_group',
              'original_group_counts': original_group_counts, 'random_permutations_requested': permutations if independent else 0,
              'seed': seed, 'min_prevalence': min_prevalence, 'min_abundance': min_abundance, 'top_n': top_n, 'p_max': p_max, 'q_max': q_max,
              'pseudocount_cpm': pseudocount, 'rho_min': rho_min, 'min_stability': min_stability, 'repeats_requested': repeats,
              'numpy': np.__version__, 'pandas': pd.__version__, 'groups': reports,
              'warning': EXPLORATION_NOTICE + (' ' + POOLED_NOTICE if use_all_samples else '') + ' 重选频率以全部子集为分母；p/q依赖置换假设。'}
    (out / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf8')
    return report
