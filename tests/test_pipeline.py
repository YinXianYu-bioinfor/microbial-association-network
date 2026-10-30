import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from network_core import analyze, select_top_features
import pandas as pd


class PipelineTests(unittest.TestCase):
    def test_pooled_global_mode_and_random_stream_independence(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            samples = [f's{i}' for i in range(10)]
            (p/'counts.tsv').write_text('id\t'+'\t'.join(samples)+'\n'+'A\t'+'\t'.join(map(str,range(1,11)))+'\n'+'B\t'+'\t'.join(map(str,range(10,0,-1)))+'\n'+'C\t'+'\t'.join(map(str,[1,3,2,5,4,7,6,9,8,10]))+'\n', encoding='utf8')
            (p/'metadata.tsv').write_text('sample-id\tgroup\n'+'\n'.join(f'{s}\t'+('g1' if i<5 else 'g2') for i,s in enumerate(samples))+'\n', encoding='utf8')
            common = dict(counts_path=p/'counts.tsv', metadata_path=p/'metadata.tsv', min_abundance=0, use_all_samples=True, permutations=20, repeats=5)
            r = analyze(output_dir=p/'confirmed', independent=True, **common)
            analyze(output_dir=p/'unconfirmed', independent=False, **common)
            self.assertEqual(r['analysis_mode'], 'pooled')
            self.assertEqual(r['original_group_counts'], {'g1': 5, 'g2': 5})
            self.assertEqual(r['groups'][0]['samples'], 10)
            self.assertEqual(r['groups'][0]['group'], 'all_samples')
            a = pd.read_csv(p/'confirmed'/'all_samples'/'edges.csv')
            b = pd.read_csv(p/'unconfirmed'/'all_samples'/'edges.csv')
            for column in ('stability_hits', 'reselected_hits', 'reselected_eligible', 'candidate'):
                self.assertTrue(a[column].equals(b[column]), column)
            grouped = analyze(p/'counts.tsv', p/'grouped', metadata_path=p/'metadata.tsv', min_abundance=0, repeats=0)
            self.assertEqual(grouped['analysis_mode'], 'per_group')
            self.assertEqual(len(grouped['groups']), 2)

    def test_optional_p_and_q_thresholds_filter_candidates_only(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'counts.tsv').write_text('id\ts1\ts2\ts3\ts4\ts5\nA\t1\t2\t3\t4\t5\nB\t5\t4\t3\t2\t1\nC\t1\t3\t2\t4\t6\n', encoding='utf8')
            report = analyze(p/'counts.tsv', p/'out', min_abundance=0, rho_min=0, min_stability=0, repeats=1, independent=True, p_max=.5, q_max=.3)
            edges = pd.read_csv(p/'out'/'all'/'edges.csv')
            selected = pd.read_csv(p/'out'/'all'/'candidate_edges.csv')
            expected = edges[(edges.p_value < .5) & (edges.q_value < .3)]
            self.assertEqual(len(edges), 3)
            self.assertEqual(selected[['source', 'target']].to_dict('records'), expected[['source', 'target']].to_dict('records'))
            self.assertEqual(report['p_max'], .5)
            self.assertEqual(report['q_max'], .3)
            with self.assertRaisesRegex(ValueError, '置换|独立'):
                analyze(p/'counts.tsv', p/'no_confirm', p_max=.05)
            self.assertFalse((p/'no_confirm').exists())
            with self.assertRaises(ValueError):
                analyze(p/'counts.tsv', p/'invalid', independent=True, q_max=float('nan'))

    def test_top_n_applied_after_prevalence_filter_per_subset(self):
        data = pd.DataFrame({'rare': [100, 0, 0, 0], 'a': [5, 5, 5, 5], 'b': [4, 4, 4, 4], 'c': [3, 3, 3, 3]})
        _, features = select_top_features(data, .75, 0, 3)
        self.assertEqual(features, ['a', 'b', 'c'])
        _, features = select_top_features(data.iloc[[0, 1, 2]], .25, 0, 3)
        self.assertIn('rare', features)

    def test_top_nodes_tie_break_and_separate_taxonomy_export(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'counts.tsv').write_text('id\ts1\ts2\ts3\ts4\nB\t20\t30\t40\t50\nA\t50\t40\t30\t20\nC\t20\t20\t20\t20\nD\t5\t5\t5\t5\nE\t5\t5\t5\t5\n', encoding='utf8')
            (p/'taxonomy.tsv').write_text('Feature ID\tTaxon\tConfidence\nA\td__Bacteria; p__A\t0.98\nC\td__Bacteria; p__C\t0.90\n', encoding='utf8')
            report = analyze(p/'counts.tsv', p/'out', taxonomy_path=p/'taxonomy.tsv', repeats=0, min_abundance=0, top_n=3)
            self.assertEqual(report['software_version'], '0.1.0')
            with (p/'out'/'all'/'selected_taxonomy.csv').open(encoding='utf8', newline='') as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([r['feature_id'] for r in rows], ['A', 'B', 'C'])
            self.assertEqual(rows[0]['taxonomy'], 'd__Bacteria; p__A')
            self.assertEqual(rows[0]['confidence'], '0.98')
            self.assertEqual(rows[1]['taxonomy'], '')
            self.assertEqual(report['top_n'], 3)
            with self.assertRaises(ValueError):
                analyze(p/'counts.tsv', p/'bad', top_n=2)

    def test_existing_directory_allows_unrelated_files_but_refuses_collisions(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            counts = p/'counts.tsv'
            counts.write_text('id\ts1\ts2\ts3\ts4\nA\t1\t2\t3\t4\nB\t4\t3\t2\t1\nC\t1\t3\t2\t4\n', encoding='utf8')
            out = p/'out'; out.mkdir()
            (out/'keep.txt').write_text('untouched', encoding='utf8')
            analyze(counts, out, existing_output=True, repeats=0, min_abundance=0)
            self.assertEqual((out/'keep.txt').read_text(), 'untouched')
            report = (out/'report.json').read_bytes()
            with self.assertRaises(FileExistsError):
                analyze(counts, out, existing_output=True, repeats=0, min_abundance=0)
            self.assertEqual((out/'report.json').read_bytes(), report)
            with self.assertRaises((FileNotFoundError, ValueError)):
                analyze(counts, p/'missing', existing_output=True)

    def test_outputs_and_existing_directory_protection(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'counts.tsv').write_text('id\ts1\ts2\ts3\ts4\ts5\nA\t1\t2\t3\t4\t5\nB\t5\t4\t3\t2\t1\nC\t1\t3\t2\t4\t6\n', encoding='utf8')
            out = p/'result'
            analyze(p/'counts.tsv', out, permutations=50, repeats=5, min_abundance=0, independent=True)
            with (out/'all'/'edges.csv').open(encoding='utf8', newline='') as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 3)
            self.assertTrue(all(r['p_value'] and r['q_value'] for r in rows))
            self.assertTrue(all('reselected_hits' in r and 'reselected_eligible' in r for r in rows))
            self.assertTrue(all(int(r['reselected_hits']) <= int(r['reselected_eligible']) <= int(r['stability_subsets']) for r in rows))
            self.assertTrue(all('loo_sign_fraction' in r and 'loo_max_delta' in r for r in rows))
            self.assertTrue(all((r['candidate'] == 'False') or float(r['reselected_frequency']) >= 0.8 for r in rows))
            self.assertTrue((out/'all'/'nodes.csv').exists())
            self.assertTrue((out/'report.json').exists())
            with self.assertRaises(FileExistsError):
                analyze(p/'counts.tsv', out)

    def test_taxonomy_annotation_on_nodes(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'counts.tsv').write_text('id\ts1\ts2\ts3\ts4\nA\t1\t2\t3\t4\nB\t4\t3\t2\t1\nC\t1\t3\t2\t4\n', encoding='utf8')
            (p/'taxonomy.tsv').write_text('Feature ID\tTaxon\tConfidence\nA\td__Bacteria; p__P\t0.8\n', encoding='utf8')
            analyze(p/'counts.tsv', p/'out', taxonomy_path=p/'taxonomy.tsv', repeats=0, min_abundance=0)
            with (p/'out'/'all'/'nodes.csv').open(encoding='utf8', newline='') as f:
                rows = {r['feature_id']: r for r in csv.DictReader(f)}
            self.assertEqual(rows['A']['taxonomy'], 'd__Bacteria; p__P')
            self.assertEqual(rows['B']['taxonomy'], '')

    def test_unconfirmed_independence_disables_p(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'counts.tsv').write_text('id\ts1\ts2\ts3\ts4\nA\t1\t2\t3\t4\nB\t4\t3\t2\t1\nC\t1\t3\t2\t4\n', encoding='utf8')
            analyze(p/'counts.tsv', p/'out', repeats=0, min_abundance=0)
            with (p/'out'/'all'/'edges.csv').open(encoding='utf8', newline='') as f:
                rows = list(csv.DictReader(f))
            self.assertTrue(all(not r['p_value'] and not r['q_value'] for r in rows))
            self.assertTrue(all(r['candidate'] == 'False' for r in rows))

    def test_small_group_uses_exact_permutations(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'counts.tsv').write_text('id\ts1\ts2\ts3\ts4\nA\t1\t2\t3\t4\nB\t4\t3\t2\t1\nC\t1\t3\t2\t4\n', encoding='utf8')
            report = analyze(p/'counts.tsv', p/'out', repeats=0, min_abundance=0, independent=True)
            self.assertEqual(report['groups'][0]['permutation_method'], 'exact')
            self.assertEqual(report['groups'][0]['permutations_used'], 24)
