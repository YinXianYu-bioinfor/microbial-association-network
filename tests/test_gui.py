import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gui import build_options

class GUIOptionsTests(unittest.TestCase):
    def test_numeric_values_reach_shared_core(self):
        options = build_options('a.tsv', 'out', '', 'sample-id', 'group', '0.5', '0.001', '0.7', '1000', '36', False)
        self.assertEqual(options['counts_path'], 'a.tsv')
        self.assertEqual(options['permutations'], 1000)
        self.assertIsNone(options['metadata_path'])
        self.assertEqual(options['top_n'], 50)
        selected = build_options('a.tsv', 'out', '', 'sample-id', 'group', '0.5', '0.001', '0.7', '1000', '36', False, top_n='35')
        self.assertEqual(selected['top_n'], 35)
        self.assertIsNone(options['p_max'])
        self.assertIsNone(options['q_max'])
        filtered = build_options('a.tsv', 'out', '', 'sample-id', 'group', '0.5', '0.001', '0.7', '1000', '36', True, p_max='0.05', q_max='0.1')
        self.assertEqual(filtered['p_max'], .05)
        self.assertEqual(filtered['q_max'], .1)
        self.assertFalse(options['use_all_samples'])
        pooled = build_options('a.tsv', 'out', '', 'sample-id', 'group', '0.5', '0.001', '0.7', '1000', '36', False, use_all_samples=True)
        self.assertTrue(pooled['use_all_samples'])
