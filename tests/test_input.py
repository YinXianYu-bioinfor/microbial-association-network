import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from network_core import load_inputs


class InputTests(unittest.TestCase):
    def test_qiime_biom_comment_and_single_group(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'table.tsv'
            p.write_text('# Constructed from biom file\n#OTU ID\ts1\ts2\nA\t1.0\t2.0\nB\t4.0\t3.0\n', encoding='utf8')
            counts, meta = load_inputs(p)
            self.assertEqual(counts.shape, (2, 2))
            self.assertEqual(meta['Group'].tolist(), ['all', 'all'])

    def test_missing_values_fail(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'table.tsv'
            p.write_text('OTUID\ts1\ts2\nA\t1\t\nB\t4\t3\n', encoding='utf8')
            with self.assertRaises(ValueError):
                load_inputs(p)
