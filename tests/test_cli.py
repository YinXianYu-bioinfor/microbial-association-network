import subprocess
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

class CLITests(unittest.TestCase):
    def test_version_option_does_not_require_inputs(self):
        result = subprocess.run([sys.executable, str(ROOT/'cli.py'), '--version'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('0.1.0', result.stdout)

    def test_global_accepts_y_and_yes_default_keeps_groups(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); ids=[f's{i}' for i in range(8)]
            (p/'data.tsv').write_text('id\t'+'\t'.join(ids)+'\nA\t1\t2\t3\t4\t5\t6\t7\t8\nB\t8\t7\t6\t5\t4\t3\t2\t1\nC\t1\t3\t2\t4\t6\t5\t8\t7\n',encoding='utf8')
            (p/'meta.tsv').write_text('sample-id\tgroup\n'+'\n'.join(f'{s}\t'+('g1' if i<4 else 'g2') for i,s in enumerate(ids))+'\n',encoding='utf8')
            for opt in ('y','yes','default'):
                cmd=[sys.executable,str(ROOT/'cli.py'),'--counts',str(p/'data.tsv'),'--metadata',str(p/'meta.tsv'),'--output',str(p/opt),'--min-abundance','0','--repeats','0']
                if opt!='default':cmd.extend(['--use-all-samples',opt])
                result=subprocess.run(cmd,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                report=json.loads((p/opt/'report.json').read_text(encoding='utf8'))
                self.assertEqual(report['analysis_mode'],'per_group' if opt=='default' else 'pooled')
                self.assertEqual(len(report['groups']),2 if opt=='default' else 1)

    def test_cli_runs_small_input(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'data.tsv').write_text('id\ts1\ts2\ts3\ts4\nA\t1\t2\t3\t4\nB\t4\t3\t2\t1\nC\t1\t3\t2\t4\n', encoding='utf8')
            run = subprocess.run([sys.executable, str(ROOT/'cli.py'), '--counts', str(p/'data.tsv'), '--output', str(p/'out'), '--min-abundance', '0', '--repeats', '0'], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertTrue((p/'out'/'report.json').exists())
