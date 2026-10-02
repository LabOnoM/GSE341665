"""Regression checks for input safety and preservation of scientific evidence."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pandas as pd

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO))
from validate import compare_frames


class ReproductionChecks(unittest.TestCase):
    def run_cli(self,*args):
        return subprocess.run([sys.executable,*map(str,args)],capture_output=True,text=True)

    def test_tiny_probability_cannot_become_zero(self):
        a=pd.DataFrame({'padj':[1e-150,None,0.0]})
        with self.assertRaisesRegex(ValueError,'zero mask'):
            compare_frames(a,pd.DataFrame({'padj':[0.0,None,0.0]}),'probabilities')

    def test_missing_probability_cannot_become_zero(self):
        with self.assertRaisesRegex(ValueError,'missing-value mask'):
            compare_frames(pd.DataFrame({'padj':[None]}),pd.DataFrame({'padj':[0.0]}),'probabilities')

    def test_corrupt_input_fails_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);data=root/'data';work=root/'work';dest=root/'results'
            first=json.loads((REPO/'data/INPUT_MANIFEST.json').read_text())[0]
            p=data/first['path'];p.parent.mkdir(parents=True);p.write_bytes(b'corrupted input')
            result=self.run_cli(REPO/'reproduce.py','--data-root',data,'--work-dir',work,'--destination',dest)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('checksum mismatch',result.stderr)
            self.assertFalse(work.exists());self.assertFalse(dest.exists())

    def test_overlapping_output_paths_fail_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);work=root/'work';dest=work/'results'
            result=self.run_cli(REPO/'reproduce.py','--data-root',root/'data','--work-dir',work,'--destination',dest)
            self.assertNotEqual(result.returncode,0);self.assertIn('must be separate',result.stderr)
            self.assertFalse(work.exists())

    def test_nonempty_workspace_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);work=root/'work';work.mkdir();p=work/'keep.txt';p.write_text('retain')
            result=self.run_cli(REPO/'reproduce.py','--data-root',root/'data','--work-dir',work,'--destination',root/'results')
            self.assertNotEqual(result.returncode,0);self.assertIn('Refusing to overwrite',result.stderr)
            self.assertEqual(p.read_text(),'retain')

    def test_partial_reference_archive_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'GSE146855_RAW.tar').write_bytes(b'incomplete download')
            result=self.run_cli(REPO/'pipeline/reference_download.py','--download','--destination',root)
            self.assertNotEqual(result.returncode,0);self.assertIn('checksum mismatch',result.stderr)
            self.assertFalse((root/'mtx').exists())

    def test_panel_requires_original_source_records(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            result=self.run_cli(REPO/'pipeline/17_fetch_ai_panel.py','--sources-dir',root/'missing','--output-dir',root/'output')
            self.assertNotEqual(result.returncode,0);self.assertFalse((root/'output').exists())


if __name__=='__main__':unittest.main()
