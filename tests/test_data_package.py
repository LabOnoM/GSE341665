"""Deposit checks that numerical reproduction alone cannot establish."""
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from verify_data_package import DIRECTORIES, verify


def fixture():
    files = {f"{directory}/input.txt": b"scientific data" for directory in DIRECTORIES}
    rows = [{"path": name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(),
             "role": "source input"} for name, data in files.items()]
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=["path", "role", "bytes", "sha256"])
    writer.writeheader(); writer.writerows(rows)
    files.update({"DATA_FILE_MANIFEST.csv": stream.getvalue().encode(),
        "INPUT_MANIFEST.json": json.dumps(rows).encode(),
        "RAW_EVIDENCE_MANIFEST.json": b"[]", "README.md": b"Study data",
        "DATA_DICTIONARY.md": b"Data definitions", "SAMPLE_METADATA.tsv": b"sample\nCON1\n",
        "CODE_VERSION.txt": ("a"*40 + "\n").encode()})
    files["SHA256SUMS.txt"] = "".join(f"{hashlib.sha256(data).hexdigest()}  {name}\n"
        for name, data in files.items()).encode()
    return files


class DataPackageChecks(unittest.TestCase):
    def check(self, files):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "data.zip"
            with zipfile.ZipFile(p, "w") as z:
                for name, data in files.items():
                    z.writestr(name, data)
            return verify(p)

    def test_complete_package(self):
        self.assertEqual(self.check(fixture())["status"], "PASS")

    def test_numbering_gap_fails(self):
        files = {name.replace("07_locus", "08_locus"): data for name, data in fixture().items()}
        with self.assertRaisesRegex(ValueError, "consecutive"):
            self.check(files)

    def test_corrupt_scientific_file_fails(self):
        files = fixture(); files["01_processed_data/input.txt"] = b"changed data"
        with self.assertRaisesRegex(ValueError, "Checksum mismatch"):
            self.check(files)

    def test_undeclared_file_fails(self):
        files = fixture(); files["01_processed_data/old_draft.txt"] = b"old manuscript"
        with self.assertRaisesRegex(ValueError, "exactly"):
            self.check(files)

    def test_hidden_metadata_fails(self):
        files = fixture(); files[".DS_Store"] = b"Finder metadata"
        with self.assertRaisesRegex(ValueError, "hidden"):
            self.check(files)

    def test_zip_traversal_fails(self):
        files = fixture(); files["../outside.txt"] = b"bad path"
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            self.check(files)

    def test_document_checksum_fails(self):
        files = fixture(); files["README.md"] = b"stale changed README"
        with self.assertRaisesRegex(ValueError, "Checksum mismatch"):
            self.check(files)

if __name__ == "__main__":
    unittest.main()
