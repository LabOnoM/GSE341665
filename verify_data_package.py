#!/usr/bin/env python3
"""Verify the complete deposited data package without extracting or modifying it."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import stat
import zipfile

DIRECTORIES = (
    "01_processed_data", "02_denovo_assembly", "03_external_panel_sources",
    "04_pipeline_objects", "05_experimental_data", "06_single_cell_reference",
    "07_locus_alignments",
)
DOCUMENTS = {
    "README.md", "DATA_DICTIONARY.md", "SAMPLE_METADATA.tsv", "CODE_VERSION.txt",
    "DATA_FILE_MANIFEST.csv", "INPUT_MANIFEST.json", "RAW_EVIDENCE_MANIFEST.json",
    "SHA256SUMS.txt",
}


def safe_path(name):
    path = PurePosixPath(name)
    if (not name or "\\" in name or path.is_absolute() or
            any(part in (".", "..") or part.startswith(".") or part == "__MACOSX"
                for part in path.parts) or str(path) != name):
        raise ValueError(f"Unsafe or hidden package path: {name}")
    return name


class Package:
    def __init__(self, path):
        self.path = Path(path)
        self.zip = None
        if self.path.is_dir():
            entries = list(self.path.rglob("*"))
            if any(p.is_symlink() for p in entries):
                raise ValueError("Symlinks are not permitted in the package")
            self.names = [safe_path(p.relative_to(self.path).as_posix())
                          for p in entries if p.is_file()]
        else:
            self.zip = zipfile.ZipFile(self.path)
            entries = self.zip.infolist()
            for entry in entries:
                safe_path(entry.filename.rstrip("/"))
                if stat.S_ISLNK(entry.external_attr >> 16):
                    raise ValueError("Symlinks are not permitted in the package")
            self.names = [e.filename for e in entries if not e.is_dir()]
            if len({e.filename for e in entries}) != len(entries):
                raise ValueError("Duplicate ZIP entries")
        if len({n.casefold() for n in self.names}) != len(self.names):
            raise ValueError("Case-insensitive path collision")

    def read(self, name):
        if self.zip:
            return self.zip.read(name)
        return (self.path / name).read_bytes()

    def close(self):
        if self.zip:
            self.zip.close()


def indexed(rows, label):
    result = {}
    for row in rows:
        name = safe_path(row["path"])
        if name in result:
            raise ValueError(f"Duplicate {label} entry: {name}")
        result[name] = row
    return result


def verify(package_path, repo=None):
    package = Package(package_path)
    try:
        names = set(package.names)
        dirs = {n.split("/", 1)[0] for n in names if "/" in n}
        if dirs != set(DIRECTORIES):
            raise ValueError(f"Expected consecutive directories 01–07; found {sorted(dirs)}")
        if {n for n in names if "/" not in n} != DOCUMENTS:
            raise ValueError("Missing or unexpected root documentation")
        manifest = indexed(csv.DictReader(io.StringIO(package.read("DATA_FILE_MANIFEST.csv").decode())), "data manifest")
        if set(manifest) != names - DOCUMENTS:
            raise ValueError("Data manifest does not cover exactly the scientific files")
        checksums = {}
        for line in package.read("SHA256SUMS.txt").decode().splitlines():
            digest, name = line.split("  ", 1)
            safe_path(name)
            if not re.fullmatch(r"[0-9a-f]{64}", digest) or name in checksums:
                raise ValueError("Invalid or duplicate checksum entry")
            checksums[name] = digest
        if set(checksums) != names - {"SHA256SUMS.txt"}:
            raise ValueError("Checksum list does not cover every other package file")
        for name, digest in checksums.items():
            data = package.read(name)
            if hashlib.sha256(data).hexdigest() != digest:
                raise ValueError(f"Checksum mismatch: {name}")
            if name in manifest:
                row = manifest[name]
                if int(row["bytes"]) != len(data) or row["sha256"] != digest:
                    raise ValueError(f"Data manifest mismatch: {name}")
        runtime = indexed(json.loads(package.read("INPUT_MANIFEST.json")), "runtime")
        evidence = indexed(json.loads(package.read("RAW_EVIDENCE_MANIFEST.json")), "evidence")
        if set(runtime) | set(evidence) != set(manifest):
            raise ValueError("Runtime and evidence manifests do not cover scientific files")
        for name, row in list(runtime.items()) + list(evidence.items()):
            if any(str(row[k]) != str(manifest[name][k]) for k in ("bytes", "sha256")):
                raise ValueError(f"Inconsistent manifest: {name}")
        if repo is not None:
            repo = Path(repo)
            for name in ("INPUT_MANIFEST.json", "RAW_EVIDENCE_MANIFEST.json"):
                if json.loads(package.read(name)) != json.loads((repo / "data" / name).read_text()):
                    raise ValueError(f"Package/code manifest mismatch: {name}")
        version = package.read("CODE_VERSION.txt").decode().splitlines()
        if not any(re.fullmatch(r"[0-9a-f]{40}", line) for line in version):
            raise ValueError("CODE_VERSION.txt must identify a complete Git commit")
        return {"status": "PASS", "files": len(names), "scientific_files": len(manifest),
                "runtime_inputs": len(runtime), "original_evidence_files": len(evidence),
                "directories": list(DIRECTORIES)}
    finally:
        package.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path, help="ZIP file or extracted data directory")
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    print(json.dumps(verify(args.package, args.repo), indent=2))
