"""Prepare the checksum-pinned original GSE146855 10x archive."""
import argparse
import hashlib
from pathlib import Path
import tarfile
import urllib.request

URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146855/suppl/GSE146855_RAW.tar'
SHA256 = 'a1c86880c26f65ad59da8b9a276541e0708d0cb8fad20b39cce0eec3364fc323'
MEMBERS = {'GSM4407907_barcodes.tsv.gz':'barcodes.tsv.gz',
           'GSM4407907_genes.tsv.gz':'features.tsv.gz',
           'GSM4407907_matrix.mtx.gz':'matrix.mtx.gz'}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument('--archive', type=Path, help='Explicit deposited original tar')
    source.add_argument('--download', action='store_true', help='Fetch the pinned original archive')
    ap.add_argument('--destination', type=Path, default=Path(__file__).resolve().parent/'output/scref/GSE146855')
    args = ap.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    archive = args.archive
    if args.download:
        archive = args.destination/'GSE146855_RAW.tar'
        if not archive.exists():
            temporary = archive.with_suffix('.tar.part')
            try:
                request = urllib.request.Request(URL,headers={'User-Agent':'Smad4-reproduction/1.0'})
                with urllib.request.urlopen(request,timeout=60) as r,temporary.open('wb') as w:
                    while block := r.read(1024*1024):w.write(block)
                if hashlib.file_digest(temporary.open('rb'),'sha256').hexdigest()!=SHA256:
                    raise ValueError('Downloaded reference checksum mismatch')
                temporary.replace(archive)
            finally:
                temporary.unlink(missing_ok=True)
    if hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()!=SHA256:
        raise ValueError(f'Reference archive checksum mismatch: {archive}')
    mtx = args.destination/'mtx';mtx.mkdir(exist_ok=True)
    with tarfile.open(archive) as tar:
        for name, output in MEMBERS.items():
            member = tar.getmember(name)
            if not member.isfile():raise ValueError(f'Invalid reference archive member: {name}')
            with tar.extractfile(member) as r,(mtx/output).open('wb') as w:
                while block:=r.read(1024*1024):w.write(block)
    print(f'Checksum-verified reference ready at {mtx}',flush=True)


if __name__=='__main__':main()
