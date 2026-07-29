#!/usr/bin/env python3
"""
20_smad4_loxp_evidence.py
=========================

Build the tidy data tables behind the Smad4 conditional-allele evidence figure.

WHAT THIS SCRIPT ESTABLISHES, FROM THE DATA, EVERY TIME IT RUNS
---------------------------------------------------------------
1. The Smad4 floxed allele is physically present in the sequencing libraries:
   RNA-seq reads soft-clip into two loxP cassettes flanking a single 184 bp
   exon (Ensembl exon ENSMUSE00000143011, chr18:73,782,797-73,782,980, minus
   strand; called "exon 8" in the Cre-lox literature and "exon 9" on the
   Ensembl canonical transcript ENSMUST00000025393 - the same 184 bp).
2. The cassette is present in BOTH genotypes.  CON is Cre-negative
   Smad4^fl/fl, so loxP reads in CON are EXPECTED and are not a control
   failure.  This must be stated in the figure legend.
3. The NEGATIVE RESULT: the Cre-mediated exon-skip junction
   chr18:73,782,125-73,790,265 has ZERO uniquely-mapped reads in all six
   samples.  PSI = 1.0 everywhere.  Smad4 mRNA in dcKO retains the floxed
   exon.  This is not hidden; it is a primary output of this script.

DISCREPANCY WITH A PRIOR PASS - READ BEFORE WRITING THE LEGEND
---------------------------------------------------------------
A prior analysis recorded that the floxed exon's coverage dip "FAILS testing
(dcKO3 > CON3 on the read-fraction metric, P = 0.20)" and that "exon 6 shows
an equal and OPPOSITE clean separation".  Neither statement reproduces here.
Under all three metrics computed or tested by this script -
  (a) fraction of per-base Smad4 exonic depth,
  (b) fraction of Smad4 exonic reads,
  (c) fraction of Smad4 exonic reads counting uniquely-mapped reads only
      (MAPQ >= 60; note these BAMs use MAPQ 60, not 255, for unique) -
the floxed exon is the ONLY exon of the twelve whose CON and dcKO values do
not overlap, giving P = 0.10 (the floor), and exon 6 never separates cleanly
(P = 0.20-0.30).  The prior metric definition could not be recovered.

This does NOT change the conclusion, and it must not be spun as if it did:
P = 0.10 is the smallest value this design can produce, so the dip is still
not significant and still must not be drawn as a demonstrated depletion.  But
the numbers "P = 0.20", "dcKO3 > CON3" and the exon-6 claim are not supported
and must not be quoted.  The script asserts what actually reproduces.

WHAT IS DERIVED RATHER THAN ASSERTED (added 2026-07-28 after audit)
-------------------------------------------------------------------
* BOTH cassette breakpoints are re-derived from the reads (derive_insertion_
  sites) and the module constants SITE_A_POS / SITE_B_POS are ASSERTED against
  them.  Previously they were the only load-bearing numbers here that nothing
  checked: setting either to a wrong value produced a confident, silently wrong
  summary JSON with exit 0.
* The ±3 bp ambiguity at site A is the spread of the per-read placements, and
  the microhomology that explains it is verified in genome.fa to be a tandem
  repeat.  Neither is quoted from memory.
* The restriction-linker composition of each cassette is derived from the base
  calls (derive_cassette_linkers), bounded by CIGAR arithmetic alone, including
  sites formed ACROSS the loxP-core/linker boundary.  The previous hard-coded
  strings named XbaI/SalI at site A but omitted EcoRV at site B, so the figure
  contradicted its own per-read decomposition.
* The floxed exon's stable-ID / exon-number ambiguity is read from the GTF
  (floxed_exon_identity).  The stable ID is NOT unambiguous: three different
  stable IDs are attached to the identical 184 bp within this one gene.
* loxP reads that are byte-identical copies of another read in the same library
  are flagged (flag_duplicate_reads), so a figure that counts evidence reads
  cannot present one fragment as two independent observations.

PROJECT RULES HONOURED HERE
---------------------------
* No data value is hard-coded into an output.  Every number is recomputed
  from the BAMs / SJ.out.tab / DEG_table.csv at run time.  The one block of
  non-measured content is STUDY_DESIGN (which Cre driver produced the sequenced
  material); it is provenance metadata, is labelled as such, and carries its
  source.
* If a source file is missing the script RAISES.  There is no fallback.
* DEG_table.csv is read with float_precision="round_trip" and header=1
  (two-row header; tiny padj values are stored in full decimal and the
  default parser flushes them to 0.0).
* Exon numbering is on the MINUS strand: exon 1 = highest coordinate.
* Key results are asserted against independently established facts.  If the
  underlying data ever change, this script breaks rather than lies.

OUTPUTS (into scripts/JDR_reanalysis/output/)
---------------------------------------------
  smad4_loxp_reads.csv      one row per loxP-containing read, all six BAMs
  smad4_exon_coverage.csv   per exon x sample depth + within-gene fraction
  smad4_exon_ratio.csv      per exon CON vs dcKO fraction, ratio, exact perm P
  smad4_junctions.csv       every junction in the locus, annotated
  smad4_locus_depth.csv     per-base depth, 10 bp bins, sample-normalised
  smad4_loxp_summary.json   scalar facts for the figure legend

Usage:  python3 scripts/JDR_reanalysis/20_smad4_loxp_evidence.py
"""

from __future__ import annotations

import itertools
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path("/mnt/Disk1/ZXY_2026")
BAM_DIR = PROJECT_ROOT / "00.RawData" / "Bulk-RNAseq"
SCRIPT_DIR = PROJECT_ROOT / "scripts" / "JDR_reanalysis"
OUT_DIR = SCRIPT_DIR / "output"
DATA_DIR = SCRIPT_DIR / "data"

GENOME_FA = BAM_DIR / "genome.fa"
GENOME_GTF = BAM_DIR / "genome.gtf"           # 3.5 GB - only scanned to build the cache
SMAD4_GTF_CACHE = DATA_DIR / "smad4_locus.gtf"  # small cached subset, committed
DEG_TABLE = OUT_DIR / "DEG_table.csv"

SAMPLES = ["CON1", "CON2", "CON3", "dcKO1", "dcKO2", "dcKO3"]
GENOTYPE = {s: ("CON" if s.startswith("CON") else "dcKO") for s in SAMPLES}


def bam_path(sample: str) -> Path:
    return BAM_DIR / f"{sample}.STAR.Aligned.sortedByCoord.out.bam"


def sj_path(sample: str) -> Path:
    return BAM_DIR / f"{sample}.STAR.SJ.out.tab"


# ---------------------------------------------------------------------------
# Locus constants.  These are COORDINATES AND SEQUENCE MOTIFS, not measurements.
# Every measured quantity below is computed from the data; these are the
# addresses at which to look, plus the 34 bp loxP core (a fixed reagent
# sequence, not a result).
# ---------------------------------------------------------------------------
CHROM = "18"
SMAD4_GENE_ID = "ENSMUSG00000024515"
SMAD4_CANONICAL_TX = "ENSMUST00000025393"
GENE_START, GENE_END = 73_769_343, 73_836_915          # 1-based inclusive
GENE_STRAND = "-"

FLOXED_EXON_ID = "ENSMUSE00000143011"
FLOXED_EXON_START, FLOXED_EXON_END = 73_782_797, 73_782_980   # 184 bp

# loxP insertion points determined from soft-clipped reads (see docstring).
SITE_A_POS = 73_782_242      # 3' of the floxed exon in genome order; +/- 3 bp (AAC microhomology)
SITE_B_POS = 73_784_709      # 5' of the floxed exon in genome order; exact
FLOXED_INTERVAL_START = SITE_A_POS + 1        # 73,782,243
FLOXED_INTERVAL_END = SITE_B_POS              # 73,784,709

LOXP_FWD = "ATAACTTCGTATAATGTATGCTATACGAAGTTAT"
LOXP_REV = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"   # reverse complement of LOXP_FWD

# Search windows
LOXP_SEARCH_START, LOXP_SEARCH_END = 73_778_000, 73_792_000
JUNCTION_WINDOW_START, JUNCTION_WINDOW_END = 73_779_000, 73_793_000
LOCUS_DEPTH_START, LOCUS_DEPTH_END = 73_780_000, 73_792_000
LOCUS_BIN_BP = 10

# Predicted Cre product: skip junction joining the two exons that flank the
# floxed exon (exon 10 end -> exon 8 start in minus-strand numbering).
SKIP_JUNCTION = (73_782_125, 73_790_265)

# Restriction recognition sequences used to ANNOTATE the polylinker that the
# reads actually carry.  This is a lookup table of reagent motifs, not a claim
# about this allele: which of them are present at each cassette is DERIVED from
# the base calls (see derive_cassette_linkers).  Palindromic sites need no
# reverse complement; the two non-palindromic ones are listed explicitly.
RESTRICTION_SITES = {
    "SacI": "GAGCTC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "XbaI": "TCTAGA",
    "SalI": "GTCGAC",
    "EcoRV": "GATATC",
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "NotI": "GCGGCCGC",
    "KpnI": "GGTACC",
    "PstI": "CTGCAG",
    "SpeI": "ACTAGT",
    "NheI": "GCTAGC",
    "ClaI": "ATCGAT",
    "BglII": "AGATCT",
}

# Study design / provenance.  Not measurements: recorded here so that every
# downstream deliverable states which animals produced the sequenced material.
# Source: CLAUDE.md section 5 (corrected 2026-07-27 against the original BBRC
# submission, whose Results and DEG-figure legend both state K14-creER for the
# sequenced material).
STUDY_DESIGN = dict(
    sequenced_material_driver="K14-creER (tamoxifen-inducible)",
    control_genotype="Cre-negative Smad4^fl/fl^ littermate",
    dcko_genotype="K14-creER Smad4^fl/fl^ (tamoxifen-induced)",
    tissue="bulk dissected tooth germ (whole tissue, not sorted epithelium)",
    caveat=(
        "The RNA-seq libraries come from K14-creER animals, so excision is "
        "tamoxifen-dependent and restricted to the K14-expressing epithelium. "
        "The libraries are bulk tissue, so any epithelium-restricted excision is "
        "diluted by mesenchyme and by epithelial cells that never expressed the "
        "driver. Figure 1 (in vivo phenotype) uses constitutive K14-cre and is a "
        "DIFFERENT experiment. A PSI measured on bulk RNA therefore bounds the "
        "population-average excision, not the excision inside a recombined cell."
    ),
    source="CLAUDE.md section 5 (genotypes differ by experiment)",
)

# A fragment flanking the loxP core is called "genome-derived" if it contains an
# exact match of >= MIN_GENOMIC_ANCHOR bp to the reference WITHIN
# ANCHOR_SEARCH_PAD bp of the cassette insertion point it belongs to.
# Choosing 14 bp: the observed genomic flank on the decisive CON1 site-B read is
# 15 bp of exact match (the read's first base is a mismatch), so 16 bp would
# discard a true positive.  Restricting the search to a +/-300 bp window around
# the breakpoint keeps the chance-match rate negligible: ~600 x 2 strands / 4^14
# = 9e-6 per candidate k-mer, ~1e-3 per 150 bp fragment.  The realised anchor
# length is reported per read (genomic_anchor_5p_bp / _3p_bp) so the threshold is
# auditable rather than load-bearing.
MIN_GENOMIC_ANCHOR = 14
ANCHOR_SEARCH_PAD = 300


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------
def require(path: Path, what: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(
            f"Required {what} is missing: {path}\n"
            "This script never substitutes a fallback value. Restore the file and re-run."
        )
    return path


def run(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\nSTDERR:\n{proc.stderr}")
    return proc.stdout


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def ref_span(cigar: str) -> int:
    return sum(n for n, op in parse_cigar(cigar) if op in REF_CONSUMING)


def lead_trail_clip(cigar: str) -> tuple[int, int]:
    ops = parse_cigar(cigar)
    lead = ops[0][0] if ops and ops[0][1] in "SH" else 0
    trail = ops[-1][0] if ops and ops[-1][1] in "SH" else 0
    return lead, trail


# ---------------------------------------------------------------------------
# Reference sequence + exon model
# ---------------------------------------------------------------------------
def fetch_reference_window(start: int, end: int) -> str:
    require(GENOME_FA, "genome FASTA")
    out = run(["samtools", "faidx", str(GENOME_FA), f"{CHROM}:{start}-{end}"])
    seq = "".join(line.strip() for line in out.splitlines() if not line.startswith(">"))
    if len(seq) != end - start + 1:
        raise RuntimeError(f"faidx returned {len(seq)} bp for {CHROM}:{start}-{end}")
    return seq.upper()


def ensure_smad4_gtf_cache() -> Path:
    """Small, committed subset of the Ensembl GTF covering only Smad4.

    Built once from genome.gtf (3.5 GB).  Never fabricated: if neither the
    cache nor genome.gtf exists, this raises.
    """
    if SMAD4_GTF_CACHE.exists() and SMAD4_GTF_CACHE.stat().st_size > 0:
        return SMAD4_GTF_CACHE
    require(GENOME_GTF, "Ensembl GTF (needed once to build the Smad4 subset cache)")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    sys.stderr.write(
        f"[20] {SMAD4_GTF_CACHE} not found; extracting {SMAD4_GENE_ID} from "
        f"{GENOME_GTF} (one-time full scan)...\n"
    )
    tmp = SMAD4_GTF_CACHE.with_suffix(".gtf.tmp")
    with open(tmp, "w") as fout:
        proc = subprocess.run(["grep", "-F", SMAD4_GENE_ID, str(GENOME_GTF)],
                              stdout=fout, stderr=subprocess.PIPE, text=True)
    if proc.returncode not in (0, 1):
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"grep on {GENOME_GTF} failed:\n{proc.stderr}")
    tmp.replace(SMAD4_GTF_CACHE)
    if SMAD4_GTF_CACHE.stat().st_size == 0:
        raise RuntimeError(f"{SMAD4_GENE_ID} not found in {GENOME_GTF}")
    return SMAD4_GTF_CACHE


ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def parse_gtf_attrs(field: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(field))


def load_exon_model() -> tuple[pd.DataFrame, list[tuple[int, int]]]:
    """Canonical 12-exon model (minus strand: exon 1 = highest coordinate),
    plus the merged exonic blocks across ALL Smad4 transcripts."""
    gtf = ensure_smad4_gtf_cache()
    canon: list[dict] = []
    all_exon_iv: list[tuple[int, int]] = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            start, end, strand = int(f[3]), int(f[4]), f[6]
            attrs = parse_gtf_attrs(f[8])
            if attrs.get("gene_id") != SMAD4_GENE_ID:
                continue
            all_exon_iv.append((start, end))
            if attrs.get("transcript_id") == SMAD4_CANONICAL_TX:
                canon.append(
                    dict(
                        exon_number=int(attrs["exon_number"]),
                        exon_id=attrs.get("exon_id", ""),
                        start=start,
                        end=end,
                        strand=strand,
                        length=end - start + 1,
                    )
                )
    if not canon:
        raise RuntimeError(f"no exons for canonical transcript {SMAD4_CANONICAL_TX} in {gtf}")
    ex = pd.DataFrame(canon).sort_values("exon_number").reset_index(drop=True)

    # --- assertions against the established gene model -------------------
    assert len(ex) == 12, f"canonical transcript should have 12 exons, got {len(ex)}"
    assert (ex["strand"] == GENE_STRAND).all(), "Smad4 must be on the minus strand"
    # minus strand: exon 1 = highest coordinate, strictly decreasing thereafter
    starts = ex.sort_values("exon_number")["start"].to_numpy()
    assert (np.diff(starts) < 0).all(), (
        "exon numbering is not monotonically decreasing in coordinate; "
        "minus-strand exon order is wrong"
    )
    hit = ex[(ex["start"] == FLOXED_EXON_START) & (ex["end"] == FLOXED_EXON_END)]
    assert len(hit) == 1, (
        f"expected exactly one canonical exon at {FLOXED_EXON_START}-{FLOXED_EXON_END}, "
        f"got {len(hit)}"
    )
    assert hit.iloc[0]["exon_id"] == FLOXED_EXON_ID, (
        f"exon at {FLOXED_EXON_START}-{FLOXED_EXON_END} is "
        f"{hit.iloc[0]['exon_id']}, expected {FLOXED_EXON_ID}"
    )
    assert int(hit.iloc[0]["length"]) == 184, "floxed exon must be 184 bp"

    # merged exonic blocks over ALL transcripts
    merged: list[list[int]] = []
    for s, e in sorted(set(all_exon_iv)):
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    blocks = [(s, e) for s, e in merged]

    # the floxed interval must contain exactly ONE exonic block, and that
    # block must be the floxed exon itself
    inside = [b for b in blocks if b[0] >= FLOXED_INTERVAL_START and b[1] <= FLOXED_INTERVAL_END]
    assert len(inside) == 1, (
        f"floxed interval {FLOXED_INTERVAL_START}-{FLOXED_INTERVAL_END} must contain exactly "
        f"one exonic block (any transcript); found {len(inside)}: {inside}"
    )
    assert inside[0] == (FLOXED_EXON_START, FLOXED_EXON_END), (
        f"the single exonic block in the floxed interval is {inside[0]}, expected "
        f"({FLOXED_EXON_START}, {FLOXED_EXON_END})"
    )
    # and no exon from any transcript may straddle either breakpoint
    for s, e in blocks:
        for bp, name in ((SITE_A_POS, "A"), (SITE_B_POS, "B")):
            assert not (s <= bp <= e), (
                f"loxP site {name} at {bp} falls inside exonic block {s}-{e}; "
                "the cassettes are supposed to be intronic"
            )
    return ex, blocks


def floxed_exon_identity() -> dict:
    """Every GTF record for the 184 bp floxed interval, DERIVED not asserted.

    The stable exon ID does NOT unambiguously identify this 184 bp: within the
    same gene, some transcripts label the identical interval with a different
    stable ID.  Querying one ID in Ensembl therefore silently loses transcripts.
    The unambiguous identifier is the coordinate interval.  This function
    reports what the annotation actually contains so that no downstream text has
    to state a transcript count, an exon-ID list or an exon_number list from
    memory.
    """
    gtf = ensure_smad4_gtf_cache()
    recs = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            if int(f[3]) != FLOXED_EXON_START or int(f[4]) != FLOXED_EXON_END:
                continue
            a = parse_gtf_attrs(f[8])
            if a.get("gene_id") != SMAD4_GENE_ID:
                continue
            recs.append(dict(transcript_id=a.get("transcript_id", ""),
                             exon_id=a.get("exon_id", ""),
                             exon_number=int(a["exon_number"])))
    if not recs:
        raise RuntimeError(
            f"no exon record at {CHROM}:{FLOXED_EXON_START}-{FLOXED_EXON_END} in {gtf}")
    df = pd.DataFrame(recs)
    by_id = df.groupby("exon_id")["transcript_id"].nunique().sort_values(ascending=False)
    primary = str(by_id.index[0])
    assert primary == FLOXED_EXON_ID, (
        f"the majority stable exon ID for {FLOXED_EXON_START}-{FLOXED_EXON_END} is now "
        f"{primary}, not {FLOXED_EXON_ID}")
    numbers = sorted(df["exon_number"].unique())
    return dict(
        interval=f"{CHROM}:{FLOXED_EXON_START}-{FLOXED_EXON_END}",
        n_transcript_records=int(len(df)),
        n_transcripts=int(df["transcript_id"].nunique()),
        primary_exon_id=primary,
        n_transcripts_with_primary_exon_id=int(by_id.iloc[0]),
        exon_ids={str(k): int(v) for k, v in by_id.items()},
        other_exon_ids=[str(k) for k in by_id.index[1:]],
        exon_number_values=[int(v) for v in numbers],
        exon_number_counts={str(k): int(v)
                            for k, v in df["exon_number"].value_counts().sort_index().items()},
        unambiguous_identifier=(
            f"chr{CHROM}:{FLOXED_EXON_START:,}-{FLOXED_EXON_END:,} — the COORDINATE "
            f"interval. {by_id.iloc[0]} of the {df['transcript_id'].nunique()} transcripts "
            f"carrying this exon label it {primary}; the remaining "
            f"{int(df['transcript_id'].nunique() - by_id.iloc[0])} label the identical "
            f"interval {', '.join(str(k) for k in by_id.index[1:])}. Querying one stable "
            "ID in Ensembl silently drops the others."
        ),
    )


# ---------------------------------------------------------------------------
# 1. loxP-containing reads
# ---------------------------------------------------------------------------
def collect_loxp_reads(ref_win_start: int, ref_win: str) -> pd.DataFrame:
    # Local reference windows around each insertion point.  Anchors must match
    # HERE, not anywhere in the locus, so a chance match cannot manufacture a
    # "cassette-spanning" read.
    local_ref: dict[str, tuple[str, str]] = {}
    for site, bp in (("A", SITE_A_POS), ("B", SITE_B_POS)):
        lo = bp - ANCHOR_SEARCH_PAD - ref_win_start
        hi = bp + ANCHOR_SEARCH_PAD - ref_win_start + 1
        if lo < 0 or hi > len(ref_win):
            raise RuntimeError(f"site {site} anchor window falls outside the fetched reference")
        w = ref_win[lo:hi]
        local_ref[site] = (w, revcomp(w))

    rows: list[dict] = []

    for sample in SAMPLES:
        bam = require(bam_path(sample), f"BAM for {sample}")
        require(Path(str(bam) + ".bai"), f"BAM index for {sample}")
        sam = run(
            ["samtools", "view", "-F", "0x900", str(bam),
             f"{CHROM}:{LOXP_SEARCH_START}-{LOXP_SEARCH_END}"]
        )
        for line in sam.splitlines():
            if not line:
                continue
            f = line.split("\t")
            qname, flag, pos, cigar, seq = f[0], int(f[1]), int(f[3]), f[5], f[9]
            for loxp_orient, motif in (("fwd", LOXP_FWD), ("rev", LOXP_REV)):
                off = seq.find(motif)
                if off < 0:
                    continue
                lead, trail = lead_trail_clip(cigar)
                aln_end = pos + ref_span(cigar) - 1
                d_a = min(abs(pos - SITE_A_POS), abs(aln_end - SITE_A_POS))
                d_b = min(abs(pos - SITE_B_POS), abs(aln_end - SITE_B_POS))
                site = "A" if d_a <= d_b else "B"

                lw, lw_rc = local_ref[site]
                prefix = seq[:off]
                suffix = seq[off + len(motif):]
                pre_gen = longest_genomic_anchor(prefix, lw, lw_rc)
                suf_gen = longest_genomic_anchor(suffix, lw, lw_rc)
                spans = bool(pre_gen >= MIN_GENOMIC_ANCHOR and suf_gen >= MIN_GENOMIC_ANCHOR)

                rows.append(
                    dict(
                        sample=sample,
                        genotype=GENOTYPE[sample],
                        qname=qname,
                        pos=pos,
                        cigar=cigar,
                        aln_end=aln_end,
                        orient="rev" if (flag & 0x10) else "fwd",
                        site=site,
                        lead_clip=lead,
                        trail_clip=trail,
                        loxp_offset_in_read=off,
                        spans_cassette=spans,
                        seq=seq,
                        loxp_orient=loxp_orient,
                        dist_to_site_A=d_a,
                        dist_to_site_B=d_b,
                        genomic_anchor_5p_bp=pre_gen,
                        genomic_anchor_3p_bp=suf_gen,
                        linkers=";".join(
                            sorted(n for n, s in RESTRICTION_SITES.items() if s in seq)
                        ),
                    )
                )
                break  # a read is counted once, under the orientation that matched

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("no loxP-containing reads found - refusing to emit an empty table")
    df = df.sort_values(["sample", "pos", "qname"]).reset_index(drop=True)
    return df


def core_side(row) -> str:
    """Which soft-clip carries the loxP core: 'lead', 'trail', or 'straddle'.

    Read coordinates: [0, lead_clip) is the leading clip, [lead_clip,
    lead_clip + aln_len) is the CIGAR M block, the rest is the trailing clip.
    """
    lead = int(row["lead_clip"])
    aln_len = int(row["aln_end"]) - int(row["pos"]) + 1
    off = int(row["loxp_offset_in_read"])
    if off + len(LOXP_FWD) <= lead:
        return "lead"
    if off >= lead + aln_len:
        return "trail"
    return "straddle"


def derive_insertion_sites(lox: pd.DataFrame, ref_win_start: int, ref_win: str) -> dict:
    """Re-derive both cassette breakpoints FROM THE READS.

    Nothing here consults SITE_A_POS / SITE_B_POS; the caller asserts the two
    agree.  Previously those two coordinates were the only load-bearing numbers
    in this script that were never checked against the data, so a wrong constant
    produced a confident, silently wrong summary JSON.

    For a read whose loxP core sits in the LEADING clip, the alignment re-enters
    the genome at `pos`, so the cassette was inserted between pos-1 and pos.
    For a read whose core sits in the TRAILING clip, the alignment leaves the
    genome after `aln_end`, so the insertion is between aln_end and aln_end+1.
    The breakpoint reported per site is the modal value across its reads; the
    spread of the cluster IS the microhomology ambiguity, and the microhomology
    is verified to be a genuine tandem repeat in genome.fa.
    """
    out: dict[str, dict] = {}
    for site, g in lox.groupby("site"):
        bps: list[int] = []
        for _, r in g.iterrows():
            sd = core_side(r)
            if sd == "lead":
                bps.append(int(r["pos"]) - 1)
            elif sd == "trail":
                bps.append(int(r["aln_end"]))
        if not bps:
            raise RuntimeError(f"site {site}: no read places the cassette breakpoint")
        counts = pd.Series(bps).value_counts()
        top = counts.max()
        bp = int(max(v for v, c in counts.items() if c == top))   # ties -> higher coordinate
        lo, hi = min(bps), max(bps)
        amb = hi - lo
        micro = ""
        if amb:
            micro = ref_win[lo + 1 - ref_win_start: hi + 1 - ref_win_start]
            nxt = ref_win[hi + 1 - ref_win_start: hi + 1 + amb - ref_win_start]
            if micro != nxt:
                raise RuntimeError(
                    f"site {site}: reads place the breakpoint over {lo}-{hi} but "
                    f"{micro!r} at {lo+1}-{hi} is not repeated at {hi+1}-{hi+amb} "
                    f"({nxt!r}); the spread is not explained by microhomology"
                )
        out[site] = dict(
            breakpoint=bp,
            n_reads_used=len(bps),
            breakpoint_min=lo,
            breakpoint_max=hi,
            ambiguity_bp=amb,
            microhomology=micro,
            per_read_breakpoints=sorted(bps),
        )
    return out


def restriction_hits(frag: str, boundary_prefix: str = "") -> list[dict]:
    """Recognition sites in `frag`, DERIVED from the base calls.

    `boundary_prefix` is sequence immediately 5' of `frag` (in practice the 3'
    end of the loxP core).  A site that starts inside the prefix straddles the
    core/linker boundary and is reported with straddles_core=True -- this is why
    panel A previously said "XbaI/SalI" while panel B's per-read decomposition of
    the same segment said only "SalI": XbaI's TCTAGA is formed by the final T of
    the loxP core plus CTAGA of the linker, so it is invisible to a scan of the
    linker alone.
    """
    ctx = boundary_prefix + frag
    hits = []
    for name, motif in RESTRICTION_SITES.items():
        for m in re.finditer(f"(?={re.escape(motif)})", ctx):
            i = m.start()
            if i + len(motif) <= len(boundary_prefix):
                continue                      # entirely inside the core
            hits.append(dict(enzyme=name, motif=motif,
                             offset_in_fragment=i - len(boundary_prefix),
                             straddles_core=i < len(boundary_prefix)))
    return sorted(hits, key=lambda h: (h["offset_in_fragment"], h["enzyme"]))


def _consensus(frags: list[str]) -> str:
    """Per-position majority consensus of equal-length fragments."""
    L = max(len(f) for f in frags)
    keep = [f for f in frags if len(f) == L]
    return "".join(pd.Series([f[i] for f in keep]).value_counts().idxmax() for i in range(L))


def derive_cassette_linkers(lox: pd.DataFrame) -> dict:
    """Derive each cassette's polylinker sequences from the reads themselves.

    Both linkers at both sites are bounded by CIGAR arithmetic alone -- no
    threshold, no anchor search:
      * a read whose core is in the LEADING clip gives the 3' linker exactly,
        seq[core_end : lead_clip] (the CIGAR M block starts at lead_clip);
      * a read whose core is in the TRAILING clip gives the 5' linker exactly,
        seq[lead_clip + aln_len : core_start].
    Both read geometries occur at both sites, so both linkers are recovered at
    both sites and the result is a measurement, not an annotation copied from a
    vector map.
    """
    out: dict[str, dict] = {}
    for site, g in lox.groupby("site"):
        five, three = [], []
        for _, r in g.iterrows():
            seq, off, lead = r["seq"], int(r["loxp_offset_in_read"]), int(r["lead_clip"])
            aln_len = int(r["aln_end"]) - int(r["pos"]) + 1
            sd = core_side(r)
            if sd == "trail":
                five.append(seq[lead + aln_len: off])
            elif sd == "lead":
                three.append(seq[off + len(LOXP_FWD): lead])
        if not five or not three:
            raise RuntimeError(
                f"site {site}: cannot derive both linkers -- {len(five)} reads place the "
                f"core in the trailing clip and {len(three)} in the leading clip; both "
                "geometries are needed and neither may be assumed"
            )
        v5, v3 = _consensus(five), _consensus(three)
        h5 = restriction_hits(v5)
        h3 = restriction_hits(v3, boundary_prefix=LOXP_FWD[-8:])

        def label(hits):
            if not hits:
                return "unnamed linker"
            return " + ".join(h["enzyme"] + ("\u2020" if h["straddles_core"] else "") for h in hits)

        out[site] = dict(
            linker_5p_sequence=v5,
            linker_3p_sequence=v3,
            linker_5p_n_reads=len(five),
            linker_3p_n_reads=len(three),
            linker_5p_sites=h5,
            linker_3p_sites=h3,
            linkers=f"{label(h5)} ({v5}) ... loxP ... {label(h3)} ({v3})",
            linkers_short=f"{label(h5)}–loxP–{label(h3)}",
            straddling_site_note=(
                "\u2020 = the recognition site is formed across the loxP-core/linker boundary "
                "and is therefore absent from a scan of the linker sequence alone."
                if any(h["straddles_core"] for h in h5 + h3) else ""
            ),
        )
    return out


def flag_duplicate_reads(lox: pd.DataFrame) -> pd.DataFrame:
    """Mark reads that are the same physical fragment sequenced twice.

    Criterion: same library, same POS, same CIGAR and byte-identical SEQ.  These
    are not flagged in the BAM (no tool marked duplicates), but a figure whose
    entire content is counting evidence reads must not present one fragment as
    two independent observations.
    """
    df = lox.copy()
    key = df["sample"].astype(str) + "|" + df["pos"].astype(str) + "|" + df["cigar"] + "|" + df["seq"]
    grp = key.map({k: i for i, k in enumerate(sorted(key.unique()))})
    size = key.map(key.value_counts())
    df["duplicate_group"] = grp
    df["duplicate_group_size"] = size
    # the first member of each group is the representative; the rest are copies
    df["is_duplicate_copy"] = df.duplicated(subset=["sample", "pos", "cigar", "seq"], keep="first")
    return df


def longest_genomic_anchor(frag: str, ref_win: str, ref_win_rc: str) -> int:
    """Longest exact substring of `frag` (>= MIN_GENOMIC_ANCHOR) found in the
    LOCAL reference window (breakpoint +/- ANCHOR_SEARCH_PAD) on either strand.
    Returns 0 if none.  Used to decide whether a fragment flanking the loxP core
    is genome-derived or cassette-derived."""
    if len(frag) < MIN_GENOMIC_ANCHOR:
        return 0
    for k in range(len(frag), MIN_GENOMIC_ANCHOR - 1, -1):
        for i in range(0, len(frag) - k + 1):
            sub = frag[i:i + k]
            if sub in ref_win or sub in ref_win_rc:
                return k
    return 0


# ---------------------------------------------------------------------------
# 2. Exon coverage + exact permutation test
# ---------------------------------------------------------------------------
def samtools_depth(sample: str, start: int, end: int) -> np.ndarray:
    """Per-base depth (all reads, no MAPQ/BAQ filter), 1-based inclusive."""
    bam = require(bam_path(sample), f"BAM for {sample}")
    out = run(
        ["samtools", "depth", "-a", "-d", "0", "-Q", "0", "-q", "0",
         "-r", f"{CHROM}:{start}-{end}", str(bam)]
    )
    n = end - start + 1
    depth = np.zeros(n, dtype=np.int64)
    for line in out.splitlines():
        if not line:
            continue
        _, p, d = line.split("\t")
        depth[int(p) - start] = int(d)
    return depth


def aligned_blocks(pos: int, cigar: str) -> list[tuple[int, int]]:
    """Reference intervals actually covered by aligned bases (M/=/X), 1-based
    inclusive. N (splice gaps) and D are not covered."""
    blocks, ref = [], pos
    for n, op in parse_cigar(cigar):
        if op in "M=X":
            blocks.append((ref, ref + n - 1))
            ref += n
        elif op in "DN":
            ref += n
    return blocks


def exon_read_counts(exons: pd.DataFrame) -> dict[str, dict[int, int]]:
    """Reads with >=1 ALIGNED base in each exon (a read spliced across an exon
    does not count toward it). This is the 'read fraction' metric, which is
    independent of the per-base depth metric."""
    counts: dict[str, dict[int, int]] = {}
    exon_iv = [(int(e["exon_number"]), int(e["start"]), int(e["end"])) for _, e in exons.iterrows()]
    for sample in SAMPLES:
        bam = require(bam_path(sample), f"BAM for {sample}")
        c = {k: 0 for k, _, _ in exon_iv}
        sam = run(["samtools", "view", "-F", "0x900", str(bam),
                   f"{CHROM}:{GENE_START}-{GENE_END}"])
        for line in sam.splitlines():
            if not line:
                continue
            f = line.split("\t")
            blocks = aligned_blocks(int(f[3]), f[5])
            for k, s, e in exon_iv:
                if any(b0 <= e and b1 >= s for b0, b1 in blocks):
                    c[k] += 1
        counts[sample] = c
    return counts


def exon_coverage(exons: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    rows: list[dict] = []
    depth_cache: dict[str, np.ndarray] = {}
    read_counts = exon_read_counts(exons)
    for sample in SAMPLES:
        depth = samtools_depth(sample, GENE_START, GENE_END)
        depth_cache[sample] = depth
        per_exon = {}
        for _, e in exons.iterrows():
            sl = depth[int(e["start"]) - GENE_START: int(e["end"]) - GENE_START + 1]
            per_exon[int(e["exon_number"])] = (sl.sum(), sl.mean(), len(sl))
        total_exonic_depth = sum(v[0] for v in per_exon.values())
        if total_exonic_depth <= 0:
            raise RuntimeError(f"{sample}: zero Smad4 exonic depth - refusing to normalise")
        total_reads = sum(read_counts[sample].values())
        if total_reads <= 0:
            raise RuntimeError(f"{sample}: zero Smad4 exonic reads - refusing to normalise")
        for _, e in exons.iterrows():
            k = int(e["exon_number"])
            s, m, ln = per_exon[k]
            rows.append(
                dict(
                    sample=sample,
                    genotype=GENOTYPE[sample],
                    exon_number=k,
                    exon_id=e["exon_id"],
                    start=int(e["start"]),
                    end=int(e["end"]),
                    exon_length=int(ln),
                    total_base_depth=int(s),
                    mean_depth=float(m),
                    frac_of_smad4_exonic_depth=float(s) / float(total_exonic_depth),
                    sample_total_smad4_exonic_depth=int(total_exonic_depth),
                    n_reads=int(read_counts[sample][k]),
                    frac_of_smad4_exonic_reads=float(read_counts[sample][k]) / float(total_reads),
                    sample_total_smad4_exonic_reads=int(total_reads),
                    is_floxed_exon=bool(
                        int(e["start"]) == FLOXED_EXON_START and int(e["end"]) == FLOXED_EXON_END
                    ),
                )
            )
    return pd.DataFrame(rows), depth_cache


def exact_perm_p(values: np.ndarray, labels: list[str]) -> tuple[float, float, int]:
    """Exact two-sided permutation test on the difference of group means over
    all C(6,3) = 20 balanced 3/3 relabellings.

    NOTE ON THE FLOOR: with n = 3 vs 3 there are only 20 splits, and the
    observed split and its complement always tie the observed statistic, so
    the smallest attainable two-sided P is 2/20 = 0.10.  No exon can ever be
    called significant at P < 0.05 by this design.  A t test would be invalid
    here (n = 3, no distributional support); this exact test is the honest
    ceiling on what the design can say.
    """
    idx = list(range(len(labels)))
    obs_dcko = np.array([i for i, l in enumerate(labels) if l == "dcKO"])
    obs_con = np.array([i for i, l in enumerate(labels) if l == "CON"])
    obs_stat = values[obs_dcko].mean() - values[obs_con].mean()
    stats = []
    for combo in itertools.combinations(idx, len(obs_dcko)):
        a = np.array(combo)
        b = np.array([i for i in idx if i not in combo])
        stats.append(values[a].mean() - values[b].mean())
    stats = np.array(stats)
    n_perm = len(stats)
    p = float((np.abs(stats) >= abs(obs_stat) - 1e-15).sum()) / n_perm
    return p, float(obs_stat), n_perm


def exon_ratio_table(cov: pd.DataFrame, exons: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    order = [s for s in SAMPLES]
    labels = [GENOTYPE[s] for s in order]
    for _, e in exons.iterrows():
        k = int(e["exon_number"])
        sub = cov[cov["exon_number"] == k].set_index("sample")
        vals = np.array([sub.loc[s, "frac_of_smad4_exonic_depth"] for s in order], dtype=float)
        con = vals[[i for i, l in enumerate(labels) if l == "CON"]]
        dck = vals[[i for i, l in enumerate(labels) if l == "dcKO"]]
        p, stat, n_perm = exact_perm_p(vals, labels)
        # Second, INDEPENDENT metric: fraction of exonic reads rather than of
        # per-base depth. Reported so nobody can pick the flattering one.
        rvals = np.array([sub.loc[s, "frac_of_smad4_exonic_reads"] for s in order], dtype=float)
        rcon = rvals[[i for i, l in enumerate(labels) if l == "CON"]]
        rdck = rvals[[i for i, l in enumerate(labels) if l == "dcKO"]]
        rp, rstat, _ = exact_perm_p(rvals, labels)
        rows.append(
            dict(
                exon_number=k,
                exon_id=e["exon_id"],
                start=int(e["start"]),
                end=int(e["end"]),
                exon_length=int(e["length"]),
                is_floxed_exon=bool(
                    int(e["start"]) == FLOXED_EXON_START and int(e["end"]) == FLOXED_EXON_END
                ),
                mean_frac_CON=float(con.mean()),
                mean_frac_dcKO=float(dck.mean()),
                sd_frac_CON=float(con.std(ddof=1)),
                sd_frac_dcKO=float(dck.std(ddof=1)),
                ratio_dcKO_over_CON=float(dck.mean() / con.mean()),
                mean_diff_dcKO_minus_CON=stat,
                perm_p_two_sided=p,
                n_permutations=n_perm,
                perm_p_floor=2.0 / n_perm,
                separates_groups=bool(dck.max() < con.min() or dck.min() > con.max()),
                # --- independent read-fraction metric ---
                mean_readfrac_CON=float(rcon.mean()),
                mean_readfrac_dcKO=float(rdck.mean()),
                ratio_readfrac_dcKO_over_CON=float(rdck.mean() / rcon.mean()),
                mean_readfrac_diff_dcKO_minus_CON=rstat,
                perm_p_two_sided_readfrac=rp,
                separates_groups_readfrac=bool(rdck.max() < rcon.min() or rdck.min() > rcon.max()),
                metrics_agree_on_separation=bool(
                    (dck.max() < con.min() or dck.min() > con.max())
                    == (rdck.max() < rcon.min() or rdck.min() > rcon.max())
                ),
            )
        )
    return pd.DataFrame(rows).sort_values("exon_number").reset_index(drop=True)


# ---------------------------------------------------------------------------
# 3. Junctions
# ---------------------------------------------------------------------------
def annotate_junction(istart: int, iend: int, exons: pd.DataFrame) -> tuple[str, str, str]:
    """Return (donor_annot, acceptor_annot, class). Minus strand: the splice
    DONOR is at the higher coordinate (intron end + 1 = exon start)."""
    left_exon = exons[exons["end"] == istart - 1]     # exon immediately 5' in coordinate space
    right_exon = exons[exons["start"] == iend + 1]    # exon immediately 3' in coordinate space
    # minus strand -> right_exon (higher coord) is the donor side
    donor = f"exon{int(right_exon.iloc[0]['exon_number'])}_end" if len(right_exon) else "novel"
    acceptor = f"exon{int(left_exon.iloc[0]['exon_number'])}_start" if len(left_exon) else "novel"
    if len(right_exon) and len(left_exon):
        dn = int(right_exon.iloc[0]["exon_number"])
        an = int(left_exon.iloc[0]["exon_number"])
        cls = "canonical_adjacent" if an == dn + 1 else "canonical_exon_skip"
    elif len(right_exon) or len(left_exon):
        cls = "one_novel_end"
    else:
        cls = "both_ends_novel"
    return donor, acceptor, cls


def junction_table(exons: pd.DataFrame) -> pd.DataFrame:
    per_sample: dict[str, dict[tuple[int, int], int]] = {}
    for sample in SAMPLES:
        p = require(sj_path(sample), f"SJ.out.tab for {sample}")
        d: dict[tuple[int, int], int] = {}
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] != CHROM:
                    continue
                s, e = int(f[1]), int(f[2])
                if s < JUNCTION_WINDOW_START or e > JUNCTION_WINDOW_END:
                    continue
                d[(s, e)] = int(f[6])   # column 7: uniquely-mapped crossing reads
        per_sample[sample] = d

    keys = sorted(set().union(*(set(d) for d in per_sample.values())) | {SKIP_JUNCTION})
    rows = []
    for (s, e) in keys:
        donor, acceptor, cls = annotate_junction(s, e, exons)
        r = dict(
            intron_start=s,
            intron_end=e,
            intron_length=e - s + 1,
            donor_site=e + 1,        # minus strand: donor is the higher-coordinate exon boundary
            acceptor_site=s - 1,
            donor_annot=donor,
            acceptor_annot=acceptor,
            junction_class=cls,
            is_predicted_cre_skip=bool((s, e) == SKIP_JUNCTION),
        )
        for sample in SAMPLES:
            r[f"{sample}_unique"] = int(per_sample[sample].get((s, e), 0))
        r["total_unique"] = sum(r[f"{x}_unique"] for x in SAMPLES)
        r["CON_unique_sum"] = sum(r[f"{x}_unique"] for x in SAMPLES if GENOTYPE[x] == "CON")
        r["dcKO_unique_sum"] = sum(r[f"{x}_unique"] for x in SAMPLES if GENOTYPE[x] == "dcKO")
        r["n_CON_samples_detected"] = sum(
            1 for x in SAMPLES if GENOTYPE[x] == "CON" and r[f"{x}_unique"] > 0)
        r["n_dcKO_samples_detected"] = sum(
            1 for x in SAMPLES if GENOTYPE[x] == "dcKO" and r[f"{x}_unique"] > 0)
        if r["total_unique"] == 0:
            spec = "not_detected"
        elif r["CON_unique_sum"] > 0 and r["dcKO_unique_sum"] > 0:
            spec = "shared"
        elif r["dcKO_unique_sum"] > 0:
            spec = "dcKO_only"
        else:
            spec = "CON_only"
        r["genotype_specificity"] = spec
        r["touches_floxed_exon"] = bool(
            r["donor_site"] == FLOXED_EXON_START or r["acceptor_site"] == FLOXED_EXON_END
            or r["donor_site"] == FLOXED_EXON_END + 1 or r["acceptor_site"] == FLOXED_EXON_START - 1
        )
        rows.append(r)
    df = pd.DataFrame(rows).sort_values(["intron_start", "intron_end"]).reset_index(drop=True)

    # For each junction, the best-supported CANONICAL junction sharing one of
    # its ends. Lets a novel junction be expressed as a percentage of the
    # normal splicing event it competes with, instead of as a bare read count.
    canon = df[df["junction_class"] == "canonical_adjacent"]
    sib, pct = [], []
    for _, r in df.iterrows():
        share = canon[(canon["intron_start"] == r["intron_start"])
                      | (canon["intron_end"] == r["intron_end"])]
        if r["junction_class"] == "canonical_adjacent" or share.empty:
            sib.append(np.nan)
            pct.append(np.nan)
        else:
            best = int(share["total_unique"].max())
            sib.append(best)
            pct.append(100.0 * r["total_unique"] / best if best else np.nan)
    df["sibling_canonical_total_unique"] = sib
    df["pct_of_sibling_canonical"] = pct
    return df


# ---------------------------------------------------------------------------
# 4. Binned locus depth
# ---------------------------------------------------------------------------
def locus_depth_table(depth_cache: dict[str, np.ndarray], cov: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n = LOCUS_DEPTH_END - LOCUS_DEPTH_START + 1
    n_bins = (n + LOCUS_BIN_BP - 1) // LOCUS_BIN_BP
    for sample in SAMPLES:
        full = depth_cache[sample]
        sl = full[LOCUS_DEPTH_START - GENE_START: LOCUS_DEPTH_END - GENE_START + 1].astype(float)
        total = float(cov.loc[cov["sample"] == sample, "sample_total_smad4_exonic_depth"].iloc[0])
        for b in range(n_bins):
            lo, hi = b * LOCUS_BIN_BP, min((b + 1) * LOCUS_BIN_BP, n)
            chunk = sl[lo:hi]
            rows.append(
                dict(
                    sample=sample,
                    genotype=GENOTYPE[sample],
                    bin_start=LOCUS_DEPTH_START + lo,
                    bin_end=LOCUS_DEPTH_START + hi - 1,
                    bin_bp=hi - lo,
                    mean_depth=float(chunk.mean()),
                    max_depth=float(chunk.max()),
                    norm_depth=float(chunk.mean()) / total,
                )
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. PSI + Smad4 DE stats
# ---------------------------------------------------------------------------
def psi_table(junc: pd.DataFrame) -> pd.DataFrame:
    """PSI = inclusion / (inclusion + skip), inclusion = mean of the two
    exon-flanking junctions."""
    # Inclusion junctions are derived from the floxed exon's own boundaries and
    # the exon model, not hard-coded: the two CANONICAL junctions that abut the
    # floxed exon (its 3' and 5' flanks in genome coordinates).
    rows = []
    skip_row = junc[junc["is_predicted_cre_skip"]].iloc[0]
    abuts = (junc["intron_end"] == FLOXED_EXON_START - 1) | (junc["intron_start"] == FLOXED_EXON_END + 1)
    canon_incl = junc[abuts & (junc["junction_class"] == "canonical_adjacent")]
    if len(canon_incl) != 2:
        raise RuntimeError(
            "expected exactly two canonical inclusion junctions flanking the floxed exon, "
            f"found {len(canon_incl)}:\n{canon_incl[['intron_start', 'intron_end']]}"
        )
    # All junctions that retain the exon, canonical or not - reported alongside
    # so nothing is silently discarded from the denominator.
    all_incl = junc[abuts]

    for sample in SAMPLES:
        skip = int(skip_row[f"{sample}_unique"])
        canon_vals = [int(v) for v in canon_incl[f"{sample}_unique"]]
        # Inclusion support = mean of the two flanking junctions (each read
        # crosses only one of them, so the mean, not the sum, is the number of
        # reads informative about this exon's inclusion).
        incl = float(np.mean(canon_vals))
        rows.append(
            dict(
                sample=sample,
                genotype=GENOTYPE[sample],
                skip_reads=skip,
                inclusion_reads_mean=incl,
                inclusion_junctions=";".join(
                    f"{int(a)}-{int(b)}"
                    for a, b in zip(canon_incl["intron_start"], canon_incl["intron_end"])
                ),
                inclusion_reads_detail=";".join(str(v) for v in canon_vals),
                all_exon_retaining_junction_reads=int(all_incl[f"{sample}_unique"].sum()),
                psi=incl / (incl + skip) if (incl + skip) > 0 else float("nan"),
            )
        )
    return pd.DataFrame(rows)


def smad4_de_stats() -> dict:
    require(DEG_TABLE, "DEG_table.csv")
    deg = pd.read_csv(DEG_TABLE, header=1, float_precision="round_trip", low_memory=False)
    hit = deg[deg["ensembl_gene_id"] == SMAD4_GENE_ID]
    if len(hit) != 1:
        raise RuntimeError(f"expected exactly one {SMAD4_GENE_ID} row in DEG_table.csv, got {len(hit)}")
    r = hit.iloc[0]
    return dict(
        ensembl_gene_id=str(r["ensembl_gene_id"]),
        gene_name=str(r["external_gene_name"]),
        baseMean=float(r["baseMean"]),
        log2FoldChange=float(r["log2FoldChange"]),
        pvalue=float(r["pvalue"]),
        padj=float(r["padj"]),
        is_sig_at_project_threshold=bool(float(r["padj"]) < 0.05 and abs(float(r["log2FoldChange"])) > 1),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
_CODON = {}
for _i, _aa in enumerate("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    _b = "TCAG"
    _CODON[_b[_i // 16] + _b[(_i // 4) % 4] + _b[_i % 4]] = _aa


def coding_consequence() -> dict:
    """Derive what excising the floxed exon does to the reading frame.

    The manuscript states that the recombined transcript is frameshifted, that its
    premature termination codon sits a few bases into the next exon, and that this
    makes it a canonical nonsense-mediated-decay substrate. Those are claims about
    sequence, so they are recomputed here from genome.fa and the annotation on every
    run rather than transcribed once and trusted. Values that reach the manuscript are
    asserted; if the reference or annotation ever changes, this raises.
    """
    def translate(seq: str) -> str:
        prot = ""
        for i in range(0, len(seq) - 2, 3):
            aa = _CODON[seq[i:i + 3]]
            prot += aa
            if aa == "*":
                break
        return prot

    gtf = ensure_smad4_gtf_cache()
    cds = sorted(
        ((int(f[3]), int(f[4])) for f in (l.rstrip("\n").split("\t") for l in open(gtf))
         if len(f) > 8 and f[2] == "CDS"
         and parse_gtf_attrs(f[8]).get("transcript_id") == SMAD4_CANONICAL_TX),
        key=lambda x: -x[0],                       # minus strand: transcription order
    )
    if (FLOXED_EXON_START, FLOXED_EXON_END) not in cds:
        raise SystemExit("[20] floxed exon is not a CDS interval of the canonical transcript")

    wt = "".join(revcomp(fetch_reference_window(a, b)) for a, b in cds)
    mut_iv = [(a, b) for a, b in cds if (a, b) != (FLOXED_EXON_START, FLOXED_EXON_END)]
    mut = "".join(revcomp(fetch_reference_window(a, b)) for a, b in mut_iv)
    wt_p, mut_p = translate(wt), translate(mut)

    exon_cds_len = FLOXED_EXON_END - FLOXED_EXON_START + 1
    upstream = sum(b - a + 1
                   for a, b in cds[:cds.index((FLOXED_EXON_START, FLOXED_EXON_END))])
    stop_codon = len(mut_p)                        # translate() stops at the first '*'
    stop_nt = (stop_codon - 1) * 3 + 1
    last_ejc = sum(b - a + 1 for a, b in mut_iv[:-1])
    ptc_to_last_ejc = last_ejc - stop_nt

    res = dict(
        wildtype_protein_aa=len(wt_p),
        floxed_exon_cds_nt=exon_cds_len,
        frameshift=exon_cds_len % 3 != 0,
        residues_spanned=[upstream // 3 + 1,
                          -(-(upstream + exon_cds_len) // 3)],   # ceil, no math import
        junction_after_cds_nt=upstream,
        first_ptc_codon=stop_codon,
        ptc_nt_into_next_exon=stop_nt - upstream - 1,
        truncated_protein_aa=len(mut_p) - 1,
        novel_residues_after_junction=mut_p[upstream // 3:stop_codon - 1],
        ptc_nt_upstream_of_last_junction=ptc_to_last_ejc,
        predicted_nmd_substrate=ptc_to_last_ejc > 55,
    )

    # Values quoted in the manuscript and Supplementary Methods B7.
    for field, expect in (("wildtype_protein_aa", 551), ("floxed_exon_cds_nt", 184),
                          ("first_ptc_codon", 321), ("ptc_nt_into_next_exon", 8),
                          ("truncated_protein_aa", 320),
                          ("ptc_nt_upstream_of_last_junction", 299)):
        if res[field] != expect:
            raise SystemExit(
                f"[20] coding consequence changed: {field} = {res[field]}, manuscript says "
                f"{expect}. Update the manuscript and Supplementary Methods B7 before "
                f"proceeding - do not silently ship the old number.")
    if not res["frameshift"] or not res["predicted_nmd_substrate"]:
        raise SystemExit("[20] excision is no longer frameshifting/NMD-predicted; text is stale")
    return res


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    exons, blocks = load_exon_model()
    ref_win = fetch_reference_window(LOXP_SEARCH_START, LOXP_SEARCH_END)

    # --- 1. loxP reads ---------------------------------------------------
    lox = collect_loxp_reads(LOXP_SEARCH_START, ref_win)
    n_lox = len(lox)
    per_sample_counts = lox.groupby("sample").size().reindex(SAMPLES, fill_value=0).to_dict()

    assert n_lox == 15, (
        f"expected 15 loxP-containing reads across the six BAMs, found {n_lox}.\n"
        f"per-sample: {per_sample_counts}\n"
        "The BAMs or the search window have changed. Do not proceed with the figure."
    )
    expected_counts = {"CON1": 4, "CON2": 1, "CON3": 4, "dcKO1": 2, "dcKO2": 1, "dcKO3": 3}
    assert per_sample_counts == expected_counts, (
        f"per-sample loxP read counts changed: got {per_sample_counts}, "
        f"expected {expected_counts}"
    )
    assert set(lox["genotype"]) == {"CON", "dcKO"}, (
        "loxP reads must be present in BOTH genotypes (CON is Cre-negative Smad4^fl/fl). "
        f"Found only: {sorted(set(lox['genotype']))}"
    )
    assert set(lox["site"]) == {"A", "B"}, (
        f"expected reads at both cassette sites A and B, got {sorted(set(lox['site']))}"
    )
    money = lox[(lox["sample"] == "CON1") & (lox["pos"] == 73_784_710) & (lox["cigar"] == "70S79M")]
    assert len(money) == 1 and bool(money.iloc[0]["spans_cassette"]), (
        "the cassette-spanning CON1 read (POS 73784710, CIGAR 70S79M) is missing or no longer "
        "carries genomic flank + loxP + genomic flank; site B is no longer resolved to bp"
    )

    # --- 1a. breakpoints and linkers, re-derived from the reads ----------
    sites = derive_insertion_sites(lox, LOXP_SEARCH_START, ref_win)
    for site, const in (("A", SITE_A_POS), ("B", SITE_B_POS)):
        d = sites[site]
        assert d["breakpoint"] == const, (
            f"loxP site {site} is placed by the reads at {d['breakpoint']} "
            f"(cluster {d['breakpoint_min']}-{d['breakpoint_max']}, "
            f"{d['n_reads_used']} reads) but the constant SITE_{site}_POS says {const}. "
            "The constant is wrong, or the alignments changed. Refusing to emit a summary."
        )
    assert sites["A"]["ambiguity_bp"] == 3 and sites["A"]["microhomology"] == "AAC", (
        f"site A ambiguity re-derived as {sites['A']['ambiguity_bp']} bp / "
        f"{sites['A']['microhomology']!r}; established value is 3 bp / 'AAC'")
    assert sites["B"]["ambiguity_bp"] == 0, (
        f"site B is no longer resolved exactly: spread {sites['B']['ambiguity_bp']} bp")

    linkers = derive_cassette_linkers(lox)

    # --- 1b. optical/PCR duplicates --------------------------------------
    lox = flag_duplicate_reads(lox)
    n_dup_copies = int(lox["is_duplicate_copy"].sum())
    n_unique_fragments = int(n_lox - n_dup_copies)
    dup_groups = (lox[lox["duplicate_group_size"] > 1]
                  .groupby(["sample", "pos", "cigar"]).size().to_dict())
    assert not lox.loc[lox["is_duplicate_copy"], "spans_cassette"].any(), (
        "a cassette-SPANNING read is now a duplicate copy of another read; the "
        "spanning-read count would double-count one fragment"
    )
    lox.to_csv(OUT_DIR / "smad4_loxp_reads.csv", index=False)

    # --- 2. exon coverage + ratio ---------------------------------------
    cov, depth_cache = exon_coverage(exons)
    cov.to_csv(OUT_DIR / "smad4_exon_coverage.csv", index=False)

    ratio = exon_ratio_table(cov, exons)
    with open(OUT_DIR / "smad4_exon_ratio.csv", "w") as fh:
        fh.write(
            "# Per-exon share of each sample's total Smad4 exonic depth, CON vs dcKO.\n"
            "# perm_p_two_sided is an EXACT permutation test over all C(6,3) = 20 balanced\n"
            "#   3/3 relabellings of the six samples. The observed split and its complement\n"
            "#   always tie the observed statistic, so the smallest attainable two-sided P is\n"
            "#   2/20 = 0.10. NO EXON CAN REACH P < 0.05 IN THIS DESIGN. A t test would be\n"
            "#   invalid at n = 3; this is the honest ceiling on what six samples can say.\n"
            "# The floxed exon's coverage dip is therefore NOT a significant depletion and\n"
            "#   must not be drawn as one.\n"
            "# TWO INDEPENDENT METRICS are reported: frac of per-base DEPTH (perm_p_two_sided)\n"
            "#   and frac of exonic READS (perm_p_two_sided_readfrac). They do not always agree\n"
            "#   on which exons separate the groups cleanly - see metrics_agree_on_separation.\n"
            "#   Both are given so that neither the figure nor the text can pick the flattering\n"
            "#   one. Neither metric can reach P < 0.05 at n = 3 vs 3.\n"
        )
        ratio.to_csv(fh, index=False)

    flox = ratio[ratio["is_floxed_exon"]].iloc[0]
    for col in ("perm_p_two_sided", "perm_p_two_sided_readfrac"):
        assert float(flox[col]) >= 0.10 - 1e-12, (
            f"{col} for the floxed exon fell below the 2/20 floor - the test is wrong"
        )
        assert float(flox[col]) > 0.05, (
            f"THE NEGATIVE RESULT HAS CHANGED: {col} for the floxed exon is "
            f"{float(flox[col]):.3f}. The exon coverage dip can no longer be described as "
            "non-significant."
        )
    assert 0.90 <= float(flox["ratio_dcKO_over_CON"]) <= 1.0, (
        f"floxed-exon dcKO/CON depth ratio is {flox['ratio_dcKO_over_CON']:.3f}; "
        "established analysis found only an ~8-10% dip. Re-examine before drawing anything."
    )
    # Established: all 12 exons lie within 0.90-1.11 on the depth metric.
    assert ratio["ratio_dcKO_over_CON"].between(0.90, 1.11).all(), (
        "some exon now falls outside the 0.90-1.11 depth-ratio band; the established "
        "'all 12 exons within 0.90-1.11' statement no longer holds:\n"
        f"{ratio.loc[~ratio['ratio_dcKO_over_CON'].between(0.90, 1.11), ['exon_number', 'ratio_dcKO_over_CON']]}"
    )
    # What the data ACTUALLY show (see the DISCREPANCY note in the module
    # docstring): the floxed exon is the ONLY exon of the twelve that separates
    # CON from dcKO cleanly, on BOTH metrics. It still cannot be called
    # significant, because P = 0.10 is the floor of this design. Both halves of
    # that sentence are asserted, so neither can be dropped silently.
    sep_depth = sorted(int(r["exon_number"]) for _, r in ratio.iterrows() if r["separates_groups"])
    sep_reads = sorted(
        int(r["exon_number"]) for _, r in ratio.iterrows() if r["separates_groups_readfrac"])
    flox_n = int(flox["exon_number"])
    assert sep_depth == [flox_n] and sep_reads == [flox_n], (
        "the set of exons that cleanly separate CON from dcKO has changed "
        f"(depth: {sep_depth}, read-fraction: {sep_reads}; expected [{flox_n}] for both). "
        "The exon-coverage panel's caption is no longer accurate."
    )
    flox_dev = abs(float(flox["ratio_dcKO_over_CON"]) - 1.0)
    other_dev = (ratio.loc[~ratio["is_floxed_exon"], "ratio_dcKO_over_CON"] - 1.0).abs()

    # Which exon deviates most from 1, on EACH metric, and in which direction.
    # The two metrics do not agree, so every comparative sentence downstream has
    # to name its metric; these fields exist so no caption has to guess.
    def _dev_summary(col: str) -> dict:
        d = (ratio[col] - 1.0).abs()
        i = int(d.idxmax())
        return dict(exon_number=int(ratio.loc[i, "exon_number"]),
                    exon_id=str(ratio.loc[i, "exon_id"]),
                    ratio=float(ratio.loc[i, col]),
                    deviation=float(d.loc[i]),
                    is_floxed_exon=bool(ratio.loc[i, "is_floxed_exon"]),
                    direction="above 1" if float(ratio.loc[i, col]) > 1 else "below 1")

    most_deviating = dict(depth=_dev_summary("ratio_dcKO_over_CON"),
                          readfrac=_dev_summary("ratio_readfrac_dcKO_over_CON"))
    # The non-floxed exon that moves furthest in the SAME direction as the floxed
    # exon on the depth metric: the single most reassuring comparison available,
    # and the one a caption is most tempted to omit.
    _same = ratio[(~ratio["is_floxed_exon"]) &
                  (np.sign(ratio["ratio_dcKO_over_CON"] - 1.0)
                   == np.sign(float(flox["ratio_dcKO_over_CON"]) - 1.0))]
    if _same.empty:
        same_dir = None
    else:
        _i = int((_same["ratio_dcKO_over_CON"] - 1.0).abs().idxmax())
        same_dir = dict(exon_number=int(ratio.loc[_i, "exon_number"]),
                        ratio=float(ratio.loc[_i, "ratio_dcKO_over_CON"]),
                        deviation=float(abs(ratio.loc[_i, "ratio_dcKO_over_CON"] - 1.0)),
                        perm_p=float(ratio.loc[_i, "perm_p_two_sided"]),
                        deviates_more_than_floxed=bool(
                            abs(ratio.loc[_i, "ratio_dcKO_over_CON"] - 1.0) > flox_dev))

    # --- 3. junctions ----------------------------------------------------
    junc = junction_table(exons)
    with open(OUT_DIR / "smad4_junctions.csv", "w") as fh:
        fh.write(
            f"# Every splice junction reported by STAR in {CHROM}:{JUNCTION_WINDOW_START}-"
            f"{JUNCTION_WINDOW_END}; counts are SJ.out.tab column 7 (UNIQUELY mapped reads).\n"
            "# A row with all-zero counts was reported by STAR from multi-mapping reads only.\n"
            "# is_predicted_cre_skip marks the exon-skip junction Cre would create. It is 0 in\n"
            "#   all six samples - that is the primary negative result, not an omission.\n"
            "# genotype_specificity flags junctions seen in only one genotype. Several exist in\n"
            "#   BOTH directions at ~0.5% of the canonical junction they compete with; they are\n"
            "#   sporadic low-count noise, not a dcKO splicing phenotype. Read\n"
            "#   pct_of_sibling_canonical before interpreting any of them.\n"
        )
        junc.to_csv(fh, index=False)

    skip = junc[junc["is_predicted_cre_skip"]]
    assert len(skip) == 1, "the predicted Cre skip junction row is missing"
    skip_counts = {s: int(skip.iloc[0][f"{s}_unique"]) for s in SAMPLES}
    assert all(v == 0 for v in skip_counts.values()), (
        "THE NEGATIVE RESULT HAS CHANGED: the predicted Cre skip junction "
        f"{SKIP_JUNCTION[0]}-{SKIP_JUNCTION[1]} is no longer zero in all six samples: "
        f"{skip_counts}"
    )
    # --- is there ANY dcKO-specific splicing change? ----------------------
    # Blunt "zero dcKO-only junctions" is not what the data show: the locus
    # carries sporadic low-count novel junctions in BOTH genotypes. The claim
    # that must hold is narrower and stronger:
    #   (i)  no dcKO-specific junction touches the floxed exon at all, and
    #   (ii) every genotype-specific junction is noise-level relative to the
    #        canonical junction it competes with, in both directions.
    dcko_only = junc[junc["genotype_specificity"] == "dcKO_only"]
    con_only = junc[junc["genotype_specificity"] == "CON_only"]

    bad = dcko_only[dcko_only["touches_floxed_exon"]]
    assert len(bad) == 0, (
        "a dcKO-SPECIFIC junction involving a boundary of the floxed exon has appeared. "
        "The 'floxed exon is spliced normally in dcKO' conclusion no longer holds:\n"
        f"{bad[['intron_start', 'intron_end', 'donor_annot', 'acceptor_annot', 'total_unique']]}"
    )
    NOISE_PCT = 2.0
    loud = dcko_only[dcko_only["pct_of_sibling_canonical"].fillna(0) >= NOISE_PCT]
    assert len(loud) == 0, (
        f"a dcKO-specific junction now exceeds {NOISE_PCT}% of its sibling canonical "
        "junction and can no longer be dismissed as sporadic noise:\n"
        f"{loud[['intron_start', 'intron_end', 'total_unique', 'pct_of_sibling_canonical']]}"
    )
    assert len(con_only) > 0, (
        "genotype-specific low-count junctions were found only in dcKO. They can no longer "
        "be presented as bidirectional background noise - re-examine before claiming "
        "'no dcKO-specific junction'."
    )
    max_dcko_only_pct = float(dcko_only["pct_of_sibling_canonical"].max()) if len(dcko_only) else 0.0
    max_con_only_pct = float(con_only["pct_of_sibling_canonical"].max()) if len(con_only) else 0.0

    # --- PSI -------------------------------------------------------------
    psi = psi_table(junc)
    assert (psi["psi"] == 1.0).all(), f"PSI is no longer 1.0 in all six samples:\n{psi}"
    n_incl_total = float(psi["inclusion_reads_mean"].sum())

    # rule of three: with 0 skip events observed out of N informative reads,
    # the one-sided 95% upper bound on the skip rate is 3/N.
    dcko_incl = float(psi.loc[psi["genotype"] == "dcKO", "inclusion_reads_mean"].sum())
    # The literal number of dcKO reads crossing either flanking junction, kept so
    # that no deliverable can present the (deliberately conservative) mean-based
    # denominator above as if it were a read count.
    dcko_junc_reads = int(sum(
        int(v) for _, row in psi[psi["genotype"] == "dcKO"].iterrows()
        for v in str(row["inclusion_reads_detail"]).split(";")))
    skip_upper_95 = 3.0 / dcko_incl
    psi_lower_95 = 1.0 - skip_upper_95
    assert psi_lower_95 >= 0.994, (
        f"the one-sided 95% lower bound on dcKO PSI is now {100 * psi_lower_95:.2f}%, below the "
        "established >=99.4%. Either coverage dropped or the inclusion junctions changed."
    )

    # --- 4. locus depth --------------------------------------------------
    locus = locus_depth_table(depth_cache, cov)
    with open(OUT_DIR / "smad4_locus_depth.csv", "w") as fh:
        fh.write(
            f"# Per-base samtools depth (-a -d 0 -Q 0 -q 0) over {CHROM}:"
            f"{LOCUS_DEPTH_START}-{LOCUS_DEPTH_END}, DOWNSAMPLED to {LOCUS_BIN_BP} bp bins\n"
            f"#   (mean and max of the per-base depths within each bin) to keep the file small.\n"
            "# norm_depth = mean_depth / (that sample's total per-base depth summed over all 12\n"
            "#   canonical Smad4 exons), so tracks are comparable across samples.\n"
        )
        locus.to_csv(fh, index=False)

    # --- 5. summary JSON -------------------------------------------------
    de = smad4_de_stats()
    assert de["log2FoldChange"] > 0, (
        "Smad4 is no longer UP in dcKO; the established +0.317 log2FC has changed"
    )

    flox_ex_row = exons[(exons["start"] == FLOXED_EXON_START) & (exons["end"] == FLOXED_EXON_END)].iloc[0]
    flox_identity = floxed_exon_identity()
    spanning = lox[lox["spans_cassette"]]

    summary = dict(
        generated_by="scripts/JDR_reanalysis/20_smad4_loxp_evidence.py",
        genome_build="GRCm39 / Ensembl release 116",
        gene=dict(
            symbol="Smad4",
            ensembl_gene_id=SMAD4_GENE_ID,
            canonical_transcript=SMAD4_CANONICAL_TX,
            chrom=CHROM,
            start=GENE_START,
            end=GENE_END,
            strand=GENE_STRAND,
            n_canonical_exons=int(len(exons)),
            note="minus strand: exon 1 is the highest coordinate",
        ),
        study_design=STUDY_DESIGN,
        floxed_exon=dict(
            ensembl_exon_id=FLOXED_EXON_ID,
            chrom=CHROM,
            start=FLOXED_EXON_START,
            end=FLOXED_EXON_END,
            length_bp=int(FLOXED_EXON_END - FLOXED_EXON_START + 1),
            canonical_exon_number=int(flox_ex_row["exon_number"]),
            literature_exon_number=8,
            identity=flox_identity,
            note=(
                f"The same 184 bp carries exon_number "
                f"{', '.join(str(v) for v in flox_identity['exon_number_values'][:-1])} or "
                f"{flox_identity['exon_number_values'][-1]} depending on the transcript; the "
                "literature says 'exon 8' and the Ensembl canonical transcript says 'exon 9'. "
                "Lead with the COORDINATE INTERVAL: the stable exon ID resolves the "
                f"exon-number ambiguity for {flox_identity['n_transcripts_with_primary_exon_id']} "
                f"of the {flox_identity['n_transcripts']} transcripts that carry this exon, but "
                f"the remaining {flox_identity['n_transcripts'] - flox_identity['n_transcripts_with_primary_exon_id']} "
                f"label the identical interval {', '.join(flox_identity['other_exon_ids'])}."
            ),
            protein="encodes aa 318-379, N-terminal MH2 domain; entirely CDS",
        ),
        loxp_cassettes=dict(
            core_sequence_34bp=LOXP_FWD,
            site_A=dict(
                insertion_between=[SITE_A_POS, SITE_A_POS + 1],
                ambiguity_bp=sites["A"]["ambiguity_bp"],
                ambiguity_reason=(
                    f"{sites['A']['microhomology']} microhomology at "
                    f"{sites['A']['breakpoint_min'] + 1:,}-{sites['A']['breakpoint_max']:,}, "
                    "tandemly repeated immediately 3' of it, so the alignment can place the "
                    "boundary anywhere in that 3 bp"
                ),
                derivation=sites["A"],
                **{k: v for k, v in linkers["A"].items()},
                location=(
                    "intronic; on the LOWER-coordinate side of the floxed exon, which on the "
                    "minus strand is DOWNSTREAM (3') of it in transcription order"
                ),
            ),
            site_B=dict(
                insertion_between=[SITE_B_POS, SITE_B_POS + 1],
                ambiguity_bp=sites["B"]["ambiguity_bp"],
                ambiguity_reason=(
                    f"resolved exactly: {sites['B']['n_reads_used']} reads place the breakpoint "
                    "at the same base, and one CON1 read (POS 73784710, CIGAR 70S79M) carries "
                    "genomic sequence on both sides of the cassette"
                ),
                derivation=sites["B"],
                **{k: v for k, v in linkers["B"].items()},
                location=(
                    "intronic; on the HIGHER-coordinate side of the floxed exon, which on the "
                    "minus strand is UPSTREAM (5') of it in transcription order"
                ),
            ),
            floxed_interval=dict(
                start=FLOXED_INTERVAL_START,
                end=FLOXED_INTERVAL_END,
                length_bp=int(FLOXED_INTERVAL_END - FLOXED_INTERVAL_START + 1),
                n_exons_contained=1,
            ),
            read_provenance=(
                "Both cassettes are INTRONIC, so every loxP-containing read in an mRNA-seq "
                "library comes from unspliced pre-mRNA or from residual genomic DNA. That is "
                "why there are only a handful of them, and it is why their number carries no "
                "information about transcript abundance or about excision efficiency."
            ),
            n_loxp_reads_total=int(n_lox),
            n_loxp_reads_per_sample=per_sample_counts,
            n_duplicate_copies=n_dup_copies,
            n_unique_fragments=n_unique_fragments,
            duplicate_groups={f"{k[0]} POS {k[1]} {k[2]}": int(v) for k, v in dup_groups.items()},
            duplicate_note=(
                f"{n_dup_copies} of the {n_lox} loxP reads is a byte-identical copy of another "
                "read in the same library (same POS, same CIGAR, same SEQ) and is one fragment "
                f"sequenced twice, so the number of INDEPENDENT fragments is {n_unique_fragments}. "
                "No cassette-spanning read is affected. Nothing in the analysis depends on the "
                "count, but a figure that draws evidence reads must say so."
                if n_dup_copies else "No duplicate loxP reads."
            ),
            n_cassette_spanning_reads=int(len(spanning)),
            n_cassette_spanning_reads_by_site=spanning.groupby("site").size().to_dict(),
            n_loxp_reads_by_site=lox.groupby("site").size().to_dict(),
            cassette_spanning_reads=[
                f"{r['sample']}:{r['qname']}@{r['pos']}({r['cigar']})"
                for _, r in spanning.iterrows()
            ],
            present_in_both_genotypes=True,
            legend_caveat=(
                "CON animals are Cre-negative Smad4^fl/fl^, so they carry the identical "
                "loxP cassettes. loxP reads in CON are expected and are NOT a control failure."
            ),
        ),
        negative_result=dict(
            predicted_cre_skip_junction=dict(
                intron_start=SKIP_JUNCTION[0],
                intron_end=SKIP_JUNCTION[1],
                uniquely_mapped_reads_per_sample=skip_counts,
            ),
            psi_per_sample={r["sample"]: float(r["psi"]) for _, r in psi.iterrows()},
            psi_all_samples=1.0,
            inclusion_junctions=str(psi.iloc[0]["inclusion_junctions"]),
            inclusion_reads_per_sample={
                r["sample"]: str(r["inclusion_reads_detail"]) for _, r in psi.iterrows()
            },
            dcko_informative_reads=dcko_incl,
            dcko_inclusion_equivalent_reads=dcko_incl,
            dcko_junction_reads_total=dcko_junc_reads,
            denominator_definition=(
                f"{dcko_incl:.0f} is NOT a count of reads. Each dcKO library contributes the "
                "MEAN of its two exon-flanking junction counts, because a 150 bp read crosses "
                "only one of the two junctions and the two counts are therefore two partial "
                f"views of the same inclusion event. The actual number of dcKO reads crossing "
                f"either flanking junction is {dcko_junc_reads}. Using the mean halves the "
                f"denominator and WIDENS the interval (bound {100 * psi_lower_95:.2f}% instead "
                f"of {100 * (1 - 3.0 / dcko_junc_reads):.2f}%), so the reported bound is the "
                "conservative one."
            ),
            rule_of_three_upper_bound_on_skip=skip_upper_95,
            one_sided_95_lower_bound_on_psi=psi_lower_95,
            statement=(
                f"Zero exon-skipping junctions were observed against {dcko_incl:.0f} "
                f"inclusion-equivalent dcKO reads ({dcko_junc_reads} junction reads, counted as "
                "the mean of the two flanking junctions per library); by the rule of three at "
                f"least {100 * psi_lower_95:.1f}% of dcKO Smad4 mRNA retains the floxed exon."
            ),
            bulk_tissue_caveat=STUDY_DESIGN["caveat"],
            junction_scan_window=[JUNCTION_WINDOW_START, JUNCTION_WINDOW_END],
            n_junctions_in_window=int(len(junc)),
            no_dcko_specific_junction_at_floxed_exon=True,
            genotype_specific_junctions=dict(
                n_dcKO_only=int(len(dcko_only)),
                n_CON_only=int(len(con_only)),
                max_dcKO_only_pct_of_sibling_canonical=max_dcko_only_pct,
                max_CON_only_pct_of_sibling_canonical=max_con_only_pct,
                dcKO_only_list=[
                    f"{int(r['intron_start'])}-{int(r['intron_end'])}"
                    f" ({int(r['total_unique'])} reads,"
                    f" {r['pct_of_sibling_canonical']:.2f}% of canonical)"
                    for _, r in dcko_only.iterrows()
                ],
                CON_only_list=[
                    f"{int(r['intron_start'])}-{int(r['intron_end'])}"
                    f" ({int(r['total_unique'])} reads,"
                    f" {r['pct_of_sibling_canonical']:.2f}% of canonical)"
                    for _, r in con_only.iterrows()
                ],
                scan_window=(
                    f"chr{CHROM}:{JUNCTION_WINDOW_START:,}-{JUNCTION_WINDOW_END:,} "
                    f"({(JUNCTION_WINDOW_END - JUNCTION_WINDOW_START) / 1000:.0f} kb around the "
                    "floxed exon). These counts are NOT gene-wide; the window must be stated "
                    "wherever they are quoted."
                ),
                interpretation=(
                    "Sporadic novel junctions occur in BOTH genotypes at well under 2% of the "
                    "canonical junction they compete with, and none of them involves the floxed "
                    "exon. They are low-count background, not a dcKO splicing phenotype. Do not "
                    "state 'no dcKO-specific junction exists' without this qualification, and "
                    "do not quote the counts without the scan window."
                ),
            ),
        ),
        exon_coverage=dict(
            depth_metric=dict(
                floxed_exon_ratio_dcKO_over_CON=float(flox["ratio_dcKO_over_CON"]),
                floxed_exon_perm_p=float(flox["perm_p_two_sided"]),
                min_exon_ratio=float(ratio["ratio_dcKO_over_CON"].min()),
                max_exon_ratio=float(ratio["ratio_dcKO_over_CON"].max()),
                exons_with_clean_group_separation=[
                    int(r["exon_number"]) for _, r in ratio.iterrows() if r["separates_groups"]
                ],
            ),
            read_fraction_metric=dict(
                floxed_exon_ratio_dcKO_over_CON=float(flox["ratio_readfrac_dcKO_over_CON"]),
                floxed_exon_perm_p=float(flox["perm_p_two_sided_readfrac"]),
                min_exon_ratio=float(ratio["ratio_readfrac_dcKO_over_CON"].min()),
                max_exon_ratio=float(ratio["ratio_readfrac_dcKO_over_CON"].max()),
                exons_with_clean_group_separation=[
                    int(r["exon_number"]) for _, r in ratio.iterrows()
                    if r["separates_groups_readfrac"]
                ],
            ),
            perm_p_floor=float(flox["perm_p_floor"]),
            n_permutations=int(flox["n_permutations"]),
            exons_where_metrics_disagree_on_separation=[
                int(r["exon_number"]) for _, r in ratio.iterrows()
                if not r["metrics_agree_on_separation"]
            ],
            floxed_exon_largest_depth_deviation_from_1=bool(flox_dev >= float(other_dev.max())),
            largest_other_exon_depth_deviation_from_1=float(other_dev.max()),
            most_deviating_exon=most_deviating,
            same_direction_comparison_exon=same_dir,
            exons_separating_groups_depth=sep_depth,
            exons_separating_groups_readfrac=sep_reads,
            cross_metric_warning=(
                "THE TWO METRICS RANK THE EXONS DIFFERENTLY. On per-base depth the largest "
                f"deviation from 1 belongs to exon {most_deviating['depth']['exon_number']} "
                f"({most_deviating['depth']['deviation']:.3f}); on read fraction it belongs to "
                f"exon {most_deviating['readfrac']['exon_number']} "
                f"({most_deviating['readfrac']['deviation']:.3f}). Any sentence of the form "
                "'exon X deviates more than the floxed exon' is therefore metric-specific and "
                "must name its metric. Both metrics are reported so that neither the figure nor "
                "the text can quote only the flattering one."
            ),
            caveat=(
                "The floxed exon is the ONLY one of the twelve whose CON and dcKO values do not "
                "overlap, on both the depth and the read-fraction metric. It nevertheless CANNOT "
                "be called significant: with n = 3 vs 3 the exact permutation test bottoms out at "
                "P = 2/20 = 0.10, so clean separation is the strongest result the design can "
                "produce and it still does not reach P < 0.05. The effect is small (~8% on depth, "
                "~11% on reads) and every exon stays within 0.90-1.11 on depth. Present this as "
                "'suggestive but underpowered', never as a demonstrated depletion, and never with "
                "a significance marker. Equally, do NOT assert the absence of a depletion: this "
                "design can neither establish nor exclude an effect of this size. The defensible "
                "wording is 'no DEMONSTRABLE depletion' / 'no depletion this design can "
                "establish', never 'there is no depletion'."
            ),
            discrepancy_with_prior_analysis=(
                "A prior pass recorded P = 0.20 with dcKO3 > CON3 on a read-fraction metric, and "
                "an 'equal and opposite clean separation' at exon 6. Neither reproduces here "
                "under per-base depth fraction, read fraction, or unique-reads-only (MAPQ>=60) "
                "read fraction: the floxed exon separates cleanly (P = 0.10) on all three and "
                "exon 6 never does (P = 0.20-0.30, overlapping groups). The prior metric "
                "definition could not be recovered. The shared conclusion - not significant, do "
                "not draw as a depletion - is unaffected, but the retracted sentence must not be "
                "quoted. Note a numerical coincidence that will look alarming: exon 6's "
                "recomputed two-sided permutation P on the depth metric is also 0.20. That 0.20 "
                "is a freshly computed value attached to an exon whose groups OVERLAP (it does "
                "not separate); it is not the retracted statistic, which claimed a clean "
                "separation at exon 6."
            ),
        ),
        smad4_differential_expression=de,
        coding_consequence=coding_consequence(),
        methods=dict(
            depth_command="samtools depth -a -d 0 -Q 0 -q 0",
            junction_source="STAR *.SJ.out.tab column 7 (uniquely mapped crossing reads)",
            loxp_search_window=[LOXP_SEARCH_START, LOXP_SEARCH_END],
            loxp_orientations_searched=[LOXP_FWD, LOXP_REV],
            genomic_anchor_min_bp=MIN_GENOMIC_ANCHOR,
            locus_depth_bin_bp=LOCUS_BIN_BP,
        ),
    )
    with open(OUT_DIR / "smad4_loxp_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=False)
        fh.write("\n")

    # --- console report --------------------------------------------------
    print(f"[20] loxP reads              : {n_lox}  {per_sample_counts}")
    print(f"[20] cassette-spanning reads : {len(spanning)}")
    print(f"[20] Cre skip junction reads : {skip_counts}")
    print(f"[20] PSI                     : {sorted(set(psi['psi']))}")
    print(f"[20] 95% lower bound on PSI  : {100 * psi_lower_95:.2f}%  "
          f"(0/{dcko_incl:.0f} dcKO informative reads)")
    print(f"[20] floxed exon dcKO/CON    : depth {flox['ratio_dcKO_over_CON']:.4f} "
          f"P={flox['perm_p_two_sided']:.2f} | reads {flox['ratio_readfrac_dcKO_over_CON']:.4f} "
          f"P={flox['perm_p_two_sided_readfrac']:.2f}  (floor {flox['perm_p_floor']:.2f})")
    print(f"[20] exons separating groups : depth {sep_depth} | reads {sep_reads}  "
          f"-> only the floxed exon, but P = floor, so NOT significant")
    print(f"[20] Smad4 DE                : log2FC {de['log2FoldChange']:+.4f}  "
          f"padj {de['padj']:.3e}")
    for f in ["smad4_loxp_reads.csv", "smad4_exon_coverage.csv", "smad4_exon_ratio.csv",
              "smad4_junctions.csv", "smad4_locus_depth.csv", "smad4_loxp_summary.json"]:
        print(f"[20] wrote {OUT_DIR / f}")


if __name__ == "__main__":
    main()
