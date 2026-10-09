
    
  
#!/usr/bin/env python3
"""Calculate 12 reference-to-read substitution rates from primary BAM alignments.

Requires Python >=3.8 and samtools on PATH, e.g.:
    conda create -n bam-errors -c conda-forge -c bioconda python=3.11 samtools
    conda activate bam-errors
    python bam_error_rates.py sample1.bam sample2.bam -o error_rates.tsv

Each BAM is converted to a temporary SAM using samtools view. The SAM is
removed after analysis; --sam-dir retains conversions in the chosen directory.
Unmapped (0x4), secondary (0x100), and supplementary (0x800) records are excluded.
MD tags must be present and consistent with CIGAR and SEQ on included records.

A_to_C = count(reference A, read C) / count(reference A, read A/C/G/T).
Only aligned positions with canonical A/C/G/T bases on both sides are counted.
Indels, reference skips, and clipping are excluded from numerator and denominator.
Rates are fractions, not percentages; a zero denominator is written as NA.
Directions use reference orientation, including for reverse-strand alignments.
These are reference mismatch rates: biological variants and alignment errors
can contribute. No expected rates or pass/fail thresholds are assumed.
No mapping/base-quality, duplicate, or QC-fail filtering is applied beyond the
primary mapped alignment requirement. Paired-end overlap is counted per read.
"""

import argparse
import collections
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

BASES = "ACGT"
SUBSTITUTIONS = [(a, b) for a in BASES for b in BASES if a != b]
EXCLUDE_FLAGS = 0x4 | 0x100 | 0x800
CIGAR_PATTERN = re.compile(r"(?:[1-9][0-9]*[MIDNSHP=X])+")
MD_PATTERN = re.compile(r"[0-9]+(?:(?:[A-Z]|\^[A-Z]+)[0-9]+)*")


def alignment_columns(sequence, cigar):
    """Yield (read base, CIGAR operation) for columns represented in MD."""
    if not CIGAR_PATTERN.fullmatch(cigar):
        raise ValueError("missing or invalid CIGAR")
    operations = [(int(n), op) for n, op in re.findall(r"([0-9]+)([MIDNSHP=X])", cigar)]
    consumed = sum(n for n, op in operations if op in "MIS=X")
    if consumed != len(sequence):
        raise ValueError("CIGAR read length does not match SEQ length")
    index = 0
    for length, op in operations:
        if op in "M=X":
            for base in sequence[index:index + length]:
                yield base, op
            index += length
        elif op in "IS":
            index += length
        elif op == "D":
            for _ in range(length):
                yield None, op
        # H, P, and N consume no SEQ bases and are absent from MD.


def count_alignment(sequence, cigar, md, denominators, substitutions):
    """Walk MD and CIGAR together and update counts; reject inconsistencies."""
    if sequence == "*":
        raise ValueError("included alignment has no SEQ")
    sequence = sequence.upper()
    if "=" in sequence:
        raise ValueError("SEQ uses '=' bases; explicit bases or a reference are needed")
    if not MD_PATTERN.fullmatch(md):
        raise ValueError("invalid MD tag")
    columns = iter(alignment_columns(sequence, cigar))

    def take():
        try:
            return next(columns)
        except StopIteration:
            raise ValueError("MD extends beyond CIGAR") from None

    for token in re.findall(r"[0-9]+|\^[A-Z]+|[A-Z]", md):
        if token.isdigit():
            for _ in range(int(token)):
                base, op = take()
                if op not in "M=":
                    raise ValueError("MD match conflicts with CIGAR")
                if base in BASES:
                    denominators[base] += 1
        elif token.startswith("^"):
            for _ in token[1:]:
                _, op = take()
                if op != "D":
                    raise ValueError("MD deletion conflicts with CIGAR")
        else:
            base, op = take()

