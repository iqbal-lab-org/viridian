#!/usr/bin/env python3
import sys
from collections import defaultdict
from intervaltree import Interval

import pysam
from primers import AmpliconSet


def score(matches):
    """Assign winning amplicon set id based on match stats
    """

    # naive: take max of all bins
    m = 0
    winner = None
    for k, v in matches.items():
        if v >= m:
            m = v
            winner = k
    return winner


def read_interval(read):
    """ determine template start and end coords for either a single read or
    paired reads
    """
    if read.is_paired:
        if not read.is_reverse:
            start = read.reference_start
            end = read.reference_start + read.template_length
            return start, end

        else:
            start = read.next_reference_start
            end = read.next_reference_start - read.template_length
            return start, end
    else:
        return read.reference_start, read.reference_end


def annotate_read(read, match):
    """Set the match details for a read 
    """
    raise NotImplementedError


def match_reads(reads, amplicon_sets):
    """given a stream of reads, yield reads with a set of matched amplicons
    """
    for read in reads:
        if read.is_unmapped:
            continue

        matches = {}
        for amplicons in amplicon_sets:
            m = amplicons.match(*read_interval(read))
            if m:
                matches[amplicons.name] = m

        yield read, matches


def detect(amplicon_sets, reads):
    """Generate amplicon match stats and identify closest set of
    matching amplicons
    """

    matches = {}
    for aset in amplicon_sets:
        matches[aset.name] = 0
        # other stats for stuff like amplicon containment and
        # ambiguous match counts

    for read, amplicon_matches in match_reads(reads, amplicon_sets):
        for a in amplicon_matches:
            # unambiguous match for the read
            if len(amplicon_matches[a]) == 1:
                matches[a] += 1

    return score(matches)


if __name__ == "__main__":
    amplicons = [AmpliconSet.from_tsv(tsv) for tsv in sys.argv[1].split(",")]
    reads = pysam.AlignmentFile(sys.argv[2], "rb")
    print(detect(amplicons, reads))
