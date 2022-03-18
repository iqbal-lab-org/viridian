from collections import defaultdict

import pysam
from viridian_workflow.primers import AmpliconSet
from viridian_workflow.readstore import PairedReads, SingleRead


def score(matches, mismatches):
    """Assign winning amplicon set id based on match stats"""
    amplicon_sets = set([*matches.keys(), *mismatches.keys()])

    m = 0
    winner = None
    for amplicon_set in amplicon_sets:
        print(amplicon_set.name, mismatches[amplicon_set], matches[amplicon_set])
        if matches[amplicon_set] >= m:
            winner = amplicon_set
            m = matches[amplicon_set]
    return winner


def match_read_to_amplicons(read, amplicon_sets):
    if read.is_unmapped:
        return None
    matches = {}
    for amplicons in amplicon_sets:
        m = amplicons.match(*read_interval(read))
        if m:
            matches[amplicons.name] = m
    return matches


def match_reads(reads, amplicon_sets):
    """given a stream of reads, yield reads with a set of matched amplicons"""
    for read in reads:
        if read.is_unmapped:
            continue

        matches = match_read_to_amplicons(read, amplicon_sets)
        yield read, matches


def amplicon_set_counts_to_naive_total_counts(scheme_counts):
    counts = defaultdict(int)
    for scheme_tuple, count in scheme_counts.items():
        for scheme in scheme_tuple:
            counts[scheme] += count
    return counts


def amplicon_set_counts_to_json_friendly(scheme_counts):
    dict_out = {}
    for k, v in scheme_counts.items():
        new_key = ";".join(sorted([str(x) for x in k]))
        dict_out[new_key] = v
    return dict_out
