#!/usr/bin/env python3
from collections import defaultdict

import pysam
from viridian_workflow.primers import AmpliconSet

from readstore import PairedReads, SingleRead


def score(matches, off_target):
    """Assign winning amplicon set id based on match stats"""

    # naive: take max of all bins
    m = 0
    winner = None
    for k, v in matches.items():
        if v >= m:
            m = v
            winner = k
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


def syncronise_fragments(reads, stats):
    infile_is_paired = None
    reads_by_name = {}

    for read in reads:
        if infile_is_paired is None:
            infile_is_paired = read.is_paired
        else:
            if infile_is_paired != read.is_paired:
                raise Exception("Mix of paired and unpaired reads.")

        if read.is_secondary or read.is_supplementary:
            continue

        stats["total_reads"] += 1
        if read.is_read1:
            stats["reads1"] += 1
        elif read.is_read2:
            stats["reads2"] += 1

        if read.is_unmapped:
            continue

        stats["read_lengths"][read.query_length] += 1
        stats["mapped"] += 1

        if not read.is_paired:
            stats["unpaired_reads"] += 1
            stats["template_lengths"][abs(read.query_length)] += 1  # TODO: check this
            yield SingleRead(read)

        if not read.is_proper_pair:
            continue

        vwf_read = Read(
            read.sequence,
            read.reference_start,
            read.reference_end,
            read.qry_start,
            read.qry_end,
            read.is_reverse,
        )

        if read.is_read1:
            stats["template_lengths"][abs(read.template_length)] += 1
            reads_by_name[read.query_name] = vwf_read

        elif read.is_read2:
            read1 = reads_by_name[read.query_name]
            yield PairedReads(read1, vwf_read)
            del reads_by_name[read.query_name]


def gather_stats_from_bam(infile, bam_out, amplicon_sets):
    open_mode_in = "r" + pysam_open_mode(infile)
    aln_file_in = pysam.AlignmentFile(infile, open_mode_in)
    if bam_out is not None:
        open_mode_out = "w" + pysam_open_mode(bam_out)
        aln_file_out = pysam.AlignmentFile(bam_out, open_mode_out, template=aln_file_in)

    match_any_amplicon = 0
    amplicon_scheme_set_matches = defaultdict(int)

    stats = {
        "unpaired_reads": 0,
        "reads1": 0,
        "reads2": 0,
        "total_reads": 0,
        "mapped": 0,
        "read_lengths": defaultdict(int),
        "template_lengths": defaultdict(int),
    }

    for fragment in syncronise_fragments(aln_file_in, stats):
        for amplicon_set in amplicon_sets:
            amplicon_matches[Amplicon_set] = amplicon_set.match(fragment)

        amplicon_matches = match_read_to_amplicons(read, amplicon_sets)

        if amplicon_matches:
            match_any_amplicon += 1
            read = set_tags(amplicon_sets, read, amplicon_matches)

            amplicon_key = tuple(sorted(list(amplicon_matches.keys())))
            amplicon_scheme_set_matches[amplicon_key] += 1
            if mate:
                mate = set_tags(amplicon_sets, mate, amplicon_matches)

        if bam_out is not None:
            aln_file_out.write(read)
            if mate:
                aln_file_out.write(mate)

    aln_file_in.close()
    if bam_out is not None:
        aln_file_out.close()

    stats["match_any_amplicon"] = match_any_amplicon
    stats["amplicon_scheme_set_matches"] = amplicon_scheme_set_matches
    stats["amplicon_scheme_simple_counts"] = amplicon_set_counts_to_naive_total_counts(
        stats["amplicon_scheme_set_matches"]
    )
    stats["chosen_amplicon_scheme"] = score(stats["amplicon_scheme_simple_counts"])
    return stats
