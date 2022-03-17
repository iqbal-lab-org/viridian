from collections import defaultdict

import pysam
from viridian_workflow.primers import AmpliconSet
from viridian_workflow.readstore import PairedReads, SingleRead, read_from_pysam


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
            yield SingleRead(read_from_pysam(read))

        if not read.is_proper_pair:
            continue

        if read.is_read1:
            stats["template_lengths"][abs(read.template_length)] += 1
            reads_by_name[read.query_name] = read_from_pysam(read)

        elif read.is_read2:
            read1 = reads_by_name[read.query_name]
            yield PairedReads(read1, read_from_pysam(read))
            del reads_by_name[read.query_name]


def gather_stats_from_bam(infile, amplicon_sets):
    # open_mode_in = "r" + pysam_open_mode(infile)
    print("gathering stats")
    aln_file_in = pysam.AlignmentFile(infile, "rb")

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

    mismatches = defaultdict(int)
    matches = defaultdict(int)

    for fragment in syncronise_fragments(aln_file_in, stats):
        for amplicon_set in amplicon_sets:
            hit = amplicon_set.match(fragment)
            if hit:
                matches[amplicon_set] += 1
            else:
                mismatches[amplicon_set] += 1

    aln_file_in.close()

    stats["match_any_amplicon"] = match_any_amplicon
    stats["amplicon_scheme_set_matches"] = amplicon_scheme_set_matches
    stats["amplicon_scheme_simple_counts"] = amplicon_set_counts_to_naive_total_counts(
        stats["amplicon_scheme_set_matches"]
    )
    stats["chosen_amplicon_scheme"] = score(matches, mismatches)
    return stats
