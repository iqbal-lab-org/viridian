#!/usr/bin/env python3
from collections import defaultdict

import pysam
from primers import AmpliconSet, set_tags, get_tags


def score(matches):
    """Assign winning amplicon set id based on match stats"""

    # naive: take max of all bins
    m = 0
    winner = None
    print(matches.items())
    for k, v in matches.items():
        if v >= m:
            m = v
            winner = k
    return winner


def read_interval(read):
    """determine template start and end coords for either a single read or
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


def pysam_open_mode(filename):
    if filename.endswith(".sam"):
        return ""
    elif filename.endswith(".bam"):
        return "b"
    else:
        raise Exception(f"Filename {filename} does not end with .sam or .bam")


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


def gather_stats_from_bam(infile, bam_out, amplicon_sets):
    open_mode_in = "r" + pysam_open_mode(infile)
    aln_file_in = pysam.AlignmentFile(infile, open_mode_in)
    if bam_out is not None:
        open_mode_out = "w" + pysam_open_mode(bam_out)
        aln_file_out = pysam.AlignmentFile(bam_out, open_mode_out, template=aln_file_in)

    stats = {
        "unpaired_reads": 0,
        "reads1": 0,
        "reads2": 0,
        "total_reads": 0,
        "mapped": 0,
        "match_any_amplicon": 0,
        "read_lengths": defaultdict(int),
        "amplicon_scheme_set_matches": defaultdict(int),
    }
    infile_is_paired = None

    for read in aln_file_in:
        if read.is_secondary or read.is_supplementary:
            continue

        if infile_is_paired is None:
            infile_is_paired = read.is_paired
        elif read.is_paired != infile_is_paired:
            raise Exception("Reads must be all paired or all unpaired")

        stats["total_reads"] += 1
        stats["read_lengths"][read.query_length] += 1
        if not read.is_unmapped:
            stats["mapped"] += 1

        if read.is_paired:
            if read.is_read1:
                if amplicon_matches is not None:
                    raise Exception(
                        "Paired reads not in expected order. Cannot continue"
                    )
                stats["reads1"] += 1
                amplicon_matches = match_read_to_amplicons(read, amplicon_sets)
                read = set_tags(amplicon_sets, read, amplicon_matches)

            else:
                if amplicon_matches is None:
                    raise Exception(
                        "Paired reads not in expected order. Cannot continue"
                    )
                stats["reads2"] += 1
                read = set_tags(amplicon_sets, read, amplicon_matches)
                amplicon_matches = None
        else:
            stats["unpaired_reads"] += 1
            amplicon_matches = match_read_to_amplicons(read, amplicon_sets)
            read = set_tags(amplicon_sets, read, amplicon_matches)

        if amplicon_matches is not None and len(amplicon_matches) > 0:
            stats["match_any_amplicon"] += 1
            key = tuple(sorted(list(amplicon_matches.keys())))
            stats["amplicon_scheme_set_matches"][key] += 1

        if bam_out is not None:
            aln_file_out.write(read)

    aln_file_in.close()
    if bam_out is not None:
        aln_file_out.close()

    stats["amplicon_scheme_simple_counts"] = amplicon_set_counts_to_naive_total_counts(
        stats["amplicon_scheme_set_matches"]
    )
    stats["chosen_amplicon_scheme"] = score(stats["amplicon_scheme_simple_counts"])
    return stats


if __name__ == "__main__":
    raise NotImplementedError
