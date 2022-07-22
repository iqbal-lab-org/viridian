from intervaltree import Interval
import os
import pytest
import random
import subprocess
from unittest import mock

import pyfastaq

from viridian_workflow import primers, readstore
from viridian_workflow.subtasks import Minimap
from viridian_workflow.readstore import PairedReads, SingleRead, Bam

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "detect_primers")


def match_read_to_amplicon_sets(read, amplicon_sets):
    matches = {}
    for amplicon_set in amplicon_sets:
        read_match = amplicon_set.match(read)
        if read_match:
            matches[amplicon_set.name] = read_match.name
    return matches


def test_read_interval_paired():
    read1 = mock.Mock()
    read1.is_reverse = False
    read1.ref_start = 42
    read1.ref_end = 142

    read2 = mock.Mock()
    read2.is_reverse = True
    read2.ref_start = 142
    read2.ref_end = 242

    fragment = PairedReads(read1, read2)
    interval = fragment.ref_end - fragment.ref_start
    assert fragment.ref_start == 42
    assert fragment.ref_end == 242
    assert interval == 200


def test_read_interval_unpaired():
    read = mock.Mock()
    read.is_reverse = False
    read.ref_start = 42
    read.ref_end = 300
    fragment = SingleRead(read)

    assert (fragment.ref_start, fragment.ref_end) == (42, 300)


def test_match_read_to_amplicons():
    tsv_files = {
        "scheme1": os.path.join(data_dir, "match_reads_to_amplicons_1.tsv"),
        "scheme2": os.path.join(data_dir, "match_reads_to_amplicons_2.tsv"),
    }
    amplicon_sets = [
        primers.AmpliconSet.from_tsv(v, name=k) for k, v in tsv_files.items()
    ]
    amplicons1 = primers.AmpliconSet.from_tsv(tsv_files["scheme1"], name="scheme1")
    amplicons2 = primers.AmpliconSet.from_tsv(tsv_files["scheme2"], name="scheme2")

    # scheme1:
    # 100-300, 290-800, 790-1000
    #
    # scheme2:
    # 100-300, 290-500, 490-700, 790-1001
    fragment = mock.Mock()
    fragment.ref_start = None
    fragment.ref_end = None
    matches = {}

    # this test is no longer valid. fragments must be aligned.
    # for amplicon_set in amplicon_sets:
    #    matches[amplicon_set.name] = amplicon_set.match(fragment)
    # assert matches == {}

    fragment.ref_start = 100
    fragment.ref_end = 290
    assert match_read_to_amplicon_sets(fragment, amplicon_sets) == {
        "scheme1": "amp1",
        "scheme2": "amp1",
    }

    fragment.ref_start = 250
    fragment.ref_end = 400
    assert match_read_to_amplicon_sets(fragment, amplicon_sets) == {}

    fragment.ref_start = 100
    fragment.ref_end = 500
    assert match_read_to_amplicon_sets(fragment, amplicon_sets) == {}

    fragment.ref_start = 400
    fragment.ref_end = 750
    assert match_read_to_amplicon_sets(fragment, amplicon_sets) == {
        "scheme1": "amp2",
    }


def _write_sim_reads(ref_seq, coords, outfile, suffix="", revcomp=False):
    with open(outfile, "w") as f:
        for i, (start, end) in enumerate(coords):
            read = pyfastaq.sequences.Fasta(f">read.{i}{suffix}", ref_seq[start:end])
            if revcomp:
                read.revcomp()
            print(read, file=f)


def test_gather_stats_from_bam():
    # Make a toy genome, two amplicon schemes, and a few reads. Map to make
    # the unsorted BAM, then we can test gather_stats_from_bam().
    random.seed(42)
    ref_seq = "".join(random.choices(["A", "C", "G", "T"], k=1100))
    ref_fasta = "tmp.gather_stats_from_bam.ref.fa"
    with open(ref_fasta, "w") as f:
        print(">ref", file=f)
        print(ref_seq, file=f)

    tsv_files = {
        "scheme1": os.path.join(data_dir, "gather_stats_from_bam.amplicons_1.tsv"),
        "scheme2": os.path.join(data_dir, "gather_stats_from_bam.amplicons_2.tsv"),
    }
    # scheme1:
    # 100-300, 290-800, 790-1000
    #
    # scheme2:
    # 100-300, 290-500, 490-700, 790-1001
    unpaired_read_coords = [(100, 290), (110, 290), (300, 800), (750, 900)]
    unpaired_reads_fa = "tmp.gather_stats_from_bam.reads.fa"
    _write_sim_reads(ref_seq, unpaired_read_coords, unpaired_reads_fa)
    unpaired_bam = "tmp.gather_stats_from_bam.unpaired.bam"
    Minimap(unpaired_bam, ref_fasta, unpaired_reads_fa, sort=False).run()

    reads1_coords = [(100, 200), (110, 210), (310, 410)]
    reads2_coords = [(200, 300), (200, 300), (900, 1000)]
    reads1_fa = "tmp.gather_stats_from_bam.reads_1.fa"
    reads2_fa = "tmp.gather_stats_from_bam.reads_2.fa"
    _write_sim_reads(ref_seq, reads1_coords, reads1_fa, suffix="/1")
    _write_sim_reads(ref_seq, reads2_coords, reads2_fa, revcomp=True, suffix="/2")
    paired_bam = "tmp.gather_stats_from_bam.paired.bam"
    Minimap(paired_bam, ref_fasta, reads1_fa, fq2=reads2_fa, sort=False).run()

    amplicon_sets = [
        primers.AmpliconSet.from_tsv(v, name=k) for k, v in tsv_files.items()
    ]

    rs = Bam(unpaired_bam, template_length_threshold=20)
    try:
        rs.detect_amplicon_set(amplicon_sets)
    except Exception as error:
        # this test case should fail to match any of the sets
        if str(error) != "failed to choose amplicon scheme":
            raise Error
    got = rs.stats
    # got = detect_primers.gather_stats_from_bam(unpaired_bam, tmp_bam_out, amplicon_sets)
    for k, v in got.items():
        print(k, v)
    assert got == {
        "total_reads": 4,
        "reads1": 0,
        "reads2": 0,
        "unpaired_reads": 4,
        "mapped": 4,
        # "match_any_amplicon": 3,
        "match_no_amplicon_sets": 1,
        "read_lengths": {190: 1, 180: 1, 500: 1, 150: 1},
        "template_lengths": {190: 1, 180: 1, 500: 1, 150: 1},  # TODO: check this
        "templates_that_were_too_short": {},
        # "amplicon_scheme_set_matches": {("scheme1",): 1, ("scheme1", "scheme2"): 2},
        "amplicon_scheme_set_matches": {"scheme1": 3, "scheme2": 2},
        # "amplicon_scheme_simple_counts": {"scheme1": 3, "scheme2": 2},
        "chosen_amplicon_scheme": "scheme1",
    }

    rs = Bam(paired_bam)
    try:
        rs.detect_amplicon_set(amplicon_sets)
    except Exception as error:
        if str(error) != "failed to choose amplicon scheme":
            raise Error
    got = rs.stats
    # got = detect_primers.gather_stats_from_bam(paired_bam, tmp_bam_out, amplicon_sets)
    assert got == {
        "total_reads": 6,
        "reads1": 3,
        "reads2": 3,
        "unpaired_reads": 0,
        "mapped": 6,
        # "match_any_amplicon": 2,
        "match_no_amplicon_sets": 1,  # confirm
        "read_lengths": {100: 6},
        "template_lengths": {200: 1, 190: 1, 690: 1},  # TODO: check this
        "templates_that_were_too_short": {},
        # "amplicon_scheme_set_matches": {("scheme1", "scheme2"): 2},
        "amplicon_scheme_set_matches": {"scheme1": 2, "scheme2": 2},
        # "amplicon_scheme_simple_counts": {"scheme1": 2, "scheme2": 2},
        "chosen_amplicon_scheme": "scheme1",
    }

    os.unlink(ref_fasta)
    os.unlink(reads1_fa)
    os.unlink(reads2_fa)
    os.unlink(unpaired_reads_fa)
    os.unlink(unpaired_bam)
    os.unlink(paired_bam)
