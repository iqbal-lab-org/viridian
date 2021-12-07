from intervaltree import Interval
import os
import pytest
import random
import subprocess
from unittest import mock

import pyfastaq

from viridian_workflow import detect_primers, minimap, primers

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "detect_primers")


def test_read_interval_paired():
    read = mock.Mock()
    read.is_paired = True
    read.is_reverse = False
    read.reference_start = 42
    read.next_reference_start = 262
    read.template_length = 250
    assert detect_primers.read_interval(read) == (42, 292)

    read.is_reverse = True
    read.reference_start = 262
    read.next_reference_start = 42
    read.template_length = -250
    assert detect_primers.read_interval(read) == (42, 292)


def test_read_interval_unpaired():
    read = mock.Mock()
    read.is_paired = False
    read.is_reverse = False
    read.reference_start = 42
    read.reference_end = 300
    assert detect_primers.read_interval(read) == (42, 300)


def test_match_read_to_amplicons():
    tsv_files = {
        "scheme1": os.path.join(data_dir, "match_reads_to_amplicons_1.tsv"),
        "scheme2": os.path.join(data_dir, "match_reads_to_amplicons_2.tsv"),
    }
    amplicon_sets = [primers.AmpliconSet(k, tsv_file=v) for k, v in tsv_files.items()]
    amplicons1 = primers.AmpliconSet.from_tsv(tsv_files["scheme1"])
    amplicons2 = primers.AmpliconSet.from_tsv(tsv_files["scheme2"])

    # scheme1:
    # 100-300, 290-800, 790-1000
    #
    # scheme2:
    # 100-300, 290-500, 490-700, 790-1001
    read = mock.Mock()
    read.is_paired = False
    read.reference_start = 100
    read.reference_end = 290
    assert detect_primers.match_read_to_amplicons(read, amplicon_sets) == {
        "scheme1": {Interval(95, 316, amplicons1["amp1"])},
        "scheme2": {Interval(95, 316, amplicons2["amp1"])},
    }

    read.reference_start = 250
    read.reference_end = 400
    assert detect_primers.match_read_to_amplicons(read, amplicon_sets) == {}

    read.reference_start = 100
    read.reference_end = 500
    assert detect_primers.match_read_to_amplicons(read, amplicon_sets) == {}

    read.reference_start = 400
    read.reference_end = 750
    assert detect_primers.match_read_to_amplicons(read, amplicon_sets) == {
        "scheme1": {Interval(285, 818, amplicons1["amp2"])},
    }


def test_pysam_open_mode():
    assert detect_primers.pysam_open_mode("foo.sam") == ""
    assert detect_primers.pysam_open_mode("foo.bam") == "b"
    with pytest.raises(Exception):
        detect_primers.pysam_open_mode("foo.bar")


def test_amplicon_set_counts_to_naive_total_counts():
    dict_in = {
        (1, 2, 3): 42,
        (1,): 11,
        (2,): 100,
    }
    assert detect_primers.amplicon_set_counts_to_naive_total_counts(dict_in) == {
        1: 53,
        2: 142,
        3: 42,
    }


def test_amplicon_set_counts_to_json_friendly():
    dict_in = {
        (3, 1, 2): 42,
        (1,): 11,
        (2,): 100,
    }
    assert detect_primers.amplicon_set_counts_to_json_friendly(dict_in) == {
        "1;2;3": 42,
        "1": 11,
        "2": 100,
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
    minimap.run(unpaired_bam, ref_fasta, unpaired_reads_fa, sort=False)

    reads1_coords = [(100, 200), (110, 210), (310, 410)]
    reads2_coords = [(200, 300), (200, 300), (900, 1000)]
    reads1_fa = "tmp.gather_stats_from_bam.reads_1.fa"
    reads2_fa = "tmp.gather_stats_from_bam.reads_2.fa"
    _write_sim_reads(ref_seq, reads1_coords, reads1_fa, suffix="/1")
    _write_sim_reads(ref_seq, reads2_coords, reads2_fa, revcomp=True, suffix="/2")
    paired_bam = "tmp.gather_stats_from_bam.paired.bam"
    minimap.run(paired_bam, ref_fasta, reads1_fa, fq2=reads2_fa, sort=False)

    amplicon_sets = [primers.AmpliconSet(k, tsv_file=v) for k, v in tsv_files.items()]
    tmp_bam_out = "tmp.bam"
    subprocess.check_output(f"rm -f {tmp_bam_out}", shell=True)
    got = detect_primers.gather_stats_from_bam(unpaired_bam, tmp_bam_out, amplicon_sets)
    assert got == {
        "total_reads": 4,
        "reads1": 0,
        "reads2": 0,
        "unpaired_reads": 4,
        "mapped": 4,
        "match_any_amplicon": 3,
        "read_lengths": {190: 1, 180: 1, 500: 1, 150: 1},
        "amplicon_scheme_set_matches": {("scheme1",): 1, ("scheme1", "scheme2"): 2},
        "amplicon_scheme_simple_counts": {"scheme1": 3, "scheme2": 2},
        "chosen_amplicon_scheme": "scheme1",
    }
    assert os.path.exists(tmp_bam_out)
    os.unlink(tmp_bam_out)

    got = detect_primers.gather_stats_from_bam(paired_bam, tmp_bam_out, amplicon_sets)
    assert got == {
        "total_reads": 6,
        "reads1": 3,
        "reads2": 3,
        "unpaired_reads": 0,
        "mapped": 6,
        "match_any_amplicon": 2,
        "read_lengths": {100: 6},
        "amplicon_scheme_set_matches": {("scheme1", "scheme2"): 2},
        "amplicon_scheme_simple_counts": {"scheme1": 2, "scheme2": 2},
        "chosen_amplicon_scheme": "scheme2",
    }
    assert os.path.exists(tmp_bam_out)
    os.unlink(tmp_bam_out)

    os.unlink(ref_fasta)
    os.unlink(reads1_fa)
    os.unlink(reads2_fa)
    os.unlink(unpaired_reads_fa)
    os.unlink(unpaired_bam)
    os.unlink(paired_bam)
