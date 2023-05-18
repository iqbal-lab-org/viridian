import copy
import intervaltree
import os
import pytest
import random
from unittest import mock

import pyfastaq
import pysam

from viridian import maptools, reads, utils


def file_to_lines(filename):
    with open(filename) as f:
        return [x.rstrip() for x in f]


def make_perfect_reads(
    ref_seq,
    start,
    end,
    out_unpaired,
    out_paired1,
    out_paired2,
    number_of_reads,
    paired_length,
    left_adapter="",
    right_adapter="",
):
    assert paired_length < end - start + 1
    read1 = pyfastaq.sequences.Fasta(
        f"{ref_seq.id}.{start}.{end}",
        left_adapter + ref_seq[start : start + paired_length],
    )
    read1_rev = copy.deepcopy(read1)
    read1_rev.revcomp()
    read2 = pyfastaq.sequences.Fasta(
        f"{ref_seq.id}.{start}.{end}",
        ref_seq[end - paired_length + 1 : end + 1] + right_adapter,
    )
    read2_rev = copy.deepcopy(read2)
    read2_rev.revcomp()
    unpaired_read = pyfastaq.sequences.Fasta(
        f"{ref_seq.id}.{start}.{end}",
        left_adapter + ref_seq[start : end - 1] + right_adapter,
    )
    unpaired_read_rev = copy.deepcopy(unpaired_read)
    unpaired_read_rev.revcomp()

    for i in range(int(number_of_reads / 2)):
        print(
            f">{unpaired_read.id}.{i}", unpaired_read.seq, sep="\n", file=out_unpaired
        )
        print(
            f">{unpaired_read_rev.id}.r.{i}",
            unpaired_read_rev.seq,
            sep="\n",
            file=out_unpaired,
        )
        print(f">{read1.id}.a.{i} /1", read1.seq, sep="\n", file=out_paired1)
        print(f">{read1.id}.a.{i} /2", read2_rev.seq, sep="\n", file=out_paired2)
        print(f">{read1.id}.b.{i} /1", read2_rev.seq, sep="\n", file=out_paired1)
        print(f">{read1.id}.b.{i} /2", read1.seq, sep="\n", file=out_paired2)


def test_to_adapter_trimmed_seq_trim_start():
    left_starts = {10: 12, 11: 12, 12: 12}
    right_ends = {}
    seq = "AGCTGCGA"
    rev_seq = "TCGCAGCT"
    read = mock.Mock()
    read.get_forward_sequence = mock.MagicMock(return_value=rev_seq)
    read.query_sequence = seq
    read.template_length = 100
    read.is_paired = True

    read.is_reverse = False
    read.cigartuples = [(pysam.CMATCH, 8)]
    read.reference_start = 9
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_start = 10
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[2:]
    read.reference_start = 11
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[1:]
    read.reference_start = 12
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_start = 13
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq

    read.is_reverse = True
    read.reference_start = 9
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_start = 10
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[:-2]
    read.reference_start = 11
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[:-1]
    read.reference_start = 12
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_start = 13
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq

    read.is_reverse = False
    read.cigartuples = [(pysam.CSOFT_CLIP, 3), (pysam.CMATCH, 5)]
    read.reference_start = 9
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_start = 10
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[5:]
    read.reference_start = 11
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[4:]
    read.reference_start = 12
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[3:]
    read.reference_start = 13
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_start = 14
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq

    read.is_reverse = True
    read.cigartuples = [(pysam.CSOFT_CLIP, 3), (pysam.CMATCH, 5)]
    read.reference_start = 9
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_start = 10
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[:-5]
    read.reference_start = 11
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[:-4]
    read.reference_start = 12
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[:-3]
    read.reference_start = 13
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_start = 14
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq


def test_to_adapter_trimmed_seq_trim_end():
    left_starts = {}
    right_ends = {20: 20, 21: 20, 22: 20}
    seq = "AGCTGCGA"
    rev_seq = "TCGCAGCT"
    read = mock.Mock()
    read.get_forward_sequence = mock.MagicMock(return_value=rev_seq)
    read.query_sequence = seq
    read.is_paired = True
    read.template_length = -100

    read.is_reverse = False
    read.cigartuples = [(pysam.CMATCH, 8)]
    # reminder: reference_end is one past the end
    read.reference_end = 24
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_end = 23
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[:-2]
    read.reference_end = 22
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[:-1]
    read.reference_end = 21
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_end = 20
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_end = 19
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq

    read.is_reverse = True
    read.reference_end = 24
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_end = 23
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[2:]
    read.reference_end = 22
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[1:]
    read.reference_end = 21
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_end = 20
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_end = 19
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq

    read.is_reverse = False
    read.cigartuples = [(pysam.CMATCH, 5), (pysam.CSOFT_CLIP, 3)]
    read.reference_end = 24
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_end = 23
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[:-5]
    read.reference_end = 22
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[:-4]
    read.reference_end = 21
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq[:-3]
    read.reference_end = 20
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq
    read.reference_end = 19
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == seq

    read.is_reverse = True
    read.cigartuples = [(pysam.CMATCH, 5), (pysam.CSOFT_CLIP, 3)]
    read.reference_end = 24
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_end = 23
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[5:]
    read.reference_end = 22
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[4:]
    read.reference_end = 21
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq[3:]
    read.reference_end = 20
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq
    read.reference_end = 19
    assert reads.to_adapter_trimmed_seq(read, left_starts, right_ends) == rev_seq


def test_read_to_fragment_coords():
    read = mock.Mock()
    read.is_paired = False
    read.reference_start = 1100
    read.template_length = 400
    read.reference_end = 1200
    assert reads.read_to_fragment_coords(read) == (1100, 1199)
    read.is_paired = True
    assert reads.read_to_fragment_coords(read) == (1100, 1499)
    read.template_length = -400
    assert reads.read_to_fragment_coords(read) == (800, 1199)


def test_write_cylon_fasta_paired():
    merged_seq = "CATCGTCACATGCTTATATAAAAGCGATCTCATAGGTACTGAGTCCTATCG"
    merged_rev = pyfastaq.sequences.Fasta("x", merged_seq)
    merged_rev.revcomp()
    merged_rev = merged_rev.seq
    read1_1 = merged_seq[:-5]
    read1_2 = merged_rev[:-5]
    read1_total = len(read1_1) + len(read1_2)
    fragments = [
        reads.Fragment(0, 0, "name1", False, {"read": read1_1, "mate": read1_2}),
        reads.Fragment(0, 0, "name2", False, {"read": "ACGTA", "mate": "GGGGG"}),
        reads.Fragment(0, 0, "name3", False, {"read": "TTTTT", "mate": "TTTAT"}),
    ]
    outfile = "tmp.write_cylon_fasta_paired.fa"
    utils.syscall(f"rm -f {outfile}")
    assert reads.write_cylon_fasta_paired(fragments, outfile, 10) == read1_total
    expect = [">1", merged_seq]
    got = file_to_lines(outfile)
    assert expect == got

    assert (
        reads.write_cylon_fasta_paired(fragments, outfile, read1_total + 1)
        == read1_total + 10
    )
    expect.extend([">2", "ACGTA", ">3", "GGGGG"])
    got = file_to_lines(outfile)
    assert expect == got

    assert (
        reads.write_cylon_fasta_paired(fragments, outfile, read1_total + 10)
        == read1_total + 10
    )
    got = file_to_lines(outfile)
    assert expect == got

    assert (
        reads.write_cylon_fasta_paired(fragments, outfile, read1_total + 11)
        == read1_total + 20
    )
    expect = [
        ">1",
        merged_seq,
        ">2",
        "ACGTA",
        ">3",
        "TTTTT",
        ">4",
        "GGGGG",
        ">5",
        "TTTAT",
    ]
    got = file_to_lines(outfile)
    assert expect == got
    os.unlink(outfile)


def test_write_cylon_fasta_unpaired():
    fragments = [
        reads.Fragment(0, 0, "name1", False, {"read": "AAAAA"}),
        reads.Fragment(0, 0, "name2", False, {"read": "ACGTA"}),
        reads.Fragment(0, 0, "name3", False, {"read": "TTTTT"}),
        reads.Fragment(0, 0, "name4", False, {"read": "TATAT"}),
    ]
    outfile = "tmp.write_cylon_fasta_unpaired.fa"
    utils.syscall(f"rm -f {outfile}")
    assert reads.write_cylon_fasta_unpaired(fragments, outfile, 5) == 0

    fragments[2] = reads.Fragment(0, 0, "name3", True, {"read": "TTTTT"})
    fragments[3] = reads.Fragment(0, 0, "name4", True, {"read": "TATAT"})
    assert reads.write_cylon_fasta_unpaired(fragments, outfile, 5) == 10
    expect = [">0", "AAAAA", ">1", "TTTTT"]
    got = file_to_lines(outfile)
    assert expect == got

    assert reads.write_cylon_fasta_unpaired(fragments, outfile, 10) == 10
    got = file_to_lines(outfile)
    assert expect == got

    assert reads.write_cylon_fasta_unpaired(fragments, outfile, 11) == 20
    expect.extend([">2", "ACGTA", ">3", "TATAT"])
    got = file_to_lines(outfile)
    assert expect == got
    os.unlink(outfile)


def test_write_qc_fastas_paired():
    outprefix = "tmp.write_qc_fastas.paired"
    utils.syscall(f"rm -rf {outprefix}*")
    fragments = [
        reads.Fragment(0, 0, "name1", False, {"read": "AAA", "mate": "CCC"}),
        reads.Fragment(0, 0, "name2", False, {"read": "ACG", "mate": "GGG"}),
        reads.Fragment(0, 1, "name3", False, {"read": "TTA", "mate": "TTT"}),
    ]

    got_files = reads.write_qc_fastas(fragments, outprefix)
    expect_files = {
        (0, 0): (f"{outprefix}.0.0.1.fa", f"{outprefix}.0.0.2.fa"),
        (0, 1): (f"{outprefix}.0.1.1.fa", f"{outprefix}.0.1.2.fa"),
    }
    assert got_files == expect_files
    assert os.path.exists(expect_files[(0, 0)][0])
    assert os.path.exists(expect_files[(0, 0)][1])
    assert os.path.exists(expect_files[(0, 1)][0])
    assert os.path.exists(expect_files[(0, 1)][1])
    assert file_to_lines(expect_files[(0, 0)][0]) == [">0/1", "AAA", ">1/1", "ACG"]
    assert file_to_lines(expect_files[(0, 0)][1]) == [">0/2", "CCC", ">1/2", "GGG"]
    assert file_to_lines(expect_files[(0, 1)][0]) == [">0/1", "TTA"]
    assert file_to_lines(expect_files[(0, 1)][1]) == [">0/2", "TTT"]
    utils.syscall(f"rm -r {outprefix}*")


def test_write_qc_fastas_unpaired():
    outprefix = "tmp.write_qc_fastas.unpaired"
    utils.syscall(f"rm -rf {outprefix}*")
    fragments = [
        reads.Fragment(0, 0, "name1", False, {"read": "AAA", "mate": ""}),
        reads.Fragment(0, 0, "name2", False, {"read": "ACG", "mate": ""}),
        reads.Fragment(0, 1, "name3", True, {"read": "TTA", "mate": ""}),
    ]

    got_files = reads.write_qc_fastas(fragments, outprefix)
    expect_files = {
        (0, 0): (f"{outprefix}.0.0.fa", None),
        (0, 1): (f"{outprefix}.0.1.fa", None),
    }
    assert got_files == expect_files
    assert os.path.exists(expect_files[(0, 0)][0])
    assert os.path.exists(expect_files[(0, 1)][0])
    assert file_to_lines(expect_files[(0, 0)][0]) == [">0", "AAA", ">1", "ACG"]
    assert file_to_lines(expect_files[(0, 1)][0]) == [">0", "TTA"]
    utils.syscall(f"rm -r {outprefix}*")


def test_read_sampler():
    # One monolothic test that runs through everything from start to end.
    # There's so much setting up needed to make input data etc that not
    # worth pulling out individual functions in most cases
    outprefix = "tmp.sample_reads_test"
    utils.syscall(f"rm -rf {outprefix}.*")
    amplicons = [
        {
            "name": "amp1",
            "start": 90,
            "end": 405,
            "primers": {
                "left": [
                    {"start": 90, "end": 110, "read_count": 50},
                    {"start": 100, "end": 120, "read_count": 50},
                ],
                "right": [{"start": 390, "end": 405, "read_count": 50}],
            },
            "excluded_primers": {},
        },
        {
            "name": "amp2",
            "start": 375,
            "end": 600,
            "primers": {
                "left": [{"start": 375, "end": 385, "read_count": 50}],
                "right": [
                    {"start": 580, "end": 600, "read_count": 50},
                    {"start": 570, "end": 590, "read_count": 50},
                ],
            },
            "excluded_primers": {},
        },
        {
            "name": "amp3",
            "start": 550,
            "end": 800,
            "primers": {
                "left": [{"start": 550, "end": 570, "read_count": 50}],
                "right": [
                    {"start": 780, "end": 800, "read_count": 50},
                ],
            },
            "excluded_primers": {},
        },
    ]

    read_sampler = reads.ReadSampler(amplicons, adapter_trim_tolerance=1)
    assert read_sampler.left_primer_starts == {
        89: 90,
        90: 90,
        99: 100,
        100: 100,
        374: 375,
        375: 375,
        549: 550,
        550: 550,
    }
    assert read_sampler.right_primer_ends == {
        405: 405,
        406: 405,
        600: 600,
        601: 600,
        590: 590,
        591: 590,
        800: 800,
        801: 800,
    }

    expect_amp_tree = intervaltree.IntervalTree()
    expect_amp_tree[89:406] = (0, 0, 0)
    expect_amp_tree[99:406] = (0, 1, 0)
    expect_amp_tree[374:591] = (1, 0, 0)
    expect_amp_tree[374:601] = (1, 0, 1)
    expect_amp_tree[549:801] = (2, 0, 0)
    assert read_sampler.amp_tree == expect_amp_tree

    ref_fasta = f"{outprefix}.ref.fa"
    random.seed(42)
    ref_bases = random.choices(["A", "C", "G", "T"], k=1000)
    ref_bases[89] = "T"  # left adapter will map to this
    ref_bases[90] = "T"
    ref_bases[400] = "G"
    ref_bases[595] = "T"
    ref_bases[700] = "A"
    ref_bases[701] = "A"
    ref_seq = pyfastaq.sequences.Fasta("ref", "".join(ref_bases))
    with open(ref_fasta, "w") as f:
        print(ref_seq, file=f)

    unpaired_fa = f"{outprefix}.reads.fa"
    paired1_fa = f"{outprefix}.reads_1.fa"
    paired2_fa = f"{outprefix}.reads_2.fa"
    f_up = open(unpaired_fa, "w")
    f_p1 = open(paired1_fa, "w")
    f_p2 = open(paired2_fa, "w")
    right_adapters = {400: "C"}
    left_adapters = {90: "ACGT"}
    for amp in amplicons:
        for left_primer in amp["primers"]["left"]:
            for right_primer in amp["primers"]["right"]:
                make_perfect_reads(
                    ref_seq,
                    left_primer["start"],
                    right_primer["end"],
                    f_up,
                    f_p1,
                    f_p2,
                    100,
                    150,
                    right_adapter=right_adapters.get(right_primer["end"], ""),
                    left_adapter=left_adapters.get(left_primer["start"], ""),
                )
    f_up.close()
    f_p1.close()
    f_p2.close()

    qc_ref_bases = copy.deepcopy(ref_bases)
    qc_ref_bases[595] = "A"
    qc_ref_bases[700] = "ACGT"
    qc_ref_fa = f"{outprefix}.qc_ref.fa"
    qc_seq = pyfastaq.sequences.Fasta("qc", "".join(qc_ref_bases))
    with open(qc_ref_fa, "w") as f:
        print(qc_seq, file=f)

    unpaired_bam = f"{outprefix}.unpaired.bam"
    maptools.map_reads(unpaired_bam, ref_fasta, unpaired_fa)
    unpaired_sample_out = f"{outprefix}.unpaired.sample"
    unpaired_pileup_out = f"{outprefix}.unpaired.pileup"
    read_sampler.sample_reads(unpaired_bam, unpaired_sample_out)
    pileup = read_sampler.pileups(qc_ref_fa, unpaired_pileup_out, debug=True)

    # The pileup dict is huge. Check a few bits of it
    assert sorted(list(pileup.keys())) == [0, 1, 2]
    assert 89 not in pileup[0][(0, 0)]
    assert 90 in pileup[0][(0, 0)]
    assert pileup[0][(0, 0)][90]["T"] == 50
    assert pileup[0][(0, 0)][90]["t"] == 49
    empty_indels = {x: {} for x in ["D", "d", "I", "i"]}
    expect_700 = {"A-3CGT": 49, "a-3cgt": 50, "D": 0, "I": 0, "d": 0, "i": 0}
    expect_700["indel"] = copy.deepcopy(empty_indels)
    assert pileup[2][(0, 0)][700] == expect_700

    paired_bam = f"{outprefix}.paired.bam"
    maptools.map_reads(paired_bam, ref_fasta, paired1_fa, reads2=paired2_fa)
    paired_sample_out = f"{outprefix}.paired.sample"
    paired_pileup_out = f"{outprefix}.paired.pileup"
    read_sampler.sample_reads(paired_bam, paired_sample_out)
    pileup = read_sampler.pileups(qc_ref_fa, paired_pileup_out, debug=True)
    assert sorted(list(pileup.keys())) == [0, 1, 2]

    assert 89 not in pileup[0][(0, 0)]
    assert 90 in pileup[0][(0, 0)]
    assert pileup[0][(0, 0)][90]["T"] == 99
    assert "t" not in pileup[0][(0, 0)][90]
    expect_700 = {"a-3cgt": 99, "D": 0, "I": 0, "d": 0, "i": 0}
    expect_700["indel"] = copy.deepcopy(empty_indels)
    assert pileup[2][(0, 0)][700] == expect_700

    utils.syscall(f"rm -r {outprefix}*")
