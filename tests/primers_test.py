import os
import pytest

from collections import defaultdict

from intervaltree import Interval
from viridian_workflow import primers

import pysam

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "primers")


def test_basic_Amplicon_methods():
    amp = primers.Amplicon("name", shortname=0)
    amp.add(primers.Primer("left_primer", "ACGTA", True, True, 100))
    amp.add(primers.Primer("right_primer", "TCT", False, False, 300))
    assert amp.start == 100
    # End position is one past the end (ie in the style of python ranges)
    assert amp.end == 303
    # primer starts at 100, and ends at 302. Length should therefore be 203
    assert len(amp) == 203
    # If we add in a new primer that's contained in the existing amplicon,
    # end positions and length should not change.
    amp.add(primers.Primer("left_primer_alt", "AAAA", True, True, 120))
    assert amp.start == 100
    assert amp.end == 303
    assert len(amp) == 203
    # Add new primers outside the current coords - then start, end, length
    # should all change
    amp.add(primers.Primer("left_primer_alt2", "CCC", True, True, 90))
    amp.add(primers.Primer("rightt_primer_alt", "G", False, False, 305))
    assert amp.start == 90
    assert amp.end == 306
    assert len(amp) == 216


def test_AmpliconSet_from_json():
    with pytest.raises(NotImplementedError):
        primers.AmpliconSet.from_json("foo.json")


def test_AmpliconSet_from_tsv():
    tsv_file = os.path.join(data_dir, "AmpliconSet_from_tsv.tsv")
    got = primers.AmpliconSet.from_tsv(tsv_file)
    primer1_l = primers.Primer("amp1_left_primer", "ACGTACGTAC", True, True, 100)
    primer1_r = primers.Primer("amp1_right_primer", "TCTCTTCTCAG", False, False, 300)
    primer2_l = primers.Primer("amp2_left_primer", "GGGCGCGTAGTC", True, True, 290)
    primer2_r = primers.Primer("amp2_right_primer", "ATGCGCGTAAGCT", False, False, 500)
    amp1 = primers.Amplicon("amp1", shortname=0)
    amp1.add(primer1_l)
    amp1.add(primer1_r)
    amp2 = primers.Amplicon("amp2", shortname=1)
    amp2.add(primer2_l)
    amp2.add(primer2_r)
    expect = {
        "amp1": amp1,
        "amp2": amp2,
    }
    assert got.amplicons == expect


def test_AmpliconSet_from_tsv_viridian_workflow_format():
    tsv_file = os.path.join(
        data_dir, "AmpliconSet_from_tsv_viridian_workflow_format.tsv"
    )
    got = primers.AmpliconSet.from_tsv(tsv_file)
    primer1_l = primers.Primer("amp1_left_primer", "ACGTACGTAC", True, True, 100)
    primer1_r = primers.Primer("amp1_right_primer", "TCTCTTCTCAG", False, False, 300)
    primer2_l = primers.Primer("amp2_left_primer", "GGGCGCGTAGTC", True, True, 290)
    primer2_r = primers.Primer("amp2_right_primer", "ATGCGCGTAAGCT", False, False, 500)
    primer2_r_alt = primers.Primer(
        "amp2_right_primer_alt", "TGCGCGTAAGCTA", False, False, 501
    )
    amp1 = primers.Amplicon("amp1", shortname=0)
    amp1.add(primer1_l)
    amp1.add(primer1_r)
    amp2 = primers.Amplicon("amp2", shortname=1)
    amp2.add(primer2_l)
    amp2.add(primer2_r)
    amp2.add(primer2_r_alt)
    expect = {
        "amp1": amp1,
        "amp2": amp2,
    }
    assert got.amplicons == expect


def test_AmpliconSet_match():
    tsv_file = os.path.join(data_dir, "AmpliconSet_match.amplicons.tsv")
    amplicons = primers.AmpliconSet.from_tsv(tsv_file)
    amplicon_set = primers.AmpliconSet.from_tsv(tsv_file, name="NAME", tolerance=5)
    f = amplicon_set.match
    assert f(0, 0) is None
    assert f(0, 10000) is None
    assert f(0, 100) is None
    assert f(90, 100) is None
    assert f(90, 150) is None
    assert f(94, 150) is None
    print("f(95, 150)", f(95, 150))
    print("amp1:", amplicons.amplicons["amp1"])
    assert f(95, 150) == [amplicons.amplicons["amp1"]]
    assert f(96, 150) == [amplicons.amplicons["amp1"]]
    assert f(96, 315) == [amplicons.amplicons["amp1"]]
    assert f(96, 316) is None
    assert f(110, 120) == [amplicons.amplicons["amp1"]]
    assert f(150, 350) is None
    assert f(300, 400) == [amplicons.amplicons["amp2"]]


def test_fragment_syncronisation():
    stats = {
        "unpaired_reads": 0,
        "reads1": 0,
        "reads2": 0,
        "total_reads": 0,
        "mapped": 0,
        "read_lengths": defaultdict(int),
        "template_lengths": defaultdict(int),
    }

    name_sorted_paired_reads = pysam.AlignmentFile(
        os.path.join(data_dir, "truncated_name_sorted_40_reads.bam")
    )
    for r1, r2 in detect_primers.syncronise_fragments(name_sorted_paired_reads, stats):
        assert r1.qname == r2.qname
        assert r1.is_reverse != r2.is_reverse
