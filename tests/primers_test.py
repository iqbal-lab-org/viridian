import os
import pytest

from collections import defaultdict

from intervaltree import Interval
from viridian_workflow import primers, readstore

import pysam

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "primers")


def test_basic_Amplicon_methods():
    amp = primers.Amplicon("name", shortname=0)
    amp.add(primers.Primer("left_primer", "ACGTA", True, True, 100, 100 + 5))
    amp.add(primers.Primer("right_primer", "TCT", False, False, 300, 300 + 3))
    assert amp.start == 100
    # End position is one past the end (ie in the style of python ranges)
    assert amp.end == 303
    # primer starts at 100, and ends at 302. Length should therefore be 203
    assert len(amp) == 203
    # If we add in a new primer that's contained in the existing amplicon,
    # end positions and length should not change.
    amp.add(primers.Primer("left_primer_alt", "AAAA", True, True, 120, 124))
    assert amp.start == 100
    assert amp.end == 303
    assert len(amp) == 203
    # Add new primers outside the current coords - then start, end, length
    # should all change
    amp.add(primers.Primer("left_primer_alt2", "CCC", True, True, 90, 93))
    amp.add(primers.Primer("rightt_primer_alt", "G", False, False, 305, 306))
    assert amp.start == 90
    assert amp.end == 306
    assert len(amp) == 216


def test_AmpliconSet_from_json():
    with pytest.raises(NotImplementedError):
        primers.AmpliconSet.from_json("foo.json")


def test_AmpliconSet_from_tsv():
    tsv_file = os.path.join(data_dir, "AmpliconSet_from_tsv.tsv")
    got = primers.AmpliconSet.from_tsv(tsv_file)

    p1 = "ACGTACGTAC"
    p2 = "TCTCTTCTCAG"
    p3 = "GGGCGCGTAGTC"
    p4 = "ATGCGCGTAAGCT"

    primer1_l = primers.Primer(
        "amp1_left_primer", p1, True, True, 100, 100 + len(p1) - 1
    )
    primer1_r = primers.Primer(
        "amp1_right_primer", p2, False, False, 300, 300 + len(p2) - 1
    )
    primer2_l = primers.Primer(
        "amp2_left_primer", p3, True, True, 290, 290 + len(p3) - 1
    )
    primer2_r = primers.Primer(
        "amp2_right_primer", p4, False, False, 500, 500 + len(p4) - 1
    )
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

    p1 = "ACGTACGTAC"
    p2 = "TCTCTTCTCAG"
    p3 = "GGGCGCGTAGTC"
    p4 = "ATGCGCGTAAGCT"
    p2a = "TGCGCGTAAGCTA"
    primer1_l = primers.Primer(
        "amp1_left_primer", p1, True, True, 100, 100 + len(p1) - 1
    )
    primer1_r = primers.Primer(
        "amp1_right_primer", p2, False, False, 300, 300 + len(p2) - 1
    )
    primer2_l = primers.Primer(
        "amp2_left_primer", p3, True, True, 290, 290 + len(p3) - 1
    )
    primer2_r = primers.Primer(
        "amp2_right_primer", p4, False, False, 500, 500 + len(p4) - 1
    )

    primer2_r_alt = primers.Primer(
        "amp2_right_primer_alt", p2a, False, False, 501, 501 + len(p2a) - 1
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
    for a in expect:
        print(expect[a], got.amplicons[a])
    assert got.amplicons == expect


def test_AmpliconSet_match():
    tsv_file = os.path.join(data_dir, "AmpliconSet_match.amplicons.tsv")
    amplicons = primers.AmpliconSet.from_tsv(tsv_file)
    amplicon_set = primers.AmpliconSet.from_tsv(tsv_file, name="NAME", tolerance=5)
    matchfn = amplicon_set.match

    for interval, truth in [
        ((0, 0), None),
        ((0, 10000), None),
        ((0, 100), None),
        ((90, 100), None),
        ((90, 150), None),
        ((94, 150), None),
        ((95, 150), [amplicons.amplicons["amp1"]]),
        ((96, 150), [amplicons.amplicons["amp1"]]),
        ((96, 314), [amplicons.amplicons["amp1"]]),
        ((96, 315), None),
        #        ((96, 315), [amplicons.amplicons["amp1"]]),
        ((96, 316), None),
        ((110, 120), [amplicons.amplicons["amp1"]]),
        ((150, 350), None),
        ((300, 400), [amplicons.amplicons["amp2"]]),
    ]:
        fragment = readstore.Fragment([])
        start, end = interval
        fragment.ref_start = start
        fragment.ref_end = end
        matches = matchfn(fragment)
        if matches is not None:
            assert matches.name == truth[0].name
        else:
            assert matches == truth  # None matches


def test_fragment_syncronisation_position_sorted():
    # this case should raise an exception now
    # actually has 38 read pairs and 1 unmated
    bam_file = data_dir + "/truncated_position_sorted_40_reads.bam"

    bam = readstore.Bam(bam_file)

    r2_for_qn = {}
    qns_for_r1s = defaultdict(set)

    for read in pysam.AlignmentFile(bam_file, "rb"):
        if read.is_read2:
            r2_for_qn[read.query_name] = read.query_sequence
        elif read.is_read1:
            qns_for_r1s[read.query_sequence].add(read.query_name)

    count = 0
    try:
        for fragment in bam.syncronise_fragments():
            count += 1
            r1, r2 = fragment.reads
            names = qns_for_r1s[r1.seq]
            assert r2.seq in map(lambda x: r2_for_qn[x], names)
    except Exception as e:
        assert str(e) == "Bam file is not sorted by name"


def test_fragment_syncronisation():
    # actually has 38 read pairs and 1 unmated
    bam_file = data_dir + "/truncated_name_sorted_40_reads.bam"

    bam = readstore.Bam(bam_file)

    r2_for_qn = {}
    qns_for_r1s = defaultdict(set)

    for read in pysam.AlignmentFile(bam_file, "rb"):
        if read.is_read2:
            r2_for_qn[read.query_name] = read.query_sequence
        elif read.is_read1:
            qns_for_r1s[read.query_sequence].add(read.query_name)

    count = 0
    for fragment in bam.syncronise_fragments():
        count += 1
        r1, r2 = fragment.reads
        names = qns_for_r1s[r1.seq]
        assert r2.seq in map(lambda x: r2_for_qn[x], names)
    assert count == 18
