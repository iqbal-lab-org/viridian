import os
import pytest

from intervaltree import Interval
from viridian_workflow import primers

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "primers")


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
    amp1 = primers.Amplicon("amp1")
    amp1.add(primer1_l)
    amp1.add(primer1_r)
    amp2 = primers.Amplicon("amp2")
    amp2.add(primer2_l)
    amp2.add(primer2_r)
    expect = {
        "amp1": amp1,
        "amp2": amp2,
    }
    assert got == expect


def test_AmpliconSet_from_tsv_viridian_workflow_format():
    tsv_file = os.path.join(
        data_dir, "AmpliconSet_from_tsv_viridian_workflow_format.tsv"
    )
    got = primers.AmpliconSet.from_tsv_viridian_workflow_format(tsv_file)
    primer1_l = primers.Primer("amp1_left_primer", "ACGTACGTAC", True, True, 100)
    primer1_r = primers.Primer("amp1_right_primer", "TCTCTTCTCAG", False, False, 300)
    primer2_l = primers.Primer("amp2_left_primer", "GGGCGCGTAGTC", True, True, 290)
    primer2_r = primers.Primer("amp2_right_primer", "ATGCGCGTAAGCT", False, False, 500)
    primer2_r_alt = primers.Primer(
        "amp2_right_primer_alt", "TGCGCGTAAGCTA", False, False, 501
    )
    amp1 = primers.Amplicon("amp1")
    amp1.add(primer1_l)
    amp1.add(primer1_r)
    amp2 = primers.Amplicon("amp2")
    amp2.add(primer2_l)
    amp2.add(primer2_r)
    amp2.add(primer2_r_alt)
    expect = {
        "amp1": amp1,
        "amp2": amp2,
    }
    assert got == expect


def test_AmpliconSet_match():
    tsv_file = os.path.join(data_dir, "AmpliconSet_match.amplicons.tsv")
    amplicons = primers.AmpliconSet.from_tsv(tsv_file)
    amplicon_set = primers.AmpliconSet("NAME", tolerance=5, tsv_file=tsv_file)
    f = amplicon_set.match
    assert f(0, 0) is None
    assert f(0, 10000) is None
    assert f(0, 100) is None
    assert f(90, 100) is None
    assert f(90, 150) is None
    assert f(94, 150) is None
    assert f(95, 150) == {Interval(95, 316, amplicons["amp1"])}
    assert f(96, 150) == {Interval(95, 316, amplicons["amp1"])}
    assert f(96, 315) == {Interval(95, 316, amplicons["amp1"])}
    assert f(96, 316) is None
    assert f(110, 120) == {Interval(95, 316, amplicons["amp1"])}
    assert f(150, 350) is None
    assert f(300, 400) == {Interval(285, 518, amplicons["amp2"])}
