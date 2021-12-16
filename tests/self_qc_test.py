import os
import pytest

from intervaltree import Interval
from viridian_workflow import self_qc, primers

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "primers")


class StatsTest:
    def __init__(self, fail):
        self.fail = fail
        self.log = []

    def check_for_failure(self):
        return fail


def test_cigar_tuple_construction():
    ref = "AAA"
    query = "AAA"
    cigar = [
        (0, 3),
    ]
    assert self_qc.cigar_to_alts(ref, query, cigar) == [(0, "A"), (1, "A"), (2, "A")]

    ref = "AAA"
    query = "ATTAA"
    cigar = [(0, 1), (1, 2), (0, 2)]
    assert self_qc.cigar_to_alts(ref, query, cigar) == [
        (0, "A"),
        (1, "TT"),
        (2, "A"),
        (3, "A"),
    ]

    ref = "ATTAA"
    query = "AAA"
    cigar = [(0, 1), (2, 2), (0, 2)]
    assert self_qc.cigar_to_alts(ref, query, cigar) == [(0, "A"), (4, "A"), (5, "A")]


def test_mappy_cigar_liftover():
    amplicon = primers.Amplicon("test_amplicon")


def test_bias_test():
    assert (self_qc.test_bias(10, 100, threshold=0.3), False)
    assert (self_qc.test_bias(90, 100, threshold=0.3), False)

    assert (self_qc.test_bias(40, 100, threshold=0.3), True)
    assert (self_qc.test_bias(60, 100, threshold=0.3), True)


def test_stat_evaluation():

    test_amplicon1 = primers.Amplicon("test_amplicon1")
    test_amplicon2 = primers.Amplicon("test_amplicon2")

    fwd = self_qc.BaseProfile(False, True, test_amplicon1)
    rev = self_qc.BaseProfile(False, False, test_amplicon2)
    # 20% alt alleles
    pileup20 = ["A", "A", "C", "T", "A", "A", "A", "A", "A", "A"]
    # 0% alt alleles
    pileup0 = ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A"]
    # 100% alt alleles
    pileup100 = ["T", "T", "T", "G", "G", "G", "T", "G", "C", "C"]

    stats = self_qc.Stats()

    for base in pileup20:
        if base != "A":
            stats.add_alt(fwd)
            stats.add_alt(rev)
        else:
            stats.add_ref(fwd)
            stats.add_ref(rev)

    assert stats.check_for_failure()


def test_masking():

    fail = StatsTest(True)
    succeed = StatsTest(False)

    sequence = "ATCATC"
    stats = {0: succeed, 4: fail}
    masked, _ = self_qc.mask_sequence(sequence, stats) == "NTCANC"
    assert masked == "ATCANC"

    sequence = "ATCATC"
    stats = {0: fail, 4: fail}
    masked, _ = self_qc.mask_sequence(sequence, stats) == "NTCANC"
    assert masked == "NTCANC"
