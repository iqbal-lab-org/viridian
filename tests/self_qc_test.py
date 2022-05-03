import os
import pytest
from collections import namedtuple

from intervaltree import Interval
from viridian_workflow import self_qc, primers

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "primers")


class StatsTest:
    def __init__(self, fail):
        self.fail = fail
        self.log = []
        self.config = self_qc.default_config

    def check_for_failure(self, **kwargs):
        return self.fail


def test_cigar_tuple_construction():

    Alignment = namedtuple("Alignment", ["r_st", "cigar", "q_st"])
    ref = "AAA"
    query = "AAA"
    cigar = [(3,0)]

    alignment = Alignment(0, cigar, 0)
    assert self_qc.parse_cigar(ref, query, alignment) == [(0, "A"), (1, "A"), (2, "A")]

    ref = "AAA"
    query = "ATTAA"
    cigar = [(1, 0), (2, 1), (2, 0)]
    assert self_qc.parse_cigar(ref, query, cigar) == [
        (0, "A"),
        #        (1, "TT"),
        (1, "A"),
        (2, "A"),
    ]

    ref = "ATTAA"
    query = "AAA"
    cigar = [(1, 0), (2, 2), (2, 0)]
    assert self_qc.parse_cigar(ref, query, cigar) == [
        (0, "A"),
        (1, "-"),
        (2, "-"),
        (3, "A"),
        (4, "A"),
    ]

    ref = "ATTAA"
    query = "AAA"
    cigar = [(0, 1), (2, 2), (0, 2)]
    assert self_qc.parse_cigar(ref, query, cigar, pysam=True) == [
        (0, "A"),
        (1, "-"),
        (2, "-"),
        (3, "A"),
        (4, "A"),
    ]

    ref = "AAAA"
    query = "GGGAAAA"
    cigar = [(3, 4), (4, 0)]
    assert self_qc.parse_cigar(ref, query, cigar, q_pos=3) == [
        (0, "A"),
        (1, "A"),
        (2, "A"),
        (3, "A"),
    ]


def test_mappy_cigar_liftover():
    amplicon = primers.Amplicon("test_amplicon")
    seq = "CTTCAGGTGATGGCACAACAAGTCCTATTTGAACATAGACTCACGAGATTGCGGTTATACTTTCGAAAATGGGAATCTGGAGTAAAAGACTAAAGTTAGATACACAGTTGCTTCACTTCAGACTATTACCAGCTGTACTCAACTCAATTGAGTACAGACACTGGTGTTGAACATGTGCCATCTTCTTCATCTACAATAAAATTGTTGATGAGCCTGAAGAACATGGTCCAATTCACACAACGACGGTTCATCCGGAGTTGTTAATCCAGTAATGGAACCAATTTATGATGAACCGACGACGACTACTAGCGTGCCTTTGTGTTACTCAAGCTGATGAGTACGAACTTATGTACTCATTCGTTTCGGGAAGAGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCGT"

    cigar = [
        (4, 32),
        (0, 29),
        (2, 2),
        (0, 7),
        (1, 1),
        (0, 4),
        (1, 1),
        (0, 8),
        (2, 1),
        (0, 11),
        (1, 3),
        (0, 1),
        (2, 1),
        (0, 26),
        (1, 1),
        (0, 8),
        (2, 1),
        (0, 76),
        (1, 2),
        (0, 46),
        (1, 1),
        (0, 4),
        (2, 1),
        (0, 11),
        (2, 1),
        (0, 77),
        (1, 2),
        (0, 5),
        (2, 1),
        (0, 40),
        (1, 1),
        (0, 54),
        (4, 70),
    ]

    self_qc.parse_cigar(seq, seq, cigar)


def test_bias_test():
    return True  # TODO resolve
    assert not self_qc.test_bias(10, 100, threshold=0.3)
    assert not self_qc.test_bias(90, 100, threshold=0.3)

    assert self_qc.test_bias(40, 100, threshold=0.3)
    assert self_qc.test_bias(60, 100, threshold=0.3)


def test_stat_evaluation():
    return True  # resolve

    fwd = self_qc.BaseProfile(False, True, "test_amplicon1")
    rev = self_qc.BaseProfile(False, False, "test_amplicon2")
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

    assert stats.check_for_failure(bias_threshold=0.3)


def test_masking():

    fail = StatsTest(True)
    succeed = StatsTest(False)

    sequence = "ATCATC"
    stats = {0: succeed, 4: fail}
    masked, _ = self_qc.mask_sequence(sequence, stats)
    assert masked == "ATCANC"

    sequence = "ATCATC"
    stats = {0: fail, 4: fail}
    masked, _ = self_qc.mask_sequence(sequence, stats)
    assert masked == "NTCANC"

def cigar_logic():
    # ref TTAGTACAACTACTAACATAGTTACACGGTG---TTTAAACCGTGTTTGTACTAATTATATGCCTTATTTCTTTACTTTATTGCTACAATTGTGTACTTTTACTAGAAGTACAAATTCTAGAATTAAAGCATCTATGCCGACTACTATAG
    # qry TTAGTACAACTACTAACATAGTTACACGGTGGTGTTTAAACCGTGTTTGTACTAATTATATGCCTTATTTCTTTACTTTATTGCTACAATTGTGTACTTTTACTAGAAGTACAAATTCTAGAATTAAAGCATCTATGCCGACTACTATAG
    pass

