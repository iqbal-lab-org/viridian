import os
from pathlib import Path
import pytest
from collections import namedtuple, defaultdict

from intervaltree import Interval
from viridian_workflow import self_qc, primers

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "self_qc")


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
    cigar = [(3, 0)]

    alignment = Alignment(0, cigar, 0)
    assert self_qc.parse_cigar(query, alignment) == [(0, "A"), (1, "A"), (2, "A")]

    ref = "AAA"
    query = "ATTAA"
    alignment = Alignment(0, [(1, 0), (2, 1), (2, 0)], 0)
    assert self_qc.parse_cigar(query, alignment) == [
        (0, "A"),
        #        (1, "TT"),
        (1, "A"),
        (2, "A"),
    ]

    ref = "ATTAA"
    query = "AAA"
    alignment = Alignment(0, [(1, 0), (2, 2), (2, 0)], 0)
    assert self_qc.parse_cigar(query, alignment) == [
        (0, "A"),
        (1, "-"),
        (2, "-"),
        (3, "A"),
        (4, "A"),
    ]

    ref = "AAAA"
    query = "GGGAAAA"
    alignment = Alignment(0, [(3, 4), (4, 0)], 3)
    assert self_qc.parse_cigar(query, alignment) == [
        (0, "A"),
        (1, "A"),
        (2, "A"),
        (3, "A"),
    ]


def test_mappy_cigar_liftover():
    # TODO: assert that we parse key positions in the alignment properly
    amplicon = primers.Amplicon("test_amplicon")
    seq = "CTTCAGGTGATGGCACAACAAGTCCTATTTGAACATAGACTCACGAGATTGCGGTTATACTTTCGAAAATGGGAATCTGGAGTAAAAGACTAAAGTTAGATACACAGTTGCTTCACTTCAGACTATTACCAGCTGTACTCAACTCAATTGAGTACAGACACTGGTGTTGAACATGTGCCATCTTCTTCATCTACAATAAAATTGTTGATGAGCCTGAAGAACATGGTCCAATTCACACAACGACGGTTCATCCGGAGTTGTTAATCCAGTAATGGAACCAATTTATGATGAACCGACGACGACTACTAGCGTGCCTTTGTGTTACTCAAGCTGATGAGTACGAACTTATGTACTCATTCGTTTCGGGAAGAGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCGT"
    ref = "CTTCAGGTGATGGCACAACAAGTCCTATTTCTGAACATGACTACCAGATTGGTGGTTATACTGAAAAATGGGAATCTGGAGTAAAAGACTGTGTTGTATTACACAGTTACTTCACTTCAGACTATTACCAGCTGTACTCAACTCAATTGAGTACAGACACTGGTGTTGAACATGTTACCTTCTTCATCTACAATAAAATTGTTGATGAGCCTGAAGAACATGTCCAAATTCACACAATCGACGGTTCATCCGGAGTTGTTAATCCAGTAATGGAACCAATTTATGATGAACCGACGACGACTACTAGCGTGCCTTTGTAAGCACAAGCTGATGAGTACGAACTTATGTACTCATTCGTTTCGGAAGAGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCGT"

    Alignment = namedtuple("Alignment", ["r_st", "cigar", "q_st"])

    # these cigar ops are written in pysam order!
    pysam_cigar = [
        #        (4, 32),
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
        #        (4, 70),
    ]
    cigar = []
    for op, count in pysam_cigar:
        cigar.append((count, op))
    alignment = Alignment(0, cigar, 0)
    a = self_qc.parse_cigar(seq, alignment)

    assert len(ref) - 1 == a[-1][0]
    assert ref[-54:] == seq[-54:]
    assert ref[-54:] == "".join(map(lambda x: x[1], a))[-54:]


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


def test_pileup_msa_mapping():
    ref = "ACTGACTATCGATCGATCGATCAG"
    msa = Path(data_dir) / "ref_first.msa"

    # ensure parse_msa doesn't throw error
    _, _, _ = self_qc.parse_msa(msa)

    # in this example the consensus sequence is in the first line
    bad_msa = Path(data_dir) / "cons_first.msa"

    with open(msa) as fd:
        line1 = fd.readline().strip()
        assert line1.replace("-", "") == ref

    with open(bad_msa) as fd:
        _ = fd.readline().strip()
        line2 = fd.readline().strip()
        assert line2.replace("-", "") == ref


def test_ref_cons_position_translation():
    ref = "ACTGACTATCGATCGATCGATCAG"
    """
                  1          2
    1     7  8    3          4
    ACTGACT--ATCGATCGATCGATCAG
    ---GACTGCAGC--TCGCACG-----
       1  4  7    1     1
                  0     6
    """
    msa = Path(data_dir) / "ref_first.msa"
    ref_to_consensus, consensus_to_ref, _ = self_qc.parse_msa(msa)

    assert ref_to_consensus[1] == 0  # None
    assert ref_to_consensus[7] == 4
    assert consensus_to_ref[4] == 7
    assert ref_to_consensus[13] == 10
    assert ref_to_consensus[8] == 7
    assert consensus_to_ref[7] == 8
    assert consensus_to_ref[10] == 13
    assert ref_to_consensus[24] == 16
    assert consensus_to_ref[16] == 19
    assert ref_to_consensus[8] == 7


def test_position_table_cons_shorter():
    """
             1      1     2
    1   5    0      5     1
    GACTGCAGCGCCCT--TCGCACG
    ----ACT--ATCGATCGATT---
        1    4      1  1
                    1  4
    """
    ref = "GACTGCAGCGCCCTTCGCACG"
    msa = Path(data_dir) / "cons_shorter.msa"
    ref_to_consensus, consensus_to_ref, _ = self_qc.parse_msa(msa)

    assert ref_to_consensus[1] == 0  # None
    assert ref_to_consensus[3] == 0  # None

    assert ref_to_consensus[5] == 1
    assert consensus_to_ref[1] == 5

    assert ref_to_consensus[10] == 4
    assert consensus_to_ref[4] == 10
    assert ref_to_consensus[15] == 11
    assert consensus_to_ref[11] == 15


def test_position_table_ref_shorter():
    """
    Note that this situation is nonsensical
                    1  1
        1    4      1  4
    ----ACT--ATCGATCGATT---
    GACTGCAGCGCCCT--TCGCACG
    1   5    1      1     2
             0      5     1
    """

    ref = "ACTATCGATCGATT"
    # ref_alignment = "----ACT--ATCGATCGATT---"
    msa = Path(data_dir) / "ref_shorter.msa"
    ref_to_consensus, consensus_to_ref, _ = self_qc.parse_msa(msa)

    assert consensus_to_ref[10] == 4
    assert ref_to_consensus[4] == 10
    assert consensus_to_ref[15] == 11
    assert ref_to_consensus[11] == 15


def test_pileup_masking():

    ref = "ACTGACTATCGATCGATCGATCAG"

    stats = [
        self_qc.EvaluatedStats(self_qc.Stats(i, base, base))
        for (i, base) in enumerate(ref)
    ]

    failing_filters = {"easy_fail": (lambda x: x.total.refs >= 0, lambda _: "")}

    masked_all_failed, _, _ = self_qc.Pileup._mask(ref, stats, failing_filters)
    assert masked_all_failed == "".join(["N" for _ in ref])

    passing_filters = {"easy_pass": (lambda x: x.total.refs < 0, lambda _: "")}

    masked, _, _ = self_qc.Pileup._mask(ref, stats, passing_filters)
    assert masked == ref


def cigar_logic():
    # ref TTAGTACAACTACTAACATAGTTACACGGTG---TTTAAACCGTGTTTGTACTAATTATATGCCTTATTTCTTTACTTTATTGCTACAATTGTGTACTTTTACTAGAAGTACAAATTCTAGAATTAAAGCATCTATGCCGACTACTATAG
    # qry TTAGTACAACTACTAACATAGTTACACGGTGGTGTTTAAACCGTGTTTGTACTAATTATATGCCTTATTTCTTTACTTTATTGCTACAATTGTGTACTTTTACTAGAAGTACAAATTCTAGAATTAAAGCATCTATGCCGACTACTATAG
    pass
