from __future__ import annotations

import sys

from collections import defaultdict
from dataclasses import dataclass
from typing import NewType, Callable, Optional, Any
from pathlib import Path

from viridian_workflow.utils import Index0, Index1
from viridian_workflow.primers import Amplicon


@dataclass
class Config:
    # TODO move this somewhere cleaner
    min_frs: float
    min_depth: int


default_config = Config(0.7, 10)


@dataclass
class BaseProfile:
    """Per-base info extracted while reads are traversed
    """

    base: str
    in_primer: bool
    forward_strand: bool
    amplicon: Amplicon


@dataclass
class BaseCounts:
    """Per-amplicon base counts, populated while the pileup is built
    """

    refs: tuple[int, int]  # matching base calls in forward/reverse reads
    alts: tuple[int, int]
    in_primer: bool


@dataclass
class EvaluatedStats:
    """Per-position base counts, evaluated after the pileup is built
    """

    total: int
    total_reads: int
    refs: int
    alts: int
    calls_by_amplicon: dict[Amplicon, tuple[int, int]]
    position_failed: bool
    failures: dict[str, str]


Filter = Callable[[EvaluatedStats], bool]
FilterMsg = Callable[[EvaluatedStats], str]


class Stats:
    def __init__(
        self,
        aux_reference_pos: Index0,  # may be non-consensus sequence index
        base: str,
        config: Config = default_config,
    ):

        self.basecounts: dict[Amplicon, BaseCounts] = {}

        self.log: list[str] = []
        self.failures: list[str] = []

        self.aux_reference_pos = aux_reference_pos
        self.base = base
        self.config = config

    def update(self, profile: BaseProfile):
        """Accumulate per-position base calling stats
        """
        if profile.amplicon not in self.basecounts:
            self.basecounts[profile.amplicon] = BaseCounts(
                (0, 0), (0, 0), profile.in_primer
            )

        if profile.base != self.base:
            # alt calls
            f_calls, r_calls = self.basecounts[profile.amplicon].alts
            if profile.forward_strand:
                self.basecounts[profile.amplicon].alts = (f_calls + 1, r_calls)
            else:
                self.basecounts[profile.amplicon].alts = (f_calls, r_calls + 1)

        else:
            # ref calls
            f_calls, r_calls = self.basecounts[profile.amplicon].refs
            if profile.forward_strand:
                self.basecounts[profile.amplicon].refs = (f_calls + 1, r_calls)
            else:
                self.basecounts[profile.amplicon].refs = (f_calls, r_calls + 1)

    def total(self) -> int:
        total = 0
        for bc in self.basecounts.values():
            total += bc.alts[0] + bc.alts[1] + bc.refs[0] + bc.refs[1]
        return total

    def evaluate(self, filters: dict[str, tuple[Filter, FilterMsg]]) -> EvaluatedStats:
        """Evaluate the accumulated positon stats
        """
        total = 0
        total_reads = self.total()
        refs = 0
        alts = 0
        calls_by_amplicon: dict[Amplicon, tuple[int, int]] = {}

        # if more than one amplicon covers this position
        # we can decide to discount primer-base calls
        consider_primers = len(self.basecounts.keys()) > 1

        for amplicon, bc in self.basecounts.items():

            if bc.in_primer and not consider_primers:
                continue

            refs += bc.refs[0] + bc.refs[1]
            alts += bc.alts[0] + bc.alts[1]
            total += refs + alts
            calls_by_amplicon[amplicon] = (
                bc.refs[0] + bc.refs[1],
                bc.alts[0] + bc.alts[1],
            )

        e = EvaluatedStats(total, total_reads, refs, alts, calls_by_amplicon, False, {})

        for filter_name, (filter_func, msg_format) in filters.items():
            # if the filter fails
            if filter_func(e):
                e.failures[filter_name] = msg_format(e)
                e.position_failed = True

        return e

    def info(self) -> str:
        """Output position stats as VCF INFO field
        """
        # if self.position_failed is None:
        #    raise Exception("pileup must be evaluated for failure first")
        # amplicon_totals = ",".join(
        #    [
        #        f"{str(self.refs_in_amplicons[amplicon])}/{str(self.alts_in_amplicons[amplicon])}"
        #        for amplicon in self.amplicon_totals
        #    ]
        # )
        # depth_symbol = ">=" if self.total >= 1000 else "="
        # info = f"primer={self.refs_in_primer}/{self.alts_in_primer};total{depth_symbol}{self.total};amplicon_overlap={len(self.amplicon_totals)};amplicon_totals={amplicon_totals}"
        return "disabled"

    def tsv_row(self) -> str:
        # this will be None if the consensus is an 'N'. may not want to skip
        # if self.position_failed is None:
        #    raise Exception("pileup must be evaluated")

        # amplicon_totals = ";".join(
        #    [
        #        f"{str(self.refs_in_amplicons[amplicon])}:{self.refs_in_forward_strands[amplicon]}:{str(self.alts_in_amplicons[amplicon])}"
        #        for amplicon in self.amplicon_totals
        #    ]
        # )
        # alts = ";".join([f"{k}:{v}" for k, v in self.alt_bases.items()])
        row = "\t".join(
            map(
                str,
                [
                    self.aux_reference_pos,
                    self.base,
                    #            alts,
                    #            self.refs_in_primer,
                    #            self.alts_in_primer,
                    self.total,
                    #           len(self.amplicon_totals),
                    #           amplicon_totals,
                ],
            )
        )

        return row

    def __str__(self) -> str:
        #        alts = " ".join([f"{alt}:{count}" for alt, count in self.alt_bases.items()])
        #        return f"{self.base}: {alts}"
        return self.base


# filter closures take an EvaluatedStats struct and return True on failure
# def test_amplicon_bias(s: EvaluatedStats) -> bool:
#    passing = []
#
#    for amplicon, (ref, alt) in s.calls_by_amplicon.items():
#        total = ref + alt
#        if total < config.min_depth:
#            continue
#        if s.amplicon_totals[amplicon] == 0:
#            continue
#        if (
#            s.calls_by_amplicon[amplicon][0] / total
#            < self.config.min_frs
#        ):
#            passing.append(False)
#        else:
#            passing.append(True)
#    if len(passing) > 1 and not all(passing[0] == e for e in passing):
#        # this is a failure
#        return True
#    return False

"""
def test_strand_bias(s: EvaluatedStats) -> bool:
    passing = []
    for amplicon, total in s.amplicon_totals.items():
        if total < self.config.min_depth:
            continue
        if s.amplicon_totals[amplicon] == 0:
            continue
        if (
            s.refs_in_forward_strands[amplicon] / s.amplicon_totals[amplicon]
            < self.config.min_frs
        ):
            passing.append(False)
        else:
            passing.append(True)
    if len(passing) > 1 and not all(passing[0] == e for e in passing):
        # True = failure
        return True
    return False
"""


class Pileup:
    """A pileup is an array of Stats objects indexed by position in a sequence
    """

    def __init__(
        self,
        consensus_seq: str,
        msa: Optional[Path] = None,
        config: Config = default_config,
    ):
        self.config: Config = config
        self.consensus_seq: str = consensus_seq
        self.seq: list[Stats] = []
        self.evaluated_sequence: list[Optional[EvaluatedStats]] = [
            None for _ in range(len(consensus_seq))
        ]

        # 1-based index translation tables
        self._ref_to_consensus: dict[Index1, Index1] = {}
        self._consensus_to_ref: dict[Index1, Index1] = {}

        if msa is not None:
            with open(msa) as msa_fd:
                con_seq = msa_fd.readline().strip()  # consensus
                ref_seq = msa_fd.readline().strip()  # reference

                if len(msa_fd.readline()) > 0:
                    raise Exception("Invalid multiple sequence alignment file")

                # this is valid for testing when the consensus is complete:
                # assert con_seq.replace("-", "") == self.consensus_seq

                ref_pos = 0  # we're using 1-based coords (vcf)
                con_pos = 0

                for con_base, ref_base in zip(con_seq, ref_seq):
                    if con_base != "-" and ref_base != "-":
                        ref_pos += 1
                        con_pos += 1
                        self._ref_to_consensus[Index1(ref_pos)] = Index1(con_pos)
                        self._consensus_to_ref[Index1(con_pos)] = Index1(ref_pos)
                    elif con_base != "-" and ref_base == "-":
                        ref_pos += 1
                        self._ref_to_consensus[Index1(ref_pos)] = Index1(con_pos)
                    elif ref_base != "-" and con_base == "-":
                        con_pos += 1
                        self._consensus_to_ref[Index1(con_pos)] = Index1(ref_pos)

        else:
            # if an msa is not provided use the consensus coordinates
            for i, _ in enumerate(self.consensus_seq):
                pos = Index1(i + 1)
                self._ref_to_consensus[pos] = pos
                self._consensus_to_ref[pos] = pos

        for i, r in enumerate(self.consensus_seq):
            p = Index1(i + 1)
            self.seq.append(Stats(Index0(self.consensus_to_ref(p) - 1), r,))

        self.filters = {
            "low_depth": (
                lambda s: s.total < self.config.min_depth,
                lambda s: f"Insufficient depth; {s.total} < {self.config.min_depth}. {s.total_reads} including primer regions.",
            ),
            "low_frs": (
                lambda s: s.refs / s.total < self.config.min_frs
                if s.total > 0
                else False,
                lambda s: f"Insufficient support of consensus base; {s.refs} / {s.total} < {self.config.min_frs}. {s.total_reads} including primer regions.",
            ),
            #            "amplicon_bias": (
            #                test_amplicon_bias,
            #                lambda s: "Per-amplicon FRS failure",
            #            ),
            #            "strand_bias": (
            #                test_strand_bias,
            #                lambda s: f"Strand-biased FRS failure"),
        }

        # initialise summary for each filter
        self.summary = defaultdict(int)
        for f in self.filters:
            self.summary[f] = 0

    def ref_to_consensus(self, p: Index1) -> Index1:
        if p in self._ref_to_consensus:
            return self._ref_to_consensus[p]
        return Index1(0)

    def consensus_to_ref(self, p: Index1) -> Index1:
        if p in self._consensus_to_ref:
            return self._consensus_to_ref[p]
        # raise Exception("Failure mapping consensus coordinates to reference sequence")
        return Index1(0)

    def __getitem__(self, pos: Index0) -> Stats:
        if pos > len(self.seq):
            raise Exception(f"position too big: {pos} {len(self.seq)}")
        return self.seq[pos]

    def __setitem__(self, pos: Index0, profile: BaseProfile) -> None:
        raise Exception("don't set Pileup positions")

    def update(self, pos: Index0, profile: BaseProfile) -> None:
        self.seq[pos].update(profile)

    def __len__(self) -> int:
        return len(self.seq)

    def mask(self) -> str:
        """Evaluate all positions and determine if they pass filters
        """
        sequence: list[str] = list(self.consensus_seq)
        self.qc: dict[str, Any] = {}

        log: list[str] = []
        failures: list[str] = []

        # filters: dict[str, tuple[Filter, FilterMsg]]
        for p, raw_stats in enumerate(self.seq):
            if p >= len(sequence):
                raise Exception(
                    f"Invalid condition: mapped position {p} greater than consensus length {len(sequence)}"
                )
            position = Index0(p)
            position_failed = False

            # if a position is already masked by an upstream process skip it
            if sequence[position] == "N":
                self.summary["already_masked"] += 1
                self.summary["total_masked"] += 1
                continue

            stats = raw_stats.evaluate(self.filters)
            self.evaluated_sequence[position] = stats

            self.summary["consensus_length"] += 1

            if stats.position_failed:
                for failure in stats.failures:
                    self.summary[failure] += 1
                self.summary["total_masked"] += 1
                sequence[position] = "N"
                self.qc[str(position)] = stats.failures
            self.qc["masking_summary"] = self.summary
        return "".join(sequence)

    def dump_tsv(self, tsv: Path) -> Path:
        fd = open(tsv, "w")
        for pos, stats in enumerate(self.seq):
            cons_pos = self.consensus_to_ref(Index1(pos + 1))
            print(f"{cons_pos}\t{pos+1}\t{stats.tsv_row()}", file=fd)
        return tsv

    def annotate_vcf(self, vcf: Path) -> tuple[list[str], Any]:
        header = []
        records = []

        # TODO: assert 'chromosome' names are the same

        for line in open(vcf):
            line = line.strip()
            if line[0] == "#":
                header.append(line)
                continue
            (
                chrom,
                pos_token,
                mut_id,
                ref,
                alt,
                qual,
                original_filters,
                info,
                fmt,
                *r,
            ) = line.split("\t")

            pos: Index1 = Index1(int(pos_token))
            cons_coord: Index1 = self.ref_to_consensus(pos)

            maybe_stats: Optional[EvaluatedStats] = self.evaluated_sequence[
                Index0(cons_coord - 1)
            ]

            stats: EvaluatedStats = maybe_stats if maybe_stats is not None else self.seq[
                Index0(cons_coord - 1)
            ].evaluate(
                self.filters
            )

            # inf_field = stats.info()
            vcf_filters = original_filters
            if stats.position_failed is None:
                print(f"Warning: attemped to evaluate N basecall", file=sys.stderr)
                continue
            if stats.position_failed:
                if original_filters == "PASS":
                    vcf_filters = ";".join(stats.failures.keys())
                else:
                    vcf_filters = ";".join([original_filters, *vcf_filters])

            info_field = ""  # TODO
            records.append(
                (chrom, pos, mut_id, ref, alt, qual, vcf_filters, info_field, fmt, *r)
            )

        return header, records


def parse_cigar(ref, query: str, alignment: Any) -> list[tuple[Index0, str]]:
    """Interpret cigar string and query sequence in reference
    coords from mappy (count, op)

    Returns a list of query basecalls per reference position

    ref: AATGG
    qry: AACT-
    (1, "A")
    (2, "AC")
    (3, "T")
    (4, "-")
    (5, "-")
    etc.
    """
    positions = []
    r_pos = alignment.r_st
    cigar = alignment.cigar
    q_pos = alignment.q_st

    for count, op in cigar:
        if op == 0:
            # match/mismatch
            for _ in range(count):
                if len(query) == q_pos:
                    # done? TODO: verify that this captures the full query sequence
                    # print(f"WARNING: invalid cigar string. {len(query)}, index {q_pos}. cigar: {cigar}", file=sys.stderr)
                    break

                positions.append((r_pos, query[q_pos]))
                q_pos += 1
                r_pos += 1

        elif op == 1:
            # insertion
            # positions.append((r_pos, query[q_pos : q_pos + count + 1]))
            q_pos += count

        elif op == 2:
            # deletion
            for _ in range(count):
                positions.append((r_pos, "-"))
                r_pos += 1

        elif op == 3:
            # ref_skip
            pass

        elif op == 4:
            # soft clip
            # q_pos += count
            # may not need to be considered if q_pos offset is set
            # TODO verify
            # may need to consider where softclipping on the 3' end
            pass

        elif op == 5:
            # hard clip
            pass

        else:
            raise Exception(f"invalid cigar op {op}")

    return positions
