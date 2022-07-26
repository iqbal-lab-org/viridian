from __future__ import __annotations__

import sys

from collections import defaultdict
from dataclasses import dataclass
from typing import NewType, Callable, Optional, Any
from pathlib import Path


@dataclass
class BaseProfile:
    base: str
    in_primer: bool
    forward_strand: bool
    amplicon_name: str


@dataclass
class Config:
    min_frs: float
    min_depth: int


default_config = Config(0.7, 10)

Index0 = NewType("Index0", int)
Index1 = NewType("Index1", int)

Filter = Callable[[Stats], bool]
FilterMsg = Callable[[Stats], str]


class Stats:
    def __init__(
        self, reference_pos, ref_base=None, cons_base=None, config=default_config,
    ):
        self.alts_in_primer = 0
        self.refs_in_primer = 0
        self.alt_bases = defaultdict(int)

        self.alts_in_amplicons = defaultdict(int)
        self.refs_in_amplicons = defaultdict(int)
        self.refs_in_forward_strands = defaultdict(int)
        self.alts_in_forward_strands = defaultdict(int)
        self.amplicon_totals = defaultdict(int)

        self.alts_forward = 0
        self.refs_forward = 0

        self.alts = 0
        self.refs = 0
        self.total = 0
        self.log = []
        self.total_reads = 0
        self.failures = None

        self.alts_matching_refs = 0
        self.reference_pos = reference_pos
        self.ref_base = ref_base
        self.cons_base = cons_base

        # self.config = config

        self.position_failed = None

    def update(self, profile: BaseProfile, alt=None):

        if profile.in_primer:
            if profile.base != self.ref_base:
                self.alts_in_primer += 1
            else:
                self.refs_in_primer += 1
        else:
            self.alt_bases[profile.base] += 1
            if profile.amplicon_name:
                self.amplicon_totals[profile.amplicon_name] += 1

            if profile.base != self.ref_base:
                self.alts += 1
                if profile.amplicon_name:
                    self.alts_in_amplicons[profile.amplicon_name] += 1
                    if profile.forward_strand:
                        self.alts_in_forward_strands[profile.amplicon_name] += 1

                if profile.forward_strand:
                    self.alts_forward += 1

            else:
                self.refs += 1
                if profile.amplicon_name:
                    self.refs_in_amplicons[profile.amplicon_name] += 1
                    if profile.forward_strand:
                        self.refs_in_forward_strands[profile.amplicon_name] += 1

                if profile.forward_strand:
                    self.refs_forward += 1

            self.total += 1

    def check_for_failure(self, filters: dict[str, tuple[Filter, FilterMsg]]) -> bool:
        """return whether a position should be masked

        mutates the pileup's masking decision log
        """

        self.position_failed = False

        for filter_name, (filter_func, msg_format) in filters.items():
            if filter_func(self):
                self.log.append(msg_format(self))
                if self.failures is None:
                    self.failures = []
                self.failures.append(filter_name)
                self.position_failed = True

        return self.position_failed

    def get_failures(self) -> list[str]:
        """return list of failed filters for VCF FILTER field
        """
        if self.position_failed is None:
            raise Exception("pileup must be evaluated for failure first")
        return self.failures

    def info(self) -> str:
        """Output position stats as VCF INFO field
        """
        if self.position_failed is None:
            raise Exception("pileup must be evaluated for failure first")
        amplicon_totals = ",".join(
            [
                f"{str(self.refs_in_amplicons[amplicon])}/{str(self.alts_in_amplicons[amplicon])}"
                for amplicon in self.amplicon_totals
            ]
        )
        depth_symbol = ">=" if self.total >= 1000 else "="
        info = f"primer={self.refs_in_primer}/{self.alts_in_primer};total{depth_symbol}{self.total};amplicon_overlap={len(self.amplicon_totals)};amplicon_totals={amplicon_totals}"
        return info

    def tsv_row(self) -> str:
        # this will be None if the consensus is an 'N'. may not want to skip
        # if self.position_failed is None:
        #    raise Exception("pileup must be evaluated")

        amplicon_totals = ";".join(
            [
                f"{str(self.refs_in_amplicons[amplicon])}:{self.refs_in_forward_strands[amplicon]}:{str(self.alts_in_amplicons[amplicon])}"
                for amplicon in self.amplicon_totals
            ]
        )
        alts = ";".join([f"{k}:{v}" for k, v in self.alt_bases.items()])
        row = "\t".join(
            map(
                str,
                [
                    self.reference_pos,
                    self.ref_base,
                    alts,
                    self.refs_in_primer,
                    self.alts_in_primer,
                    self.total,
                    len(self.amplicon_totals),
                    amplicon_totals,
                ],
            )
        )

        return row

    def __str__(self) -> str:
        alts = " ".join([f"{alt}:{count}" for alt, count in self.alt_bases.items()])
        return f"{self.ref_base}: {alts}"


class Pileup:
    """A pileup is an array of Stats objects indexed by position in a reference
    """

    def __init__(
        self, refseq: str, msa: Optional[Path] = None, config: Config = default_config,
    ):
        self.config: Config = config
        self.ref: str = refseq
        self.seq: list[Stats] = []

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
                # assert seq1.replace("-", "") == self.ref

                ref_pos = 0  # we're using 1-based coords (vcf)
                con_pos = 0

                for con_base, ref_base in zip(con_seq, ref_seq):
                    if con_base != "-":
                        ref_pos += 1
                    else:
                        self._ref_to_consensus[Index1(ref_pos)] = Index1(con_pos)

                    if ref_base != "-":
                        con_pos += 1
                    else:
                        self._consensus_to_ref[Index1(con_pos)] = Index1(ref_pos)
        else:
            for i in range(1, len(self.ref) + 1):
                pos = Index1(i)
                self._ref_to_consensus[pos] = pos
                self._consensus_to_ref[pos] = pos

        for i, r in enumerate(refseq):
            p = Index1(i + 1)
            self.seq.append(Stats(self.consensus_to_ref(p), ref_base=r))

        # define the filters
        # filter closures take a Stats struct and return True on failure
        def test_amplicon_bias(s: Stats) -> bool:
            passing = []
            # if len(s.amplicon_totals) > 2:
            #    for _aa in s.amplicon_totals:
            #        print(
            #            "\t\t-->",
            #            _aa.name,
            #            _aa.start,
            #            _aa.end,
            #            s.amplicon_totals[_aa],
            #            file=sys.stderr,
            #        )

            for amplicon, total in s.amplicon_totals.items():
                if total < self.config.min_depth:
                    continue
                if s.amplicon_totals[amplicon] == 0:
                    continue
                if (
                    s.refs_in_amplicons[amplicon] / s.amplicon_totals[amplicon]
                    < self.config.min_frs
                ):
                    passing.append(False)
                else:
                    passing.append(True)
            if len(passing) > 1 and not all(passing[0] == e for e in passing):
                #                print(f"{s.reference_pos}\t{passing}")
                #                for amplicon, total in s.amplicon_totals.items():
                #                    print(f"\t{amplicon.name}\t{total}\t{s.refs_in_amplicons[amplicon]} / {s.amplicon_totals[amplicon]}")

                # this is a failure
                return True
            return False

        def test_strand_bias(s: Stats) -> bool:
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

    def ref_to_consensus(self, p: Index1) -> Index1:  # -> Optional[Index1]:
        if p in self._ref_to_consensus:
            return self._ref_to_consensus[p]
        raise Exception

    def consensus_to_ref(self, p: Index1) -> Index1:  # Optional[Index1]:
        if p in self._consensus_to_ref:
            return self._consensus_to_ref[p]
        raise Exception
        # return None

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
        sequence: list[str] = list(self.ref)
        self.qc = {}

        for p, stats in enumerate(self.seq):
            if p >= len(sequence):
                raise Exception(
                    f"Invalid condition: mapped position {p} greater than consensus length {len(sequence)}"
                )
            position = Index0(p)

            self.summary["consensus_length"] += 1
            # print(position, self.consensus_to_ref[position + 1], file=sys.stderr)
            if sequence[position] == "N":
                # if a position is already masked by an upstream process skip it
                self.summary["already_masked"] += 1
                self.summary["total_masked"] += 1
                continue
            elif stats.check_for_failure(self.filters):
                for failure in stats.get_failures():
                    self.summary[failure] += 1
                self.summary["total_masked"] += 1
                sequence[position] = "N"
                self.qc[str(position)] = stats.log
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

            stats: Stats = self.seq[Index0(cons_coord - 1)]
            info_field = stats.info()
            vcf_filters = original_filters
            if stats.position_failed is None:
                print(f"Warning: attemped to evaluate N basecall", file=sys.stderr)
                continue
            if stats.position_failed:
                if original_filters == "PASS":
                    vcf_filters = ";".join(stats.get_failures())
                else:
                    vcf_filters = ";".join([original_filters, *vcf_filters])

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
