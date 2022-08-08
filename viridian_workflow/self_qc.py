from __future__ import annotations

import sys

from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, Optional, Any
from pathlib import Path

import mappy as mp  # type: ignore

from viridian_workflow.utils import Index0, Index1, in_range
from viridian_workflow.primers import Amplicon
from viridian_workflow.readstore import ReadStore


@dataclass
class Config:
    # TODO move this somewhere cleaner
    min_frs: float
    min_depth: int


default_config = Config(0.7, 10)


@dataclass(frozen=True)
class BaseProfile:
    """Per-base info extracted while reads are traversed"""

    base: str
    in_primer: bool
    forward_strand: bool
    amplicon: Amplicon


@dataclass
class Calls:
    refs: int
    alts: int


class EvaluatedStats:
    """Per-position base counts, evaluated after the pileup is built"""

    def __init__(self, stats):
        self.base: str = stats.base
        self.aux_reference_pos: Index0 = stats.aux_reference_pos

        self.depth: int = 0
        self.total: Calls = Calls(0, 0)
        self.primer_calls: Calls = Calls(0, 0)
        self.primer_calls_ignored: Calls = Calls(0, 0)
        self.calls_by_amplicon: dict[Amplicon, Calls] = {}
        self.multiple_amplicon_support = len(stats.baseprofiles) > 1

        self.alt_bases: defaultdict[str, int] = defaultdict(int)

        for amplicon, profiles in stats.baseprofiles.items():
            if amplicon not in self.calls_by_amplicon:
                self.calls_by_amplicon[amplicon] = Calls(0, 0)
            for profile in profiles:
                is_ref = profile.base == self.base
                self.depth += 1
                self.alt_bases[profile.base] += 1

                if profile.in_primer:
                    if is_ref:
                        self.primer_calls.refs += 1
                        if self.multiple_amplicon_support:
                            self.primer_calls_ignored.refs += 1
                        else:
                            self.total.refs += 1
                            self.calls_by_amplicon[amplicon].refs += 1

                    else:
                        self.primer_calls.alts += 1
                        if self.multiple_amplicon_support:
                            self.primer_calls_ignored.alts += 1
                        else:
                            self.total.alts += 1
                            self.calls_by_amplicon[amplicon].alts += 1
                else:
                    if is_ref:
                        self.total.refs += 1
                        self.calls_by_amplicon[amplicon].refs += 1
                    else:
                        self.total.alts += 1
                        self.calls_by_amplicon[amplicon].alts += 1

    def evaluate(
        self, filters: dict[str, tuple[Filter, FilterMsg]]
    ) -> tuple[bool, dict[str, str]]:
        """Returns True if any filter fails"""
        failures: dict[str, str] = {}
        fail = False
        for filter_name, (f, msg) in filters.items():
            if f(self):
                fail = True
                failures[filter_name] = msg(self)
        return fail, failures

    def info(self) -> str:
        """Output position stats as VCF INFO field"""
        amplicon_totals = ",".join(
            [
                f"{str(self.calls_by_amplicon[amplicon].refs)}/{str(self.calls_by_amplicon[amplicon].alts)}"
                for amplicon in self.calls_by_amplicon
            ]
        )
        amplicon_names = ",".join(
            [amplicon.name for amplicon in self.calls_by_amplicon]
        )
        info_fields = [
            f"primer_calls_ignored={self.primer_calls_ignored.refs}/{self.primer_calls_ignored.alts}",
            f"total_primer_bases={self.primer_calls.refs}/{self.primer_calls.alts}",
            f"unfiltered_depth={self.depth}",
            f"total={self.total.refs}/{self.total.alts}",
            f"amplicon_overlap={len(self.calls_by_amplicon.keys())}",
            f"amplicon_totals={amplicon_totals}",
            f"amplicon_names={amplicon_names}",
        ]
        return ";".join(info_fields)

    def tsv_row(self) -> str:
        alts = ";".join([f"{k}:{v}" for k, v in self.alt_bases.items()])
        row = "\t".join(
            map(
                str,
                [
                    self.aux_reference_pos,
                    self.base,
                    self.total.refs,
                    self.total.alts,
                    self.depth,
                    self.primer_calls_ignored.refs,
                    self.primer_calls_ignored.alts,
                    self.primer_calls.refs,
                    self.primer_calls.alts,
                    len(self.calls_by_amplicon),
                    alts,
                    self.info(),
                ],
            )
        )

        return row

    def __str__(self) -> str:
        alts = ",".join([f"{alt}:{count}" for alt, count in self.alt_bases.items()])
        return f"{self.base}:{self.total.refs};{alts}"


class Stats:
    def __init__(
        self,
        aux_reference_pos: Index0,  # may be non-consensus sequence index
        base: str,
    ):

        self.baseprofiles: dict[Amplicon, list[BaseProfile]] = {}

        self.aux_reference_pos: Index0 = aux_reference_pos
        self.base: str = base  # reference base

    def update(self, profile: BaseProfile):
        """Accumulate per-position base calling stats"""
        if profile.amplicon not in self.baseprofiles:
            self.baseprofiles[profile.amplicon] = []
        self.baseprofiles[profile.amplicon].append(profile)


Filter = Callable[[EvaluatedStats], bool]
FilterMsg = Callable[[EvaluatedStats], str]


def parse_msa(msa: Path) -> tuple[dict[Index1, Index1], dict[Index1, Index1]]:

    # 1-based index translation tables
    ref_to_consensus: dict[Index1, Index1] = {}
    consensus_to_ref: dict[Index1, Index1] = {}

    with open(msa) as msa_fd:
        ref_seq = msa_fd.readline().strip()  # reference
        con_seq = msa_fd.readline().strip()  # consensus

        if len(msa_fd.readline()) > 0:
            raise Exception("Invalid multiple sequence alignment file")

        if len(con_seq) != len(ref_seq):
            raise Exception("Both sequences in MSA must be same length")

        # this is valid for testing when the consensus is complete
        # but some of the unittests break this assumption
        # assert con_seq.replace("-", "") == self.consensus_seq

        ref_pos = 0
        con_pos = 0

        for con_base, ref_base in zip(con_seq, ref_seq):
            if con_base != "-" and ref_base != "-":
                ref_pos += 1
                con_pos += 1
                ref_to_consensus[Index1(ref_pos)] = Index1(con_pos)
                consensus_to_ref[Index1(con_pos)] = Index1(ref_pos)
            elif con_base == "-" and ref_base != "-":
                ref_pos += 1
                ref_to_consensus[Index1(ref_pos)] = Index1(con_pos)
            elif ref_base == "-" and con_base != "-":
                con_pos += 1
                consensus_to_ref[Index1(con_pos)] = Index1(ref_pos)

    return ref_to_consensus, consensus_to_ref


class Pileup:
    """A pileup is an array of Stats objects indexed by position in a sequence"""

    def __init__(
        self,
        consensus_fasta: Path,
        readstore: ReadStore,
        msa: Optional[Path] = None,
        config: Config = default_config,
        minimap_presets: Optional[str] = None,
        seq: Optional[str] = None,  # Only for legacy tests
    ):
        self.config: Config = config
        self.seq: list[EvaluatedStats] = []

        # remap readstore to consensus sequence

        aligner = mp.Aligner(str(consensus_fasta), preset=minimap_presets, n_threads=1)
        if seq is None:
            if len(aligner.seq_names) != 1:
                Exception(
                    f"Consensus fasta {consensus_fasta} has more than one sequence"
                )
            self.consensus_seq: str = aligner.seq(aligner.seq_names[0])
        else:
            print(
                f"Warning: overriding sequence. Ignoring {consensus_fasta}. Use only for testing.",
                file=sys.stderr,
            )
            self.consensus_seq = seq

        if msa is None:
            raise Exception("Building pileup without MSA is not supported")
        rtoc, ctor = parse_msa(msa)

        # cannot destructure with type annotations?
        self._ref_to_consensus: dict[Index1, Index1] = rtoc
        self._consensus_to_ref: dict[Index1, Index1] = ctor

        _pileup: list[Stats] = []

        for i, base in enumerate(self.consensus_seq):
            p = Index1(i + 1)
            _pileup.append(
                Stats(
                    Index0(self.consensus_to_ref(p) - 1),
                    base,
                )
            )

        self.filters: dict[str, tuple[Filter, FilterMsg]] = {
            "low_depth": (
                lambda s: (s.total.refs + s.total.alts) < self.config.min_depth,
                lambda s: f"Insufficient depth; {s.total.refs + s.total.alts} < {self.config.min_depth}. {s.depth} including primer regions.",
            ),
            "low_frs": (
                lambda s: s.total.refs / (s.total.refs + s.total.alts)
                < self.config.min_frs
                if (s.total.refs + s.total.alts) > 0
                else False,
                lambda s: f"Insufficient support of consensus base; {s.total.refs} / {s.total.refs + s.total.alts} < {self.config.min_frs}. {s.depth} including primer regions.",
            ),
        }

        # initialise summary for each filter
        self.summary = defaultdict(int)
        for f in self.filters:
            self.summary[f] = 0

        for amplicon, fragments in readstore.amplicons.items():
            for fragment in fragments:
                l_primer, r_primer = amplicon.match_primers(fragment)
                primers = [
                    primer for primer in [l_primer, r_primer] if primer is not None
                ]
                for read in fragment.reads:
                    alns = aligner.map(read.seq)  # remap to consensus
                    alignment = None
                    for x in alns:
                        # test that the re-alignment is still within the
                        # original amplicon call
                        if x.is_primary and in_range(  # this is always true with mappy
                            (Index0(amplicon.start - 10), Index0(amplicon.end + 10)),
                            Index0(self.consensus_to_ref(Index1(x.r_st + 1)) - 1),
                        ):
                            alignment = x

                    if alignment is None:
                        continue

                    assert alignment.q_en > alignment.q_st
                    assert alignment.r_en > alignment.r_st

                    aln: list[tuple[Index0, str]] = parse_cigar(read.seq, alignment)

                    # testing for indel conditions:
                    # ex = "".join(map(lambda x: x[1] if len(x[1]) == 1 else x[1], aln))
                    # c = consensus_seq[alignment.r_st : alignment.r_en]

                    for consensus_pos, call in aln:
                        reference_pos = Index0(
                            self.consensus_to_ref(Index1(consensus_pos + 1)) - 1
                        )
                        if consensus_pos >= len(self.consensus_seq):
                            print(
                                f"consensus pos out of bounds: {consensus_pos} >= {len(self.consensus_seq)}",
                                file=sys.stderr,
                            )
                            continue

                        in_primer = any(
                            [
                                in_range(
                                    (primer.ref_start, primer.ref_end), reference_pos
                                )
                                for primer in primers
                            ]
                        )

                        _pileup[consensus_pos].update(
                            BaseProfile(
                                call,
                                in_primer,
                                read.is_reverse,
                                amplicon,
                            )
                        )

        # Finalise the pileup object by evaluating bases

        for stat in _pileup:
            self.seq.append(EvaluatedStats(stat))

    def ref_to_consensus(self, p: Index1) -> Index1:
        if p in self._ref_to_consensus:
            return self._ref_to_consensus[p]
        return Index1(0)

    def consensus_to_ref(self, p: Index1) -> Index1:
        if p in self._consensus_to_ref:
            return self._consensus_to_ref[p]
        return Index1(0)

    def __getitem__(self, pos: Index0) -> EvaluatedStats:
        if pos > len(self.seq):
            raise Exception(f"position too big: {pos} {len(self.seq)}")
        return self.seq[pos]

    def __setitem__(self, pos: Index0, profile: BaseProfile) -> None:
        raise Exception("don't set Pileup positions")

    def __len__(self) -> int:
        return len(self.seq)

    def mask(self) -> str:
        sequence, qc, summary = self._mask(self.consensus_seq, self.seq, self.filters)
        self.qc = qc
        for key, value in summary.items():
            self.summary[key] = value
        return sequence

    @staticmethod
    def _mask(
        consensus_seq: str, stats_seq: list[EvaluatedStats], filters
    ) -> tuple[str, dict[str, Any], defaultdict[str, Any]]:
        """Evaluate all positions and determine if they pass filters"""
        sequence: list[str] = list(consensus_seq)
        qc: dict[str, Any] = {}
        summary: defaultdict[str, int] = defaultdict(int)

        failures: dict[str, str] = {}

        # filters: dict[str, tuple[Filter, FilterMsg]]
        for p, stats in enumerate(stats_seq):
            if p >= len(sequence):
                raise Exception(
                    f"Invalid condition: mapped position {p} greater than consensus length {len(sequence)}"
                )
            position = Index0(p)

            # if a position is already masked by an upstream process skip it
            if sequence[position] == "N":
                summary["already_masked"] += 1
                summary["total_masked"] += 1
                continue

            position_failed, failures = stats.evaluate(filters)
            assert sequence[position] == stats.base

            summary["consensus_length"] += 1

            if position_failed:
                for failure in failures:
                    summary[failure] += 1
                summary["total_masked"] += 1
                sequence[position] = "N"
                qc[str(position)] = failures
            qc["masking_summary"] = summary
        return "".join(sequence), qc, summary

    def dump_tsv(self, tsv: Path) -> Path:
        fd = open(tsv, "w")
        for pos, stats in enumerate(self.seq):
            ref_pos = self.consensus_to_ref(Index1(pos + 1))
            print(f"{ref_pos}\t{pos}\t{stats.tsv_row()}", file=fd)
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
                _,
                fmt,
                *r,
            ) = line.split("\t")

            pos: Index1 = Index1(int(pos_token))
            cons_coord: Index1 = self.ref_to_consensus(pos)

            stats = self.seq[Index0(cons_coord - 1)]

            info_field = stats.info()
            vcf_filters = original_filters

            position_failed, failures = stats.evaluate(self.filters)
            if position_failed:
                if original_filters == "PASS":
                    vcf_filters = ";".join(failures.keys())
                else:
                    vcf_filters = ";".join([original_filters, *vcf_filters])

            records.append(
                (chrom, pos, mut_id, ref, alt, qual, vcf_filters, info_field, fmt, *r)
            )

        return header, records


def parse_cigar(query: str, alignment: Any) -> list[tuple[Index0, str]]:
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
