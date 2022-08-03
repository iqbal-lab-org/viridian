from __future__ import annotations

import sys

from collections import defaultdict
from dataclasses import dataclass
from typing import NewType, Callable, Optional, Any
from pathlib import Path

import mappy as mp  # type: ignore

from viridian_workflow.utils import Index0, Index1
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


class EvaluatedStats:
    """Per-position base counts, evaluated after the pileup is built
    """

    def __init__(
        self,
        base: str,
        aux_reference_pos: Index0,
        ref_alt_total: tuple[int, int],  # refs, alts
        total_reads: int,
        # (ref, alt) calls by amplicon
        calls_by_amplicon: dict[Amplicon, tuple[int, int]],
        position_failed: bool,
        failures: dict[str, str],
        primer_bases_considered: tuple[int, int],
        primer_bases_total: int,
        alt_bases: defaultdict[str, int],
    ):
        self.base = base
        self.total = (
            ref_alt_total[0] + ref_alt_total[1]
        )  # number of base calls considered
        self.total_reads = total_reads  # total bases covering this position
        self.refs: int = ref_alt_total[0]
        self.alts: int = ref_alt_total[1]
        self.calls_by_amplicon: dict[Amplicon, tuple[int, int]] = calls_by_amplicon
        self.position_failed: bool = position_failed
        self.failures: dict[str, str] = failures
        self.primer_bases_considered: tuple[int, int] = primer_bases_considered
        self.primer_bases_total: int = primer_bases_total
        self.alt_bases: dict[str, int] = alt_bases
        self.aux_reference_pos: Index0 = aux_reference_pos

    def info(self) -> str:
        """Output position stats as VCF INFO field
        """
        amplicon_totals = ",".join(
            [
                f"{amplicon.name}:{str(self.calls_by_amplicon[amplicon][0])}/{str(self.calls_by_amplicon[amplicon][1])}"
                for amplicon in self.calls_by_amplicon
            ]
        )
        info_fields = [
            f"primer={self.primer_bases_considered[0]}/{self.primer_bases_considered[1]}",
            f"total_primer_bases={self.primer_bases_total}",
            f"total{self.total}",
            f"amplicon_overlap={len(self.calls_by_amplicon)}",
            f"amplicon_totals={amplicon_totals}",
        ]
        return ";".join(info_fields)

    def tsv_row(self) -> str:
        amplicon_totals = ",".join(
            [
                f"{amplicon.name}:{str(self.calls_by_amplicon[amplicon][0])}/{str(self.calls_by_amplicon[amplicon][1])}"
                for amplicon in self.calls_by_amplicon
            ]
        )
        alts = ";".join([f"{k}:{v}" for k, v in self.alt_bases.items()])
        row = "\t".join(
            map(
                str,
                [
                    self.aux_reference_pos,
                    self.base,
                    alts,
                    self.refs,
                    self.alts,
                    self.primer_bases_considered[0],
                    self.primer_bases_considered[1],
                    self.primer_bases_total,
                    len(self.calls_by_amplicon),
                    amplicon_totals,
                ],
            )
        )

        return row

    def __str__(self) -> str:
        alts = " ".join([f"{alt}:{count}" for alt, count in self.alt_bases.items()])
        return f"{self.base}:{self.refs};{alts}"


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

        self.aux_reference_pos: Index0 = aux_reference_pos
        self.base: str = base
        self.config = config

        self.alt_bases: defaultdict[str, int] = defaultdict(int)

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

            self.alt_bases[profile.base] += 1

        else:
            # ref calls
            f_calls, r_calls = self.basecounts[profile.amplicon].refs
            if profile.forward_strand:
                self.basecounts[profile.amplicon].refs = (f_calls + 1, r_calls)
            else:
                self.basecounts[profile.amplicon].refs = (f_calls, r_calls + 1)

    def evaluate(self, filters: dict[str, tuple[Filter, FilterMsg]]) -> EvaluatedStats:
        """Evaluate the accumulated positon stats
        """
        self.total: int = 0
        self.refs: int = 0
        self.alts: int = 0
        self.refs_in_primer: int = 0
        self.alts_in_primer: int = 0
        self.total_reads: int = 0

        # if more than one amplicon covers this position
        # we can decide to discount primer-base calls
        consider_primers: bool = len(self.basecounts.keys()) > 1
        self.primer_bases_considered: int = 0

        calls_by_amplicon: dict[Amplicon, tuple[int, int]] = {}
        for amplicon, bc in self.basecounts.items():
            if bc.in_primer:
                self.refs_in_primer += bc.refs[0] + bc.refs[1]
                self.alts_in_primer += bc.alts[0] + bc.alts[1]

                if consider_primers:
                    self.primer_bases_considered += (
                        self.refs_in_primer + self.alts_in_primer
                    )
                else:
                    continue

            self.refs += bc.refs[0] + bc.refs[1]
            self.alts += bc.alts[0] + bc.alts[1]
            self.total += self.refs + self.alts
            calls_by_amplicon[amplicon] = (
                bc.refs[0] + bc.refs[1],
                bc.alts[0] + bc.alts[1],
            )

        failures: dict[str, str] = {}
        position_failed: bool = False

        for filter_name, (filter_func, msg_format) in filters.items():
            # if the filter fails
            if filter_func(self):
                failures[filter_name] = msg_format(self)
                position_failed = True

        return EvaluatedStats(
            self.base,
            self.aux_reference_pos,
            (self.refs, self.alts),
            self.total_reads,
            calls_by_amplicon,
            position_failed,
            failures,
            (self.refs_in_primer, self.alts_in_primer),
            self.primer_bases_considered,
            self.alt_bases,
        )


Filter = Callable[[Stats], bool]
FilterMsg = Callable[[Stats], str]


# filter closures take a Stats struct and return True on failure
# def test_amplicon_bias(s: Stats) -> bool:
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

# def test_strand_bias(s: EvaluatedStats) -> bool:
#    passing = []
#    for amplicon, total in s.amplicon_totals.items():
#        if total < self.config.min_depth:
#            continue
#        if s.amplicon_totals[amplicon] == 0:
#            continue
#        if (
#            s.refs_in_forward_strands[amplicon] / s.amplicon_totals[amplicon]
#            < self.config.min_frs
#        ):
#            passing.append(False)
#        else:
#            passing.append(True)
#    if len(passing) > 1 and not all(passing[0] == e for e in passing):
#        # True = failure
#        return True
#    return False


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
    """A pileup is an array of Stats objects indexed by position in a sequence
    """

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

        intermediate_sequence: list[Stats] = []

        for i, base in enumerate(self.consensus_seq):
            p = Index1(i + 1)
            intermediate_sequence.append(
                Stats(Index0(self.consensus_to_ref(p) - 1), base,)
            )

        self.filters: dict[str, tuple[Filter, FilterMsg]] = {
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
        }

        # initialise summary for each filter
        self.summary = defaultdict(int)
        for f in self.filters:
            self.summary[f] = 0

        for amplicon, fragments in readstore.amplicons.items():
            for fragment in fragments:
                for read in fragment.reads:
                    alns = aligner.map(read.seq)  # remap to consensus
                    alignment = None
                    for x in alns:
                        # test that the re-alignment is still within the
                        # original amplicon call
                        # if (
                        #    self.consensus_to_ref(Index1(x.r_st + 1)) is None
                        #    or self.consensus_to_ref(Index1(x.r_en + 1)) is None
                        # ):
                        #    continue
                        if (
                            x.is_primary  # this is always true with mappy
                            and Index0(self.consensus_to_ref(Index1(x.r_st + 1)) - 1)
                            > (amplicon.start - 10)
                            and Index0(self.consensus_to_ref(Index1(x.r_en + 1)) - 1)
                            < (amplicon.end + 10)
                        ):
                            alignment = x

                    if alignment is None:
                        continue

                    assert alignment.q_en > alignment.q_st
                    assert alignment.r_en > alignment.r_st

                    aln: list[tuple[Index0, str]] = parse_cigar(
                        self.consensus_seq, read.seq, alignment
                    )

                    # testing for indel conditions:
                    # ex = "".join(map(lambda x: x[1] if len(x[1]) == 1 else x[1], aln))
                    # c = consensus_seq[alignment.r_st : alignment.r_en]

                    # TODO: this is now made redundant. We could store the results
                    # of match_primers from push_fragments
                    primers = amplicon.match_primers(fragment)

                    for consensus_pos, call in aln:
                        reference_pos = Index0(
                            self.consensus_to_ref(Index1(consensus_pos + 1)) - 1
                        )
                        if consensus_pos >= len(self.consensus_seq):
                            print(
                                f"consensus pos out of bounds: {consensus_pos}/{len(self.consensus_seq)}",
                                file=sys.stderr,
                            )
                            continue
                        reference_pos = Index0(
                            self.consensus_to_ref(Index1(consensus_pos + 1)) - 1
                        )

                        in_primer = amplicon.position_in_primer(reference_pos)

                        profile = BaseProfile(
                            call, in_primer, read.is_reverse, amplicon,
                        )
                        intermediate_sequence[consensus_pos].update(profile)

        # Finalise the pileup object by evaluating bases

        for stat in intermediate_sequence:
            self.seq.append(stat.evaluate(self.filters))

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
        sequence, qc, summary = self._mask(self.consensus_seq, self.seq)
        self.qc = qc
        for key, value in summary.items():
            self.summary[key] = value
        return sequence

    @staticmethod
    def _mask(
        consensus_seq: str, stats_seq: list[EvaluatedStats]
    ) -> tuple[str, dict[str, Any], defaultdict[str, Any]]:
        """Evaluate all positions and determine if they pass filters
        """
        sequence: list[str] = list(consensus_seq)
        qc: dict[str, Any] = {}
        summary: defaultdict[str, int] = defaultdict(int)

        log: list[str] = []
        failures: list[str] = []

        # filters: dict[str, tuple[Filter, FilterMsg]]
        for p, stats in enumerate(stats_seq):
            if p >= len(sequence):
                raise Exception(
                    f"Invalid condition: mapped position {p} greater than consensus length {len(sequence)}"
                )
            position = Index0(p)
            position_failed = False

            # if a position is already masked by an upstream process skip it
            if sequence[position] == "N":
                summary["already_masked"] += 1
                summary["total_masked"] += 1
                continue

            summary["consensus_length"] += 1

            if stats.position_failed:
                for failure in stats.failures:
                    summary[failure] += 1
                summary["total_masked"] += 1
                sequence[position] = "N"
                qc[str(position)] = stats.failures
            qc["masking_summary"] = summary
        return "".join(sequence), qc, summary

    def dump_tsv(self, tsv: Path) -> Path:
        fd = open(tsv, "w")
        for pos, stats in enumerate(self.seq):
            cons_pos = self.consensus_to_ref(Index1(pos + 1))
            print(f"{cons_pos}\t{stats.aux_reference_pos}\t{stats.tsv_row()}", file=fd)
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

            stats = self.seq[Index0(cons_coord - 1)]

            info_field = stats.info()
            vcf_filters = original_filters

            if stats.position_failed:
                if original_filters == "PASS":
                    vcf_filters = ";".join(stats.failures.keys())
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
