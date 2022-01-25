import sys

import mappy as mp
import pysam

from viridian_workflow.primers import AmpliconSet, get_tags, load_amplicon_schemes

from collections import namedtuple, defaultdict

BaseProfile = namedtuple(
    "BaseProfile",
    [
        "base",
        "cons_base",
        "ref_base",
        "in_primer",
        "forward_strand",
        "amplicon_name",
        "reference_pos",
    ],
)

Config = namedtuple("Config", ["min_frs", "min_depth", "trim_5prime", "log_liftover"])

default_config = Config(
    min_frs=0.7, min_depth=10, trim_5prime=False, log_liftover=False
)


def mask_sequence(sequence, position_stats, config=default_config):
    sequence = list(sequence)
    summary = defaultdict(int)
    qc = {}
    for position, stats in position_stats.items():
        if position >= len(sequence):
            print(
                f"Invalid condition: mapped position {position} greater than consensus length {len(sequence)}",
                file=sys.stderr,
            )
            continue
        summary["consensus_length"] += 1
        if sequence[position] == "N":
            # if a position is already masked by an upstream process skip it
            summary["already_masked"] += 1
            summary["total_masked"] += 1
            continue
        elif stats.check_for_failure(summary=summary):
            summary["total_masked"] += 1
            sequence[position] = "N"
            qc[position] = stats.log
        qc["masking_summary"] = summary
    return "".join(sequence), qc


def test_bias(n, trials, threshold=0.3):
    """Test whether a number of boolean trials fit the binomial
    distribution (returns true if unbiased)
    """
    if trials == 0:
        return False

    # TODO: decide whether to include scipy to use binom_test
    # or a pure python implementation.
    bias = abs(0.5 - (n / trials))

    return bias >= threshold


class Stats:
    def __init__(
        self,
        ref_base=None,
        cons_base=None,
        reference_position=None,
        config=default_config,
    ):
        self.alts_in_primer = 0
        self.refs_in_primer = 0

        self.alts_in_amplicons = defaultdict(int)
        self.refs_in_amplicons = defaultdict(int)
        self.amplicon_totals = defaultdict(int)

        self.alts_forward = 0
        self.refs_forward = 0

        self.alts = 0
        self.refs = 0
        self.total = 0
        self.log = []
        self.total_reads = 0

        self.alts_matching_refs = 0
        self.reference_pos = reference_position
        self.ref_base = ref_base
        self.cons_base = cons_base

        self.config = config

    def add_alt(self, profile, alt=None):
        if not self.reference_pos:
            self.reference_pos = profile.reference_pos

        if profile.in_primer:
            self.alts_in_primer += 1
            return

        self.alts += 1
        if profile.base == profile.ref_base:
            self.alts_matching_refs += 1

        if profile.amplicon_name:
            # when unambiguous amplicon call cannot be made, do not
            # consider (or consider differently)
            self.alts_in_amplicons[profile.amplicon_name] += 1
            self.amplicon_totals[profile.amplicon_name] += 1

        if profile.forward_strand:
            self.alts_forward += 1

        self.total += 1

    def add_ref(self, profile):
        if not self.reference_pos:
            self.reference_pos = profile.reference_pos

        if profile.in_primer:
            self.refs_in_primer += 1
            return

        self.refs += 1

        if profile.amplicon_name:
            self.refs_in_amplicons[profile.amplicon_name] += 1
            self.amplicon_totals[profile.amplicon_name] += 1

        if profile.forward_strand:
            self.refs_forward += 1

        self.total += 1

    def check_for_failure(self, summary=None):
        """return whether a position should be masked

        optionally accepts a mutable reference to a dictionary that summarises masking decisions
        """

        position_failed = False
        bias_threshold = 0.3

        if self.total < self.config.min_depth:
            self.log.append(
                f"Insufficient depth to evaluate consensus; {self.total} < {self.config.min_depth}. {self.total_reads} including primer regions."
            )
            if summary:
                summary["insufficient_depth"] += 1
            return True  # position failed

        # test total percentage of bases supporting consensus
        if self.refs / self.total < self.config.min_frs:
            self.log.append(
                f"Insufficient support of consensus base; {self.refs} / {self.total} < {self.config.min_frs}. {self.total_reads} including primer regions."
            )
            if summary:
                summary["low_frs"] += 1
            return True

        # look for overrepresentation of alt alleles in regions covered
        # by primer sequences. This is reported but not as a failure

        if not self.alts_in_primer + self.refs_in_primer == 0:

            if not test_bias(
                self.refs_in_primer / (self.alts_in_primer + self.refs_in_primer),
                self.refs / (self.alts + self.refs),
                threshold=bias_threshold,
            ):
                self.log.append(
                    f"Consensus base calls are biased in primer region; {self.refs_in_primer}/{self.alts_in_primer}+{self.refs_in_primer} vs. {self.refs}/{self.alts}+{self.refs}"
                )
                # position_failed = True

        # strand bias in alt calls
        if not test_bias(self.refs_forward, self.refs, threshold=bias_threshold):
            self.log.append(
                f"Strand bias for reads with consensus alleles; {self.refs_forward} / {self.refs}"
            )
            # position_failed = True

        # amplicon bias
        for amplicon, total in self.amplicon_totals.items():
            amplicon_frs = self.refs_in_amplicons[amplicon] / total
            if amplicon_frs < self.config.min_frs:
                self.log.append(
                    f"Per-amplicon FRS failure, amplicon {amplicon}: {self.refs_in_amplicons[amplicon]} / {total}"
                )
                if summary:
                    summary["low_amplicon_specific_frs"] += 1
                position_failed = True
        return position_failed

    def __str__(self):
        f = []
        if len(self.amps_total) > 1:
            return "-".join([f"{k}:{v}" for k, v in self.amps_total.items()])
        if self.alts / self.total > 0.2:
            return f"{self.alts}/{self.total}"
        return "-"


def cigar_to_alts(ref, query, cigar, q_pos=0, pysam=False):
    """Interpret cigar string and query sequence in reference
    coords from mappy (count, op) or pysam (op, count)
    """
    positions = []
    r_pos = 0

    for count, op in cigar:
        if pysam:
            count, op = op, count  # this makes me sad
        if op == 0:
            # match/mismatch
            for i in range(count):
                positions.append((r_pos + i, query[q_pos + i]))
            q_pos += count
            r_pos += count

        elif op == 1:
            pass
            # insertion
            #            positions.append((q_pos, query[q_pos : q_pos + count]))
            q_pos += count
            r_pos += 0

        elif op == 2:
            # deletion
            for n in range(count):
                positions.append((r_pos + n, "-"))
            r_pos += count

        elif op == 3:
            # ref_skip
            pass

        elif op == 4:
            # soft clip
            # q_pos += count
            # may not need to be considered if q_pos offset is set
            # TODO verify
            pass

        elif op == 5:
            # hard clip
            pass

        else:
            raise Exception(f"invalid cigar op {op}")

    return positions


def remap(
    ref_genome,
    consensus_fasta,
    minimap_presets,
    amplicon_set,
    tagged_bam,
    config=default_config,
):
    stats = {}
    cons = mp.Aligner(consensus_fasta, preset=minimap_presets)
    if len(cons.seq_names) != 1:
        Exception(f"Consensus fasta {consensus_fasta} has more than one sequence")
    consensus_seq = cons.seq(cons.seq_names[0])

    ref = mp.Aligner(ref_genome)
    reference_seq = ref.seq(ref.seq_names[0])

    multi_amplicons = 0
    no_amplicons = 0
    tagged = 0

    for r in pysam.AlignmentFile(tagged_bam):
        a = cons.map(r.seq)  # remap to consensus
        alignment = None
        for x in a:
            if x.is_primary:
                alignment = x

        if not alignment:
            continue

        assert alignment.q_en > alignment.q_st
        assert alignment.r_en > alignment.r_st

        amplicons = get_tags(amplicon_set, r)

        amplicon = None

        if len(amplicons) == 1:
            amplicon = amplicons[0]
            tagged += 1
        elif len(amplicons) > 1:
            multi_amplicons += 1
            continue
        else:
            no_amplicons += 1
            continue

        strand = False  # strand is forward
        if alignment.strand == 0:
            # should be error
            raise Exception()
        elif alignment.strand == 1:
            strand = True

        # problem: if the consensus is shorter than ref (like if a primer sequence is ellided) then we can't zip the cigar_to_alt lists
        # slice_start = max(alignment.r_st, r.reference_start)
        # slice_end = min(alignment.r_en, r.reference_end)
        # but the start and end coords in both cases are in different reference frames

        cons_alts = cigar_to_alts(
            consensus_seq[alignment.r_st : alignment.r_en],
            r.seq,
            alignment.cigar,
            q_pos=alignment.q_st,
        )

        ref_alts = cigar_to_alts(
            reference_seq[r.reference_start : r.reference_end],
            r.query_alignment_sequence,
            r.cigar,
            q_pos=r.query_alignment_start,
            pysam=True,
        )

        # TODO test softclip

        for (
            read_pos,
            ((read_cons_pos, read_cons_base), (read_ref_pos, read_ref_base)),
        ) in enumerate(zip(cons_alts, ref_alts)):
            consensus_position = read_cons_pos + alignment.r_st
            reference_position = read_ref_pos + r.reference_start

            if consensus_position >= len(consensus_seq):
                print(
                    f"Warning: Position {read_pos}+{alignment.r_st} extends beyond length of consensus ({len(consensus_seq)}) {consensus_position}:{reference_position} ({read_pos}: {read_cons_base})",
                    file=sys.stderr,
                )
                continue

            cons_base = consensus_seq[consensus_position]
            ref_base = reference_seq[reference_position]

            # TODO resolve assumption: if there is an ambiguous amplicon id, in_primer is false
            # in_primer = amplicon.position_in_primer(reference_position)
            in_primer = False

            # the trim 5' option assumes that all bases max_primer_len from the 5' end of the read are
            # guaranteed to be inside of a primer
            if config.trim_5prime and read_pos < amplicon.max_length:
                in_primer = True

            base_profile = BaseProfile(
                read_cons_pos,
                cons_base,
                ref_base,
                in_primer,
                strand,
                amplicon.name,
                reference_position,
            )

            if consensus_position not in stats:
                stats[consensus_position] = Stats(
                    ref_base=ref_base,
                    cons_base=cons_base,
                    reference_position=reference_position,
                    config=config,
                )
            stats[consensus_position].total_reads += 1

            if (
                config.log_liftover
                and reference_position != stats[consensus_position].reference_pos
            ):
                print(
                    f"{config} {config.min_frs} {config.trim_5prime} {config.log_liftover} WARNING: liftover error; Ref:{reference_position}:{read_ref_base} (expected {stats[consensus_position].reference_pos})\t{r.cigarstring}\tCons:{consensus_position}:{read_cons_base}\t{alignment.cigar_str}\t{ref_base}/{cons_base}\t{r.seq}",
                    file=sys.stderr,
                )

            if read_cons_base != cons_base:
                stats[consensus_position].add_alt(base_profile)
            else:
                stats[consensus_position].add_ref(base_profile)

    print(
        f"reads tagged with more than one amplicon: {multi_amplicons}, with zero: {no_amplicons}. tagged: {tagged}.",
        file=sys.stderr,
    )
    return stats


def mask(fasta, stats, outpath=None, name=None, config=default_config):
    if not outpath:
        output = f"{name}.masked.fasta"
    ref = mp.Aligner(fasta)
    if len(ref.seq_names) != 1:
        Exception(f"Consensus fasta {fasta} has more than one sequence")
    sequence = ref.seq(ref.seq_names[0])
    if not name:
        name = ref.seq_names[0]

    masked, log = mask_sequence(sequence, stats, config=default_config)

    # write masked fasta
    with open(outpath, "w") as maskfd:
        print(f">{name}\n{masked}", file=maskfd, end="")
    return outpath, log
