import sys

import mappy as mp
import pysam

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

Config = namedtuple(
    "Config",
    ["min_frs", "min_depth", "trim_5prime", "log_liftover", "test_amplicon_frs"],
)

default_config = Config(
    min_frs=0.7,
    min_depth=10,
    trim_5prime=False,
    log_liftover=False,
    test_amplicon_frs=False,
)


class Pileup:
    """A pileup is an array of Stats objects indexed by position in a reference
    """

    def __init__(self, ref):
        self.ref = ref
        for r in ref:
            self.seq.append(Stats(ref=r))

    def __getitem__(self, pos):
        return self.seq[pos]

    def __setitem__(self, pos, profile):
        self.seq[pos].update(profile)

    def mask():
        sequence = list(self.ref)
        summary = defaultdict(int)
        self.qc = {}
        for position, stats in enumerate(self.seq):
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
                self.qc[position] = stats.log
            self.qc["masking_summary"] = summary
        return "".join(sequence)

    def annotate_vcf(vcf, msa=None):
        header = []
        records = []

        coords = {}

        if msa:
            with open(msa) as msa_fd:
                seq1 = msa_fd.readline().strip()
                seq2 = msa_fd.readline().strip()
                ref = 0  # 1-based coords in vcf
                con = 0

                for a, b in zip(seq1, seq2):
                    if a != "-":
                        ref += 1
                    if b != "-":
                        con += 1

                    coords[ref] = con
                    stats[con] = f"{ref}:{a}-{con}{b}"
        else:
            for i in range(len(self.ref)):
                coords[i] = i

        # TODO: assert 'chromosome' names are the same

        for line in open(vcf):
            line = line.strip()
            if line[0] == "#":
                header.append(line)
                continue

            #        MN908947.3	12781	6	C	T	.	PASS	.	GT	1/1
            _, pos, _, x, y, *r = line.split("\t")
            pos = int(pos)

            if pos in stats:
                records.append((*line.split("\t"), stats[coords[pos]]))
            else:
                print(f"test state {pos} -- check indel logic")
        return records


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
        if False:
            if not test_bias(self.refs_forward, self.refs, threshold=bias_threshold):
                self.log.append(
                    f"Strand bias for reads with consensus alleles; {self.refs_forward} / {self.refs}"
                )
                # position_failed = True

        # amplicon bias
        if self.config.test_amplicon_frs:
            for amplicon, total in self.amplicon_totals.items():
                if (
                    total < self.config.min_depth
                ):  # only evaluate per amplicon frs if there're enough reads
                    continue
                amplicon_frs = self.refs_in_amplicons[amplicon] / total
                if amplicon_frs < self.config.min_frs:
                    self.log.append(
                        f"Per-amplicon FRS failure, amplicon {amplicon}: {self.refs_in_amplicons[amplicon]} / {total}"
                    )
                    if summary:
                        summary["low_amplicon_specific_frs"] += 1
                    position_failed = True
                    break  # to prevent over-counting if both amplicons fail
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
                if len(query) <= q_pos + i:
                    # TODO this condition is being hit when there's soft-clipping at the end of the read (I think)
                    # print(f"WARNING: invalid cigar string. {query}, index {q_pos + i}. cigar: {cigar}", file=sys.stderr)
                    continue
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
