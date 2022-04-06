import sys

import mappy as mp
import pysam

from collections import namedtuple, defaultdict

BaseProfile = namedtuple(
    "BaseProfile", ["base", "in_primer", "forward_strand", "amplicon_name",],
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
        self.seq = []
        for r in ref:
            self.seq.append(Stats(ref_base=r))

    def __getitem__(self, pos):
        if pos >= len(self.seq):
            print(f"position too big: {pos} {len(self.seq)}")
            return None
        return self.seq[pos]

    def __setitem__(self, pos, profile):
        raise Exception(f"don't set Pileup positions")

    def update(self, pos, profile):
        self.seq[pos].update(profile)

    def __len__(self):
        return len(self.seq)

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

    def annotate_vcf(self, vcf, msa=None):
        header = []
        records = []

        # 1-based coord translation table
        ref_to_consensus = {}

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

                    ref_to_consensus[ref] = con
        else:
            for i in range(len(self.ref)):
                ref_to_consensus[i] = i

        # TODO: assert 'chromosome' names are the same

        for line in open(vcf):
            line = line.strip()
            if line[0] == "#":
                header.append(line)
                continue

            #        MN908947.3	12781	6	C	T	.	PASS	.	GT	1/1
            _, pos, _, x, y, *r = line.split("\t")
            pos = int(pos)

            cons_coord = ref_to_consensus[pos]
            records.append(
                (
                    *line.split("\t"),
                    str(pos),
                    str(cons_coord),
                    str(self.seq[cons_coord - 1]),
                )
            )

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
        self.alt_bases = defaultdict(int)

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

    def update(self, profile, alt=None):
        # TODO: check if alt
        self.alt_bases[profile.base] += 1

        if profile.in_primer:
            self.alts_in_primer += 1
            return

        self.alts += 1
        # if profile.base == profile.ref_base:
        #    self.alts_matching_refs += 1

        if profile.amplicon_name:
            # when unambiguous amplicon call cannot be made, do not
            # consider (or consider differently)
            self.alts_in_amplicons[profile.amplicon_name] += 1
            self.amplicon_totals[profile.amplicon_name] += 1

        if profile.forward_strand:
            self.alts_forward += 1

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
        alts = " ".join([f"{alt}:{count}" for alt, count in self.alt_bases.items()])
        return f"{self.ref_base}: {alts}"


def parse_cigar(ref, query, alignment):
    """Interpret cigar string and query sequence in reference
    coords from mappy (count, op) or pysam (op, count)

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
                if len(query) == q_pos + 1:
                    # done? TODO: verify that this captures the full query sequence
                    # print(f"WARNING: invalid cigar string. {len(query)}, index {q_pos}. cigar: {cigar}", file=sys.stderr)
                    break

                positions.append((r_pos, query[q_pos]))
                q_pos += 1
                r_pos += 1

        elif op == 1:
            # insertion
            positions.append((r_pos, query[q_pos : q_pos + count + 1]))
            q_pos += count + 1
            r_pos += 1

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
