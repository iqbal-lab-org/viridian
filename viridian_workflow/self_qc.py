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


def mask_sequence(sequence, position_stats):
    sequence = list(sequence)
    qc = {}
    for position, stats in position_stats.items():
        if position >= len(sequence):
            print(
                f"Invalid condition: mapped position {position} greater than consensus length {len(sequence)}",
                file=sys.stderr,
            )
            continue

        if sequence[position] == "N":
            # if a position is already masked by an upstream process skip it
            continue
        elif stats.check_for_failure():
            sequence[position] = "N"
            qc[position] = stats.log
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
    def __init__(self):
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
        self.reference_pos = None
        self.ref_base = None
        self.cons_base = None

    def add_alt(self, profile, alt=None):
        if not self.reference_pos:
            self.reference_pos = profile.reference_pos
        assert self.reference_pos == profile.reference_pos

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
        if self.reference_pos != profile.reference_pos:
            print(
                f"expected to be the same {self.reference_pos} {profile.reference_pos}"
            )
            assert False

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

    def check_for_failure(self, minimum_depth=10, minimum_frs=0.7, bias_threshold=0.3):
        """return whether a position should be masked
        """

        position_failed = False

        if self.total < minimum_depth:
            self.log.append(
                f"Insufficient depth to evaluate consensus; {self.total} < {minimum_depth}. {self.total_reads} including primer regions."
            )
            return True  # position failed

        # test total percentage of bases supporting consensus
        if self.refs / self.total < minimum_frs:
            self.log.append(
                f"Insufficient support of consensus base; {self.refs} / {self.total} < {minimum_frs}. {self.total_reads} including primer regions."
            )
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
            if not test_bias(
                self.refs_in_amplicons[amplicon], total, threshold=bias_threshold
            ):
                self.log.append(
                    f"Amplicon bias in consensus allele calls, amplicon {amplicon}: {self.refs_in_amplicons[amplicon]} / {total}"
                )
                # position_failed = True
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
            q_pos += count
            pass

        elif op == 5:
            # hard clip
            pass

        else:
            raise Exception(f"invalid cigar op {op}")

    return positions


def remap(ref_genome, consensus_fasta, minimap_presets, amplicon_set, tagged_bam):
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

        # print(consensus_seq[alignment.r_st : alignment.r_en], file=sys.stderr)
        # print(r.seq, file=sys.stderr)
        # print(alignment.cigar, file=sys.stderr)

        for (read_pos, base), (ref_pos, _) in zip(cons_alts, ref_alts):
            consensus_position = read_pos + alignment.r_st
            reference_position = ref_pos + r.reference_start

            if consensus_position >= len(consensus_seq):
                print(
                    f"Position {read_pos}+{alignment.r_st} extends beyond length of consensus ({len(consensus_seq)}) {consensus_position}:{reference_position} ({read_pos}: {base})",
                    file=sys.stderr,
                )
                continue

            cons_base = consensus_seq[consensus_position]
            ref_base = reference_seq[reference_position]

            # TODO resolve assumption: if there is an ambiguous amplicon id, in_primer is false
            # in_primer = amplicon.position_in_primer(reference_position)
            in_primer = False

            base_profile = BaseProfile(
                base,
                cons_base,
                ref_base,
                in_primer,
                strand,
                amplicon.name,
                reference_position,
            )

            if consensus_position not in stats:
                stats[consensus_position] = Stats()
            stats[consensus_position].total_reads += 1

            if base != cons_base:
                stats[consensus_position].add_alt(base_profile)
            else:
                stats[consensus_position].add_ref(base_profile)

    print(
        f"reads tagged with more than one amplicon: {multi_amplicons}, with zero: {no_amplicons}. tagged: {tagged}.",
        file=sys.stderr,
    )
    return stats


def mask(fasta, stats, outpath=None, name=None):
    if not outpath:
        output = f"{name}.masked.fasta"
    ref = mp.Aligner(fasta)
    if len(ref.seq_names) != 1:
        Exception(f"Consensus fasta {fasta} has more than one sequence")
    sequence = ref.seq(ref.seq_names[0])
    if not name:
        name = ref.seq_names[0]

    masked, log = mask_sequence(sequence, stats)

    # write masked fasta
    with open(outpath, "w") as maskfd:
        print(f">{name}\n{masked}", file=maskfd, end="")
    return outpath, log


if __name__ == "__main__":
    amplicon_sets = load_amplicon_sets(sys.argv[1])

    shortname = None
    ref_seq = None
    for s in mp.fastx_read(sys.argv[1]):
        ref_seq = s[1]

    amplicon_set = sys.argv[2]

    amplicons = {}
    for a, aset in amplicon_sets.items():
        if aset.name == sys.argv[2]:
            shortname = aset.shortname
            for amplicon in aset.tree:
                amplicon = amplicon.data
                amplicons[amplicon.shortname] = amplicon

    stats = remap(sys.argv[3], amplicon_set, sys.argv[4])

    for p in sorted(stats.keys()):
        if stats[p].total < 5:
            continue
        if str(stats[p]) == "-":
            continue
        print(p, stats[p])
