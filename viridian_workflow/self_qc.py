import sys

import mappy as mp
import pysam

from primers import AmpliconSet, get_tags, load_amplicon_sets

from collections import namedtuple, defaultdict

Stats = namedtuple("Stats", ["alt_in_primer", "ref", "alts", "total"])
BaseProfile = namedtuple("BaseProfile", ["in_primer", "forward_strand", "amplicon"])


def mask_sequence(sequence, position_stats):
    qc = {}
    for position, stats in position_stats.items():
        if stats.check_for_failure():
            sequence[position] = "N"
            qc[position] = stats.log
    return sequence, qc


class Stats:
    def __init__(self):
        self.alts_in_primer = 0
        self.refs_in_primer = 0

        self.alts_in_amplicons = defaultdict(int)
        self.amplicon_totals = defaultdict(int)

        self.alts_forward = 0
        self.refs_forward = 0

        self.total = 0
        self.log = []

    def add_alt(self, profile, alt=None):
        self.alts += 1

        if profile.amplicon:
            # when unambiguous amplicon call cannot be made, do not
            # consider (or consider differently)
            self.alts_in_amplicons[profile.amplicon] += 1
            self.amplicon_totals[profile.amplicon] += 1

        if profile.forward_strand:
            self.alts_forward += 1

        if profile.in_primer:
            self.alts_in_primer += 1

        self.total += 1

    def add_ref(self, profile):
        if profile.amplicon:
            self.amplicon_totals[profile.amplicon] += 1

        if profile.forward_strand:
            self.refs_forward += 1

        if profile.in_primer:
            self.refs_in_primer += 1

        self.total += 1

    def check_for_failure(self):
        """return whether a position should be masked
        """

        position_failed = False

        # look for overrepresentation of alt alleles in regions covered
        # by primer sequences
        if self.alts_in_primer / self.alts < 0.3:
            self.log.append("alts in primer")
            position_failed = True

        # strand bias in alt calls
        if self.alts_forward / self.alts < 0.3:
            self.log.append("strand bias in alt calls")
            position_failed = True

        # amplicon bias
        for amplicon, total in amplicon_totals:
            if self.alts_in_amplicons[amplicon] / self.alts < 0.3:
                self.log.append(f"amplicon bias in alt cals, amplicon {amplicon}")
                position_failed = True

        return position_failed

    def __str__(self):
        f = []
        if len(self.amps_total) > 1:
            return "-".join([f"{k}:{v}" for k, v in self.amps_total.items()])
        if self.alts / self.total > 0.2:
            return f"{self.alts}/{self.total}"
        return "-"


def cigar_to_alts(ref, query, cigar):
    """Interpret cigar string and query sequence in reference
    coords
    """
    positions = []
    q_pos = 0
    for op, count in cigar:
        if op == 0:
            # match/mismatch
            for i in range(q_pos, q_pos + count):
                positions.append((q_pos + i, query[i]))
            q_pos += count

        elif op == 1:
            # insertion
            positions.append((q_pos, query[q_pos : q_pos + count]))
            q_pos += count

        elif op == 2:
            # deletion
            pass

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


def remap(ref, minimap_presets, amplicon_set, tagged_bam):
    stats = {}
    ref = mp.Aligner(sys.argv[1], preset=sys.argv[3])
    for r in pysam.AlignmentFile(sys.argv[4]):
        a = ref.map(r.seq)
        alignment = None
        for x in a:
            if x.is_primary:
                alignment = x

        if not alignment:
            continue

        tags = get_tags(None, r)[amplicon_set.shortname]

        amplicon = None
        if len(tags) == 1:
            amplicon = amplicons[tags[0]]

        strand = False  # strand is forward
        if alignment.strand == 0:
            # should be error
            raise Exception()
        elif alignment.strand == 1:
            strand = True

        base_profile = BaseProfile(in_primer, strand, amplicon)
        alts = cigar_to_alts(ref_seq[alignment.r_st : alignment.r_en], r.seq, r.cigar)

        for read_pos, base in alts:
            position = read_pos + alignment.r_st

            if position not in stats:
                stats[position] = Stats()
            if base != ref_seq[position]:
                stats[position].add_alt(base_profile)
            else:
                stats[position].add_ref(base_profile)
    return stats


def mask(fasta, stats, name=None, prefix=None):
    seq = ""
    masked, log = mask_sequence(seq, stats)

    # write masked fasta
    with open(outpath, "w") as maskfd:
        print(">{name}\n{seq}", file=maskfd, end="")
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
