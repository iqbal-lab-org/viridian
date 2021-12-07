import sys

import mappy as mp
import pysam

from primers import AmpliconSet, get_tags

from collections import namedtuple, defaultdict

ref = mp.Aligner(sys.argv[1], preset="sr")
ref_seq = None
for s in mp.fastx_read(sys.argv[1]):
    ref_seq = s[1]

amplicon_set = sys.argv[2]

Stats = namedtuple("Stats", ["alt_in_primer", "ref", "alts", "total"])


class Stats:
    def __init__(self):
        self.alt_in_primer = 0
        self.alts = 0
        self.amps = defaultdict(int)
        self.amps_total = defaultdict(int)
        self.total = 0

    def add_alt(self, alt, amplicon):
        self.alts += 1
        self.amps[amplicon] += 1
        self.amps_total[amplicon] += 1
        self.total += 1

    def add_ref(self, amplicon):
        self.total += 1
        self.amps_total[amplicon] += 1

    def score(self):
        """return whether a position should be masked
        """
        return True

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
            print(f"invalid cigar op {op}")

    return positions


def read_error_profile(read):
    pass


if __name__ == "__main__":
    amplicon_tsvs = [
        "../data/midnight-1200.qcovid.tsv",
        "../data/artic-v3.qcovid.tsv",
        "../data/artic-v4.qcovid.tsv",
    ]
    amplicon_sets = dict(
        [
            (s, AmpliconSet(tsv, tsv_file=tsv, shortname=s))
            for tsv, s in zip(amplicon_tsvs, ["a", "b", "c"])
        ]
    )
    sn = None
    amplicons = {}
    for a, aset in amplicon_sets.items():
        if aset.name == sys.argv[2]:
            sn = aset.shortname
            for amplicon in aset.tree:
                amplicon = amplicon.data
                amplicons[amplicon.shortname] = amplicon

    stats = {}
    for r in pysam.AlignmentFile(sys.argv[3]):
        a = ref.map(r.seq)
        tags = get_tags(None, r)[sn]
        if len(tags) < 1:
            amplicon = "none"
        else:
            amplicon = amplicons[tags[0]].name

        alignment = None
        for x in a:
            if x.is_primary:
                alignment = x

        if not alignment:
            continue

        #            if alignment.strand != 1:

        alts = cigar_to_alts(ref_seq[alignment.r_st : alignment.r_en], r.seq, r.cigar)
        for read_pos, base in alts:
            position = read_pos + alignment.r_st

            if position not in stats:
                stats[position] = Stats()
            if base != ref_seq[position]:
                stats[position].add_alt(None, amplicon)
            else:
                stats[position].add_ref(amplicon)

    for p in sorted(stats.keys()):
        if stats[p].total < 5:
            continue
        if str(stats[p]) == "-":
            continue
        print(p, stats[p])
