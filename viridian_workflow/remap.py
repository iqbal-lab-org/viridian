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
            continue
        amplicon = amplicons[tags[0]]
        for x in a:
            if not x.is_primary:
                continue

            if x.strand != 1:
                continue

            if x.q_en - x.q_st != x.r_en - x.r_st:
                continue

            for ref_base, q, p in zip(
                ref_seq[x.r_st : x.r_en], r.seq, range(x.r_st, x.r_en)
            ):
                if p not in stats:
                    stats[p] = Stats(0, ref_base, 0, 0)
                stats[p]._replace(total=stats[p].total + 1)
                if ref_base != q:
                    if p < amplicon.start + 25 or p > amplicon.end - 25:
                        stats[p]._replace(alt_in_primer=stats[p].alt_in_primer + 1)
                    else:
                        stats[p]._replace(alts=stats[p].alts + 1)
    #            print("\t", x.r_st, x.r_en, tags, in_primer)

    for p in stats:
        print(stats[p])
