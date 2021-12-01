import sys
from collections import namedtuple, defaultdict
from intervaltree import IntervalTree

Primer = namedtuple("Primer", ["name", "seq", "left", "forward", "pos"])


class Amplicon:
    def __init__(self, name):
        self.name = name
        self.start = int.max()
        self.end = 0
        self.left = []
        self.right = []

    def add(self, primer):
        length = len(primer.seq)
        if primer.left:
            self.left.append(primer)
        else:
            self.right.append(primer)

        if left and primer.pos < self.start:
            self.start = primer.pos

        if not left and (primer.pos + length) > self.end:
            self.end = primer.pos + length


class AmpliconSet:
    def __init__(self, amplicons, tolerance):
        """AmpliconSet supports various membership operations
        """
        self.tree = IntervalTree()
        self.name = fn
        self.seqs = {}

        for amplicon in amplicons:
            for primer in amplicon.left:
                self.seqs[primer.seq] = amplicon
            for primer in amplicon.right:
                # note: you may want to reverse complement this
                self.seqs[primer.seq] = amplicon

            # interval containment tolerance
            start = amplicon.start - self.tolerance
            end = amplicon.end + self.tolerance
            self.tree[start:end] = amplicon

        # the internal sequences table allows lookup by primer sequence
        for k, v in sequences.items():
            self.seqs[k[: self.min]] = v

    def from_json(self, fn):
        raise NotImplementedError

    def from_tsv(self, fn, tolerance=5):
        """Import primer set from tsv (QCovid style)
        """

        amplicons = defaultdict(lambda name: Amplicon(name))

        for l in open(fn):
            amplicon_name, name, seq, left, forward, pos = l.strip().split("\t")
            pos = int(pos)
            forward = forward.lower() in ["t", "+", "forward", "true"]
            left = left.lower() in ["left", "true", "t"]

            primer = Primer(name, seq, left, forward, pos)
            amplicons[amplicon_name].add(primer, left)
            # exact matching: also store the reverse complement of the primer

        return AmpliconSet(amplicons, tolerance=tolerance)

    def match(self, start, end):
        """Identify a template's mapped interval based on the start and end
        positions

        returns a set of matching amplicons
        """

        # amplicons which contain the start and end
        hits = self.tree[start].intersection(self.tree[end])

        # amplicons contained by the start and end
        # this should never happen in tiled amplicons
        enveloped = self.tree.envelops(start, end)

        if enveloped:
            return None

        if len(hits) == 0:
            return None
        elif len(hits) > 2:
            # there should not be any more than 2 ambiguous matches under any
            # known primer set. The interval tree can confirm this at the time
            # the primer set is parsed
            raise Exception
        else:
            return hits
