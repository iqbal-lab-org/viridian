import sys
from collections import namedtuple, defaultdict
from intervaltree import IntervalTree

Primer = namedtuple("Primer", ["name", "seq", "left", "forward", "pos"])


class Amplicon:
    def __init__(self, name):
        self.name = name
        self.start = None
        self.end = None
        self.left = []
        self.right = []

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def add(self, primer):
        length = len(primer.seq)
        if primer.left:
            self.left.append(primer)
        else:
            self.right.append(primer)

        if self.start is None and self.end is None:
            self.start = primer.pos
            self.end = primer.pos + length

        if primer.left and primer.pos < self.start:
            self.start = primer.pos

        if not primer.left and (primer.pos + length) > self.end:
            self.end = primer.pos + length


class AmpliconSet:
    def __init__(self, name, tolerance=5, amplicons=None, tsv_file=None, json_file=None):
        """AmpliconSet supports various membership operations"""
        self.tree = IntervalTree()
        self.name = name
        self.seqs = {}
        if not len([x for x in (amplicons, tsv_file, json_file) if x is not None]) == 1:
            raise Exception(
                "Must provide exactly one of amplicons, tsv_file, json_file"
            )
        if tsv_file is not None:
            amplicons = AmpliconSet.from_tsv(tsv_file)
        elif json_file is not None:
            amplicons = AmpliconSet.from_json(json_file)
        assert amplicons is not None

        primer_lengths = set()
        sequences = {}
        for amplicon_name in amplicons:
            amplicon = amplicons[amplicon_name]

            for primer in amplicon.left:
                sequences[primer.seq] = amplicon
                primer_lengths.add(len(primer.seq))
            for primer in amplicon.right:
                # note: you may want to reverse complement this
                sequences[primer.seq] = amplicon
                primer_lengths.add(len(primer.seq))

            # interval containment tolerance
            start = amplicon.start - tolerance
            end = amplicon.end + tolerance
            self.tree[start:end] = amplicon

        self.min_primer_length = min(primer_lengths)
        # the internal sequences table allows lookup by primer sequence
        for k, v in sequences.items():
            self.seqs[k[: self.min_primer_length]] = v

    @classmethod
    def from_json(cls, fn):
        raise NotImplementedError

    @classmethod
    def from_tsv(cls, fn, tolerance=5):
        """Import primer set from tsv (QCovid style)"""
        amplicons = {}

        for l in open(fn):
            amplicon_name, name, seq, left, forward, pos = l.strip().split("\t")
            pos = int(pos)
            forward = forward.lower() in ["t", "+", "forward", "true"]
            left = left.lower() in ["left", "true", "t"]

            if amplicon_name not in amplicons:
                amplicons[amplicon_name] = Amplicon(amplicon_name)

            primer = Primer(name, seq, left, forward, pos)
            amplicons[amplicon_name].add(primer)
            # exact matching: also store the reverse complement of the primer

        return amplicons

    def match(self, start, end):
        """Identify a template's mapped interval based on the start and end
        positions

        returns a set of matching amplicons
        """

        # amplicons which contain the start and end
        hits = self.tree[start].intersection(self.tree[end])

        # amplicons contained by the start and end
        # this should never happen in tiled amplicons
        enveloped = self.tree.envelop(start, end)

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
