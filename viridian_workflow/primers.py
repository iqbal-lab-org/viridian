import csv
import sys
from collections import namedtuple, defaultdict
from intervaltree import IntervalTree

Primer = namedtuple("Primer", ["name", "seq", "left", "forward", "pos"])


def load_amplicon_schemes(amplicon_tsvs):
    return dict(
        [
            (s, AmpliconSet.from_tsv(tsv, shortname=s))
            for tsv, s in zip(amplicon_tsvs, ["a", "b", "c"])
        ]
    )


def set_tags(amplicon_sets, read, matches):
    tags = []
    for amplicon_set in amplicon_sets:
        if amplicon_set.name not in matches:
            continue
        shortname = amplicon_set.shortname
        for amplicon in matches[amplicon_set.name]:
            tags.append((f"Z{shortname}", amplicon.shortname, "i"))
    read.set_tags(tags)
    return read


def get_tags(read, shortname):
    matches = []
    for tag, value in read.get_tags():
        if tag[0] == "Z" and tag[1] == shortname:
            matches.append(value)
    return matches


class Amplicon:
    def __init__(self, name, shortname=0):
        self.shortname = shortname
        self.name = name
        self.start = None
        self.end = None
        self.left = []
        self.right = []

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __str__(self):
        return ", ".join(
            [self.name, str(self.start), str(self.end), str(self.left), str(self.right)]
        )

    def position_in_primer(self, position):
        """Test wheter a reference position falls inside the primer
        """
        return position > self.start and position < self.end

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
    def __init__(
        self, name, amplicons, tolerance=5, shortname=None,
    ):
        """AmpliconSet supports various membership operations"""
        if not shortname and name:
            self.shortname = chr((sum(map(ord, name)) - ord("A")) % 54)
        self.tree = IntervalTree()
        self.name = name
        self.seqs = {}
        self.amplicons = amplicons

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

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    @classmethod
    def from_json(cls, fn, tolerance=5):
        raise NotImplementedError

    @classmethod
    def from_tsv(cls, fn, name=None, **kwargs):
        amplicons = {}
        required_cols = {
            "Amplicon_name",
            "Primer_name",
            "Left_or_right",
            "Sequence",
            "Position",
        }
        n = 0
        with open(fn) as f:
            reader = csv.DictReader(f, delimiter="\t")
            missing_cols = required_cols.difference(set(reader.fieldnames))
            if len(missing_cols) > 0:
                missing_cols = ",".join(sorted(list(missing_cols)))
                raise Exception(
                    f"Amplicon scheme TSV missing these columns: {missing_cols}. Got these columns: {reader.fieldnames}"
                )

            for d in reader:
                if d["Amplicon_name"] not in amplicons:
                    amplicons[d["Amplicon_name"]] = Amplicon(
                        d["Amplicon_name"], shortname=n
                    )
                    n += 1

                left = d["Left_or_right"].lower() == "left"
                # We assume that primer is always left+forward, or right+reverse
                forward = left
                pos = int(d["Position"])
                primer = Primer(d["Primer_name"], d["Sequence"], left, forward, pos)
                amplicons[d["Amplicon_name"]].add(primer)

        name = fn if not "name" in kwargs else kwargs["name"]
        return cls(name, amplicons, **kwargs)

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
            return [hit.data for hit in hits]
