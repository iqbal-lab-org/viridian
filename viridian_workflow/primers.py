import csv
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

    def __str__(self):
        return ", ".join(
            [self.name, str(self.start), str(self.end), str(self.left), str(self.right)]
        )

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
        self,
        name,
        tolerance=5,
        amplicons=None,
        vwf_tsv_file=None,
        tsv_file=None,
        json_file=None,
    ):
        """AmpliconSet supports various membership operations"""
        self.tree = IntervalTree()
        self.name = name
        self.seqs = {}
        if (
            not len(
                [
                    x
                    for x in (amplicons, vwf_tsv_file, tsv_file, json_file)
                    if x is not None
                ]
            )
            == 1
        ):
            raise Exception(
                "Must provide exactly one of amplicons, vwf_tsv_file, tsv_file, json_file"
            )
        if vwf_tsv_file is not None:
            amplicons = AmpliconSet.from_tsv_viridian_workflow_format(
                vwf_tsv_file, tolerance=tolerance
            )
        elif tsv_file is not None:
            amplicons = AmpliconSet.from_tsv(tsv_file, tolerance=tolerance)
        elif json_file is not None:
            amplicons = AmpliconSet.from_json(json_file, tolerance=tolerance)
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

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    @classmethod
    def from_json(cls, fn, tolerance=5):
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

    @classmethod
    def from_tsv_viridian_workflow_format(cs, fn, tolerance=5):
        amplicons = {}
        required_cols = {
            "Amplicon_name",
            "Primer_name",
            "Left_or_right",
            "Sequence",
            "Position",
        }
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
                    amplicons[d["Amplicon_name"]] = Amplicon(d["Amplicon_name"])

                left = d["Left_or_right"].lower() == "left"
                # We assume that primer is always left+forward, or right+reverse
                forward = left
                pos = int(d["Position"])
                primer = Primer(d["Primer_name"], d["Sequence"], left, forward, pos)
                amplicons[d["Amplicon_name"]].add(primer)
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
