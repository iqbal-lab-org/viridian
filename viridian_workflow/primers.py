import csv
from collections import namedtuple, defaultdict
from intervaltree import IntervalTree
from viridian_workflow.readstore import PairedReads, SingleRead, ReadStore

Primer = namedtuple(
    "Primer", ["name", "seq", "left", "forward", "ref_start", "ref_end"]
)


def in_range(interval, position):
    start, end = interval
    return position < end and position > start


class Amplicon:
    def __init__(self, name, shortname=0):
        self.shortname = shortname
        self.name = name
        self.start = None
        self.end = None
        self.left = []
        self.right = []
        self.left_primer_region = None
        self.right_primer_region = None
        self.max_length = 0

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __str__(self):
        return ", ".join(
            [self.name, str(self.start), str(self.end), str(self.left), str(self.right)]
        )

    def __hash__(self):
        return hash(self.name)

    def __len__(self):
        """Returns the 'length' of the amplicon, by using the outermost
        coordinates (ie distance from start of leftmost primer to the end
        of rightmost primer)"""
        return self.end - self.start

    def position_in_primer(self, position):
        """Test whether a reference position falls inside the primer"""
        return in_range(self.left_primer_region, position) or in_range(
            self.right_primer_region, position
        )

    def match_primers(self, fragment):
        """Attempt to match either end of a fragment against the amplicon's primers
        """
        p1, p2 = None, None

        min_dist = config.primer_match_threshold
        # closest leftmost position
        for primer in self.left:
            dist = fragment.ref_start - primer.ref_start
            if dist >= 0 and dist < config.primer_match_threshold:
                if dist <= min_dist:
                    min_dist = dist
                    p1 = primer

        min_dist = config.primer_match_threshold
        # closest leftmost position
        for primer in self.right:
            dist = primer.ref_end - fragment.ref_end
            if dist >= 0 and dist < config.primer_match_threshold:
                if dist <= min_dist:
                    min_dist = dist
                    p2 = primer

        return p1, p2

    def add(self, primer):
        if len(primer.seq) > self.max_length:
            self.max_length = len(primer.seq)

        if primer.left:
            self.left.append(primer)
            if not self.left_primer_region:
                self.left_primer_region = (primer.ref_start, primer.ref_end)
                self.start = primer.ref_start
            else:
                left_start, left_end = self.left_primer_region
                self.left_primer_region = (
                    min(primer.ref_start, left_start),
                    max(primer.ref_end, left_end),
                )
                self.start = min(self.start, primer.ref_start)
        else:
            self.right.append(primer)
            if not self.right_primer_region:
                self.right_primer_region = (primer.ref_start, primer.ref_end)
                self.end = primer.ref_end
            else:
                right_start, right_end = self.right_primer_region
                self.right_primer_region = (
                    min(primer.ref_start, right_start),
                    max(primer.ref_end, right_end),
                )
                self.end = max(self.end, primer.ref_end)


class AmpliconSet:
    def __init__(
        self, name, amplicons, tolerance=5, shortname=None,
    ):
        """AmpliconSet supports various membership operations"""
        if not shortname:
            # base-54 hash
            self.shortname = chr(((sum(map(ord, name)) - ord("A")) % 54) + 65)
        self.tree = IntervalTree()
        self.name = name
        self.seqs = {}
        self.amplicons = amplicons
        self.amplicon_ids = {}

        primer_lengths = set()
        sequences = {}
        for amplicon_name in amplicons:
            amplicon = amplicons[amplicon_name]
            self.amplicon_ids[amplicon.shortname] = amplicon

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

    def __hash__(self):
        return hash(self.name)

    def __iter__(self):
        for amplicon in self.amplicons.values():
            yield amplicon

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
                primer = Primer(
                    d["Primer_name"],
                    d["Sequence"],
                    left,
                    forward,
                    pos,
                    pos + len(d["Sequence"]),
                )
                amplicons[d["Amplicon_name"]].add(primer)

        name = fn if not name else name
        return cls(name, amplicons, **kwargs)

    def score(self, readstore):
        for fragment in readstore:
            self.match(fragment)

    def match(self, fragment):
        """Identify a template's mapped interval based on the start and end
        positions

        The return value is None if the fragment does not belong to the
        AmpliconSet or if the fragment matches more than one Amplicon.
        """

        # amplicons which contain the start and end
        hits = self.tree[fragment.ref_start].intersection(self.tree[fragment.ref_end])

        # amplicons contained by the start and end
        # this should never happen in tiled amplicons and may be a strong
        # signal for negatively identifying an amplicon set
        enveloped = self.tree.envelop(fragment.ref_start, fragment.ref_end)

        if enveloped:
            return None

        if len(hits) == 1:
            return hits[0].data

        return None
