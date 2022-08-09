"""
AmpliconSets are collections of Amplicons.

Amplicons may have multiple associated primers.
"""
from __future__ import annotations

from typing import Optional
import csv
from pathlib import Path
from dataclasses import dataclass
from intervaltree import IntervalTree  # type: ignore
from viridian_workflow.utils import Index0, in_range
from viridian_workflow.reads import Fragment


@dataclass(frozen=True)
class Primer:
    """A short sequence corresponding to part of the reference. Immutable."""

    name: str
    seq: str
    left: bool
    forward: bool
    ref_start: Index0
    ref_end: Index0


class Amplicon:
    """A target region of the reference to be amplified by PCR"""

    def __init__(self, name: str, shortname: int = 0):
        self.shortname: int = shortname
        self.name: str = name
        self.start: Index0
        self.end: Index0
        self.left: list[Primer] = []
        self.right: list[Primer] = []
        self.left_primer_region: Optional[tuple[Index0, Index0]] = None
        self.right_primer_region: Optional[tuple[Index0, Index0]] = None
        self.max_length: int = 0

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __str__(self):
        return ", ".join(
            map(str, [self.name, self.start, self.end, self.left, self.right])
        )

    def __hash__(self):
        return hash(self.name)

    def __len__(self):
        """Returns the 'length' of the amplicon, by using the outermost
        coordinates (ie distance from start of leftmost primer to the end
        of rightmost primer)"""
        return self.end - self.start

    def position_in_primer(self, position: Index0) -> bool:
        """Test whether a reference position falls inside any of the primers
        associated with this amplicon
        """

        for primer in self.left + self.right:
            if in_range((primer.ref_start, primer.ref_end), position):
                return True
        return False

    def match_primers(
        self, fragment: Fragment, primer_match_threshold: int = 5
    ) -> tuple[Optional[Primer], Optional[Primer]]:
        """Attempt to match either end of a fragment against the amplicon's primers"""
        p1, p2 = None, None

        min_dist = primer_match_threshold
        # closest leftmost position
        for primer in self.left:
            dist = abs(fragment.ref_start - primer.ref_start)
            if dist < primer_match_threshold:
                if dist <= min_dist:
                    min_dist = dist
                    p1 = primer

        min_dist = primer_match_threshold
        # closest rightmost position
        for primer in self.right:
            dist = abs(primer.ref_end - fragment.ref_end)
            if dist < primer_match_threshold:
                if dist <= min_dist:
                    min_dist = dist
                    p2 = primer

        return p1, p2

    def add(self, primer: Primer):
        """Push a primer into the collection of primers"""
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
    """A set of amplicons which are amplified at the same time"""

    def __init__(
        self,
        name: str,
        amplicons: dict[str, Amplicon],
        tolerance: int = 5,
        shortname: str = None,
        fn: Path = None,
    ):
        """AmpliconSet supports various membership operations"""
        if not shortname:
            # base-54 hash for packing annotations into bam file fields
            # (not used)
            self.shortname = chr(((sum(map(ord, name)) - ord("A")) % 54) + 65)
        self.tree = IntervalTree()
        self.name: str = name
        self.seqs = {}
        self.amplicons = amplicons
        self.amplicon_ids = {}
        if fn:
            self.fn: Path = Path(fn)

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
    def from_json(cls, fn: Path, tolerance=5):
        raise NotImplementedError

    @classmethod
    def from_tsv(cls, fn: Path, name: Optional[str] = None, **kwargs):
        amplicons: dict[str, Amplicon] = {}
        required_cols = {
            "Amplicon_name",
            "Primer_name",
            "Left_or_right",
            "Sequence",
            "Position",
        }
        n = 0
        with open(fn, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            assert reader.fieldnames is not None
            missing_cols_set = required_cols.difference(set(reader.fieldnames))
            if len(missing_cols_set) > 0:
                missing_cols: str = ",".join(sorted(list(missing_cols_set)))
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
                    Index0(pos),
                    # off-by-one hazard: end position is index of last pos
                    Index0(pos + len(d["Sequence"]) - 1),
                )
                amplicons[d["Amplicon_name"]].add(primer)

        name = str(fn) if name is None else name
        return cls(name, amplicons, fn=fn, **kwargs)

    def match(self, fragment: Fragment) -> Optional[Amplicon]:
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
            return list(hits)[0].data

        return None
