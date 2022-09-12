"""
Read objects

A fragment with either a pair of reads or a single read. 0-based indexing
is always assumed unless otherwise indicated
"""
from __future__ import annotations

from dataclasses import dataclass
from viridian_workflow.utils import Index0


@dataclass(frozen=True)
class Read:
    """Reads are sequences of dna that have been mapped to a reference
    coordinates (ref_start, ref_end)

    "seq" is the read sequence in the direction of the reference genome, ie what
    you get in a BAM file.
    We will enforce that  ref_start < ref_end, and qry_start < qry_end. Then we
    can use is_reverse to resolve the direction of the read.
    qry_end and ref_end one past the position, so slicing and subtracting
    coords follow the python string convention.
    """

    seq: str
    ref_start: Index0
    ref_end: Index0
    qry_start: Index0
    qry_end: Index0
    is_reverse: bool


class Fragment:
    """A fragment can be either two reads or a single read depending
    on the sequencing platform
    """

    def __init__(self, reads: list[Read]):
        """fragment ref bounds ignore softclipping"""
        self.reads: list[Read] = reads

        self.ref_start: Index0
        self.ref_end: Index0

    def total_mapped_bases(self) -> int:
        """The total number of bases that were sequenced for this fragment.
        This counts positions that were sequenced twice in overlapping
        fragments.
        """
        return sum(r.qry_end - r.qry_start for r in self.reads)


class PairedReads(Fragment):
    """A pair of reads are sequenced from the same template"""

    def __init__(self, read1: Read, read2: Read):
        super().__init__([read1, read2])
        (ref_start, ref_end) = (
            (read1.ref_start, read2.ref_end)
            if read1.ref_start < read2.ref_start
            else (read2.ref_start, read1.ref_end)
        )
        self.ref_start: Index0 = Index0(ref_start)
        self.ref_end: Index0 = Index0(ref_end)

        if read1.is_reverse and not read2.is_reverse:
            strand = False
        elif not read1.is_reverse and read2.is_reverse:
            strand = True
        else:
            raise Exception("Read pair is in invalid orientation F1F2/R1R2")
        self.strand: bool = strand


class SingleRead(Fragment):
    """Nanopore sequences a single read per template. The direction may still
    be either forward or reverse wrt the reference.
    """

    def __init__(self, read: Read):
        super().__init__([read])
        self.ref_start: Index0 = Index0(read.ref_start)
        self.ref_end: Index0 = Index0(read.ref_end)
        self.strand: bool = not read.is_reverse
