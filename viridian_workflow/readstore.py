from collections import namedtuple
import random

# "seq" is the read sequence in the direction of the reference genome, ie what
# you get in a BAM file.
# We will enforce that  ref_start < ref_end, and qry_start < qry_end. Then we
# can use is_reverse to resolve the direcrtion of the read.
# qry_end and ref_end one past the position, so slicing/subtracting
# coords works easily.
Read = namedtuple(
    "Read",
    [
        "seq",
        "ref_start",
        "ref_end",
        "qry_start",
        "qry_end",
        "is_reverse",
    ],
)


class Fragment:
    def __init__(self, reads):
        self.ref_start = None
        self.ref_end = None
        self.reads = reads

    def total_mapped_bases(self):
        return sum([r.qry_end - r.qry_start for r in self.reads])


class PairedReads(Fragment):
    def __init__(self, read1, read2):
        super().__init__([read1, read2])


class SingleRead:
    def __init__(self, read):
        super().__init__([read])


class ReadStore:
    def __init__(self, name, amplicon_set, shortname=0):
        self.amplicons = {}
        self.amplicon_set = amplicon_set
        self.reads_all_paired = None

    def __eq__(self, other):
        pass

    def __str__(self):
        pass

    def __iter__(self):
        pass

    def __getitem__(self, amplicon):
        """Given an amplicon, returns list of Fragments"""
        return self.amplicons[amplicon]

    def fetch(self, start=0, end=None):
        pass

    @classmethod
    def from_bam(cls, amplicon_set, bam):
        # construct read namedtuples
        # note - can set this from bam contents: self.reads_all_paired
        pass

    @staticmethod
    def sample_paired_reads(fragments, outfile, target_depth):
        # FIXME
        pass

    @staticmethod
    def sample_unpaired_reads(fragments, outfile, target_depth):
        # FIXME
        pass

    def make_dir_of_reads_for_viridian(self, outdir, target_depth):
        for amplicon in self.amplicon_set:
            fragments = self[amplicon]
            random.shuffle(fragments)
            # FIXME - write sampled reads to file
