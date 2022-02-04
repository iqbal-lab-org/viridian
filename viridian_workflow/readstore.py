from collections import namedtuple

Read = namedtuple("Read", ["ref start", "primer", "end", "is_reverse"])


class Fragment:
    def __init__(self, read):
        self.ref_start = None
        self.ref_end = None
        self.reads = []

    @classmethod
    def from_pair(cls, r1, r2):
        return


class PairedReads(Fragment):
    def __init__(self, read1, read2):
        self.reads = [read1, read2]


class Single:
    def __init__(self, read):
        self.reads = [read]


class ReadStore:
    def __init__(self, name, shortname=0):
        self.amplicons = {}

    def __eq__(self, other):
        pass

    def __str__(self):
        pass

    def __iter__(self):
        pass

    def fetch(self, start=0, end=None):
        pass

    def fetch_amplicon(self, amplicon):
        for fragment in self.amplcions[amplicon]:
            for read in fragment.reads:
                yield read

    @classmethod
    def from_bam(cls, amplicon_set, bam):
        # construct read namedtuples
        pass
