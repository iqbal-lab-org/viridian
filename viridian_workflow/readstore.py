from collections import namedtuple

Read = namedtuple("Read", [])


class Fragment:
    def __init__(self, read):
        self.ref_start = None
        self.ref_end = None

    @classmethod
    def from_pair(cls, r1, r2):
        return


class ReadStore:
    def __init__(self, name, shortname=0):
        pass

    def __eq__(self, other):
        pass

    def __str__(self):
        pass

    def __iter__(self):
        pass

    def fetch(self, start=0, end=None):
        pass

    @classmethod
    def from_bam(cls, amplicon_set, bam):
        pass
