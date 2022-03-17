from collections import namedtuple
import json
import os
import random

from viridian_workflow import utils

# "seq" is the read sequence in the direction of the reference genome, ie what
# you get in a BAM file.
# We will enforce that  ref_start < ref_end, and qry_start < qry_end. Then we
# can use is_reverse to resolve the direcrtion of the read.
# qry_end and ref_end one past the position, so slicing/subtracting
# coords follows the python string convention.
Read = namedtuple(
    "Read", ["seq", "ref_start", "ref_end", "qry_start", "qry_end", "is_reverse",],
)


def read_from_pysam(read):
    return Read(
        read.query_sequence,
        read.reference_start,
        read.reference_end,
        read.query_alignment_start,
        read.query_alignment_end,
        read.is_reverse,
    )


class Fragment:
    def __init__(self, reads):
        """fragment ref bounds ignore softclipping
        """
        self.ref_start = None
        self.ref_end = None
        self.reads = reads

    def total_mapped_bases(self):
        return sum([r.qry_end - r.qry_start for r in self.reads])


class PairedReads(Fragment):
    def __init__(self, read1, read2):
        super().__init__([read1, read2])
        (self.ref_start, self.ref_end) = (
            (read1.ref_start, read2.ref_end)
            if read1.ref_start < read2.ref_start
            else (read2.ref_start, read1.ref_end)
        )


class SingleRead(Fragment):
    def __init__(self, read):
        super().__init__([read])
        self.ref_start, self.ref_end = read.ref_start, read.ref_end


class ReadStore:
    def __init__(self, amplicon_set, bam):
        self.amplicons = {}
        self.amplicon_set = amplicon_set
        self.reads_all_paired = None
        self.unmatched_reads = 0

        for fragment in syncronise_fragments(bam):
            self.push_fragment(fragment)

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

    def push_fragment(self, fragment):
        amplicon = self.amplicon_set.match(fragment)
        if amplicon:
            self.amplicons[amplicon].append(fragment)
        else:
            self.unmatched_reads += 1

    @staticmethod
    def sample_paired_reads(fragments, outfile, target_bases):
        if len(fragments) == 0:
            return 0
        bases_out = 0
        with open(outfile, "w") as f:
            for i, fragment in enumerate(fragments):
                for j, read in enumerate(fragment.reads):
                    print(
                        f">{i}.{j}",
                        utils.revcomp(read.seq) if read.is_reverse else read.seq,
                        sep="\n",
                        file=f,
                    )
                bases_out += fragment.total_mapped_bases()
                if bases_out >= target_bases:
                    break

        return bases_out

    @staticmethod
    def sample_unpaired_reads(fragments, outfile, target_bases):
        rev_indexes = []
        fwd_indexes = []
        for i, fragment in enumerate(fragments):
            if fragment.reads[0].is_reverse:
                rev_indexes.append(i)
            else:
                fwd_indexes.append(i)

        if len(fwd_indexes) == 0 or len(rev_indexes) == 0:
            return 0

        bases_out = 0
        with open(outfile, "w") as f:
            for i, (fwd_i, rev_i) in enumerate(zip(fwd_indexes, rev_indexes)):
                fwd_frag = fragments[fwd_i]
                rev_frag = fragments[rev_i]
                print(
                    f">f{i}",
                    fwd_frag.reads[0].seq,
                    f">r{i}",
                    utils.revcomp(rev_frag.reads[0].seq),
                    sep="\n",
                    file=f,
                )
                bases_out += (
                    fwd_frag.total_mapped_bases() + rev_frag.total_mapped_bases()
                )
                if bases_out >= target_bases:
                    break

        return bases_out

    def make_reads_dir_for_viridian(self, outdir, target_depth):
        """Makes a directory of reads for each amplicon, in the format required
        by `viridian assemble --reads_per_amp_dir`. Returns a set of amplicon
        names that should be failed because they had no reads"""
        random.seed(42)
        os.mkdir(outdir)
        manifest_data = {}
        failed_amplicons = set()

        for amplicon in self.amplicon_set:
            outname = f"{len(manifest_data)}.fa"
            outfile = os.path.join(outdir, outname)
            fragments = self[amplicon]
            random.shuffle(fragments)
            target_bases = target_depth * len(amplicon)
            if self.reads_all_paired:
                bases_out = self.sample_paired_reads(fragments, outfile, target_bases)
            else:
                bases_out = self.sample_unpaired_reads(fragments, outfile, target_bases)
            if bases_out == 0:
                failed_amplicons.add(amplicon)
            else:
                manifest_data[amplicon.name] = outname

        with open(os.path.join(outdir, "manifest.json"), "w") as f:
            json.dump(manifest_data, f, indent=2)

        return failed_amplicons
