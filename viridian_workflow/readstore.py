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

class Bam:
    def __init__(self, bam, infile_is_paired=None):
        self.stats = None
        self.infile_is_paired = infile_is_paired
        assert os.fileexists(bam)
        self.bam = bam

    @classmethod
    def from_pe_fastqs(cls, fq1, fq2):
        pass

    def from_se_fastq(cls, fq):
        pass

    def syncronise_fragments(self):
        reads_by_name = {}

        reads = pysam
        for read in reads:
            if self.infile_is_paired is None:
                self.infile_is_paired = read.is_paired
            else:
                if self.infile_is_paired != read.is_paired:
                    raise Exception("Mix of paired and unpaired reads.")

            if read.is_secondary or read.is_supplementary:
                continue

            self.stats["total_reads"] += 1
            if read.is_read1:
                self.stats["reads1"] += 1
            elif read.is_read2:
                self.stats["reads2"] += 1

            if read.is_unmapped:
                continue

            self.stats["read_lengths"][read.query_length] += 1
            self.stats["mapped"] += 1

            if not read.is_paired:
                self.stats["unpaired_reads"] += 1
                self.stats["template_lengths"][abs(read.query_length)] += 1  # TODO: check this
                yield SingleRead(read_from_pysam(read))

            if not read.is_proper_pair:
                continue

            if read.is_read1:
                stats["template_lengths"][abs(read.template_length)] += 1
                reads_by_name[read.query_name] = read_from_pysam(read)

            elif read.is_read2:
                read1 = reads_by_name[read.query_name]
                yield PairedReads(read1, read_from_pysam(read))
                del reads_by_name[read.query_name]



    def detect_amplicon_set(self, amplicon_sets):
        """return inferred amplicon set from list
        """
        #aln_file_in = pysam.AlignmentFile(self.bam, "rb")

        match_any_amplicon = 0

        self.stats = {
            "unpaired_reads": 0,
            "reads1": 0,
            "reads2": 0,
            "total_reads": 0,
            "mapped": 0,
            "read_lengths": defaultdict(int),
            "template_lengths": defaultdict(int),
        }

        mismatches = defaultdict(int)
        matches = defaultdict(int)

        for fragment in self.syncronise_fragments():
            for amplicon_set in amplicon_sets:
                hit = amplicon_set.match(fragment)
                if hit:
                    matches[amplicon_set] += 1
                else:
                    mismatches[amplicon_set] += 1

        aln_file_in.close()

        self.stats["match_any_amplicon"] = match_any_amplicon
        self.stats["amplicon_scheme_set_matches"] = matches
        self.stats["amplicon_scheme_simple_counts"] = amplicon_set_counts_to_naive_total_counts(
            self.stats["amplicon_scheme_set_matches"]
        )
        chosen_scheme = score(matches, mismatches)
        self.stats["chosen_amplicon_scheme"] = chosen_scheme.name
        return chosen_scheme


    def stats(self):
        """return pre-computed stats (or compute)
        """
        pass

    def readstore(self, amplicon_set):
        """construct readstore object
        """
        pass

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

        for fragment in detect_primers.syncronise_fragments(bam):
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
