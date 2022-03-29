from collections import namedtuple, defaultdict
from pathlib import Path
import json
import os
import random

from viridian_workflow import utils

import pysam

# "seq" is the read sequence in the direction of the reference genome, ie what
# you get in a BAM file.
# We will enforce that  ref_start < ref_end, and qry_start < qry_end. Then we
# can use is_reverse to resolve the direcrtion of the read.
# qry_end and ref_end one past the position, so slicing/subtracting
# coords follows the python string convention.
Read = namedtuple(
    "Read", ["seq", "ref_start", "ref_end", "qry_start", "qry_end", "is_reverse",],
)


def score(matches, mismatches):
    """Assign winning amplicon set id based on match stats"""
    amplicon_sets = set([*matches.keys(), *mismatches.keys()])

    m = 0
    winner = None
    for amplicon_set in amplicon_sets:
        print(amplicon_set.name, mismatches[amplicon_set], matches[amplicon_set])
        if matches[amplicon_set] >= m:
            winner = amplicon_set
            m = matches[amplicon_set]
    return winner


def match_read_to_amplicons(read, amplicon_sets):
    if read.is_unmapped:
        return None
    matches = {}
    for amplicons in amplicon_sets:
        m = amplicons.match(*read_interval(read))
        if m:
            matches[amplicons.name] = m
    return matches


def match_reads(reads, amplicon_sets):
    """given a stream of reads, yield reads with a set of matched amplicons"""
    for read in reads:
        if read.is_unmapped:
            continue

        matches = match_read_to_amplicons(read, amplicon_sets)
        yield read, matches


def amplicon_set_counts_to_naive_total_counts(scheme_counts):
    counts = defaultdict(int)
    for scheme_tuple, count in scheme_counts.items():
        for scheme in scheme_tuple:
            counts[scheme] += count
    return counts


def amplicon_set_counts_to_json_friendly(scheme_counts):
    dict_out = {}
    for k, v in scheme_counts.items():
        new_key = ";".join(sorted([str(x) for x in k]))
        dict_out[new_key] = v
    return dict_out


class Bam:
    def __init__(self, bam, infile_is_paired=None):
        self.stats = None
        self.infile_is_paired = infile_is_paired
        if not Path(bam).is_file():
            raise Exception(f"bam file {bam} does not exist")
        self.bam = bam

    @staticmethod
    def read_from_pysam(read):
        return Read(
            read.query_sequence,
            read.reference_start,
            read.reference_end,
            read.query_alignment_start,
            read.query_alignment_end,
            read.is_reverse,
        )

    @classmethod
    def from_pe_fastqs(cls, fq1, fq2):
        pass

    @classmethod
    def from_se_fastq(cls, fq):
        pass

    def syncronise_fragments(self):
        reads_by_name = {}

        reads = pysam.AlignmentFile(self.bam, "rb")

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
                self.stats["template_lengths"][
                    abs(read.query_length)
                ] += 1  # TODO: check this
                yield SingleRead(Bam.read_from_pysam(read))

            if not read.is_proper_pair:
                continue

            if read.is_read1:
                self.stats["template_lengths"][abs(read.template_length)] += 1
                reads_by_name[read.query_name] = Bam.read_from_pysam(read)

            elif read.is_read2:
                read1 = reads_by_name[read.query_name]
                yield PairedReads(read1, Bam.read_from_pysam(read))
                del reads_by_name[read.query_name]

    def detect_amplicon_set(self, amplicon_sets):
        """return inferred amplicon set from list
        """
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

        self.stats["match_any_amplicon"] = match_any_amplicon
        self.stats["amplicon_scheme_set_matches"] = matches
        self.stats[
            "amplicon_scheme_simple_counts"
        ] = amplicon_set_counts_to_naive_total_counts(
            self.stats["amplicon_scheme_set_matches"]
        )
        chosen_scheme = score(matches, mismatches)
        if chosen_scheme:
            self.stats["chosen_amplicon_scheme"] = chosen_scheme.name
        else:
            # TODO: decide on behaviour when no appropriate scheme is chosen
            # current policy: abort
            raise Exception("failed to choose amplicon scheme")
        return chosen_scheme

    def stats(self):
        """return pre-computed stats (or compute)
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
        self.amplicons = defaultdict(list)
        self.amplicon_set = amplicon_set
        self.reads_all_paired = None
        self.unmatched_reads = 0

        for fragment in bam.syncronise_fragments():
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
        self.failed_amplicons = set()

        summary = {}

        for amplicon in self.amplicon_set:
            summary[amplicon.name] = {
                "start": amplicon.start + 1,
                "end": amplicon.end + 1,
                "total_mapped_bases": 0,
                "total_depth": 0,
                "sampled_bases": 0,
                "sampled_depth": 0,
                "pass": False,
            }

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
                self.failed_amplicons.add(amplicon)
            else:
                manifest_data[amplicon.name] = outname

        return manifest_data
