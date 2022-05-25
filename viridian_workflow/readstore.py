from collections import namedtuple, defaultdict
import sys
from pathlib import Path
import os
import random

from viridian_workflow import utils, self_qc

import pysam
import mappy as mp

# "seq" is the read sequence in the direction of the reference genome, ie what
# you get in a BAM file.
# We will enforce that  ref_start < ref_end, and qry_start < qry_end. Then we
# can use is_reverse to resolve the direcrtion of the read.
# qry_end and ref_end one past the position, so slicing/subtracting
# coords follows the python string convention.
Read = namedtuple(
    "Read", ["seq", "ref_start", "ref_end", "qry_start", "qry_end", "is_reverse"],
)


def in_range(interval, position):
    start, end = interval
    return position < end and position > start


def score(matches, mismatches):
    """Assign winning amplicon set id based on match stats"""
    amplicon_sets = set([*matches.keys(), *mismatches.keys()])

    m = 0
    winner = None
    for amplicon_set in amplicon_sets:
        total = matches[amplicon_set] + mismatches[amplicon_set]
        if mismatches[amplicon_set] / total >= 0.1:
            # if more than 10% of reads break the amplicon boundaries
            # disqualify this amplicon set
            print(
                amplicon_set.name,
                mismatches[amplicon_set],
                matches[amplicon_set],
                "[disqualified]",
            )
            continue
        else:
            print(amplicon_set.name, mismatches[amplicon_set], matches[amplicon_set])
        if matches[amplicon_set] >= m:
            winner = amplicon_set
            m = matches[amplicon_set]
    return winner


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
        improper_pairs = 0

        self.stats = {
            "unpaired_reads": 0,
            "reads1": 0,
            "reads2": 0,
            "total_reads": 0,
            "mapped": 0,
            "read_lengths": defaultdict(int),
            "template_lengths": defaultdict(int),
            "match_no_amplicon_sets": 0,
        }

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
                improper_pairs += 1
                continue

            if read.is_read1:
                self.stats["template_lengths"][abs(read.template_length)] += 1
                reads_by_name[read.query_name] = Bam.read_from_pysam(read)

            elif read.is_read2:
                if read.query_name not in reads_by_name:
                    raise Exception("Bam file is not sorted by name")
                read1 = reads_by_name[read.query_name]
                yield PairedReads(read1, Bam.read_from_pysam(read))
                del reads_by_name[read.query_name]
        print(f"{improper_pairs} improper pairs", file=sys.stderr)

    def detect_amplicon_set(self, amplicon_sets):
        """return inferred amplicon set from list
        """

        mismatches = defaultdict(int)
        matches = defaultdict(int)

        for fragment in self.syncronise_fragments():
            match_any = False
            for amplicon_set in amplicon_sets:
                hit = amplicon_set.match(fragment)
                if hit:
                    match_any = True
                    matches[amplicon_set] += 1
                else:
                    mismatches[amplicon_set] += 1
            if not match_any:
                self.stats["match_no_amplicon_sets"] += 1

        #        self.stats["match_any_amplicon"] = match_any_amplicon
        self.stats["amplicon_scheme_set_matches"] = {}
        for match in matches:
            self.stats["amplicon_scheme_set_matches"][match.name] = matches[match]

        #        self.stats[
        #            "amplicon_scheme_simple_counts"
        #        ] = amplicon_set_counts_to_naive_total_counts(
        #            self.stats["amplicon_scheme_set_matches"]
        #        )
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
        self.strand = None

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
        if read1.is_reverse and not read2.is_reverse:
            self.strand = False
        elif not read1.is_reverse and read2.is_reverse:
            self.strand = True
        else:
            raise Exception(f"Read pair is in invalid orientation F1F2/R1R2")


class SingleRead(Fragment):
    def __init__(self, read):
        super().__init__([read])
        self.ref_start, self.ref_end = read.ref_start, read.ref_end
        self.strand = not read.is_reverse


class ReadStore:
    def __init__(self, amplicon_set, bam, target_depth=1000):
        self.amplicons = defaultdict(list)
        self.reads_per_amplicon = defaultdict(int)
        self.amplicon_set = amplicon_set
        self.reads_all_paired = bam.infile_is_paired
        self.unmatched_reads = 0

        self.target_depth = target_depth
        # TODO find a home for this magic number
        self.cylon_target_depth_factor = 200

        self.start_pos = None
        self.end_pos = None
        self.amplicon_stats = {}

        self.summary = {}
        self.cylon_json = {
            "name": amplicon_set.name,
            "source": amplicon_set.fn,
            "amplicons": {},
        }

        for _, amplicon in amplicon_set.amplicons.items():
            self.summary[amplicon.name] = {
                "start": amplicon.start,
                "end": amplicon.end,
                "total_mapped_bases": 0,
                "total_depth": 0,
                "sampled_bases": 0,
                "sampled_depth": 0,
                "pass": False,
            }

            # store the global start and end position for the entire
            # primer scheme. This is used to help varifier.
            if not self.start_pos:
                self.start_pos = amplicon.start
            if not self.end_pos:
                self.end_pos = amplicon.end

            if amplicon.start < self.start_pos:
                self.start_pos = amplicon.start
            if amplicon.end > self.end_pos:
                self.end_pos = amplicon.end

            left_start, left_end = amplicon.left_primer_region
            right_start, right_end = amplicon.right_primer_region

            self.cylon_json["amplicons"][amplicon.name] = {
                "start": left_start,
                "end": right_end,
                "left_primer_end": left_end,
                "right_primer_start": right_start,
            }

        for fragment in bam.syncronise_fragments():
            self.count_fragment(fragment)

        random.seed(42)
        for fragment in bam.syncronise_fragments():
            # truncate number of reads to target count per amplicon
            self.push_fragment(fragment)

        for amplicon in self.amplicons:
            # we still want to randomise the order of the downsampled
            # amplicons. Cylon will further downsample from these
            # lists
            random.shuffle(self.amplicons[amplicon])

        self.summarise_amplicons()

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

    def count_fragment(self, fragment):
        amplicon = self.amplicon_set.match(fragment)
        if not amplicon:
            self.unmatched_reads += 1
            return

        self.reads_per_amplicon[amplicon] += 1

        self.summary[amplicon.name][
            "total_mapped_bases"
        ] += fragment.total_mapped_bases()
        self.summary[amplicon.name]["total_depth"] += 1

    def push_fragment(self, fragment):
        amplicon = self.amplicon_set.match(fragment)
        if not amplicon:
            return
        #        if len(self.amplicons[amplicon]) >= self.target_depth:
        #            return

        frags = self.reads_per_amplicon[amplicon]
        sample_rate = self.target_depth / frags

        if frags < self.target_depth or random.random() < sample_rate:
            # TODO count the observed primer extrema
            self.amplicons[amplicon].append(fragment)
            self.summary[amplicon.name][
                "sampled_bases"
            ] += fragment.total_mapped_bases()
            self.summary[amplicon.name]["sampled_depth"] += 1

    def summarise_amplicons(self):
        # normalise the bases per amplicons and such
        for amplicon in self.amplicons:
            self.amplicon_stats[amplicon] = {}
            self.amplicon_stats[amplicon][False] = 0
            self.amplicon_stats[amplicon][True] = 0
            for fragment in self.amplicons[amplicon]:
                self.amplicon_stats[amplicon][fragment.strand] += 1

            #print(
            #    f"amplicon strand bias:\t{amplicon.name}\t{self.amplicon_stats[amplicon][False]}/{self.amplicon_stats[amplicon][True]}",
            #    file=sys.stderr,
            #)

    def pileup(self, fasta, msa=None, minimap_presets=None):
        """remap reads to consensus
        """

        cons = mp.Aligner(str(fasta), preset=minimap_presets)
        if len(cons.seq_names) != 1:
            Exception(f"Consensus fasta {fasta} has more than one sequence")
        consensus_seq = cons.seq(cons.seq_names[0])
        pileup = self_qc.Pileup(consensus_seq, msa=msa)
        conspos_oob = 0  # consensus positions out-of-bounds
        for amplicon in self.amplicons:
            fragments = (
                self.amplicons[amplicon]
                if len(self.amplicons[amplicon]) >= self.target_depth
                else self.amplicons[amplicon]
            )
            for fragment in fragments:
                for read in fragment.reads:
                    a = cons.map(read.seq)  # remap to consensus
                    alignment = None
                    for x in a:
                        if x.is_primary:
                            alignment = x

                    if not alignment:
                        continue

                    assert alignment.q_en > alignment.q_st
                    assert alignment.r_en > alignment.r_st

                    aln = self_qc.parse_cigar(consensus_seq, read.seq, alignment)

                    # testing for indel conditions:
                    # ex = "".join(map(lambda x: x[1] if len(x[1]) == 1 else x[1], aln))
                    # c = consensus_seq[alignment.r_st : alignment.r_en]

                    primers = amplicon.match_primers(fragment)

                    for consensus_pos, call in self_qc.parse_cigar(
                        consensus_seq, read.seq, alignment
                    ):
                        if consensus_pos >= len(pileup):
                            # print(f"consensus pos out of bounds: {consensus_pos}/{len(pileup)}", file=sys.stderr)
                            conspos_oob += 1
                            continue

                        in_primer = False
                        for primer in primers:
                            # primers can be None
                            if primer and in_range(
                                (primer.ref_start, primer.ref_end),
                                # consensus_to_ref is 1-indexed!!
                                pileup.consensus_to_ref[consensus_pos + 1],
                            ):
                                in_primer = True
                                break

                        profile = self_qc.BaseProfile(
                            call, in_primer, read.is_reverse, amplicon,
                        )
                        pileup[consensus_pos].update(profile)
        print(f"consensus positions out of bounds: {conspos_oob}", file=sys.stderr)
        return pileup

    def reads_to_fastas(self, amplicon, outfile, target_bases):
        bases_out = 0
        with open(outfile, "w") as f:
            for i, fragment in enumerate(self[amplicon]):
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

    def make_reads_dir_for_cylon(self, outdir):
        """Makes a directory of reads for each amplicon, in the format required
        by `cylon assemble --reads_per_amp_dir`. Returns a set of amplicon
        names that should be failed because they had no reads"""
        os.mkdir(outdir)
        manifest_data = {}
        self.failed_amplicons = set()

        fasta_number = 0  # let's find another way
        for amplicon in self.amplicon_set:
            if len(self[amplicon]) == 0:
                self.failed_amplicons.add(amplicon)
                manifest_data[amplicon.name] = None
                continue
            outname = f"{fasta_number}.fa"
            outfile = os.path.join(outdir, outname)
            target_bases = self.cylon_target_depth_factor * len(amplicon)
            bases_out = self.reads_to_fastas(amplicon, outfile, target_bases)
            print(
                f"writing out {amplicon.name} reads {len(self[amplicon])}, {len(manifest_data)}.fa",
                file=sys.stderr,
            )
            fasta_number += 1

            # TODO: define failure
            # if bases_out < target_bases:
            # manifest_data[amplicon.name] = None
            #    self.failed_amplicons.add(amplicon)

            # TODO: check if we should output failed but not empty amplicon fastas
            manifest_data[amplicon.name] = outname

        return manifest_data
