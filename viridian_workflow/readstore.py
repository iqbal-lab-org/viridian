from __future__ import annotations

from typing import Optional, Any
from collections import namedtuple, defaultdict
import sys
from pathlib import Path
import os
import random

from viridian_workflow import utils
from viridian_workflow.primers import Amplicon, AmpliconSet, Primer
from viridian_workflow.reads import Read, Fragment, PairedReads, SingleRead

import pysam  # type: ignore
import mappy as mp  # type: ignore


def score(
    matches: defaultdict[AmpliconSet, int],
    mismatches: defaultdict[AmpliconSet, int],
    disqualification_threshold: float = 0.35,
) -> Optional[AmpliconSet]:
    """Assign winning amplicon set id based on match stats"""
    amplicon_sets = set([*matches.keys(), *mismatches.keys()])

    m = 0
    winner = None
    for amplicon_set in amplicon_sets:
        total = matches[amplicon_set] + mismatches[amplicon_set]
        mismatch_proportion = mismatches[amplicon_set] / total
        if mismatch_proportion > disqualification_threshold:
            # if more than 3% of reads break the amplicon boundaries
            # disqualify this amplicon set
            print(
                amplicon_set.name,
                mismatches[amplicon_set],
                matches[amplicon_set],
                f" disqualified: {mismatch_proportion * 100}% ({disqualification_threshold * 100}% threshold)",
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
    def __init__(
        self,
        bam: Path,
        infile_is_paired: Optional[bool] = None,
        template_length_threshold: int = 150,
    ):
        self.infile_is_paired: Optional[bool] = infile_is_paired
        if not Path(bam).is_file():
            raise Exception(f"bam file {bam} does not exist")
        self.bam: Path = bam
        self.template_length_threshold: int = template_length_threshold
        self.stats: dict[str, Any] = {}

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
            "templates_that_were_too_short": defaultdict(int),
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
                single_read: Fragment = SingleRead(Bam.read_from_pysam(read))
                tlen = single_read.ref_end - single_read.ref_start
                self.stats["template_lengths"][tlen] += 1
                if tlen < self.template_length_threshold:
                    self.stats["templates_that_were_too_short"][tlen] += 1
                    continue
                yield single_read

            if not read.is_proper_pair:
                improper_pairs += 1
                continue

            if read.is_read1:
                reads_by_name[read.query_name] = Bam.read_from_pysam(read)

            elif read.is_read2:
                if read.query_name not in reads_by_name:
                    raise Exception("Bam file is not sorted by name")
                read1 = reads_by_name[read.query_name]
                paired_reads: Fragment = PairedReads(read1, Bam.read_from_pysam(read))
                tlen = paired_reads.ref_end - paired_reads.ref_start
                self.stats["template_lengths"][tlen] += 1
                if tlen < self.template_length_threshold:
                    self.stats["templates_that_were_too_short"][tlen] += 1
                    continue
                yield paired_reads
                del reads_by_name[read.query_name]
        print(f"{improper_pairs} improper pairs", file=sys.stderr)

    def detect_amplicon_set(
        self, amplicon_sets: list[AmpliconSet], disqualification_threshold: float = 0.35
    ) -> AmpliconSet:
        """return inferred amplicon set from list of amplicon sets
        """

        mismatches: defaultdict[AmpliconSet, int] = defaultdict(int)
        matches: defaultdict[AmpliconSet, int] = defaultdict(int)

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
        chosen_scheme = score(
            matches, mismatches, disqualification_threshold=disqualification_threshold
        )
        if chosen_scheme:
            self.stats["chosen_amplicon_scheme"] = chosen_scheme.name
        else:
            # TODO: decide on behaviour when no appropriate scheme is chosen
            # current policy: abort
            raise Exception("failed to choose amplicon scheme")
        return chosen_scheme


class ReadStore:
    def __init__(self, amplicon_set: AmpliconSet, bam: Bam, target_depth: int = 1000):
        self.amplicons: defaultdict[Amplicon, list[Fragment]] = defaultdict(list)
        self.reads_per_amplicon: defaultdict[Amplicon, int] = defaultdict(int)
        self.amplicon_set: AmpliconSet = amplicon_set
        self.reads_all_paired: Optional[bool] = bam.infile_is_paired
        self.unmatched_reads: int = 0

        # Index of positions (0-based, wrt Reference) where different amplicons
        # have contributed base calls
        self.multiple_amplicon_support: list[bool] = []

        self.target_depth: int = target_depth

        # TODO find a home for this magic number
        self.cylon_target_depth_factor = 200

        self.start_pos = None
        self.end_pos = None
        self.amplicon_stats: dict[Amplicon, dict[bool, int]] = {}

        self.summary = {}
        self.cylon_json: dict[str, Any] = {
            "name": amplicon_set.name,
            "source": amplicon_set.fn,
            "amplicons": {},
        }

        self.primer_histogram: dict[Amplicon, dict[str, defaultdict[Primer, int]]] = {}
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
            if amplicon not in self.primer_histogram:
                self.primer_histogram[amplicon] = {
                    "left": defaultdict(int),
                    "right": defaultdict(int),
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

        for fragment in bam.syncronise_fragments():
            self.count_fragment(fragment)

        random.seed(42)
        for fragment in bam.syncronise_fragments():
            # truncate number of reads to target count per amplicon
            self.push_fragment(fragment)

        for amplicon in self.amplicon_set:

            # decide if threshold for primers is met
            p1_min, p2_max = None, None
            if amplicon in self.primer_histogram:
                p1_min, p2_max = self.filter_primer_counts(
                    self.primer_histogram[amplicon],
                    start=self.start_pos,
                    end=self.end_pos,
                )
            # default to the inner-most primer coords if an extrema isn't found
            if (
                amplicon.left_primer_region is None
                or amplicon.right_primer_region is None
            ):
                raise Exception(
                    "Attempted to test primer regions before initialising amplicon"
                )
            p1_start = (
                amplicon.left_primer_region[0] if p1_min is None else p1_min.ref_start
            )
            p2_end = (
                amplicon.right_primer_region[1] if p2_max is None else p2_max.ref_end
            )

            p1_end = (
                amplicon.left_primer_region[1] if p1_min is None else p1_min.ref_end
            )
            p2_start = (
                amplicon.right_primer_region[0] if p2_max is None else p2_max.ref_start
            )

            self.cylon_json["amplicons"][amplicon.name] = {
                "start": p1_start,
                "end": p2_end,
                "left_primer_end": p1_end,
                "right_primer_start": p2_start,
            }

        for amplicon in self.amplicons:
            # we still want to randomise the order of the downsampled
            # amplicons. Cylon will further downsample from these
            # lists
            random.shuffle(self.amplicons[amplicon])
        self.summarise_amplicons()

    @staticmethod
    def filter_primer_counts(primer_counts, threshold=100, start=0, end=sys.maxsize):
        p1 = None
        p2 = None
        p1_min = end
        p2_max = start
        for primer, count in primer_counts["left"].items():
            if count < threshold:  # primer occurence threshold
                # exclude this primer
                continue
            if primer.ref_start <= p1_min:
                p1_min = primer.ref_start
                p1 = primer
        for primer, count in primer_counts["right"].items():
            if count < threshold:
                continue
            if primer.ref_end >= p2_max:
                p2_max = primer.ref_end
                p2 = primer
        return p1, p2

    def __eq__(self, other):
        raise NotImplementedError

    def __str__(self):
        raise NotImplementedError

    def __iter__(self):
        raise NotImplementedError

    def __getitem__(self, amplicon: Amplicon) -> list[Fragment]:
        """Given an amplicon, returns list of Fragments"""
        return self.amplicons[amplicon]

    def fetch(self, start=0, end=None):
        raise NotImplementedError

    def count_fragment(self, fragment: Fragment):
        amplicon = self.amplicon_set.match(fragment)
        if not amplicon:
            self.unmatched_reads += 1
            return

        self.reads_per_amplicon[amplicon] += 1

        self.summary[amplicon.name][
            "total_mapped_bases"
        ] += fragment.total_mapped_bases()
        self.summary[amplicon.name]["total_depth"] += 1

    def push_fragment(self, fragment: Fragment):
        amplicon: Optional[Amplicon] = self.amplicon_set.match(fragment)
        if amplicon is None:
            return

        frags = self.reads_per_amplicon[amplicon]
        sample_rate = self.target_depth / frags
        if frags < self.target_depth or random.random() < sample_rate:
            p1, p2 = amplicon.match_primers(fragment)
            if p1 is not None:
                self.primer_histogram[amplicon]["left"][p1] += 1
            if p2 is not None:
                self.primer_histogram[amplicon]["right"][p2] += 1

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

    def reads_to_fastas(
        self, amplicon: Amplicon, outfile: Path, target_bases: int
    ) -> int:
        bases_out: int = 0
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
                # manifest_data[amplicon.name] = None
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
