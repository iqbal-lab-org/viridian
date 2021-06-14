import json
import logging
import os
import pysam
import random
import subprocess

from viridian_workflow import utils

random.seed(42)


class ReadSampler:
    def __init__(
        self,
        ref_fasta,
        bam_in,
        outprefix,
        amplicons_bed,
        target_depth,
        min_read_overlap_proportion=0.5,
    ):
        self.ref_genome = utils.load_single_seq_fasta(ref_fasta)
        self.amplicons = utils.load_amplicons_bed_file(amplicons_bed)
        self.bam_in = bam_in
        self.aln_file_in = None
        self.aln_file_out = None
        self.outprefix = outprefix
        self.bam_out = self.outprefix + ".bam"
        self.json_out = self.outprefix + ".json"
        self.fq_out = self.outprefix + ".fastq"
        self.fq_out1 = self.outprefix + ".reads_1.fastq"
        self.fq_out2 = self.outprefix + ".reads_2.fastq"
        self.reads_are_paired = None
        self.target_depth = target_depth
        self.mates_wanted = None
        self.used_reads = None
        self.min_read_overlap_proportion = min_read_overlap_proportion
        self.read_cache = {}
        self.output_json_data = {}

    def bam_is_paired_reads(self, aln_file):
        for read in aln_file.fetch():
            return read.is_paired

    def read_overlap_length_with_amplicon(self, read, amplicon):
        if amplicon is None:
            return 0
        try:
            end = min(amplicon.end, read.reference_end)
            start = max(amplicon.start, read.reference_start)
        except:
            return 0
        overlap = max(0, end - start + 1)
        if overlap / read.query_length > self.min_read_overlap_proportion:
            return overlap
        else:
            return 0

    def get_reads_for_amplicon(self, amplicon_index, amplicon, next_amplicon):
        total_bases_mapped = 0
        reads = {}
        previous_index = None if amplicon_index == 0 else amplicon_index - 1

        for read in self.aln_file_in.fetch(
            self.ref_genome.id, amplicon.start, amplicon.end + 1
        ):
            if read.is_supplementary or read.is_secondary:
                continue

            if self.reads_are_paired != read.is_paired:
                raise Exception(
                    "Mix of paired and unpaired reads in BAM file {self.bam_in}. Cannot continue"
                )

            if read.query_name in self.mates_wanted:
                mate = self.mates_wanted[read.query_name]
                if read.is_read1 != mate.is_read1:
                    if read.query_name not in self.used_reads:
                        self.aln_file_out.write(read)
                        self.aln_file_out.write(self.mates_wanted[read.query_name])
                        total_bases_mapped += self.read_overlap_length_with_amplicon(
                            read, amplicon
                        )
                        self.used_reads.add(read.query_name)
                    del self.mates_wanted[read.query_name]
                continue

            if (
                read.query_name in self.used_reads
                or read.mate_is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_qcfail
            ):
                continue

            overlap_len = self.read_overlap_length_with_amplicon(read, amplicon)
            next_overlap_len = self.read_overlap_length_with_amplicon(
                read, next_amplicon
            )
            total_bases_mapped += overlap_len

            if overlap_len == 0 or overlap_len < next_overlap_len:
                if self.reads_are_paired:
                    self.read_cache[amplicon_index][
                        (read.query_name, read.is_read1)
                    ] = read
                continue

            if read.query_name not in reads:
                reads[read.query_name] = {}

            reads[read.query_name][read.is_read1] = (overlap_len, read)
            if self.reads_are_paired and previous_index is not None:
                self.read_cache[previous_index].pop(
                    (read.query_name, read.is_read1), None
                )

        return reads, total_bases_mapped

    def write_reads_for_one_amplicon(self, reads, amplicon, amplicon_index):
        read_names = sorted(list(reads.keys()))
        random.shuffle(read_names)
        total_overlap_bases = 0
        amplicon_length = amplicon.end - amplicon.start + 1
        target_output_bases = self.target_depth * amplicon_length
        previous_index = None if amplicon_index == 0 else amplicon_index - 1

        for read_name in read_names:
            read_dict = reads[read_name]
            if self.reads_are_paired:
                overlap1, read1 = read_dict.get(True, (0, None))
                if read1 is None and previous_index is not None:
                    read1 = self.read_cache[previous_index].get((read_name, True), None)
                    overlap1 = self.read_overlap_length_with_amplicon(read1, amplicon)
            else:
                overlap1 = 0
                read1 = None

            overlap2, read2 = read_dict.get(False, (0, None))
            if read2 is None and previous_index is not None:
                read2 = self.read_cache[previous_index].get((read_name, False), None)
                overlap2 = self.read_overlap_length_with_amplicon(read2, amplicon)

            if read1 is None and read2 is None:
                continue

            if self.reads_are_paired:
                if read1 is not None and read2 is not None:
                    self.aln_file_out.write(read1)
                    self.aln_file_out.write(read2)
                    self.used_reads.add(read1.query_name)
                elif read1 is not None and read2 is None:
                    self.mates_wanted[read1.query_name] = read1
                elif read1 is None and read2 is not None:
                    self.mates_wanted[read2.query_name] = read2
            elif read2 is not None:
                self.aln_file_out.write(read2)
                self.used_reads.add(read2)
            else:
                continue

            total_overlap_bases += overlap1 + overlap2
            if total_overlap_bases >= target_output_bases:
                break

        return total_overlap_bases

    def initialise_output_json_data(self):
        for amplicon in self.amplicons:
            self.output_json_data[amplicon.name] = {
                "start": amplicon.start + 1,
                "end": amplicon.end + 1,
                "total_mapped_bases": 0,
                "total_depth": 0,
                "sampled_bases": 0,
                "sampled_depth": 0,
            }

    def finalise_output_json_data(self):
        for name, data in self.output_json_data.items():
            amp_len = data["end"] - data["start"] + 1
            data["total_depth"] = round(data["total_mapped_bases"] / amp_len, 2)
            data["sampled_depth"] = round(data["sampled_bases"] / amp_len, 2)

    def write_fastq(self):
        if self.reads_are_paired:
            tmp_bam = self.bam_out + ".tmp.sort_by_name.bam"
            pysam.sort("-n", "-o", tmp_bam, self.bam_out)
            pysam.fastq("-N", "-1", self.fq_out1, "-2", self.fq_out2, tmp_bam)
            os.unlink(tmp_bam)
        else:
            pysam.fastq("-0", self.fq_out, self.bam_out)

    def run(self):
        logging.info(
            f"Start sampling reads from {self.bam_in}, target depth is {self.target_depth}"
        )
        self.aln_file_in = pysam.AlignmentFile(self.bam_in, "rb")
        self.reads_are_paired = self.bam_is_paired_reads(self.aln_file_in)
        unsorted_bam = self.bam_out + ".tmp.unsorted.bam"
        self.aln_file_out = pysam.AlignmentFile(
            unsorted_bam, "wb", template=self.aln_file_in
        )
        self.mates_wanted = {}
        self.used_reads = set()
        self.initialise_output_json_data()
        self.read_cache = {i: {} for i in range(len(self.amplicons))}

        for i, amplicon in enumerate(self.amplicons):
            logging.info(
                f"Start processing amplicon {amplicon.name} ({i+1}/{len(self.amplicons)})"
            )
            if i <= len(self.amplicons) - 2:
                next_amplicon = self.amplicons[i + 1]
            else:
                next_amplicon = None

            amplicon_reads, bases_mapped = self.get_reads_for_amplicon(
                i, amplicon, next_amplicon
            )
            sampled_bases = self.write_reads_for_one_amplicon(
                amplicon_reads, amplicon, i
            )
            if i > 0:
                del self.read_cache[i - 1]

            self.output_json_data[amplicon.name]["total_mapped_bases"] = bases_mapped
            self.output_json_data[amplicon.name]["sampled_bases"] = sampled_bases
            logging.info(
                f"Finish processing amplicon {amplicon.name} ({i+1}/{len(self.amplicons)})"
            )

        self.finalise_output_json_data()
        self.aln_file_out.close()
        logging.info(f"Sorting sampled BAM file")
        pysam.sort("-o", self.bam_out, unsorted_bam)
        os.unlink(unsorted_bam)
        logging.info(f"Indexing sampled BAM file")
        pysam.index(self.bam_out)
        logging.info(f"Writing fastq file(s)")
        self.write_fastq()

        logging.info(f"Writing summary JSON file {self.json_out}")
        with open(self.json_out, "w") as f:
            json.dump(self.output_json_data, f, indent=2)

        logging.info(f"Finished sampling reads. Final BAM file is {self.bam_out}")


def sample_reads(
    ref_fasta,
    bam_in,
    outprefix,
    amplicons_bed,
    target_depth=1000,
    min_read_overlap_proportion=0.5,
):
    sampler = ReadSampler(
        ref_fasta,
        bam_in,
        outprefix,
        amplicons_bed,
        target_depth,
        min_read_overlap_proportion=min_read_overlap_proportion,
    )
    sampler.run()
    return sampler
