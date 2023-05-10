import copy
from collections import defaultdict, namedtuple
import itertools
import logging
from operator import itemgetter
import os
import random

import pyfastaq
import intervaltree
import pysam

from viridian import maptools, utils

random.seed(42)

Fragment = namedtuple(
    "Fragment", ["left_primer", "right_primer", "name", "is_reverse", "seqs"]
)


def adapter_left_to_trim(read, left_primer_starts):
    ref_trim_pos = left_primer_starts.get(read.reference_start, None)
    if ref_trim_pos is None:
        return 0
    sclip = read.cigartuples[0][1] if read.cigartuples[0][0] == pysam.CSOFT_CLIP else 0
    return max(sclip + ref_trim_pos - read.reference_start, 0)


def adapter_right_to_trim(read, right_primer_ends):
    end = read.reference_end - 1
    ref_trim_pos = right_primer_ends.get(end, None)
    if ref_trim_pos is None:
        return 0
    sclip = (
        read.cigartuples[-1][1] if read.cigartuples[-1][0] == pysam.CSOFT_CLIP else 0
    )
    return sclip + end - ref_trim_pos


def to_adapter_trimmed_seq(read, left_primer_starts, right_primer_ends):
    if read.is_paired:
        if read.template_length > 0:
            trim_left = adapter_left_to_trim(read, left_primer_starts)
            trim_right = 0
        else:
            trim_left = 0
            trim_right = adapter_right_to_trim(read, right_primer_ends)
    else:
        trim_left = adapter_left_to_trim(read, left_primer_starts)
        trim_right = adapter_right_to_trim(read, right_primer_ends)

    if read.is_reverse:
        trim_left = None if trim_left == 0 else -trim_left
        return read.get_forward_sequence()[trim_right:trim_left]
    else:
        trim_right = None if trim_right == 0 else -trim_right
        return read.query_sequence[trim_left:trim_right]


def read_to_fragment_coords(read):
    if read.is_paired:
        if read.template_length > 0:
            return read.reference_start, read.reference_start + read.template_length - 1
        else:
            return read.reference_end + read.template_length, read.reference_end - 1
    else:
        return read.reference_start, read.reference_end - 1


def run_ngmerge(fq1, fq2, outfile):
    unmerged = f"{outfile}.tmp.unmerged"
    merged = f"{outfile}.tmp.merged"
    unmerged1 = f"{unmerged}_1.fastq"
    unmerged2 = f"{unmerged}_2.fastq"
    silent = logging.getLogger().level > logging.DEBUG
    utils.syscall(f"NGmerge -1 {fq1} -2 {fq2} -o {merged} -f {unmerged}", silent=silent)
    read_count = 1
    with open(outfile, "w") as f_out:
        for infile in (merged, unmerged1, unmerged2):
            reader = pyfastaq.sequences.file_reader(infile)
            for seq in reader:
                print(f">{read_count}", seq.seq, sep="\n", file=f_out)
                read_count += 1
    utils.syscall(f"rm -f {merged} {unmerged1} {unmerged2}", silent=True)


def write_cylon_fasta_paired(fragments, outfile, target_bases):
    total_bases = 0
    tmp_out1 = f"{outfile}.tmp.1.fq"
    tmp_out2 = f"{outfile}.tmp.2.fq"
    with open(tmp_out1, "w") as f1, open(tmp_out2, "w") as f2:
        for i, frag in enumerate(fragments):
            print(
                f"@{i} /1",
                frag.seqs["read"],
                "+",
                "I" * len(frag.seqs["read"]),
                sep="\n",
                file=f1,
            )
            print(
                f"@{i} /2",
                frag.seqs["mate"],
                "+",
                "I" * len(frag.seqs["mate"]),
                sep="\n",
                file=f2,
            )
            total_bases += sum(map(len, frag.seqs.values()))
            if total_bases >= target_bases:
                break

    run_ngmerge(tmp_out1, tmp_out2, outfile)
    os.unlink(tmp_out1)
    os.unlink(tmp_out2)
    return total_bases


def write_cylon_fasta_unpaired(fragments, outfile, target_bases):
    total_bases = 0
    fwd_reads = [f.seqs["read"] for f in fragments if not f.is_reverse]
    rev_reads = [f.seqs["read"] for f in fragments if f.is_reverse]
    min_reads = min(len(fwd_reads), len(rev_reads))
    if min_reads == 0:
        return 0

    with open(outfile, "w") as f:
        for i in range(min_reads):
            print(f">{2*i}", fwd_reads[i], sep="\n", file=f)
            print(f">{2*i+1}", rev_reads[i], sep="\n", file=f)
            total_bases += len(fwd_reads[i]) + len(rev_reads[i])
            if total_bases >= target_bases:
                break
    return total_bases


def write_qc_fastas(fragments, outprefix):
    by_primers = defaultdict(list)
    for frag in fragments:
        by_primers[(frag.left_primer, frag.right_primer)].append(frag.seqs)
    paired = fragments[0].seqs["mate"] != ""
    files_written = {}

    for (left_i, right_i), seqs in by_primers.items():
        out = f"{outprefix}.{left_i}.{right_i}"

        if paired:
            out1 = f"{out}.1.fa"
            out2 = f"{out}.2.fa"
            with open(out1, "w") as f1, open(out2, "w") as f2:
                for i, seq_dict in enumerate(seqs):
                    print(f">{i}/1", seq_dict["read"], sep="\n", file=f1)
                    print(f">{i}/2", seq_dict["mate"], sep="\n", file=f2)
        else:
            out1 = f"{out}.fa"
            out2 = None
            with open(out1, "w") as f:
                for i, seq_dict in enumerate(seqs):
                    print(f">{i}", seq_dict["read"], sep="\n", file=f)

        files_written[(left_i, right_i)] = (out1, out2)

    return files_written


class ReadSampler:
    def __init__(
        self,
        amplicons,
        qc_depth=1000,
        cylon_depth=200,
        adapter_trim_tolerance=3,
    ):
        self.bam_file = None
        self.amplicons = copy.deepcopy(amplicons)
        self.amplicons.sort(key=itemgetter("start"))
        self.qc_depth = qc_depth
        self.cylon_depth = cylon_depth
        self.adapter_trim_tolerance = adapter_trim_tolerance
        self.outdir = None
        self.qc_dir = None
        self.cylon_dir = None
        self.sample_stats = {}
        self._init_lookups()
        self.scheme_start = self.amplicons[0]["start"]
        self.scheme_end = self.amplicons[-1]["end"]
        self.read_counts = {
            "unpaired_reads": 0,
            "reads1": 0,
            "reads2": 0,
            "total_reads": 0,
            "mapped": 0,
            "match_any_amplicon": 0,
            "read_lengths": {},
        }

    def _init_lookups(self):
        self.left_primer_starts = {}
        self.right_primer_ends = {}
        self.amp_tree = intervaltree.IntervalTree()
        for amp_index, amp in enumerate(self.amplicons):
            amp["primers"]["left"].sort(key=itemgetter("start"))
            amp["primers"]["right"].sort(key=itemgetter("start"))
            starts = {}
            ends = {}
            for primer_index, d in enumerate(amp["primers"]["left"]):
                starts[d["start"] - self.adapter_trim_tolerance] = primer_index
                for i in range(self.adapter_trim_tolerance + 1):
                    self.left_primer_starts[d["start"] - i] = d["start"]

            for primer_index, d in enumerate(amp["primers"]["right"]):
                ends[d["end"] + self.adapter_trim_tolerance] = primer_index
                for i in range(self.adapter_trim_tolerance + 1):
                    self.right_primer_ends[d["end"] + i] = d["end"]

            for start, end in itertools.product(starts, ends):
                self.amp_tree[start:end] = (amp_index, starts[start], ends[end])

    def reject_read_and_update_read_counts(self, read):
        if read.is_secondary or read.is_supplementary:
            return True

        self.read_counts["total_reads"] += 1
        self.read_counts["read_lengths"][read.query_length] = (
            self.read_counts["read_lengths"].get(read.query_length, 0) + 1
        )

        if read.is_paired:
            if read.is_read1:
                self.read_counts["reads1"] += 1
            else:
                self.read_counts["reads2"] += 1
        else:
            self.read_counts["unpaired_reads"] += 1

        if read.is_unmapped:
            return True

        self.read_counts["mapped"] += 1
        return read.is_paired and not read.is_proper_pair

    def overhangs_first_or_last_amplicon(self, start, end):
        if end <= self.amplicons[0]["end"]:
            return self.best_amp_containing_interval(self.scheme_start, end)
        elif start >= self.amplicons[-1]["start"]:
            return self.best_amp_containing_interval(start, self.scheme_end)
        else:
            return None

    def best_amp_containing_interval(self, start, end):
        hits = [
            x[2] for x in self.amp_tree[start:end] if x.begin <= start < end <= x.end
        ]
        if len(hits) == 0:
            return None
        elif len(hits) == 1:
            return hits[0]

        # if more than one hit, choose an amplicon at random, then return
        # the shortest match from that amplicon
        random_amp = random.choice(list(set(x[0] for x in hits)))
        return sorted([x for x in hits if x[0] == random_amp], key=lambda x: len(x))[0]

    def write_reads_files_one_amp(self, amplicon_index, fragments):
        amplicon = self.amplicons[amplicon_index]
        logging.debug(f"Making sampled reads files for amplicon {amplicon['name']}")
        stats = {
            "total_reads": len(fragments),
            "total_reads_fwd_strand": 0,
            "total_reads_rev_strand": 0,
            "qc_reads": 0,
            "qc_bases": 0,
            "qc_depth": 0,
            "assemble_bases": 0,
        }
        self.sample_stats[amplicon["name"]] = stats
        if len(fragments) == 0:
            return
        random.shuffle(fragments)
        amp_length = amplicon["end"] - amplicon["start"] + 1
        target_qc_bases = amp_length * self.qc_depth

        for frag in fragments:
            if frag.seqs["mate"] == "":
                if frag.is_reverse:
                    stats["total_reads_rev_strand"] += 1
                else:
                    stats["total_reads_fwd_strand"] += 1
            else:
                stats["total_reads_rev_strand"] += 1
                stats["total_reads_fwd_strand"] += 1

        for i, frag in enumerate(fragments):
            if stats["qc_bases"] > target_qc_bases:
                break
            stats["qc_bases"] += sum(map(len, frag.seqs.values()))

        fragments = fragments[:i]
        if len(fragments) == 0:
            return
        self.qc_reads_files[amplicon_index] = write_qc_fastas(
            fragments, os.path.join(self.qc_dir, f"{amplicon_index}.reads")
        )

        cylon_fa = f"{len(self.cylon_fa_manifest)}.fa"
        self.cylon_fa_manifest[amplicon["name"]] = cylon_fa
        cylon_fa = os.path.join(self.cylon_dir, cylon_fa)
        cylon_target_bases = amp_length * self.cylon_depth

        stats["qc_reads"] = len(fragments)
        paired = fragments[0].seqs["mate"] != ""
        if paired:
            stats["qc_reads"] *= 2
            stats["total_reads"] *= 2
            stats["assemble_bases"] = write_cylon_fasta_paired(
                fragments, cylon_fa, cylon_target_bases
            )
        else:
            stats["assemble_bases"] = write_cylon_fasta_unpaired(
                fragments, cylon_fa, cylon_target_bases
            )
            if stats["total_reads"] > 0:
                fwd_percent = round(100 * stats["total_reads_fwd_strand"] / stats["total_reads"])
                if not 0.10 <= fwd_percent <= 0.90:
                    logging.warning(f"Reads for amplicon {amplicon['name']} mostly all mapped to one strand. {fwd_percent}% mapped to forwards strand")
        stats["qc_depth"] = round(stats["qc_bases"] / amp_length, 2)
        if stats["assemble_bases"] == 0:
            del self.cylon_fa_manifest[amplicon["name"]]

        logging.debug(
            f"Finished making sampled reads files for amplicon {amplicon['name']}"
        )

    def sample_reads(self, bam_file, outdir):
        self.outdir = os.path.join(outdir)
        self.qc_dir = os.path.join(outdir, "qc")
        self.cylon_dir = os.path.join(outdir, "cylon")
        os.mkdir(self.outdir)
        os.mkdir(self.qc_dir)
        os.mkdir(self.cylon_dir)
        aln_file = pysam.AlignmentFile(bam_file)
        wait_for_mate = {}
        next_to_finish_index = 0
        next_to_finish = self.amplicons[next_to_finish_index]
        reads = {next_to_finish_index: []}
        self.cylon_fa_manifest = {}
        self.sample_stats = {}
        self.qc_reads_files = {}
        at_final_reads = False

        for read in aln_file.fetch(until_eof=True):
            if self.reject_read_and_update_read_counts(read) or at_final_reads:
                continue

            while read.pos > next_to_finish["end"]:
                self.write_reads_files_one_amp(
                    next_to_finish_index, reads.get(next_to_finish_index, [])
                )
                try:
                    del reads[next_to_finish_index]
                except:
                    pass
                next_to_finish_index += 1
                if next_to_finish_index >= len(self.amplicons):
                    at_final_reads = True
                    break
                else:
                    next_to_finish = self.amplicons[next_to_finish_index]

            if at_final_reads or next_to_finish_index >= len(self.amplicons):
                at_final_reads = True
                continue

            if read.is_paired and read.query_name in wait_for_mate:
                self.read_counts["match_any_amplicon"] += 1
                amp_index, fragment = wait_for_mate.pop(read.query_name)
                key = "read" if read.is_read1 else "mate"
                fragment.seqs[key] = to_adapter_trimmed_seq(
                    read, self.left_primer_starts, self.right_primer_ends
                )
                assert "" not in fragment.seqs.values()
                if amp_index not in reads:
                    reads[amp_index] = []
                reads[amp_index].append(fragment)
                continue

            frag_start, frag_end = read_to_fragment_coords(read)
            assert frag_start < frag_end
            containing_amp = self.best_amp_containing_interval(frag_start, frag_end)
            if containing_amp is None:
                containing_amp = self.overhangs_first_or_last_amplicon(
                    frag_start, frag_end
                )

            if containing_amp is None:
                continue

            seq = to_adapter_trimmed_seq(
                read, self.left_primer_starts, self.right_primer_ends
            )

            if read.is_paired:
                if read.is_read1:
                    seqs = {"read": seq, "mate": ""}
                else:
                    seqs = {"read": "", "mate": seq}

                fragment = Fragment(
                    containing_amp[1],
                    containing_amp[2],
                    read.query_name,
                    read.is_reverse,
                    seqs,
                )
                wait_for_mate[read.query_name] = (containing_amp[0], fragment)
            else:
                self.read_counts["match_any_amplicon"] += 1
                if containing_amp[0] not in reads:
                    reads[containing_amp[0]] = []
                reads[containing_amp[0]].append(
                    Fragment(
                        containing_amp[1],
                        containing_amp[2],
                        read.query_name,
                        read.is_reverse,
                        {"read": seq, "mate": ""},
                    )
                )

        for i in range(next_to_finish_index, len(self.amplicons)):
            self.write_reads_files_one_amp(i, reads.get(i, []))

        utils.write_json(
            os.path.join(self.cylon_dir, "manifest.json"), self.cylon_fa_manifest
        )

    def pileups(self, consensus_fasta, tmp_outdir, debug=False):
        pileups = {}
        os.mkdir(tmp_outdir)
        for amplicon, reads_dict in self.qc_reads_files.items():
            pileups[amplicon] = {}
            for (left_i, right_i), (reads1, reads2) in reads_dict.items():
                logging.debug(
                    f"Getting pileup for amplicon/primers {amplicon}/{left_i}.{right_i}"
                )
                tmp_bam = os.path.join(tmp_outdir, f"{amplicon}.{left_i}.{right_i}.bam")
                pileups[amplicon][(left_i, right_i)] = maptools.pileup(
                    tmp_bam,
                    consensus_fasta,
                    reads1,
                    reads2,
                    debug=debug,
                )

            if len(pileups) % 10 == 0:
                logging.info(
                    f"(pileup) Got pileup for {len(pileups)}/{len(self.amplicons)} amplicons"
                )

        logging.info("(pileup) Got pileup data for all amplicons")
        return pileups
