from collections import Counter
import os
import pysam
import re

from viridian import utils


indel_re = re.compile(
    r"^(?P<ref>[a-zA-Z]+)(?P<sign>[-+])(?P<length>\d+)(?P<alt>[a-zA-Z]+)$"
)


def map_reads(
    bam_out,
    ref_fasta,
    reads1,
    reads2=None,
    threads=1,
    sample_name="sample",
    minimap_x_opt=None,
    debug=False,
    silent=False,
):
    if reads2 is None:
        if minimap_x_opt is None:
            minimap_x_opt = "-x map-ont"
        reads_opt = reads1
    else:
        if minimap_x_opt is None:
            minimap_x_opt = "-x sr"
        reads_opt = f"{reads1} {reads2}"

    r_opt = f"-R '@RG\\tID:1\\tSM:{sample_name}'"
    command = f"minimap2 {r_opt} -a -t {threads} {minimap_x_opt} {ref_fasta} {reads_opt} | samtools sort -O BAM -o {bam_out}"
    utils.syscall(command, quiet=not debug, silent=silent)
    utils.syscall(f"samtools index {bam_out}", quiet=not debug, silent=silent)


def pileup_indel_parse(indel):
    match = indel_re.search(indel)
    return match if match is None else match.groupdict()


def indels_from_pileup_counts(counts):
    indels = {x: {} for x in ["D", "d", "I", "i"]}
    for k, count in counts.items():
        match = pileup_indel_parse(k)
        if match is not None:
            match["length"] = int(match["length"])
            if match["sign"] == "-":
                key = "D" if match["ref"].isupper() else "d"
                key2 = match["length"]
            elif match["sign"] == "+":
                key = "I" if match["ref"].isupper() else "i"
                key2 = match["alt"]
            else:
                raise NotImplementedError()
            indels[key][key2] = indels[key].get(key2, 0) + count

    return indels


def pileup(
    tmp_bam,
    ref_fasta,
    reads1,
    reads2=None,
    threads=1,
    minimap_x_opt=None,
    debug=False,
):
    map_reads(
        tmp_bam,
        ref_fasta,
        reads1,
        reads2=reads2,
        threads=threads,
        minimap_x_opt=minimap_x_opt,
        debug=debug,
        silent=not debug,
    )
    fa = pysam.FastaFile(ref_fasta)
    aln_file = pysam.AlignmentFile(tmp_bam, "rb")
    pileup_counts = {}
    indels = {x: {} for x in ["D", "d", "I", "i"]}

    for p in aln_file.pileup(
        fastafile=fa, stepper="samtools", ignore_overlaps=False, compute_baq=False
    ):
        pileup_counts[p.pos] = dict(Counter(p.get_query_sequences(add_indels=True)))
        try:
            del pileup_counts[p.pos]["*"]
        except:
            pass
        pileup_counts[p.pos]["D"] = 0
        pileup_counts[p.pos]["d"] = 0
        pileup_counts[p.pos]["indel"] = {
            "D": indels["D"],
            "d": indels["d"],
        }
        indels = indels_from_pileup_counts(pileup_counts[p.pos])
        pileup_counts[p.pos]["indel"]["I"] = indels["I"]
        pileup_counts[p.pos]["indel"]["i"] = indels["i"]
        pileup_counts[p.pos]["I"] = sum(indels["I"].values())
        pileup_counts[p.pos]["i"] = sum(indels["i"].values())

    for pos, counts in pileup_counts.items():
        for d_or_D in "d", "D":
            for del_length, count in counts["indel"][d_or_D].items():
                for i in range(pos, pos + del_length):
                    pileup_counts[i][d_or_D] += count

    if not debug:
        os.unlink(tmp_bam)
        os.unlink(tmp_bam + ".bai")
    return pileup_counts
