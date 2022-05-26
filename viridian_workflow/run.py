import sys
import shutil
from pathlib import Path

# import tempfile
import json

from viridian_workflow import primers, readstore
from viridian_workflow.subtasks import Cylon, Minimap, Varifier


def run_pipeline(
    work_dir,
    platform,
    fqs,
    amplicon_sets,
    ref="../covid/MN908947.fasta",
    force_amp_scheme=None,
    keep_intermediate=False,
    keep_bam=False,
    sample_name="sample",
    frs_threshold=0.1,
    self_qc_depth=20,
    consensus_max_n_percent=50,
    max_percent_amps_fail=50,
    command_line_args={},
):

    log = {}
    work_dir = Path(work_dir)
    if work_dir.exists():
        raise Exception
    work_dir.mkdir()

    # generate name-sorted bam from fastqs
    if platform == "illumina":
        fq1, fq2 = fqs
        minimap = Minimap(work_dir / "name_sorted.bam", ref, fq1, fq2=fq2, sort=False)
    elif platform == "ont":
        fq = fqs[0]
        minimap = Minimap(work_dir / "name_sorted.bam", ref, fq, sort=False)
    elif platform == "iontorrent":
        raise NotImplementedError
    else:
        print(f"Platform {platform} is not supported.", file=sys.stderr)
        exit(1)

    unsorted_bam = minimap.run()
    # add minimap task log to result log
    # log["minimap"] = minimap.log

    # pre-process input bam
    bam = readstore.Bam(unsorted_bam)

    # detect amplicon set
    amplicon_set = bam.detect_amplicon_set(amplicon_sets)
    # log["amplicons"] = bam.stats

    # construct readstore
    # this subsamples the reads
    rs = (
        readstore.ReadStore(amplicon_set, bam)
        if force_amp_scheme is None
        else readstore.ReadStore(force_amp_scheme, bam)
    )

    log["amplicons"] = rs.summary

    # save reads for cylon assembly
    amp_dir = work_dir / "amplicons"
    manifest_data = rs.make_reads_dir_for_cylon(amp_dir)

    # run cylon
    cylon = Cylon(work_dir, platform, ref, amp_dir, manifest_data, rs.cylon_json)
    consensus = cylon.run()
    log["initial_assembly"] = cylon.log

    # varifier
    varifier = Varifier(
        work_dir / "varifier",
        ref,
        consensus,
        min_coord=rs.start_pos,
        max_coord=rs.end_pos,
    )
    vcf, msa, varifier_consensus = varifier.run()
    log["varifier"] = varifier.log

    # self qc: remap reads to consensus
    pileup = rs.pileup(varifier_consensus, msa=msa)

    # mask output
    masked_fasta = pileup.mask()
    # log["self_qc"] = pileup.log
    log["qc"] = pileup.summary

    with open(work_dir / "masked.fasta", "w") as fasta_out:
        print(masked_fasta, file=fasta_out)

    # annotate vcf
    annotated_vcf = pileup.annotate_vcf(vcf)

    with open(work_dir / "final.vcf", "w") as vcf_out:
        header, records = annotated_vcf
        for h in header:
            print(h, file=vcf_out)
        for rec in records:
            print("\t".join(map(str, rec)), file=vcf_out)

    return log
