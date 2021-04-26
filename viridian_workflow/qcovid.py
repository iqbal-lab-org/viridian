"""qcovid wrapper
"""
import logging
import subprocess
import os
from viridian_workflow.utils import run_process, check_file
from viridian_workflow.config import qcovid as cfg


def bin_amplicons(outdir, ref_genome, amplicon_bed, bam):
    mask = os.path.join(outdir, "mask")
    run_process(
        f"bin_amplicons.py --min_coverage {cfg.min_coverage} --min_template_match_75 {cfg.min_template_match_75} --mask {mask} {amplicon_bed} {bam}"
    )
    check_file(mask)
    return mask


def self_qc(outdir, assembly, bam):
    masked_fasta = os.path.join(outdir, "masked.fa")
    run_process(
        f"self_qc.py {assembly} {bam} --freq_threshold {cfg.variant_freq} --fasta",
        stdout=masked_fasta,
    )
    check_file(masked_fasta)
    return masked_fasta
