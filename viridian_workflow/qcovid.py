"""qcovid wrapper
"""
import logging
import subprocess
import os
from viridian_workflow.utils import run_process, check_file


def bin_amplicons(outdir, ref_genome, amplicon_bed, bam):
    mask = os.path.join(outdir, "mask")
    run_process(f"bin_amplicons.py --mask {mask} {amplicon_bed} {bam}")
    check_file(mask)
    return mask


def self_qc(outdir, assembly, bam):
    masked_fasta = os.path.join(outdir, "masked.fa")
    run_process(f"self_qc.py {bam} {assembly} --fasta", stdout=masked_fasta)
    check_file(masked_fasta)
    return masked_fasta
