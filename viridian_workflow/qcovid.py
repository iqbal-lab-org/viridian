"""qcovid wrapper
"""
import logging
import subprocess
import os
import json
from viridian_workflow.utils import (
    run_process,
    run_process_stdout,
    check_file,
    PrimerError,
)
from viridian_workflow.config import qcovid as cfg


def detect_primers_pe(outdir, ref_genome, fq1, fq2, primers=None):
    if not primers:
        primers = cfg.primers

    primer_string = ",".join(primers)

    primer_set = run_process_stdout(
        ["detect_primers.py", "--json", ref_genome, primer_string, fq1, fq2,]
    )
    logging.error(f"primer set: {primer_set}")
    primer_set = json.loads(primer_set)
    if primer_set is None:
        raise PrimerError()

    if primer_set["status"] != "success":
        logging.info("Failed to detect primers: {primer_set['status']}")
        raise PrimerError()

    return primer_set["primer_set"]


def detect_primers_se(outdir, ref_genome, fq, primers=None):
    if not primers:
        primers = cfg.primers

    primer_string = ",".join(primers)

    primer_set = json.loads(
        run_process_stdout(
            ["detect_primers.py", "--json", ref_genome, primer_string, fq,]
        )
    )
    if primer_set is None:
        raise Exception("primer detection failed totally")

    if primer_set["status"] != "success":
        logging.info("Failed to detect primers: {primer_set['status']}")
        raise Exception("could not infer detect primer set")

    return primer_set["primer_set"]


def bin_amplicons(outdir, ref_genome, amplicon_bed, bam):
    mask = os.path.join(outdir, "mask")
    run_process(
        f"bin_amplicons.py --min_coverage {cfg.min_coverage} --min_template_match_75 {cfg.min_template_match_75} --mask {mask} {amplicon_bed} {bam}"
    )
    check_file(mask)
    return mask


def bin_amplicons_se(outdir, ref_genome, amplicon_bed, bam):
    mask = os.path.join(outdir, "mask")
    run_process(
        f"bin_amplicons.py --se --min_coverage {cfg.min_coverage} --min_template_match_75 {cfg.min_template_match_75} --mask {mask} {amplicon_bed} {bam}"
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
