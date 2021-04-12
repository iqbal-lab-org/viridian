"""varifier wrapper
"""
import logging
import subprocess
import os
from viridian_workflow.utils import run_process, check_file


def run(outdir, ref_genome, assembly):
    vcf = os.path.join(outdir, "04.truth.vcf")
    run_process(f"varifier make_truth_vcf {assembly} {ref_genome} {outdir}")
    check_file(vcf)
    return vcf
