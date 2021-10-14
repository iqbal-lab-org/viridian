"""varifier wrapper
"""
import logging
import subprocess
import os
from viridian_workflow.utils import run_process, check_file


def run(outdir, ref_genome, assembly, min_coord=None, max_coord=None):
    vcf = os.path.join(outdir, "04.truth.vcf")
    options = ["--global_align"]
    if min_coord is not None:
        options.append(f"--global_align_min_coord {min_coord + 1}")
    if max_coord is not None:
        options.append(f"--global_align_max_coord {max_coord + 1}")
    options = " ".join(options)
    run_process(f"varifier make_truth_vcf {options} {assembly} {ref_genome} {outdir}")
    check_file(vcf)
    return vcf
