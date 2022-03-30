"""varifier wrapper
"""
import os
from viridian_workflow.utils import run_process, check_file


def run(outdir, ref_genome, assembly, min_coord=None, max_coord=None):
    vcf = os.path.join(outdir, "04.truth.vcf")
    options = ["--global_align"]
    if min_coord is not None:
        options += ["--global_align_min_coord", min_coord + 1]
    if max_coord is not None:
        options += ["--global_align_max_coord", max_coord + 1]

    run_process(
        [
            "varifier",
            "make_truth_vcf",
            "--sanitise_truth_gaps",
            *options,
            assembly,
            ref_genome,
            outdir,
        ],
        ignore_error=True,
    )
    check_file(vcf)
    return vcf
