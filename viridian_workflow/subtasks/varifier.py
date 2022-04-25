"""Varifier wrapper
"""
import os
from viridian_workflow.utils import run_process, check_file
from .task import Task


class Varifier(Task):
    def __init__(self, outdir, ref, consensus, min_coord=0, max_coord=None):
        vcf = os.path.join(outdir, "04.truth.vcf")
        msa = os.path.join(outdir, "04.msa")
        self.output = (vcf, msa, consensus)

        self.options = ["--global_align"]
        if min_coord is not None:
            self.options += ["--global_align_min_coord", min_coord + 1]
        if max_coord is not None:
            self.options += ["--global_align_max_coord", max_coord + 1]
        self.cmd = [
            "varifier",
            "make_truth_vcf",
            "--sanitise_truth_gaps",
            *self.options,
            assembly,
            ref,
            outdir,
        ]
