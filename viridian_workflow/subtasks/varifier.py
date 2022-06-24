"""Varifier wrapper
"""
import os
from .task import Task


class Varifier(Task):
    def __init__(self, outdir, ref, consensus, min_coord=0, max_coord=None, sanitise_gaps=True):
        vcf = os.path.join(outdir, "04.truth.vcf")
        msa = os.path.join(outdir, "04.msa")
        consensus_out = os.path.join(outdir, "04.qry_sanitised_gaps.fa")
        self.output = (vcf, msa, consensus_out)

        self.options = ["--global_align"]
        if min_coord is not None:
            self.options += ["--global_align_min_coord", min_coord + 1]
        if max_coord is not None:
            self.options += ["--global_align_max_coord", max_coord + 1]
        self.cmd = [
            "varifier",
            "make_truth_vcf",
            ]

        if sanitise_gaps:
            self.cmd.append("--sanitise_truth_gaps")

        self.cmd += [
            *self.options,
            consensus,
            ref,
            outdir,
        ]
