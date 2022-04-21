"""Varifier wrapper
"""
import os
from viridian_workflow.utils import run_process, check_file


class Varifier(Task):
    def __init__(
        self, outdir, ref, consensus, min_coord=rs.start_pos, max_coord=rs.end_pos
    ):
        vcf = os.path.join(outdir, "04.truth.vcf")
        msa = os.path.join(outdir, "04.msa")
        self.output = (vcf, msa)
        self.outdir = outdir
        self.ref = ref
        self.consensus = consensus

        self.options = ["--global_align"]
        if min_coord is not None:
            self.options += ["--global_align_min_coord", min_coord + 1]
        if max_coord is not None:
            self.options += ["--global_align_max_coord", max_coord + 1]

    def run():

        run_process(
            [
                "varifier",
                "make_truth_vcf",
                "--sanitise_truth_gaps",
                *self.options,
                assembly,
                self.ref,
                self.outdir,
            ],
        )
        vcf, msa = self.output
        check_file(vcf)
        check_file(msa)
        return self.output
