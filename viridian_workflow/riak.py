"""Read it and keep wrapper
"""
import subprocess
from viridian_workflow.utils import run_process_safe, check_file


def riak_se(outdir, ref_genome, fq):
    pass


def riak_pe(outdir, ref_genome, fq1, fq2):
    sample = fq1.stem()
    riak = "readItAndKeep"
    run_process_safe(
        [
            riak,
            "--ref_fasta",
            ref_genome,
            "--reads1",
            fq1,
            "--reads2",
            fq2,
            "--outprefix",
            outdir / sample,
        ]
    )
    check_file(outdir / f"{sample}.clean.fq.gz")
    return mask
