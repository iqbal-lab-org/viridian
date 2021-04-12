"""minimap2 and samtools utilities
"""
import os
import subprocess
from viridian_workflow.utils import run_process, check_file


def run(outdir, ref_genome, fq1, fq2):
    bam = os.path.join(outdir, "reference_mapped.bam")
    minimap_cmd = ["minimap2", "-ax", "sr", ref_genome, fq1, fq2]
    sort_cmd = ["samtools", "sort", "-o", bam]

    map_proc = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
    sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
    sort_proc.wait()
    check_file(bam)

    run_process(f"samtools index {bam}")
    check_file(bam + ".bai")
    return bam
