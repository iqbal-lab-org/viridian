"""minimap2 and samtools utilities
"""
import os
import subprocess
from viridian_workflow.utils import run_process, check_file


def run_se(outdir, ref_genome, fq, prefix=None, threads=1, sample_name="sample"):
    """map single fastq ONT reads with minimap2
    """
    bam = os.path.join(outdir, "reference_mapped.bam")
    if prefix:
        bam = os.path.join(outdir, f"{prefix}-reference_mapped.bam")

    # Note: was getting issues with escaping the tab characters in the -R
    # option. In the end got it to work by giving Popen a string instead of
    # a list and setting shell=True (using a list and shell=False failed using
    # exactly the same minimap_cmd contents).
    minimap_cmd = ["minimap2", f"-R '@RG\\tID:1\\tSM:{sample_name}'", "-t", str(threads), "-ax", "map-ont", ref_genome, fq]
    sort_cmd = ["samtools", "sort", "-o", bam]

    map_proc = subprocess.Popen(" ".join(minimap_cmd), shell=True, stdout=subprocess.PIPE)
    sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
    sort_proc.wait()
    check_file(bam)

    run_process(f"samtools index {bam}")
    check_file(bam + ".bai")
    return bam


def run(outdir, ref_genome, fq1, fq2, prefix=None, threads=1, sample_name="sample"):
    bam = os.path.join(outdir, "reference_mapped.bam")
    if prefix:
        bam = os.path.join(outdir, f"{prefix}-reference_mapped.bam")

    minimap_cmd = ["minimap2", f"-R '@RG\\tID:1\\tSM:{sample_name}'", "-t", str(threads), "-ax", "sr", ref_genome, fq1, fq2]
    sort_cmd = ["samtools", "sort", "-o", bam]

    map_proc = subprocess.Popen(" ".join(minimap_cmd), shell=True, stdout=subprocess.PIPE)
    sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
    sort_proc.wait()
    check_file(bam)

    run_process(f"samtools index {bam}")
    check_file(bam + ".bai")
    return bam
