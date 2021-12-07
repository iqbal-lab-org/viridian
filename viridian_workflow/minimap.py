"""minimap2 and samtools utilities
"""
import os
import subprocess
from viridian_workflow.utils import run_process, check_file


def run(
    bam,
    ref_genome,
    fq1,
    fq2=None,
    threads=1,
    sample_name="sample",
    minimap_x_opt=None,
    sort=True,
):
    # Note: was getting issues with escaping the tab characters in the -R
    # option. In the end got it to work by giving Popen a string instead of
    # a list and setting shell=True (using a list and shell=False failed using
    # exactly the same minimap_cmd contents).
    minimap_cmd = [
        "minimap2",
        f"-R '@RG\\tID:1\\tSM:{sample_name}'",
        "-t",
        str(threads),
        "-a",
    ]
    if fq2 is None:
        if minimap_x_opt is None:
            minimap_x_opt = "-x map-ont"
        minimap_cmd.extend(minimap_x_opt.split())
        reads_list = [fq1]
    else:
        if minimap_x_opt is None:
            minimap_x_opt = "-x sr"
        minimap_cmd.extend(minimap_x_opt.split())
        reads_list = [fq1, fq2]
    minimap_cmd.append(ref_genome)
    minimap_cmd.extend(reads_list)

    if sort:
        sort_cmd = ["samtools", "sort", "-o", bam]
        map_proc = subprocess.Popen(
            " ".join(minimap_cmd), shell=True, stdout=subprocess.PIPE
        )
        sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
        sort_proc.wait()
    else:
        minimap_cmd.extend([">", bam])
        subprocess.check_output(" ".join(minimap_cmd), shell=True)

    check_file(bam)

    if sort:
        run_process(f"samtools index {bam}")
        check_file(bam + ".bai")

    return bam
