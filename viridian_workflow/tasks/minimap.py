"""minimap2 and samtools utilities
"""
import subprocess
from viridian_workflow.utils import run_process, check_file


class Minimap(Task):
    def __init__(self):
        pass

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
        minimap_cmd = [
            "minimap2",
            "-R",
            fr"@RG\tID:1\tSM:{sample_name}",
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
            if sort_proc.returncode:
                raise Exception("minimap2 subprocess failed")

        else:
            with open(bam, "w") as out_fd:
                map_proc = subprocess.Popen(minimap_cmd, stdout=out_fd)
                map_proc.wait()
                if map_proc.returncode:
                    raise Exception("minimap2 subprocess failed")

        check_file(bam)

        if sort:
            run_process(["samtools", "index", bam])
            check_file(bam + ".bai")

        return bam
