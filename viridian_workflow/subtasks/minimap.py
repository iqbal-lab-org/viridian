"""minimap2 and samtools utilities
"""
import sys
import subprocess
from viridian_workflow.utils import run_process, check_file
from .task import Task


class Minimap(Task):
    def __init__(
        self,
        bam,
        ref_genome,
        fq1,
        fq2=None,
        threads=1,
        sample_name="sample",
        minimap_x_opt=None,
        sort=True,
    ):

        self.output = bam
        self.sort = sort
        self.cmd = [
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
            self.cmd.extend(minimap_x_opt.split())
            reads_list = [fq1]
        else:
            if minimap_x_opt is None:
                minimap_x_opt = "-x sr"
            self.cmd.extend(minimap_x_opt.split())
            reads_list = [fq1, fq2]
        self.cmd.append(ref_genome)
        self.cmd.extend(reads_list)

    def run(self):
        if self.sort:
            sort_cmd = ["samtools", "sort", "-o", self.output]
            map_proc = subprocess.Popen(
                " ".join(self.cmd), shell=True, stdout=subprocess.PIPE
            )
            sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
            sort_proc.wait()
            if sort_proc.returncode:
                raise Exception("minimap2 subprocess failed")

        else:
            with open(self.output, "w") as out_fd:
                print(f"running: {' '.join([str(c) for c in self.cmd])}", sys.stderr)
                map_proc = subprocess.Popen(self.cmd, stdout=out_fd)
                map_proc.wait()
                if map_proc.returncode:
                    raise Exception("minimap2 subprocess failed")

        check_file(self.output)

        if self.sort:
            run_process(["samtools", "index", self.output])
            check_file(self.output + ".bai")

        return self.output
