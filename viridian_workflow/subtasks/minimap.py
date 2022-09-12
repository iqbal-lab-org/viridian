"""minimap2 and samtools utilities
"""
from __future__ import annotations

import sys
import subprocess
from typing import Optional
from pathlib import Path
from .task import Task


class Minimap(Task):
    def __init__(
        self,
        bam: Path,
        ref_genome: Path,
        fq1: Path,
        fq2: Path = None,
        threads: int = 1,
        sample_name: str = "sample",
        minimap_x_opt: Optional[str] = None,
        sort: bool = True,
    ):

        self.output: Path = Path(bam)
        self.sort: bool = sort
        self.cmd: list[str] = [
            "minimap2",
            "-R",
            rf"@RG\tID:1\tSM:{sample_name}",
            "-t",
            str(threads),
            "-a",
        ]
        reads_list: list[str] = []
        if fq2 is None:
            if minimap_x_opt is None:
                minimap_x_opt = "-x map-ont"
            self.cmd.extend(minimap_x_opt.split())
            reads_list = [str(fq1)]
        else:
            if minimap_x_opt is None:
                minimap_x_opt = "-x sr"
            self.cmd.extend(minimap_x_opt.split())
            reads_list = [str(fq1), str(fq2)]
        self.cmd.append(str(ref_genome))
        self.cmd.extend(reads_list)
        super(Minimap, self).__init__(name="minimap")

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
            subprocess.Popen(["samtools", "index", self.output]).wait()

        else:
            with open(self.output, "w") as out_fd:
                print(f"running: {' '.join([str(c) for c in self.cmd])}", sys.stderr)
                map_proc = subprocess.Popen(self.cmd, stdout=out_fd)
                map_proc.wait()
                if map_proc.returncode:
                    raise Exception("minimap2 subprocess failed")

        self.check_output()
        self.log["Success"] = True
        return self.output
