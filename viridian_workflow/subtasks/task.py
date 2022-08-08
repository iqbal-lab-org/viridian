"""
External pipeline process wrapper
"""
from __future__ import annotations

import subprocess
import time
from typing import Union
from pathlib import Path


class Task:
    """A prototype Task

    Tasks are associated with external process invocations and return
    a list of output files
    """

    def __init__(self):
        self.cmd: list[str]
        self.log: dict[str, str]
        self.output: Union[Path, list[Path]]

    def check_output(self):
        """Test if outfile files exist"""
        if isinstance(self.output, Path):
            if not self.output.exists():
                raise Exception(f"Failed to generate output file:\t{self.output}")
        else:
            for fn in self.output:
                if not fn.exists():
                    raise Exception(f"Failed to generate output file:\t{fn}")

    def run(self, ignore_error=False, stdout=None):
        """Launch a pipeline task"""
        self.log = {}
        cmd = [str(c) for c in self.cmd]
        self.log["subprocess"] = " ".join(cmd)

        stdout_fd = subprocess.PIPE
        if stdout:
            stdout_fd = open(stdout, "w", encoding="utf-8")
        start_time = time.time()
        result = subprocess.run(
            cmd,
            shell=False,
            stderr=subprocess.PIPE,
            stdout=stdout_fd,
            universal_newlines=True,
        )
        if stdout:
            stdout_fd.close()
        time_elapsed = time.time() - start_time
        self.log["duration"] = str(time_elapsed)

        if not ignore_error and result.returncode != 0:
            raise Exception(f"Process returned {result.returncode}: {result.stderr}")

        # this is awkward and probably not useful
        if not stdout:
            print(result.stdout)

        self.check_output()
        return self.output
