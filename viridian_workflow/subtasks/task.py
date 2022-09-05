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

    def __init__(self, name=None):
        self.cmd: list[str]
        self.output: Union[Path, list[Path]]

        if name is None:
            self.name = self.cmd[0]
        else:
            self.name = name

        self.start_time = time.time()
        self.log: dict[str, str] = {
            "Task": self.name,
            "Success": False,
            "start": self.start_time,
            "error": None,
        }

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
        cmd = [str(c) for c in self.cmd]

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
            check=True,
        )
        if stdout:
            stdout_fd.close()

        self.log["end"] = time.time()

        if not ignore_error and result.returncode != 0:
            self.log["error"] = result.stderr
            raise Exception(f"Process returned {result.returncode}: {result.stderr}")

        # this is awkward and probably not useful
        if not stdout:
            print(result.stdout)

        self.check_output()
        self.log["Success"] = True
        return self.output
