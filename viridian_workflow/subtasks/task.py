import subprocess
import time
import os


class Task:
    def __init__(self):
        pass

    @staticmethod
    def _check_file(fn):
        if not os.path.isfile(fn):
            raise Exception(f"Failed to generate output file:\t{fn}")

    def run(self, ignore_error=False, stdout=None):
        self.log = {}
        cmd = [str(c) for c in self.cmd]
        self.log["subprocess"] = " ".join(cmd)

        stdout_fd = subprocess.PIPE
        if stdout:
            stdout_fd = open(stdout, "w")
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

        if isinstance(self.output, tuple):
            for i in self.output:
                Task._check_file(i)
        else:
            Task._check_file(self.output)
        return self.output
