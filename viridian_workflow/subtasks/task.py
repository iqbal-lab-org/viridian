import sys

from viridian_workflow.utils import run_process, check_file


class Task:
    def __init__(self):
        pass

    def run(self):
        self.log = {}
        #        print(
        #            f"Running subprocess: {' '.join([str(c) for c in self.cmd])}",
        #            file=sys.stderr,
        #        )
        self.log["subprocess"] = " ".join([str(c) for c in self.cmd])
        run_process(self.cmd)
        if type(self.output) == type(tuple()):
            for i in self.output:
                check_file(i)
        else:
            check_file(self.output)
        return self.output
