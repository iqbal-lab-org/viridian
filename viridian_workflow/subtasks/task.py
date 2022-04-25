import sys


class Task:
    log = {}
    output = None

    def __init__(self):
        pass

    def run(self):
        print(
            f"Running subprocess: {' '.join([str(c) for c in self.cmd])}",
            file=sys.stderr,
        )
        utils.run_process(self.cmd)
        if type(self.output) == type(tuple()):
            for i in self.output:
                utils.check_file(i)
        else:
            utils.check_file(self.output)
        return self.output
