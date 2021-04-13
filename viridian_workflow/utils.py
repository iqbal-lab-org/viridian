import logging
import os
import subprocess


class OutputFileError(Exception):
    pass


class PipelineProcessError(Exception):
    pass


def check_file(fn):
    if not os.path.isfile(fn):
        raise OutputFileError(os.path.abspath(fn))


def rm(fn):
    subprocess.run(f"rm {fn}")


def run_process(cmd, ignore_error=False):
    logging.info(f"Running: {cmd}")
    result = subprocess.run(
        cmd,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    logging.info(result.stdout)
    if not ignore_error and result.returncode != 0:
        raise PipelineProcessError(f"Process returned {result.returncode}")
        logging.error(result.stderr)
