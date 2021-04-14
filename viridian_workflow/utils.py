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
    logging.info(f"Deleting file: {os.path.abspath(fn)}")
    subprocess.run(f"rm {os.path.abspath(fn)}")


def run_process(cmd, ignore_error=False, stdout=None):
    logging.info(f"Running: {cmd}")
    stdout_fd = subprocess.PIPE
    if stdout:
        stdout_fd = open(stdout, "w")
    result = subprocess.run(
        cmd,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=stdout_fd,
        universal_newlines=True,
    )
    if not stdout:
        logging.info(result.stdout)
    if not ignore_error and result.returncode != 0:
        raise PipelineProcessError(f"Process returned {result.returncode}")
        logging.error(result.stderr)
