from collections import namedtuple
import json
import logging
from operator import itemgetter
import os
import subprocess
import time

import pyfastaq


class OutputFileError(Exception):
    pass


class PipelineProcessError(Exception):
    pass


def check_file(fn):
    if not os.path.isfile(fn):
        raise OutputFileError(os.path.abspath(fn))


def rm(fn):
    logging.info(f"Deleting file: {os.path.abspath(fn)}")
    os.remove(os.path.abspath(fn))


def run_process(cmd, ignore_error=False, stdout=None):
    logging.info(f"Running: {cmd}")
    stdout_fd = subprocess.PIPE
    if stdout:
        stdout_fd = open(stdout, "w")
    start_time = time.time()
    result = subprocess.run(
        cmd,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=stdout_fd,
        universal_newlines=True,
    )
    time_elapsed = time.time() - start_time
    logging.info(f"Process ({cmd}) completed in {time_elapsed} seconds.")
    if not stdout:
        logging.info(result.stdout)
    if not ignore_error and result.returncode != 0:
        raise PipelineProcessError(f"Process returned {result.returncode}")
        logging.error(result.stderr)


def amplicons_json_to_bed(infile, outfile):
    with open(infile) as f:
        data = json.load(f)

    data_out = []
    for amplicon, d in data["amplicons"].items():
        data_out.append((amplicon, d["start"], d["end"] + 1))
    data_out.sort(key=itemgetter(1))
    with open(outfile, "w") as f:
        for t in data_out:
            print(*t, sep="\t", file=f)


def load_amplicons_bed_file(infile):
    Amplicon = namedtuple("Amplicon", ("name", "start", "end"))
    amplicons = []

    with open(infile) as f:
        for line in f:
            if line.startswith("#"):
                continue
            name, start, end = line.rstrip().split("\t")
            amplicons.append(Amplicon(name, int(start), int(end) - 1))

    return amplicons


def load_single_seq_fasta(infile):
    d = {}
    pyfastaq.tasks.file_to_dict(infile, d)
    if len(d) != 1:
        raise Exception(
            f"Expected exatcly 1 sequence in {infile} but got {len(d)} sequences"
        )
    ref = list(d.values())[0]
    ref.id = ref.id.split()[0]
    return ref
