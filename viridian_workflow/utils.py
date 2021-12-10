from collections import namedtuple
import json
import logging
from operator import itemgetter
import os
import subprocess
import time

import pyfastaq


class PrimerError(Exception):
    pass


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


def load_json(infile):
    with open(infile) as f:
        return json.load(f)


def write_json(outfile, data):
    with open(outfile, "w") as f:
        json.dump(data, f, indent=2)


def run_process_stdout(cmd, ignore_error=False):
    logging.info(f"Running: {cmd}")
    print(" ".join(cmd))

    start_time = time.time()
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )

    time_elapsed = time.time() - start_time
    logging.info(f"Process ({cmd}) completed in {time_elapsed} seconds.")
    logging.info(result.stdout)

    if not ignore_error and result.returncode != 0:
        raise PipelineProcessError(f"Process returned {result.returncode}")
        logging.error(result.stderr)

    return result.stdout


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
    if stdout:
        stdout_fd.close()
    time_elapsed = time.time() - start_time
    logging.info(f"Process ({cmd}) completed in {time_elapsed} seconds.")

    if not stdout:
        logging.info(result.stdout)
        return result.stdout

    if not ignore_error and result.returncode != 0:
        raise PipelineProcessError(f"Process returned {result.returncode}")
        logging.error(result.stderr)


def amplicons_json_to_bed_and_range(infile, outfile):
    """Converts the amplicons JSON file to a BED file, which is used in
    various stages of the pipeline. Returns the 0-based inclusive coordinates
    of the start of the first amplicon and end of the last amplicon, as
    a tuple (start, end)"""
    data = load_json(infile)
    start = float("inf")
    end = -1

    data_out = []
    for amplicon, d in data["amplicons"].items():
        data_out.append((amplicon, d["start"], d["end"] + 1))
        start = min(start, d["start"])
        end = max(end, d["end"])
    data_out.sort(key=itemgetter(1))
    with open(outfile, "w") as f:
        for t in data_out:
            print(*t, sep="\t", file=f)
    return start, end


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


def set_sample_name_in_vcf_file(infile, outfile, sample_name):
    with open(infile) as f_in, open(outfile, "w") as f_out:
        for line in f_in:
            if line.startswith("#CHROM\tPOS"):
                fields = line.rstrip().split("\t")
                fields[-1] = sample_name
                print(*fields, sep="\t", file=f_out)
            else:
                print(line, end="", file=f_out)


def set_seq_name_in_fasta_file(infile, outfile, new_name):
    """Changes name in FASTA file. Assumes that there is only one sequence
    in the file, raises Exception if it finds >1 sequence"""
    seq = load_single_seq_fasta(infile)
    seq.id = new_name
    with open(outfile, "w") as f:
        print(seq, file=f)


def check_tech_and_reads_opts_and_get_reads(options):
    if options.tech == "ont":
        if options.reads1 is not None or options.reads2 is not None:
            raise Exception(
                "Tech is 'ont'. Cannot use reads1 or reads2 options, use reads option instead"
            )
        if options.reads is None:
            raise Exception("Tech is 'ont'. Must use reads option")
        fq1 = options.reads
        fq2 = None
    elif options.tech == "illumina":
        if options.reads is not None:
            raise Exception(
                "Tech is 'illumina'. Cannot use reads option, use both reads1 and reads2 instead"
            )
        if options.reads1 is None or options.reads2 is None:
            raise Exception(
                "Tech is 'illumina'. Must use both reads1 and reads2 instead"
            )
        fq1 = options.reads1
        fq2 = options.reads2
    return fq1, fq2
