import gzip
import json
import logging
import os
import subprocess

import pyfastaq

from viridian import constants


def syscall(command, stdouterr=None, quiet=False, silent=False, cwd=None):
    if not silent:
        logging.info(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        cwd=cwd,
    )
    completed_process.stdout = completed_process.stdout.rstrip()
    completed_process.stderr = completed_process.stderr.rstrip()

    if completed_process.returncode == 0:
        log_out = logging.info
    else:
        logging.error(f"Error running this command: {command}")
        logging.error(f"cwd: {os.getcwd()}")
        logging.error(f"Return code: {completed_process.returncode}")
        log_out = logging.error

    if completed_process.returncode != 0 or not (silent or quiet):
        if completed_process.stdout != "":
            log_out(f"stdout:\n{completed_process.stdout.rstrip()}")
        if completed_process.stderr != "":
            log_out(f"stderr:\n{completed_process.stderr.rstrip()}")

    if completed_process.returncode != 0:
        raise RuntimeError("Error in system call. Cannot continue")

    if not silent:
        logging.info(f"Return code {completed_process.returncode} from: {command}")

    if stdouterr is not None:
        with open(f"{stdouterr}.out", "w") as f:
            print(completed_process.stdout.rstrip(), file=f)
        with open(f"{stdouterr}.err", "w") as f:
            print(completed_process.stderr.rstrip(), file=f)

    return completed_process


def load_json(infile):
    if infile.endswith(".gz"):
        with gzip.open(infile, "rt") as f:
            return json.load(f)
    else:
        with open(infile) as f:
            return json.load(f)


def write_json(outfile, data):
    if outfile.endswith(".gz"):
        with gzip.open(outfile, "wt") as f:
            json.dump(data, f, indent=2)
    else:
        with open(outfile, "w") as f:
            json.dump(data, f, indent=2)


def load_single_seq_fasta(infile):
    d = {}
    pyfastaq.tasks.file_to_dict(infile, d)
    if len(d) != 1:
        raise Exception(
            f"Expected exactly 1 sequence in {infile} but got {len(d)} sequences"
        )
    ref = list(d.values())[0]
    ref.id = ref.id.split()[0]
    return ref


def seq_length_of_single_seq_fasta(infile):
    return len(load_single_seq_fasta(infile))


def check_tech_and_reads_options(args):
    if args.run_accession is not None:
        return True

    if getattr(args, "tech") is None:
        raise Exception("Must use --tech option. Cannot continue")

    if args.tech not in constants.ALLOWED_TECH:
        raise Exception(
            f"Sequencing tech '{args.tech}' not recognised. Cannot continue"
        )

    reads_bam = getattr(args, "reads_bam") is not None
    reads1 = getattr(args, "reads1") is not None
    reads2 = getattr(args, "reads2") is not None
    reads = getattr(args, "reads") is not None

    if reads_bam:
        if reads or reads1 or reads2:
            raise Exception(
                "--reads_bam provided, which is incompatible with --reads, --reads1, --reads2"
            )
        if args.decontam:
            raise Exception("--decontam is not compatible with --reads_bam")
        return True

    if args.tech == "ont":
        if not (reads and not reads1 and not reads2):
            raise Exception(
                "For ont tech, must provide --reads, and provide --reads1, --reads2"
            )
    elif args.tech in ["illumina", "iontorrent"]:
        ok1 = reads and not (reads1 or reads2)
        ok2 = reads1 and reads2 and not reads
        if not (ok1 or ok2):
            raise Exception(
                "For illumina/iontorrent tech, must provide either --reads, or alternatively provide both --reads1, --reads2"
            )
    else:
        raise NotImplementedError()

    return True


def write_fasta(name, seq_string, outfile):
    seq = pyfastaq.sequences.Fasta(name, seq_string)
    if outfile.endswith(".gz"):
        with gzip.open(outfile, "wt") as f:
            print(seq, file=f)
    else:
        with open(outfile, "w") as f:
            print(seq, file=f)


def gunzip_file(infile, outfile):
    with gzip.open(infile, "rt") as f_in, open(outfile, "w") as f_out:
        for line in f_in:
            print(line, end="", file=f_out)
