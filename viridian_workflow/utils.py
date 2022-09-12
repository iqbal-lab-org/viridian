"""
Helper functions shared across modules
"""
from __future__ import annotations

import sys
from typing import NewType, Any, Optional
from collections import namedtuple
import json
import logging
from operator import itemgetter
from pathlib import Path

import pyfastaq  # type: ignore


class PrimerError(Exception):
    """Exception raised by malformed Primer record"""


class OutputFileError(Exception):
    """Output file cannot be written to or doesn't exist"""


class PipelineProcessError(Exception):
    """Pipeline subprocess error"""


Index0 = NewType("Index0", int)
Index1 = NewType("Index1", int)

TRANSLATE_TABLE = str.maketrans("ATCGatcg", "TAGCtagc")


def revcomp(seq: str) -> str:
    """Reverse complement sequence"""
    return seq.translate(TRANSLATE_TABLE)[::-1]


def check_file(filename: Path) -> bool:
    """Check if file exists"""
    if not filename.exists():
        raise OutputFileError(filename)
    return True


def in_range(interval: tuple[Index0, Index0], position: Index0) -> bool:
    """Test whether a position is inside a range"""
    start, end = interval
    return start <= position < end


def rm(filename: Path):
    """File removal wrapper"""
    filename = filename.resolve()
    logging.info("Deleting file %s", filename)
    filename.unlink()


def amplicons_json_to_bed_and_range(
    infile: Path, outfile: Path
) -> tuple[Index0, Index0]:
    """Converts the amplicons JSON file to a BED file, which is used in
    various stages of the pipeline. Returns the 0-based inclusive coordinates
    of the start of the first amplicon and end of the last amplicon, as
    a tuple (start, end)"""
    data = json.load(open(infile, encoding="utf-8"))
    start = Index0(sys.maxsize)
    end = Index0(0)

    data_out = []
    for amplicon, d in data["amplicons"].items():
        data_out.append((amplicon, d["start"], d["end"] + 1))
        start = Index0(min(start, d["start"]))
        end = Index0(max(end, d["end"]))
    data_out.sort(key=itemgetter(1))
    with open(outfile, "w", encoding="utf-8") as f:
        for t in data_out:
            print(*t, sep="\t", file=f)
    return start, end


Amplicon = namedtuple("Amplicon", ("name", "start", "end"))


def load_amplicons_bed_file(infile: Path) -> list[Amplicon]:
    """Load 0-based amplicon coords from bedfile"""
    amplicons: list[Amplicon] = []

    with open(infile, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            name, start, end = line.rstrip().split("\t")
            amplicons.append(Amplicon(name, Index0(int(start)), Index0(int(end) - 1)))

    return amplicons


def load_single_seq_fasta(infile: Path) -> Any:
    """Load single fastaq sequence from a fasta"""
    seq_dict: Any = {}
    pyfastaq.tasks.file_to_dict(infile, seq_dict)
    if len(seq_dict) != 1:
        raise Exception(
            f"Expected exactly 1 sequence in {infile} but got {len(seq_dict)} sequences"
        )
    ref = list(seq_dict.values())[0]
    ref.id = ref.id.split()[0]
    return ref


def set_sample_name_in_vcf_file(infile: Path, outfile: Path, sample_name: str):
    """Rename sample in vcf and write out to new file"""
    with open(infile, encoding="utf-8") as f_in, open(
        outfile, "w", encoding="utf-8"
    ) as f_out:
        for line in f_in:
            if line.startswith("#CHROM\tPOS"):
                fields = line.rstrip().split("\t")
                fields[-1] = sample_name
                print(*fields, sep="\t", file=f_out)
            else:
                print(line, end="", file=f_out)


def set_seq_name_in_fasta_file(infile: Path, outfile: Path, new_name: str):
    """Changes name in FASTA file. Assumes that there is only one sequence
    in the file, raises Exception if it finds >1 sequence"""
    seq = load_single_seq_fasta(infile)
    seq.id = new_name
    with open(outfile, "w", encoding="utf-8") as f:
        print(seq, file=f)


def check_tech_and_reads_opts_and_get_reads(
    options: Any,
) -> tuple[Path, Optional[Path]]:
    """Normalise reads passed through options"""
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
