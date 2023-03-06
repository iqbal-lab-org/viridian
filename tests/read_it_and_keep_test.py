import os
import pytest
import random

import pyfastaq

from viridian import read_it_and_keep, utils


def test_read_it_and_keep():
    outprefix = "tmp.riak"
    utils.syscall(f"rm -f {outprefix}*")
    ref_fasta = f"{outprefix}.ref.fa"
    reads1 = f"{outprefix}.reads.1.fa"
    reads2 = f"{outprefix}.reads.2.fa"
    random.seed(42)
    ref_bases = "".join(random.choices(["A", "C", "G", "T"], k=1000))
    with open(ref_fasta, "w") as f:
        print(">ref", ref_bases, sep="\n", file=f)
    with open(reads1, "w") as f:
        print(">read1/1", ref_bases[100:300], sep="\n", file=f)
        print(">read2/1", "A" * 200, sep="\n", file=f)
    with open(reads2, "w") as f:
        read = pyfastaq.sequences.Fasta("read1/2", ref_bases[300:500])
        read.revcomp()
        print(read, file=f)
        print(">read2/2", "A" * 200, sep="\n", file=f)

    expect_counts = {
        "Input reads file 1": 2,
        "Input reads file 2": 0,
        "Kept reads 1": 1,
        "Kept reads 2": 0,
    }
    expect_files = {
        "reads_1": f"{outprefix}.reads.fasta.gz",
        "reads_2": None,
    }
    got_counts, got_files, got_error = read_it_and_keep.run_riak(
        outprefix, ref_fasta, reads1, reads2=None
    )
    assert got_error is None
    assert got_counts == expect_counts
    assert got_files == expect_files
    os.unlink(expect_files["reads_1"])

    expect_files = {
        "reads_1": f"{outprefix}.reads_1.fasta.gz",
        "reads_2": f"{outprefix}.reads_2.fasta.gz",
    }
    expect_counts["Input reads file 2"] = 2
    expect_counts["Kept reads 2"] = 1
    got_counts, got_files, got_error = read_it_and_keep.run_riak(
        outprefix, ref_fasta, reads1, reads2=reads2
    )
    assert got_counts == expect_counts
    assert got_files == expect_files
    assert got_error is None
    utils.syscall(f"rm {outprefix}*")

    got_counts, got_files, got_error = read_it_and_keep.run_riak(
        outprefix, "not_a_file", "also_not_a_file"
    )
    assert got_counts == {}
    assert got_files == {}
    assert got_error == "Error running readItAndKeep"
