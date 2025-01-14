import os
import pytest
from unittest import mock

import pyfastaq

from viridian import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_load_write_json():
    json_file = "tmp.load_write.json"
    utils.syscall(f"rm -f {json_file}")
    d = {"a": "b"}
    utils.write_json(json_file, d)
    assert utils.load_json(json_file) == d
    os.unlink(json_file)


def test_load_write_json_gz():
    json_file = "tmp.load_write_gz.json.gz"
    utils.syscall(f"rm -f {json_file}")
    d = {"a": "b"}
    utils.write_json(json_file, d)
    assert utils.load_json(json_file) == d
    os.unlink(json_file)


def test_load_single_seq_fasta():
    expect = pyfastaq.sequences.Fasta("seq", "ACGT")
    infile = os.path.join(data_dir, "load_single_seq_fasta.ok.fa")
    assert expect == utils.load_single_seq_fasta(infile)
    infile = os.path.join(data_dir, "load_single_seq_fasta.bad1.fa")
    with pytest.raises(Exception):
        utils.load_single_seq_fasta(infile)
    infile = os.path.join(data_dir, "load_single_seq_fasta.bad2.fa")
    with pytest.raises(Exception):
        utils.load_single_seq_fasta(infile)


def test_seq_length_of_single_seq_fasta():
    infile = os.path.join(data_dir, "load_single_seq_fasta.ok.fa")
    assert utils.seq_length_of_single_seq_fasta(infile) == 4


def test_check_tech_and_reads_options():
    f = utils.check_tech_and_reads_options
    options = mock.Mock()
    options.reads_bam = None
    options.reads = None
    options.reads1 = None
    options.reads2 = None
    options.tech = None
    options.run_accession = None
    with pytest.raises(Exception):
        f(options)

    options.run_accession = "ERR12345"
    assert f(options)
    options.run_accession = None

    options.tech = "ont"
    with pytest.raises(Exception):
        f(options)
    options.reads1 = "r1.fq"
    with pytest.raises(Exception):
        f(options)
    options.reads1 = None
    options.reads2 = "r2.fq"
    with pytest.raises(Exception):
        f(options)
    options.reads = "r.fq"
    with pytest.raises(Exception):
        f(options)
    options.reads2 = None
    assert f(options)

    for tech in ["illumina", "iontorrent"]:
        options = mock.Mock()
        options.reads_bam = None
        options.reads = None
        options.reads1 = None
        options.reads2 = None
        options.tech = None
        options.run_accession = None
        with pytest.raises(Exception):
            f(options)
        options.tech = tech
        with pytest.raises(Exception):
            f(options)
        options.reads = "r.fq"
        assert f(options)
        options.reads1 = "r1.fq"
        with pytest.raises(Exception):
            f(options)
        options.reads = None
        with pytest.raises(Exception):
            f(options)
        options.reads2 = "r2.fq"
        assert f(options)
