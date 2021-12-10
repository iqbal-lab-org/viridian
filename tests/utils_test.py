from collections import namedtuple
import filecmp
import os
import pytest
import subprocess
from unittest import mock

import pyfastaq

from viridian_workflow import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


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


def test_amplicons_json_to_bed_and_range():
    json_in = os.path.join(data_dir, "amplicons_json_to_bed.json")
    expect_bed = os.path.join(data_dir, "amplicons_json_to_bed.bed")
    tmp_out = "tmp.amplicons_json_to_bed.bed"
    subprocess.check_output(f"rm -f {tmp_out}", shell=True)
    got_start, got_end = utils.amplicons_json_to_bed_and_range(json_in, tmp_out)
    assert filecmp.cmp(tmp_out, expect_bed, shallow=False)
    assert got_start == 100
    assert got_end == 250
    os.unlink(tmp_out)


def test_load_amplicons_bed_file():
    Amplicon = namedtuple("Amplicon", ("name", "start", "end"))
    expect = [
        Amplicon("name1", 42, 99),
        Amplicon("name2", 85, 150),
    ]
    infile = os.path.join(data_dir, "load_amplicons_bed_file.bed")
    assert expect == utils.load_amplicons_bed_file(infile)


def test_set_sample_name_in_vcf_file():
    infile = os.path.join(data_dir, "set_sample_name_in_vcf_file.in.vcf")
    tmp_out = "tmp.set_sample_name_in_vcf_file.vcf"
    subprocess.check_output(f"rm -f {tmp_out}", shell=True)
    utils.set_sample_name_in_vcf_file(infile, tmp_out, "NEW_NAME")
    expect = os.path.join(data_dir, "set_sample_name_in_vcf_file.expect.vcf")
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)


def test_set_seq_name_in_fasta_file():
    infile = os.path.join(data_dir, "set_seq_name_in_fasta_file.in.fasta")
    tmp_out = "tmp.set_seq_name_in_fasta_file.fasta"
    subprocess.check_output(f"rm -f {tmp_out}", shell=True)
    utils.set_seq_name_in_fasta_file(infile, tmp_out, "NEW_NAME")
    expect = os.path.join(data_dir, "set_seq_name_in_fasta_file.expect.fasta")
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)


def test_check_tech_and_reads_opts_and_get_reads():
    f = utils.check_tech_and_reads_opts_and_get_reads
    options = mock.Mock()
    options.reads = options.reads1 = options.reads2 = options.tech = None
    with pytest.raises(Exception):
        f(options)
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
    assert f(options) == ("r.fq", None)

    options = mock.Mock()
    options.reads = options.reads1 = options.reads2 = options.tech = None
    with pytest.raises(Exception):
        f(options)
    options.tech = "illumina"
    with pytest.raises(Exception):
        f(options)
    options.reads = "r.fq"
    with pytest.raises(Exception):
        f(options)
    options.reads1 = "r1.fq"
    with pytest.raises(Exception):
        f(options)
    options.reads = None
    with pytest.raises(Exception):
        f(options)
    options.reads2 = "r2.fq"
    assert f(options) == ("r1.fq", "r2.fq")
