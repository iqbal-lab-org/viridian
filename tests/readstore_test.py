import filecmp
import json
import os
import pytest
import subprocess
from unittest import mock

from viridian_workflow import primers, readstore

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "readstore")


class Bam:
    def __init__(self, infile_is_paired=None):
        self.infile_is_paired = infile_is_paired

    def syncronise_fragments(self):
        for _ in []:
            yield


def test_make_reads_dir_for_cylon_unpaired():
    # This amplicon set is 3 amplicons, each has length 200.
    # We'll make a fwd and reverse read, each 200 long., and just make copies of
    # them for the lists of fragments from which to sample reads
    amplicons_tsv = os.path.join(data_dir, "make_reads_dir_for_cylon.amplicons.tsv")
    amplicon_set = primers.AmpliconSet.from_tsv(amplicons_tsv)
    amplicons = list(amplicon_set)
    assert len(amplicons) == 3
    read_store = readstore.ReadStore(amplicon_set, Bam())
    read_store.reads_all_paired = False
    read_fwd = readstore.Read("A" * 100, 100, 199, 0, 99, False)
    read_rev = readstore.Read("C" * 100, 101, 200, 1, 98, True)
    frag_fwd = readstore.SingleRead(read_fwd)
    frag_rev = readstore.SingleRead(read_rev)
    read_store.amplicons = {
        amplicons[0]: [frag_fwd, frag_rev],
        amplicons[1]: [],
        amplicons[2]: [frag_fwd, frag_rev] * 100,
    }
    outdir = "tmp.make_reads_dir_for_cylon_unpaired"
    subprocess.check_output(f"rm -rf {outdir}", shell=True)
    manifest = read_store.make_reads_dir_for_cylon(outdir)
    got_failed = read_store.failed_amplicons
    assert got_failed == {amplicons[1]}
    expect_manifest = {
        "amp1": "0.fa",
        "amp2": None,
        "amp3": "1.fa",
    }
    assert manifest == expect_manifest
    got_amp1 = os.path.join(outdir, "0.fa")
    got_amp3 = os.path.join(outdir, "1.fa")
    expect_amp1 = os.path.join(data_dir, "make_reads_dir_for_cylon.unpaired.0.fa")
    expect_amp3 = os.path.join(data_dir, "make_reads_dir_for_cylon.unpaired.1.fa")
    assert filecmp.cmp(got_amp1, expect_amp1, shallow=False)
    # TODO isn't this supposed to fail?
    #    assert filecmp.cmp(got_amp3, expect_amp3, shallow=False)
    subprocess.check_output(f"rm -rf {outdir}", shell=True)


def test_make_reads_dir_for_cylon_paired():
    amplicons_tsv = os.path.join(data_dir, "make_reads_dir_for_cylon.amplicons.tsv")
    amplicon_set = primers.AmpliconSet.from_tsv(amplicons_tsv)
    amplicons = list(amplicon_set)
    amplicon_set_fn_dict = {}
    assert len(amplicons) == 3
    read_store = readstore.ReadStore(amplicon_set, Bam())
    read_store.cylon_target_depth_factor = 50
    read_fwd = readstore.Read("A" * 100, 100, 199, 0, 99, False)
    read_rev = readstore.Read("C" * 100, 101, 200, 1, 98, True)
    frag = readstore.PairedReads(read_fwd, read_rev)
    read_store.reads_all_paired = True
    read_store.amplicons = {
        amplicons[0]: [frag] * 1,
        amplicons[1]: [],
        amplicons[2]: [frag] * 100,
    }
    outdir = "tmp.make_reads_dir_for_cylon_paired"
    subprocess.check_output(f"rm -rf {outdir}", shell=True)
    manifest = read_store.make_reads_dir_for_cylon(outdir)
    expected_failures = set(["amp2",])
    assert set(map(lambda a: a.name, read_store.failed_amplicons)) == expected_failures
    expect_manifest = {
        #        "amp1": "0.fa", # ?????
        "amp3": "1.fa",
    }
    for k in expect_manifest:
        assert manifest[k] == expect_manifest[k]

    got_amp1 = os.path.join(outdir, "0.fa")
    got_amp3 = os.path.join(outdir, "1.fa")
    expect_amp1 = os.path.join(data_dir, "make_reads_dir_for_cylon.paired.0.fa")
    expect_amp3 = os.path.join(data_dir, "make_reads_dir_for_cylon.paired.1.fa")

    assert filecmp.cmp(got_amp1, expect_amp1, shallow=False)
    # TODO what is this testing?
    #    assert filecmp.cmp(got_amp3, expect_amp3, shallow=False)

    subprocess.check_output(f"rm -rf {outdir}", shell=True)
