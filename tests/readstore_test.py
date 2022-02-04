import filecmp
import json
import os
import pytest
import subprocess
from unittest import mock

from viridian_workflow import primers, readstore

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "readstore")


def test_make_reads_dir_for_viridian_unpaired():
    # This amplicon set is 3 amplicons, each has length 200.
    # We'll make a fwd and reverse read, each 200 long., and just make copies of
    # them for the lists of fragments from which to sample reads
    amplicons_tsv = os.path.join(data_dir, "make_reads_dir_for_viridian.amplicons.tsv")
    amplicon_set = primers.AmpliconSet.from_tsv(amplicons_tsv)
    amplicons = list(amplicon_set)
    assert len(amplicons) == 3
    read_store = readstore.ReadStore("name", amplicon_set)
    read_store.reads_all_paired = False
    read_fwd = readstore.Read("A"*100, 100, 199, 0, 99, False)
    read_rev = readstore.Read("C"*100, 101, 200, 1, 98, True)
    frag_fwd = readstore.SingleRead(read_fwd)
    frag_rev = readstore.SingleRead(read_rev)
    read_store.amplicons = {
        amplicons[0]: [frag_fwd, frag_rev],
        amplicons[1]: [],
        amplicons[2]: [frag_fwd, frag_rev]*100,
    }
    outdir = "tmp.make_reads_dir_for_viridian_unpaired"
    subprocess.check_output(f"rm -rf {outdir}", shell=True)
    got_failed = read_store.make_reads_dir_for_viridian(outdir, 3)
    assert got_failed == {amplicons[1]}
    expect_manifest = {
        "amp1": "0.fa",
        "amp3": "1.fa",
    }
    manifest_json = os.path.join(outdir, "manifest.json")
    assert os.path.exists(manifest_json)
    with open(manifest_json) as f:
        got_manifest = json.load(f)
    assert got_manifest == expect_manifest
    got_amp1 = os.path.join(outdir, "0.fa")
    got_amp3 = os.path.join(outdir, "1.fa")
    expect_amp1 = os.path.join(data_dir, "make_reads_dir_for_viridian.unpaired.0.fa")
    expect_amp3 = os.path.join(data_dir, "make_reads_dir_for_viridian.unpaired.1.fa")
    assert filecmp.cmp(got_amp1, expect_amp1, shallow=False)
    assert filecmp.cmp(got_amp3, expect_amp3, shallow=False)
    subprocess.check_output(f"rm -rf {outdir}", shell=True)


def test_make_reads_dir_for_viridian_paired():
    amplicons_tsv = os.path.join(data_dir, "make_reads_dir_for_viridian.amplicons.tsv")
    amplicon_set = primers.AmpliconSet.from_tsv(amplicons_tsv)
    amplicons = list(amplicon_set)
    assert len(amplicons) == 3
    read_store = readstore.ReadStore("name", amplicon_set)
    read_fwd = readstore.Read("A"*100, 100, 199, 0, 99, False)
    read_rev = readstore.Read("C"*100, 101, 200, 1, 98, True)
    frag = readstore.PairedReads(read_fwd, read_rev)
    read_store.reads_all_paired = True
    read_store.amplicons = {
        amplicons[0]: [frag],
        amplicons[1]: [],
        amplicons[2]: [frag]*100,
    }
    outdir = "tmp.make_reads_dir_for_viridian_paired"
    subprocess.check_output(f"rm -rf {outdir}", shell=True)
    got_failed = read_store.make_reads_dir_for_viridian(outdir, 3)
    assert got_failed == {amplicons[1]}
    expect_manifest = {
        "amp1": "0.fa",
        "amp3": "1.fa",
    }
    manifest_json = os.path.join(outdir, "manifest.json")
    assert os.path.exists(manifest_json)
    with open(manifest_json) as f:
        got_manifest = json.load(f)
    assert got_manifest == expect_manifest
    got_amp1 = os.path.join(outdir, "0.fa")
    got_amp3 = os.path.join(outdir, "1.fa")
    expect_amp1 = os.path.join(data_dir, "make_reads_dir_for_viridian.paired.0.fa")
    expect_amp3 = os.path.join(data_dir, "make_reads_dir_for_viridian.paired.1.fa")
    assert filecmp.cmp(got_amp1, expect_amp1, shallow=False)
    assert filecmp.cmp(got_amp3, expect_amp3, shallow=False)
    subprocess.check_output(f"rm -rf {outdir}", shell=True)
