import os
import pytest
import subprocess

import pysam

from viridian_workflow import minimap

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "minimap")


def bam_is_sorted_and_indexed(bam):
    if not os.path.exists(f"{bam}.bai"):
        return False
    aln_file = pysam.AlignmentFile(bam, "rb")
    header = aln_file.header.to_dict()
    # If it's sorted, the a line like this is present:
    # @HD VN:1.6 SO:coordinate
    if "HD" not in header:
        return False
    return "HD" in header and header["HD"].get("SO", None) == "coordinate"

def test_paired_unsorted():
    ref = os.path.join(data_dir, "ref.fa")
    reads1 = os.path.join(data_dir, "reads_1.fq")
    reads2 = os.path.join(data_dir, "reads_2.fq")
    bam = "tmp.minimap_paired_unsorted.bam"
    subprocess.check_output(f"rm -f {bam}", shell=True)
    minimap.run(bam, ref, reads1, fq2=reads2, sort=False)
    assert os.path.exists(bam)
    assert not bam_is_sorted_and_indexed(bam)
    os.unlink(bam)

def test_paired_sorted():
    ref = os.path.join(data_dir, "ref.fa")
    reads1 = os.path.join(data_dir, "reads_1.fq")
    reads2 = os.path.join(data_dir, "reads_2.fq")
    bam = "tmp.minimap_paired_sorted.bam"
    subprocess.check_output(f"rm -f {bam}", shell=True)
    minimap.run(bam, ref, reads1, fq2=reads2, sort=True)
    assert os.path.exists(bam)
    assert bam_is_sorted_and_indexed(bam)
    os.unlink(bam)
    os.unlink(f"{bam}.bai")

def test_unpaired_unsorted():
    ref = os.path.join(data_dir, "ref.fa")
    reads = os.path.join(data_dir, "reads.fq")
    bam = "tmp.minimap_unpaired_unsorted.bam"
    subprocess.check_output(f"rm -f {bam}", shell=True)
    minimap.run(bam, ref, reads, sort=False)
    assert os.path.exists(bam)
    assert not bam_is_sorted_and_indexed(bam)
    os.unlink(bam)

def test_unpaired_sorted():
    ref = os.path.join(data_dir, "ref.fa")
    reads = os.path.join(data_dir, "reads.fq")
    bam = "tmp.minimap_unpaired_sorted.bam"
    subprocess.check_output(f"rm -f {bam}", shell=True)
    minimap.run(bam, ref, reads, sort=True)
    assert os.path.exists(bam)
    assert bam_is_sorted_and_indexed(bam)
    os.unlink(bam)
    os.unlink(f"{bam}.bai")
