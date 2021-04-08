"""qcovid wrapper
"""
import logging
import subprocess
import os


def bin_amplicons(outdir, ref_genome, amplicon_bed, bam):
    mask = os.path.join(outdir, "mask")
    bin_command = ["bin_amplicons.py", "--mask", mask, amplicon_bed, bam]
    subprocess.Popen(bin_command).wait()
    assert os.path.isfile(mask)
    return mask


def self_qc(outdir, assembly, bam):
    raise NotImplementedError
