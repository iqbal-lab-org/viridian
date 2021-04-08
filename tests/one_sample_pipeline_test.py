import os
import pytest
import tempfile

from viridian_workflow import one_sample_pipeline

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "one_sample_pipeline")

def test_one_sample():
    pass
    # TODO: decide on minimised test data (don't want to pollute git)
    #with tempfile.TemporaryDirectory() as outdir:
    #    ref_genome = os.path.join(data_dir, "MN908947.fasta")
    #    fq1 = os.path.join(data_dir, "ERR-trunc1.fq.gz")
    #    fq2 = os.path.join(data_dir, "ERR-trunc2.fq.gz")
    #    amplicon_bed = os.path.join(data_dir, "ARTIC3-truncated.bed")
    #    one_sample_pipeline.run_one_sample(outdir, ref_genome, amplicon_bed, fq1, fq2)
