import os
import pytest

from viridian_workflow import one_sample_pipeline

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "one_sample_pipeline")

def test_foo():
    assert False

