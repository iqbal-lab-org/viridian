import os
import pytest

from viridian_workflow import amplicon_schemes


def test_get_built_in_schemes():
    found_schemes = amplicon_schemes.get_built_in_schemes()
    assert len(found_schemes) > 0
    for filename in found_schemes.values():
        assert os.path.exists(filename)
