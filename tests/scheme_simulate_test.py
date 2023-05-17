import json
import os
import pytest

from viridian import amplicon_schemes, scheme_simulate, utils


def test_simulate_all_schemes():
    # just a basic test - run the code, check we got sane output. Should
    # call the expected scheme each time
    outdir = "tmp.simulate_all_schemes"
    utils.syscall(f"rm -rf {outdir}")
    scheme_simulate.simulate_all_schemes(
        amplicon_schemes.REF_FASTA,
        outdir,
    )

    summary_json = os.path.join(outdir, "summary.json")
    assert os.path.exists(summary_json)
    with open(summary_json) as f:
        results = json.load(f)

    for test_scheme in results:
        if test_scheme == "wgs":
            continue
        scheme_name, read_type = test_scheme.split(":", maxsplit=1)
        assert read_type in ["full_length", "fragmented"]
        assert results[test_scheme]["best_schemes"] == [scheme_name]

    utils.syscall(f"rm -r {outdir}")
