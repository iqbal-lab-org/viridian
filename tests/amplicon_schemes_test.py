import json
import os
import pytest
import subprocess

from viridian_workflow import amplicon_schemes, primers

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicon_schemes")


def test_get_built_in_schemes():
    found_schemes = amplicon_schemes.get_built_in_schemes()
    assert len(found_schemes) > 0
    for filename in found_schemes.values():
        assert os.path.exists(filename)
        foo = primers.AmpliconSet.from_tsv_viridian_workflow_format(filename)


def test_convert_tsv_to_viridian_json():
    tsv_in = os.path.join(data_dir, "convert_tsv_to_viridian_json.tsv")
    expect_json = os.path.join(data_dir, "convert_tsv_to_viridian_json.json")
    tmp_json = "tmp.convert_tsv_to_viridian_json.json"
    subprocess.check_output(f"rm -f {tmp_json}", shell=True)
    amplicon_schemes.convert_tsv_to_viridian_json(
        tsv_in, tmp_json, scheme_name="test_name"
    )
    with open(tmp_json) as f:
        got = json.load(f)
    with open(expect_json) as f:
        expect = json.load(f)
    # We can't expect source filename to be the same because is abs path.
    # Change it to basename
    got["source_file"] = os.path.basename(got["source_file"])
    expect["source_file"] = os.path.basename(expect["source_file"])
    assert got == expect
    os.unlink(tmp_json)


def test_load_list_of_amplicon_sets():
    with pytest.raises(Exception):
        amplicon_schemes.load_list_of_amplicon_sets(
            built_in_names_to_use=None, tsv_others_to_use=None
        )
    scheme1_tsv = os.path.join(data_dir, "load_list_of_amplicon_sets.scheme.tsv")
    tmp_tsv = "tmp.load_list_of_amplicon_sets.tsv"
    subprocess.check_output(f"rm -f {tmp_tsv}", shell=True)
    with pytest.raises(Exception):
        amplicon_schemes.load_list_of_amplicon_sets(tsv_others_to_use=tmp_tsv)
    with open(tmp_tsv, "w") as f:
        print("Name", "File", sep="\t", file=f)
        print("Scheme1", scheme1_tsv, sep="\t", file=f)
    expect = [primers.AmpliconSet("Scheme1", vwf_tsv_file=scheme1_tsv)]
    got = amplicon_schemes.load_list_of_amplicon_sets(tsv_others_to_use=tmp_tsv)
    assert got == expect

    with pytest.raises(Exception):
        amplicon_schemes.load_list_of_amplicon_sets(
            built_in_names_to_use={"does not exist"}, tsv_others_to_use=tmp_tsv
        )

    built_in_schemes = amplicon_schemes.get_built_in_schemes()
    expect = [
        primers.AmpliconSet(
            "COVID-ARTIC-V4", vwf_tsv_file=built_in_schemes["COVID-ARTIC-V4"]
        ),
        primers.AmpliconSet("Scheme1", vwf_tsv_file=scheme1_tsv),
    ]
    got = amplicon_schemes.load_list_of_amplicon_sets(
        built_in_names_to_use={"COVID-ARTIC-V4"}, tsv_others_to_use=tmp_tsv
    )
    assert got == expect
    os.unlink(tmp_tsv)
