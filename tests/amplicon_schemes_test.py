import json
import os
import pytest
import subprocess
from unittest import mock

from viridian_workflow import primers, readstore, amplicon_schemes

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicon_schemes")


# def test_get_built_in_schemes():
#    found_schemes = amplicon_schemes.get_built_in_schemes()
#    assert len(found_schemes) > 0
#    for filename in found_schemes.values():
#        assert os.path.exists(filename)
#        primers.AmpliconSet.from_tsv(filename)


class MockBam:
    def __init__(self, l):
        self.l = l

    def syncronise_fragments(self):
        yield self.l


def test_convert_tsv_to_cylon_json():
    tsv_in = os.path.join(data_dir, "convert_tsv_to_cylon_json.tsv")
    expect_json = os.path.join(data_dir, "convert_tsv_to_cylon_json.json")
    tmp_json = "tmp.convert_tsv_to_cylon_json.json"
    subprocess.check_output(f"rm -f {tmp_json}", shell=True)
    amplicon_set = primers.AmpliconSet.from_tsv(tsv_in, name="test_name")
    bam = mock.Mock()
    bam.syncronise_fragments = list
    got = readstore.ReadStore(amplicon_set, bam).cylon_json["amplicons"]
    with open(expect_json) as f:
        expect = json.load(f)
    expect = expect["amplicons"]

    # Only test the subset of old fields
    for amplicon in got:
        for field in got[amplicon]:
            print(expect[amplicon][field], got[amplicon][field])
            assert expect[amplicon][field] == got[amplicon][field]


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
    expect_list = [primers.AmpliconSet.from_tsv(scheme1_tsv, name="Scheme1")]
    expect_dict = {"Scheme1": scheme1_tsv}

    got_dict = {}
    got_list = []

    got_dict, got_list = amplicon_schemes.load_list_of_amplicon_sets(
        tsv_others_to_use=tmp_tsv
    )
    assert got_list == expect_list
    assert got_dict == expect_dict

    with pytest.raises(Exception):
        amplicon_schemes.load_list_of_amplicon_sets(
            built_in_names_to_use={"does not exist"}, tsv_others_to_use=tmp_tsv
        )

    built_in_schemes = amplicon_schemes.get_built_in_schemes()
    expect_dict["COVID-ARTIC-V4"] = built_in_schemes["COVID-ARTIC-V4"]
    expect_list = [
        primers.AmpliconSet.from_tsv(
            built_in_schemes["COVID-ARTIC-V4"], name="COVID-ARTIC-V4"
        ),
        primers.AmpliconSet.from_tsv(scheme1_tsv, name="Scheme1"),
    ]
    got_dict, got_list = amplicon_schemes.load_list_of_amplicon_sets(
        built_in_names_to_use={"COVID-ARTIC-V4"}, tsv_others_to_use=tmp_tsv
    )
    assert got_list == expect_list
    assert got_dict == expect_dict
    os.unlink(tmp_tsv)
