import os
import pytest

from viridian import amplicon_schemes, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicon_schemes")


def test_get_built_in_schemes():
    found_schemes = amplicon_schemes.get_built_in_schemes()
    assert len(found_schemes) > 0
    for filename in found_schemes.values():
        assert os.path.exists(filename)


def test_load_list_of_amplicon_sets():
    with pytest.raises(Exception):
        amplicon_schemes.load_list_of_amplicon_sets(
            built_in_names_to_use=None, tsv_others_to_use=None
        )
    scheme1_tsv = os.path.join(data_dir, "load_list_of_amplicon_sets.scheme.tsv")
    tmp_tsv = "tmp.load_list_of_amplicon_sets.tsv"
    utils.syscall(f"rm -f {tmp_tsv}")
    with pytest.raises(Exception):
        amplicon_schemes.load_list_of_amplicon_sets(tsv_others_to_use=tmp_tsv)
    with open(tmp_tsv, "w") as f:
        print("Name", "File", sep="\t", file=f)
        print("Scheme1", scheme1_tsv, sep="\t", file=f)
    expect_dict = {"Scheme1": scheme1_tsv}
    got_dict = amplicon_schemes.load_list_of_amplicon_sets(tsv_others_to_use=tmp_tsv)
    assert got_dict == expect_dict

    with pytest.raises(Exception):
        amplicon_schemes.load_list_of_amplicon_sets(
            built_in_names_to_use={"does not exist"}, tsv_others_to_use=tmp_tsv
        )

    built_in_schemes = amplicon_schemes.get_built_in_schemes()
    expect_dict["COVID-ARTIC-V4.1"] = built_in_schemes["COVID-ARTIC-V4.1"]
    got_dict = amplicon_schemes.load_list_of_amplicon_sets(
        built_in_names_to_use={"COVID-ARTIC-V4.1"}, tsv_others_to_use=tmp_tsv
    )
    assert got_dict == expect_dict
    os.unlink(tmp_tsv)
