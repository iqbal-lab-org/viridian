import copy
import os
import pytest
import random

from viridian import varifier, utils


def test_run_varifier():
    ref_bases = random.choices(["A", "C", "G", "T"], k=1000)
    ref_bases[11] = "A"
    ref_bases[501] = "A"
    ref_bases[951] = "A"
    cons_bases = copy.deepcopy(ref_bases)
    cons_bases[11] = "T"
    cons_bases[501] = "T"
    cons_bases[951] = "T"
    root_dir = "tmp.varifier_test"
    utils.syscall(f"rm -rf {root_dir}")
    os.mkdir(root_dir)
    ref_fasta = os.path.join(root_dir, "ref.fa")
    cons_fasta = os.path.join(root_dir, "cons.fa")
    with open(ref_fasta, "w") as f:
        print(">ref", "".join(ref_bases), sep="\n", file=f)
    with open(cons_fasta, "w") as f:
        print(">cons", "".join(cons_bases), sep="\n", file=f)
    varifier_out = os.path.join(root_dir, "varifier")
    got_seqs, got_vcf_header, got_vcf_records, got_error = varifier.run_varifier(
        cons_fasta,
        ref_fasta,
        varifier_out,
        50,
        900,
    )
    assert len(got_vcf_records) == 1
    assert got_vcf_records[0].CHROM == "ref"
    assert got_vcf_records[0].POS == 501
    assert got_vcf_records[0].REF == "A"
    assert got_vcf_records[0].ALT == ["T"]

    utils.syscall(f"rm -r {varifier_out}")
    got_seqs, got_vcf_header, got_vcf_records, got_error = varifier.run_varifier(
        cons_fasta,
        "not_a_file",
        varifier_out,
        50,
        900,
        force_use_mafft=True,
    )
    assert set(got_seqs.values()) == {None}
    assert got_vcf_header is None
    assert got_vcf_records is None
    assert got_error == "Error running varifier"
    utils.syscall(f"rm -r {root_dir}")
