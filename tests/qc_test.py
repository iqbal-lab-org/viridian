import copy
import pytest

from cluster_vcf_records import vcf_record

from viridian import qc, utils


def test_unique_amp_coords():
    amplicons = [
        {"start": 10, "end": 50},
        {"start": 40, "end": 90},
        {"start": 70, "end": 120},
    ]
    expect = [(10, 39), (51, 69), (91, 120)]
    assert qc.unique_amp_coords(amplicons) == expect


def test_get_dropped_amplicons():
    amplicons = [
        {"start": 3, "end": 8},
        {"start": 7, "end": 13},
        {"start": 11, "end": 16},
    ]
    ref_msa = "0123456789012--3456"
    con_msa = "AAAAAAAAAAAAAAAA--A"
    assert qc.get_dropped_amplicons(amplicons, ref_msa, con_msa) == set()
    con_msa = "AAANAAAAAAAAAAAA--A"
    assert qc.get_dropped_amplicons(amplicons, ref_msa, con_msa) == set()
    con_msa = "AAANNAAAAAAAAAAA--A"
    assert qc.get_dropped_amplicons(amplicons, ref_msa, con_msa) == {0}
    con_msa = "AAANNAAAANNAAAAA--A"
    assert qc.get_dropped_amplicons(amplicons, ref_msa, con_msa) == {0, 1}


def test_base_counts_to_iupac():
    counts = {x: 0 for x in ["a", "A", "c", "C", "g", "G", "t", "T"]}
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == None
    counts["A"] = 9
    counts["a"] = 9
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == None
    counts["g"] = 5
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == None
    counts["G"] = 5
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == "R"

    counts = {x: 0 for x in ["a", "A", "c", "C", "g", "G", "t", "T"]}
    counts["A"] = 9
    counts["G"] = 5
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == "R"

    counts = {x: 0 for x in ["a", "A", "c", "C", "g", "G", "t", "T"]}
    counts["G"] = counts["g"] = 75
    counts["t"] = counts["T"] = 25
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == "K"

    counts = {x: 0 for x in ["a", "A", "c", "C", "g", "G", "t", "T"]}
    counts["c"] = 3
    counts["C"] = 17
    counts["T"] = 3
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == None

    counts["T"] = 4
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == None
    counts["T"] = 5
    assert qc.base_counts_to_iupac(counts, min_total=5, min_pc=10.0) == "Y"


def test_base_counts_to_insertion():
    indel_counts = {"I": {}, "i": {}}
    total_depth = 100
    min_prop = 0.1
    assert qc.base_counts_to_insertion(indel_counts, total_depth, min_prop) == False
    indel_counts["I"]["A"] = 8
    assert qc.base_counts_to_insertion(indel_counts, total_depth, min_prop) == False
    indel_counts["i"]["g"] = 8
    assert qc.base_counts_to_insertion(indel_counts, total_depth, min_prop) == False
    indel_counts["i"]["a"] = 1
    assert qc.base_counts_to_insertion(indel_counts, total_depth, min_prop) == False
    indel_counts["i"]["a"] = 2
    assert qc.base_counts_to_insertion(indel_counts, total_depth, min_prop) == True


def test_qc():
    amplicons = [
        {
            "name": "amp1",
            "start": 10,
            "end": 20,
            "primers": {
                "left": [
                    {"start": 10, "end": 12},
                    {"start": 13, "end": 15},
                ],
                "right": [{"start": 19, "end": 20}],
            },
        },
        {
            "name": "amp2",
            "start": 17,
            "end": 28,
            "primers": {
                "left": [{"start": 17, "end": 18}],
                "right": [{"start": 26, "end": 28}],
            },
        },
    ]
    default_indels = {x: 0 for x in ["D", "d", "I", "i"]}
    default_indels["indel"] = {x: {} for x in ["D", "d", "I", "i"]}

    pileups = {
        0: {
            (0, 0): {
                9: {"G": 1, "g": 1},
                10: {"T": 2, "t": 3},
                11: {"A": 11, "a": 11},
                12: {"G": 11, "G": 11},
                13: {"T": 11, "G": 11},
                14: {"A": 11, "G": 11},
                15: {"C": 11, "G": 11},
                16: {"G": 11, "G": 11},
                17: {"T": 11, "G": 11},
                18: {"G": 11, "G": 11},
                19: {"C": 10, "c": 10},
            },
            (1, 0): {
                13: {"T": 40, "t": 41},
                14: {"A": 40, "a": 41},
                15: {"C": 40, "c": 41},
                16: {"G": 40, "g": 41},
                17: {"T": 40, "t": 41},
                18: {"G": 40, "g": 41},
                19: {"A": 40, "a": 41},
            },
        },
        1: {
            (0, 0): {
                19: {"A": 200, "a": 200},
                20: {"C": 200, "c": 201},
                21: {"G": 201, "g": 199},
                22: {"A": 200, "a": 201, "C": 1, "c": 2},
                23: {"A": 200, "a": 199},
                24: {"C": 201, "c": 199},
                25: {"G": 199, "g": 200},
                26: {"T": 200, "t": 201},
            },
        },
    }
    for amp in pileups:
        for primers in pileups[amp]:
            for d in pileups[amp][primers].values():
                d.update(copy.deepcopy(default_indels))

    ref_msa = "ACGTACGTACGTACGTACGT-ACGTACGTAC"
    con_msa = "-CGTACGTACGTA-GTANGTGACGAACGTA-"
    # con coord 012345678901 2345678901234567
    tsv_out = "tmp.qc.make_qc_tsv.tsv"
    vcf_out = "tmp.qc.annotated_vcf_file.vcf"
    utils.syscall(f"rm -f {tsv_out} {vcf_out}")
    q = qc.Qc(amplicons, ref_msa, con_msa, max_amp_n_percent=50)
    q.make_pileup(pileups)
    q.make_tsv_lines_and_masked_cons()
    q.write_qc_tsv_and_make_masked_cons_msa(tsv_out)
    vcf_header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=ref,length=500>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
    ]
    vcf_records = [
        vcf_record.VcfRecord("ref\t24\t.\tT\tA\t.\tPASS\t.\tGT\t1/1"),
    ]
    q.annotated_vcf_file(vcf_header, vcf_records, vcf_out)
    # FIXME add more tests (check tsv_out is same as expected tsv_out).
    # - can do this when decided on filtering methods and indel handling
    # FIXME check VCF file - need to add in a couple more variants first, and
    # decide on how to handle indels
    assert q.masked_cons_msa == "-NNNNNNNNNNNN-GKANGNNACGAACNNN-"
    assert q.masked_cons_msa_indel_as_ref == "NNNNNNNNNNNNNNGKANGNACGAACNNNN"
    assert q.masked_cons_msa_indel_as_N == "NNNNNNNNNNNNNNGKANGNACGAACNNNN"
    utils.syscall(f"rm {tsv_out} {vcf_out}")
