import copy
import os
import pytest
import random
import statistics

import pyfastaq

from viridian import maptools, scheme_id, utils

random.seed(42)

SCHEME_TSV_HEADER_FIELDS = [
    "Amplicon_name",
    "Primer_name",
    "Left_or_right",
    "Sequence",
    "Position",
]


@pytest.fixture()
def amp_scheme_data():
    amps = [
        {
            "name": "amp1",
            "start": 15,
            "end": 52,
            "primers": {
                "left": [[15, 17], [20, 23]],
                "right": [[50, 52]],
            },
        },
        {
            "name": "amp2",
            "start": 40,
            "end": 103,
            "primers": {
                "left": [[40, 40]],
                "right": [[90, 94], [100, 103]],
            },
        },
    ]
    return {
        "amplicons": amps,
        "left_starts": {15: 0, 20: 0, 40: 1},
        "right_ends": {52: 0, 94: 1, 103: 1},
        "amplicon_name_indexes": {"amp1": 0, "amp2": 1},
        "mean_amp_length": statistics.mean(
            [52 - 15 + 1, 40 - 15, 103 - 52, 103 - 40 + 1]
        ),
    }


def test_scheme_init_from_tsv(amp_scheme_data):
    tmp_tsv = "test_scheme_init_from_tsv.tsv"
    with open(tmp_tsv, "w") as f:
        print(*SCHEME_TSV_HEADER_FIELDS, sep="\t", file=f)
        print("amp1", "amp1_left", "left", "ACGT", "20", sep="\t", file=f)
        print("amp1", "amp1_left_alt", "left", "CAG", "15", sep="\t", file=f)
        print("amp1", "amp1_left_alt2", "left", "CA", "15", sep="\t", file=f)
        print("amp1", "amp1_right", "right", "AAA", "50", sep="\t", file=f)
        print("amp2", "amp2_left", "left", "A", "40", sep="\t", file=f)
        print("amp2", "amp2_right", "right", "ATGTT", "90", sep="\t", file=f)
        print("amp2", "amp2_right_alt", "right", "GTA", "101", sep="\t", file=f)
        print("amp2", "amp2_right_alt2", "right", "GGTA", "100", sep="\t", file=f)

    scheme = scheme_id.Scheme(amp_scheme_file=tmp_tsv)
    assert scheme.amplicons == amp_scheme_data["amplicons"]
    assert scheme.left_starts == amp_scheme_data["left_starts"]
    assert scheme.right_ends == amp_scheme_data["right_ends"]
    assert scheme.amplicon_name_indexes == amp_scheme_data["amplicon_name_indexes"]
    assert scheme.mean_amp_length == amp_scheme_data["mean_amp_length"]
    os.unlink(tmp_tsv)


def test_scheme_init_from_bed(amp_scheme_data):
    tmp_bed = "test_scheme_init_from_tsv.bed"
    with open(tmp_bed, "w") as f:
        print("REF", 20, 24, "amp1_LEFT_1", "1", "+", "ACGT", sep="\t", file=f)
        print("REF", 15, 18, "amp1_LEFT_alt", "1", "+", "CAG", sep="\t", file=f)
        print("REF", 15, 17, "amp1_LEFT_alt2", "1", "+", "CA", sep="\t", file=f)
        print("REF", 50, 53, "amp1_RIGHT_1", "1", "+", "AAA", sep="\t", file=f)
        print("REF", 40, 41, "amp2_LEFT_1", "1", "+", "A", sep="\t", file=f)
        print("REF", 90, 95, "amp2_RIGHT_1", "1", "+", "ATGTT", sep="\t", file=f)
        print("REF", 101, 104, "amp2_RIGHT_1", "1", "+", "GTA", sep="\t", file=f)
        print("REF", 100, 104, "amp2_RIGHT_2", "1", "+", "GGTA", sep="\t", file=f)

    scheme = scheme_id.Scheme(amp_scheme_file=tmp_bed)
    assert scheme.amplicons == amp_scheme_data["amplicons"]
    assert scheme.left_starts == amp_scheme_data["left_starts"]
    assert scheme.right_ends == amp_scheme_data["right_ends"]
    assert scheme.amplicon_name_indexes == amp_scheme_data["amplicon_name_indexes"]
    assert scheme.mean_amp_length == amp_scheme_data["mean_amp_length"]
    os.unlink(tmp_bed)


def test_scheme_init_distance_lists():
    ref_length = 14
    scheme = scheme_id.Scheme(end_tolerance=0)
    # for this test, values of the start positions don't matter, use None
    scheme.left_starts = {3: None, 7: None, 8: None, 11: None}
    scheme.right_ends = {5: None, 9: None, 10: None, 11: None}
    scheme.init_distance_lists(ref_length)
    # indexes:                   0    1  2   3  4  5  6  7  8  9 10 11 12 13
    assert scheme.left_dists == [-3, -2, -1, 0, 1, 2, 3, 0, 0, 1, 2, 0, 1, 2]
    # indexes:                    0  1  2  3  4  5  6  7  8  9 10 11  12  13
    assert scheme.right_dists == [5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 0, 0, -1, -2]

    scheme.end_tolerance = 1
    scheme.init_distance_lists(ref_length)
    # indexes:                   0    1  2  3  4  5  6  7  8  9 10 11 12 13
    assert scheme.left_dists == [-3, -2, 0, 0, 1, 2, 0, 0, 0, 1, 0, 0, 1, 2]
    # indexes:                    0  1  2  3  4  5  6  7  8  9 10 11  12  13
    assert scheme.right_dists == [5, 4, 3, 2, 1, 0, 0, 2, 1, 0, 0, 0, 0, -2]


def test_make_primer_distance_hist():
    left_hits = [3, 100, 1, 0, 1, 80, 4, 3]
    right_hits = [0, 1, 60, 2, 1, 20, 50, 3]
    scheme = scheme_id.Scheme()
    scheme.left_dists = [-1, 0, 1, 2, 3, 0, 1, 2]
    scheme.right_dists = [2, 1, 0, 3, 2, 1, 0, -1]
    scheme.make_primer_distance_hist(left_hits, right_hits)
    assert scheme.primer_dist_hist == {
        -1: 3 + 3,
        0: 100 + 80 + 60 + 50,
        1: 1 + 4 + 1 + 20,
        2: 3 + 1,
        3: 1 + 2,
    }


def test_make_normalised_distance_hists():
    scheme = scheme_id.Scheme()
    scheme.mean_amp_length = 10
    scheme.primer_dist_hist = {
        -1: 5,
        0: 70,
        1: 10,
        2: 15,
    }
    scheme.make_normalised_distance_hists()
    assert scheme.norm_primer_dist_hist == {-10: 5, 0: 70, 10: 10, 20: 15}
    assert (
        scheme.norm_cumul_primer_dists == [70] * 10 + [80] * 10 + [95] * 70 + [100] * 11
    )
    assert scheme.score == 4200


def test_count_primer_hits():
    tmp_tsv = "test_count_primer_hits.tsv"
    with open(tmp_tsv, "w") as f:
        print(*SCHEME_TSV_HEADER_FIELDS, sep="\t", file=f)
        print("amp1", "amp1_left", "left", "ACGTA", "15", sep="\t", file=f)
        print("amp1", "amp1_left_alt1", "left", "CGTACGTA", "16", sep="\t", file=f)
        print("amp1", "amp1_left_alt2", "left", "AAAAAA", "30", sep="\t", file=f)
        print("amp1", "amp1_right", "right", "ACGTACG", "50", sep="\t", file=f)
        print("amp1", "amp1_right_alt1", "right", "ACGT", "54", sep="\t", file=f)
    scheme = scheme_id.Scheme(amp_scheme_file=tmp_tsv)
    scheme.init_distance_lists(65)
    expect_amplicons = copy.deepcopy(scheme.amplicons)
    left_hits = [0] * 65
    left_hits[11] = 32
    left_hits[12] = 1
    left_hits[13] = 2
    left_hits[14] = 4
    left_hits[15] = 100
    left_hits[16] = 8
    left_hits[20] = 16
    left_hits[21] = 1000
    right_hits = [0] * 65
    right_hits[49] = 2000
    right_hits[50] = 1000
    right_hits[51] = 1
    right_hits[52] = 2
    right_hits[53] = 4
    right_hits[54] = 8
    right_hits[55] = 16
    right_hits[56] = 32
    right_hits[57] = 64
    right_hits[58] = 128
    right_hits[59] = 256
    right_hits[60] = 512
    right_hits[61] = 1024
    scheme.count_primer_hits(left_hits, right_hits, max_dist=5)

    expect_amplicons[0]["primer_counts"] = {
        "left": [1, 130, 0],
        "right": [511, 512],
    }
    assert scheme.amplicons == expect_amplicons
    os.unlink(tmp_tsv)


def test_to_cleaned_amp_list_and_cylon_dict():
    scheme = scheme_id.Scheme()
    scheme.amplicons = [
        {
            "name": "amp2",
            "start": 40,
            "end": 103,
            "primers": {
                "left": [[40, 40]],
                "right": [[90, 94], [100, 103]],
            },
            "primer_counts": {
                "left": [50],
                "right": [49, 49],
            },
        },
        {
            "name": "amp1",
            "start": 15,
            "end": 52,
            "primers": {
                "left": [[15, 17], [20, 23]],
                "right": [[50, 52]],
            },
            "primer_counts": {
                "left": [49, 50],
                "right": [49],
            },
        },
    ]
    expect_list = [
        {
            "name": "amp1",
            "start": 20,
            "end": 52,
            "primers": {
                "left": [{"start": 20, "end": 23, "read_count": 50}],
                "right": [{"start": 50, "end": 52, "read_count": 49}],
            },
            "excluded_primers": {
                "left": [{"start": 15, "end": 17, "read_count": 49}],
                "right": [],
            },
        },
        {
            "name": "amp2",
            "start": 40,
            "end": 103,
            "primers": {
                "left": [{"start": 40, "end": 40, "read_count": 50}],
                "right": [
                    {"start": 90, "end": 94, "read_count": 49},
                    {"start": 100, "end": 103, "read_count": 49},
                ],
            },
            "excluded_primers": {
                "left": [],
                "right": [],
            },
        },
    ]
    expect_dict = {
        "amp1": {
            "start": 20,
            "end": 52,
            "left_primer_end": 23,
            "right_primer_start": 50,
        },
        "amp2": {
            "start": 40,
            "end": 103,
            "left_primer_end": 40,
            "right_primer_start": 100,
        },
    }
    got_list, got_dict = scheme.to_cleaned_amp_list_and_cylon_dict(50)
    assert got_list == expect_list
    assert got_dict == expect_dict


def test_depth_per_position_plot():
    outfile = "tmp.depth_per_position_plot.pdf"
    utils.syscall(f"rm -f {outfile}")
    depths = [0] * 10 + [100] * 90 + [200] * 10 + [150] * 100 + [0] * 20
    coords = [(10, 110), (100, 210)]
    scheme_id.depth_per_position_plot(depths, outfile, amp_coords=coords, title="TEST")
    assert os.path.exists(outfile)
    os.unlink(outfile)


def test_cumulative_score_plot():
    outfile = "tmp.cumulative_score_plot.pdf"
    utils.syscall(f"rm -f {outfile}")
    scheme1 = scheme_id.Scheme()
    scheme1.norm_cumul_primer_dists = [50] * 50 + [100] * 51
    scheme2 = scheme_id.Scheme()
    scheme2.norm_cumul_primer_dists = [30] * 50 + [80] * 40 + [100] * 11
    schemes = {"s1": scheme1, "s2": scheme2}
    scheme_id.cumulative_score_plot(schemes, outfile, title="TEST")
    assert os.path.exists(outfile)
    os.unlink(outfile)


def test_get_scores_from_schemes():
    scheme1 = scheme_id.Scheme()
    scheme1.score = None
    scheme2 = scheme_id.Scheme()
    scheme2.score = 100
    scheme3 = scheme_id.Scheme()
    scheme3.score = 200
    schemes = {"s1": scheme1}
    expect = {
        "best_scheme": None,
        "best_schemes": [],
        "best_score": None,
        "scores": {"s1": None},
        "score_ratio": None,
    }
    assert scheme_id.get_scores_from_schemes(schemes) == expect

    schemes["s2"] = scheme2
    expect = {
        "best_scheme": "s2",
        "best_schemes": ["s2"],
        "best_score": 100,
        "scores": {"s1": None, "s2": 100},
        "score_ratio": None,
    }
    assert scheme_id.get_scores_from_schemes(schemes) == expect

    schemes["s2.2"] = scheme2
    expect = {
        "best_scheme": "s2",
        "best_schemes": ["s2", "s2.2"],
        "best_score": 100,
        "scores": {"s1": None, "s2": 100, "s2.2": 100},
        "score_ratio": 1.0,
    }
    assert scheme_id.get_scores_from_schemes(schemes) == expect

    schemes["s3"] = scheme3
    expect = {
        "best_scheme": "s3",
        "best_schemes": ["s3"],
        "best_score": 200,
        "scores": {"s1": None, "s2": 100, "s2.2": 100, "s3": 200},
        "score_ratio": 0.5,
    }
    assert scheme_id.get_scores_from_schemes(schemes) == expect


def test_parse_bam():
    ref_length = 1000
    bam = "tmp.parse_bam.bam"
    ref_fa = "tmp.parse_bam.ref.fa"
    ref_seq = pyfastaq.sequences.Fasta(
        "ref", "".join(random.choices(["A", "C", "G", "T"], k=ref_length))
    )
    with open(ref_fa, "w") as f:
        print(ref_seq, file=f)

    # ------------------- unpaired reads ----------------------------------------
    reads_fa = "tmp.parse_bam.reads.fa"
    with open(reads_fa, "w") as f:
        print(">r1", ref_seq[100:200], sep="\n", file=f)
        print(">r2", ref_seq[30:250], sep="\n", file=f)
        print(">r3", "A" * 100, sep="\n", file=f)
    maptools.map_reads(bam, ref_fa, reads_fa)
    os.unlink(reads_fa)
    got_ref_length, got_lefts, got_rights, got_pileup, got_counts = scheme_id.parse_bam(
        bam
    )
    assert got_ref_length == ref_length
    expect_lefts = [0] * ref_length
    expect_lefts[30] = 1
    expect_lefts[100] = 1
    assert got_lefts == expect_lefts
    expect_rights = [0] * ref_length
    expect_rights[199] = 1
    expect_rights[249] = 1
    assert got_rights == expect_rights
    expect_pileup = {
        "depth_hist": {0: 780, 1: 120, 2: 100},
        "depth_per_position": [0] * 30 + [1] * 70 + [2] * 100 + [1] * 50 + [0] * 750,
        "depth_at_least": {1: 220, 2: 100, 5: 0, 10: 0, 15: 0, 20: 0, 50: 0, 100: 0},
        "percent_at_least_x_depth": {
            1: 22.0,
            2: 10.0,
            5: 0.0,
            10: 0.0,
            15: 0.0,
            20: 0.0,
            50: 0.0,
            100: 0.0,
        },
        "mean_depth": 0.32,
        "mode_depth": 0,
        "median_depth": 0.0,
    }
    assert got_pileup == expect_pileup
    expect_counts = {
        "mapped": 2,
        "reads1": 0,
        "reads2": 0,
        "total_reads": 3,
        "unpaired_reads": 3,
    }
    assert got_counts == expect_counts
    os.unlink(bam)
    os.unlink(bam + ".bai")

    # --------------------- paired reads ----------------------------------------
    reads_1_fa = "tmp.parse_bam.reads_1.fa"
    reads_2_fa = "tmp.parse_bam.reads_2.fa"
    with open(reads_1_fa, "w") as f:
        print(">r1 /1", ref_seq[100:150], sep="\n", file=f)
    with open(reads_2_fa, "w") as f:
        seq = pyfastaq.sequences.Fasta("r1 /2", ref_seq[160:210])
        seq.revcomp()
        print(seq, file=f)
    maptools.map_reads(bam, ref_fa, reads_1_fa, reads2=reads_2_fa)
    os.unlink(reads_1_fa)
    os.unlink(reads_2_fa)
    got_ref_length, got_lefts, got_rights, got_pileup, got_counts = scheme_id.parse_bam(
        bam
    )
    assert got_ref_length == ref_length
    expect_lefts = [0] * ref_length
    expect_lefts[100] = 1
    assert got_lefts == expect_lefts
    expect_rights = [0] * ref_length
    expect_rights[209] = 1
    assert got_rights == expect_rights
    expect_pileup = {
        "depth_hist": {0: 900, 1: 100},
        "depth_per_position": [0] * 100 + [1] * 50 + [0] * 10 + [1] * 50 + [0] * 790,
        "depth_at_least": {1: 100, 2: 0, 5: 0, 10: 0, 15: 0, 20: 0, 50: 0, 100: 0},
        "percent_at_least_x_depth": {
            1: 10.0,
            2: 0.0,
            5: 0.0,
            10: 0.0,
            15: 0.0,
            20: 0.0,
            50: 0.0,
            100: 0.0,
        },
        "mean_depth": 0.1,
        "mode_depth": 0,
        "median_depth": 0.0,
    }
    assert got_pileup == expect_pileup
    expect_counts = {
        "mapped": 2,
        "reads1": 1,
        "reads2": 1,
        "total_reads": 2,
        "unpaired_reads": 0,
    }
    assert got_counts == expect_counts
    os.unlink(bam)
    os.unlink(bam + ".bai")
    os.unlink(ref_fa)


def test_analyse_bam():
    root_outdir = "tmp.analyse_bam"
    utils.syscall(f"rm -rf {root_outdir}")
    os.mkdir(root_outdir)

    ref_length = 500
    ref_fa = os.path.join(root_outdir, "ref.fa")
    ref_seq = pyfastaq.sequences.Fasta(
        "ref", "".join(random.choices(["A", "C", "G", "T"], k=ref_length))
    )
    with open(ref_fa, "w") as f:
        print(ref_seq, file=f)

    scheme_tsvs = {
        "s1": os.path.join(root_outdir, "scheme.1.tsv"),
        "s2": os.path.join(root_outdir, "scheme.2.tsv"),
    }
    with open(scheme_tsvs["s1"], "w") as f:
        print(*SCHEME_TSV_HEADER_FIELDS, sep="\t", file=f)
        print("amp1", "amp1_left", "left", "ACGT", "20", sep="\t", file=f)
        print("amp1", "amp1_left_alt", "left", "CAG", "15", sep="\t", file=f)
        print("amp1", "amp1_right", "right", "AAA", "200", sep="\t", file=f)
        print("amp2", "amp2_left", "left", "A", "170", sep="\t", file=f)
        print("amp2", "amp2_right", "right", "ATGTT", "450", sep="\t", file=f)

    with open(scheme_tsvs["s2"], "w") as f:
        print(*SCHEME_TSV_HEADER_FIELDS, sep="\t", file=f)
        print("amp1", "amp1_left", "left", "ACGT", "100", sep="\t", file=f)
        print("amp1", "amp1_right", "right", "AAA", "300", sep="\t", file=f)
        print("amp2", "amp2_left", "left", "A", "280", sep="\t", file=f)
        print("amp2", "amp2_right", "right", "ATGTT", "400", sep="\t", file=f)

    reads_fa = os.path.join(root_outdir, "reads.fa")
    read1 = pyfastaq.sequences.Fasta("read.1", ref_seq[15:199])
    read1_rev = copy.deepcopy(read1)
    read1_rev.revcomp()
    read2 = pyfastaq.sequences.Fasta("read.2", ref_seq[170:448])
    read2_rev = copy.deepcopy(read2)
    read2_rev.revcomp()
    with open(reads_fa, "w") as f:
        for i in range(10):
            print(f">a.{i}", read1.seq, sep="\n", file=f)
            print(f">a.rev.{i}", read1_rev.seq, sep="\n", file=f)
            print(f">b.{i}", read2.seq, sep="\n", file=f)
            print(f">b.rev.{i}", read2_rev.seq, sep="\n", file=f)

    bam = os.path.join(root_outdir, "map.bam")
    maptools.map_reads(bam, ref_fa, reads_fa)
    outdir = os.path.join(root_outdir, "analyse_bam")
    got_json, got_err = scheme_id.analyse_bam(
        bam,
        scheme_tsvs,
        outdir,
        end_tolerance=3,
        sample_name="TEST",
        debug=True,
        max_primer_dist=10,
    )
    assert len(got_json["amplicons"]) == 2
    assert got_json["scheme_choice"]["best_scheme"] == "s1"
    assert os.path.exists(os.path.join(outdir, "depth_across_genome.pdf"))
    assert os.path.exists(os.path.join(outdir, "score_plot.pdf"))
    utils.syscall(f"rm -r {root_outdir}")
