import os
import pytest
import random
import subprocess
import tempfile

import pyfastaq

from viridian_workflow import one_sample_pipeline

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "one_sample_pipeline")


def perfect_reads_from_sublist(seq, start, end, read_length, read_pairs, out1, out2):
    assert read_length < end - start
    qual_string = "I" * read_length
    with open(out1, "w") as f1, open(out2, "w") as f2:
        for i in range(read_pairs):
            read1 = pyfastaq.sequences.Fastq(
                f"{i}/1", "".join(seq[start : start + read_length]), qual_string
            )
            read2 = pyfastaq.sequences.Fastq(
                f"{i}/2", "".join(seq[end - read_length : end]), qual_string
            )
            read2.revcomp()
            print(read1, file=f1)
            print(read2, file=f2)


def make_amplicons(outfile):
    amplicons = [
        ("amplicon1", 30, 410),
        ("amplicon2", 320, 726),
        ("amplicon3", 642, 1028),
        ("amplicon4", 943, 1337),
        ("amplicon5", 1242, 1651),
    ]

    with open(outfile, "w") as f:
        for amplicon in amplicons:
            print(*amplicon, sep="\t", file=f)

    return amplicons


def make_fwd_rev_reads_for_each_amplicon(ref_seq, amplicons, outprefix):
    filenames = {}
    for name, start, end in amplicons:
        out1 = f"{outprefix}.{name}.1.fastq"
        out2 = f"{outprefix}.{name}.2.fastq"
        filenames[name] = {"fwd": out1, "rev": out2}
        perfect_reads_from_sublist(ref_seq, start, end, 250, 200, out1, out2)
    return filenames


def make_catted_reads_for_amplicon_set(test_data, wanted_amplicons, out1, out2):
    fwd_files = []
    rev_files = []
    import pprint

    pprint.pprint(test_data)
    for name, _, _ in test_data["amplicons"]:
        if name in wanted_amplicons:
            fwd_files.append(test_data["reads"][name]["fwd"])
            rev_files.append(test_data["reads"][name]["rev"])

    subprocess.check_output("rm -rf {out1} {out2}", shell=True)
    subprocess.check_output(f"cat {' ' .join(fwd_files)} > {out1}", shell=True)
    subprocess.check_output(f"cat {' ' .join(rev_files)} > {out2}", shell=True)


@pytest.fixture(scope="session")
def test_data():
    random.seed(42)
    outdir = "tmp.one_sample_pipeline_test_data"
    subprocess.check_output(f"rm -rf {outdir}", shell=True)
    os.mkdir(outdir)
    data = {
        "dirname": outdir,
        "amplicons_bed": os.path.join(outdir, "amplicons.bed"),
        "ref_fasta": os.path.join(outdir, "ref.fasta"),
    }
    data["amplicons"] = make_amplicons(data["amplicons_bed"])
    ref_seq = ["A"] * 20 + random.choices(["A", "C", "G", "T"], k=1680) + ["A"] * 20

    reads_prefix = os.path.join(outdir, "reads")
    data["reads"] = make_fwd_rev_reads_for_each_amplicon(
        ref_seq, data["amplicons"], reads_prefix
    )

    ref_fa = pyfastaq.sequences.Fasta("ref", "".join(ref_seq))
    with open(data["ref_fasta"], "w") as f:
        print(ref_fa, file=f)

    yield data


def test_complete_assembly_from_all_good_amplicons(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.complete_assembly_from_all_good_amplicons"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    all_amplicon_names = set([x[0] for x in test_data["amplicons"]])
    make_catted_reads_for_amplicon_set(test_data, all_amplicon_names, fq1, fq2)
    outdir = f"{pre_out}.out"
    # TODO: get this so it runs, then check the output looks "correct",
    # for whatever we decide "correct" means
    one_sample_pipeline.run_one_sample(
        outdir, test_data["ref_fasta"], test_data["amplicons_bed"], fq1, fq2
    )

    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
