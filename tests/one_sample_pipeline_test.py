import copy
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


def tiling_reads(seq, read_length, frag_length, out1, out2, step=2):
    qual_string = "I" * read_length
    with open(out1, "w") as f1, open(out2, "w") as f2:
        for i in range(0, len(seq) - frag_length, step):
            read1 = pyfastaq.sequences.Fastq(
                f"{i}/1", "".join(seq[i : i + read_length]), qual_string
            )
            read2 = pyfastaq.sequences.Fastq(
                f"{i}/2",
                "".join(seq[i + frag_length - read_length : i + frag_length]),
                qual_string,
            )
            read2.revcomp()
            print(read1, file=f1)
            print(read2, file=f2)


def make_amplicons(outfile, amplicons=None):
    if amplicons is None:
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


def nucleotides_list_to_fasta_file(nuc_list, seq_name, outfile):
    fa = pyfastaq.sequences.Fasta(seq_name, "".join(nuc_list))
    with open(outfile, "w") as f:
        print(fa, file=f)


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
    data["ref_seq"] = (
        ["A"] * 20 + random.choices(["A", "C", "G", "T"], k=1680) + ["A"] * 20
    )

    reads_prefix = os.path.join(outdir, "reads")
    data["reads"] = make_fwd_rev_reads_for_each_amplicon(
        data["ref_seq"], data["amplicons"], reads_prefix
    )

    nucleotides_list_to_fasta_file(data["ref_seq"], "ref", data["ref_fasta"])
    yield data
    subprocess.check_output(f"rm -rf {outdir}", shell=True)


def test_complete_assembly_no_reads_map(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.no_reads_map"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    # make some garbage reads so they don't map
    with open(fq1, "w") as f1, open(fq2, "w") as f2:
        print("@read1/1", "A" * 100, "+", "I" * 100, sep="\n", file=f1)
        print("@read1/2", "A" * 100, "+", "I" * 100, sep="\n", file=f2)
    outdir = f"{pre_out}.out"
    one_sample_pipeline.run_one_sample(
        outdir, test_data["ref_fasta"], test_data["amplicons_bed"], fq1, fq2
    )
    # TODO: check that we got the expected output
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)


def test_complete_assembly_from_all_good_amplicons(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.complete_assembly_from_all_good_amplicons"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    all_amplicon_names = set([x[0] for x in test_data["amplicons"]])
    make_catted_reads_for_amplicon_set(test_data, all_amplicon_names, fq1, fq2)
    outdir = f"{pre_out}.out"
    one_sample_pipeline.run_one_sample(
        outdir, test_data["ref_fasta"], test_data["amplicons_bed"], fq1, fq2
    )
    # TODO: check that we got the expected output
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)


def test_assembly_amplicon_3_no_reads(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.assembly_amplicon_3_no_reads"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    amplicon_names = {"amplicon1", "amplicon2", "amplicon4", "amplicon5"}
    make_catted_reads_for_amplicon_set(test_data, amplicon_names, fq1, fq2)
    outdir = f"{pre_out}.out"
    one_sample_pipeline.run_one_sample(
        outdir, test_data["ref_fasta"], test_data["amplicons_bed"], fq1, fq2
    )
    # TODO: check that we got the expected output
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)


def test_complete_assembly_with_snps_and_indels(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.complete_assembly_with_snps_and_indels"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    all_amplicon_names = set([x[0] for x in test_data["amplicons"]])
    make_catted_reads_for_amplicon_set(test_data, all_amplicon_names, fq1, fq2)
    outdir = f"{pre_out}.out"
    ref_fasta = f"{pre_out}.ref.fa"
    ref_seq = copy.copy(test_data["ref_seq"])
    ref_seq[500] = "A" if ref_seq[500] != "A" else "G"
    ref_seq.insert(1100, "A")
    ref_seq.pop(1200)
    nucleotides_list_to_fasta_file(ref_seq, "ref", ref_fasta)
    one_sample_pipeline.run_one_sample(
        outdir, ref_fasta, test_data["amplicons_bed"], fq1, fq2
    )
    # TODO: check that we got the expected output
    # Should be something like this in the VCF (which doesn't yet exist):
    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
    # ref	501	0	A	C	.	PASS	.	GT	1/1
    # ref	1099	1	GA	G	.	PASS	.	GT	1/1
    # ref	1199	2	C	CG	.	PASS	.	GT	1/1
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)


def test_reads_are_wgs_not_amplicon(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.reads_are_wgs_not_amplicon"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    tiling_reads(test_data["ref_seq"], 150, 350, fq1, fq2, step=2)
    outdir = f"{pre_out}.out"
    one_sample_pipeline.run_one_sample(
        outdir, test_data["ref_fasta"], test_data["amplicons_bed"], fq1, fq2
    )
    # TODO: check that we got the expected output
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)


def test_not_expected_amplicons(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.not_expected_amplicons"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    all_amplicon_names = set([x[0] for x in test_data["amplicons"]])
    make_catted_reads_for_amplicon_set(test_data, all_amplicon_names, fq1, fq2)
    amplicons = [
        ("amplicon1", 50, 900),
        ("amplicon2", 850, 1200),
        ("amplicon3", 1165, 1700),
    ]
    amplicons_bed = f"{pre_out}.amplicons.bed"
    make_amplicons(amplicons_bed, amplicons=amplicons)
    outdir = f"{pre_out}.out"
    one_sample_pipeline.run_one_sample(
        outdir, test_data["ref_fasta"], amplicons_bed, fq1, fq2
    )
    # TODO: check that we got the expected output
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
