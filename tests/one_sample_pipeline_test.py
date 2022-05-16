import copy
import os
import pytest
import random
import subprocess

import pyfastaq

from viridian_workflow import utils, run, primers

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "one_sample_pipeline")


def run_one_sample(
    platform,
    outdir,
    ref,
    fq1,
    fq2=None,
    tsv_of_amp_schemes=None,
    keep_intermediate=None,
):
    amplicon_sets = []
    for line in open(tsv_of_amp_schemes):
        name, path = line.strip().split("\t")
        if name == "Name":
            continue
        amplicon_sets.append(primers.AmpliconSet.from_tsv(path, name=name))
    fqs = [fq1] if fq2 is None else [fq1, fq2]
    run.run_pipeline(outdir, "illumina", fqs, amplicon_sets, ref=ref)


def perfect_paired_reads_from_sublist(
    seq, start, end, read_length, read_pairs, out1, out2
):
    assert read_length < end - start
    qual_string = "I" * read_length
    with open(out1, "w") as f1, open(out2, "w") as f2:
        for i in range(read_pairs):
            read1 = pyfastaq.sequences.Fastq(
                f"{start}-{end}.{i}/1",
                "".join(seq[start : start + read_length]),
                qual_string,
            )
            read2 = pyfastaq.sequences.Fastq(
                f"{start}-{end}.{i}/2",
                "".join(seq[end - read_length : end]),
                qual_string,
            )
            read2.revcomp()
            print(read1, file=f1)
            print(read2, file=f2)


def perfect_unpaired_reads_from_sublist(seq, start, end, reads, out):
    qual_string = "I" * (end - start)
    with open(out, "w") as f:
        for i in range(reads):
            read = pyfastaq.sequences.Fastq(
                f"{start}-{end}.{i}", "".join(seq[start:end]), qual_string,
            )
            if i % 2 == 0:
                read.revcomp()
            print(read, file=f)


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


def make_amplicons(amps_outfile, schemes_outfile, ref_seq, amplicons=None):
    if amplicons is None:
        amplicons = [
            ("amplicon1", 30, 410),
            ("amplicon2", 320, 726),
            ("amplicon3", 642, 1028),
            ("amplicon4", 943, 1337),
            ("amplicon5", 1242, 1651),
        ]

    primer_length = 10
    with open(amps_outfile, "w") as f:
        print(
            "Amplicon_name",
            "Primer_name",
            "Left_or_right",
            "Sequence",
            "Position",
            sep="\t",
            file=f,
        )
        for name, start, end in amplicons:
            left_seq = "".join(ref_seq[start : start + primer_length])
            right_seq = pyfastaq.sequences.Fasta(
                "x", "".join(ref_seq[end - primer_length : end + 1])
            )
            right_seq.revcomp()
            print(name, name + "_left", "left", left_seq, start, sep="\t", file=f)
            print(
                name,
                name + "_right",
                "right",
                right_seq.seq,
                end - primer_length,
                sep="\t",
                file=f,
            )

    with open(schemes_outfile, "w") as f:
        print("Name\tFile", file=f)
        print("scheme1", os.path.abspath(amps_outfile), sep="\t", file=f)

    return amplicons


def make_fwd_rev_reads_for_each_amplicon(ref_seq, amplicons, outprefix):
    filenames = {}
    for name, start, end in amplicons:
        out1 = f"{outprefix}.{name}.1.fastq"
        out2 = f"{outprefix}.{name}.2.fastq"
        filenames[name] = {"fwd": out1, "rev": out2}
        perfect_paired_reads_from_sublist(ref_seq, start, end, 250, 200, out1, out2)
    return filenames


def make_unpaired_reads_for_each_amplicon(ref_seq, amplicons, outprefix):
    filenames = {}
    for name, start, end in amplicons:
        out = f"{outprefix}.{name}.fastq"
        perfect_unpaired_reads_from_sublist(ref_seq, start, end, 200, out)
        filenames[name] = out
    return filenames


def make_catted_paired_reads_for_amplicon_set(test_data, wanted_amplicons, out1, out2):
    fwd_files = []
    rev_files = []
    for name, _, _ in test_data["amplicons"]:
        if name in wanted_amplicons:
            fwd_files.append(test_data["reads"][name]["fwd"])
            rev_files.append(test_data["reads"][name]["rev"])

    subprocess.check_output("rm -rf {out1} {out2}", shell=True)
    subprocess.check_output(f"cat {' ' .join(fwd_files)} > {out1}", shell=True)
    subprocess.check_output(f"cat {' ' .join(rev_files)} > {out2}", shell=True)


def make_catted_unpaired_reads_for_amplicon_set(test_data, wanted_amplicons, out):
    files = []
    for name, _, _ in test_data["amplicons"]:
        if name in wanted_amplicons:
            files.append(test_data["unpaired_reads"][name])

    subprocess.check_output("rm -rf {out}", shell=True)
    subprocess.check_output(f"cat {' ' .join(files)} > {out}", shell=True)


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
        "amplicons_tsv": os.path.join(outdir, "amplicons.tsv"),
        "schemes_tsv": os.path.join(outdir, "amplicon_schemes.tsv"),
        "ref_fasta": os.path.join(outdir, "ref.fasta"),
    }
    data["ref_seq"] = (
        ["A"] * 20 + random.choices(["A", "C", "G", "T"], k=1680) + ["A"] * 20
    )
    data["amplicons"] = make_amplicons(
        data["amplicons_tsv"], data["schemes_tsv"], data["ref_seq"]
    )

    reads_prefix = os.path.join(outdir, "reads")
    data["reads"] = make_fwd_rev_reads_for_each_amplicon(
        data["ref_seq"], data["amplicons"], reads_prefix
    )
    data["unpaired_reads"] = make_unpaired_reads_for_each_amplicon(
        data["ref_seq"], data["amplicons"], reads_prefix
    )

    nucleotides_list_to_fasta_file(data["ref_seq"], "ref", data["ref_fasta"])
    yield data
    subprocess.check_output(f"rm -rf {outdir}", shell=True)


def _test_complete_assembly_no_reads_map(test_data):
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

    amplicon_sets = []
    ref = test_data["ref_fasta"]
    for line in open(test_data["amplicons_tsv"]):
        print(line)
        assert False
        name, path = line.strip().split("\t")
        amplicon_sets.append((name, path))

    try:
        run_one_sample(
            "illumina",
            outdir,
            test_data["ref_fasta"],
            fq1,
            fq2=fq2,
            amplicon_json=test_data["amplicons_tsv"],
        )
        # This test should fail on viridian, producing no consensus
        # TODO specify that it was the consensus file that's missing
    except Exception as error:
        if str(error) != "failed to choose amplicon scheme":
            raise error

    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)


def test_complete_assembly_from_all_good_amplicons(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.complete_assembly_from_all_good_amplicons"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    all_amplicon_names = set([x[0] for x in test_data["amplicons"]])
    make_catted_paired_reads_for_amplicon_set(test_data, all_amplicon_names, fq1, fq2)
    outdir = f"{pre_out}.out"
    run_one_sample(
        "illumina",
        outdir,
        test_data["ref_fasta"],
        fq1,
        fq2=fq2,
        tsv_of_amp_schemes=test_data["schemes_tsv"],
        keep_intermediate=True,
    )
    # TODO: check that we got the expected output
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)


def test_complete_assembly_from_all_good_amplicons_unpaired(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.complete_assembly_from_all_good_amplicons_unpaired"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq = f"{pre_out}.fq"
    all_amplicon_names = set([x[0] for x in test_data["amplicons"]])
    make_catted_unpaired_reads_for_amplicon_set(test_data, all_amplicon_names, fq)
    outdir = f"{pre_out}.out"
    run_one_sample(
        "ont",
        outdir,
        test_data["ref_fasta"],
        fq,
        tsv_of_amp_schemes=test_data["schemes_tsv"],
        keep_intermediate=True,
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
    make_catted_paired_reads_for_amplicon_set(test_data, amplicon_names, fq1, fq2)
    outdir = f"{pre_out}.out"
    run_one_sample(
        "illumina",
        outdir,
        test_data["ref_fasta"],
        fq1,
        fq2=fq2,
        tsv_of_amp_schemes=test_data["schemes_tsv"],
        keep_intermediate=True,
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
    make_catted_paired_reads_for_amplicon_set(test_data, all_amplicon_names, fq1, fq2)
    outdir = f"{pre_out}.out"
    ref_fasta = f"{pre_out}.ref.fa"
    ref_seq = copy.copy(test_data["ref_seq"])
    ref_seq[500] = "A" if ref_seq[500] != "A" else "G"
    ref_seq.insert(1100, "A")
    ref_seq.pop(1200)
    nucleotides_list_to_fasta_file(ref_seq, "ref", ref_fasta)
    run_one_sample(
        "illumina",
        outdir,
        ref_fasta,
        fq1,
        fq2=fq2,
        tsv_of_amp_schemes=test_data["schemes_tsv"],
        keep_intermediate=True,
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
    try:
        run_one_sample(
            "illumina",
            outdir,
            test_data["ref_fasta"],
            fq1,
            fq2=fq2,
            tsv_of_amp_schemes=test_data["schemes_tsv"],
            keep_intermediate=True,
        )
    # we should throw fail on detecting the amplicon scheme
    except Exception as error:
        if str(error) != "failed to choose amplicon scheme":
            raise error
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)


# TODO: at the time of writing, this test fails because viridian makes no
# output, because all the amplicons are failed. It hits this error:
# viridian_workflow.utils.OutputFileError: /viridian_workflow/tmp.not_expected_amplicons.out/viridian/consensus.final_assembly.fa
def _test_not_expected_amplicons(test_data):
    assert os.path.exists(test_data["dirname"])
    pre_out = "tmp.not_expected_amplicons"
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
    fq1 = f"{pre_out}.1.fq"
    fq2 = f"{pre_out}.2.fq"
    all_amplicon_names = set([x[0] for x in test_data["amplicons"]])
    make_catted_paired_reads_for_amplicon_set(test_data, all_amplicon_names, fq1, fq2)
    amplicons = [
        ("amplicon1", 50, 900),
        ("amplicon2", 850, 1200),
        ("amplicon3", 1165, 1700),
    ]
    amplicons_json = f"{pre_out}.amplicons.json"
    make_amplicons(amplicons_json, amplicons=amplicons)
    outdir = f"{pre_out}.out"
    try:
        run_one_sample(
            "illumina",
            outdir,
            test_data["ref_fasta"],
            fq1,
            fq2=fq2,
            amplicon_json=amplicons_json,
            keep_intermediate=True,
        )
    except:
        pass
    # TODO: check that we got the expected output
    subprocess.check_output(f"rm -rf {pre_out}*", shell=True)
