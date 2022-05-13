import os
import pytest
import random
import subprocess

import pyfastaq

from viridian_workflow import readstore, primers

# Make a 1kb reference genome with 4 amplicons, and then some read pairs
# that do different things (both mapped, one mapped, mapped to different
# amplicons ... etc)
def make_test_ref_genome(outfile):
    ref = "".join(random.choices(["A", "C", "G", "T"], k=1000))
    with open(outfile, "w") as f:
        print(">ref foo", file=f)
        print(ref, file=f)
    return ref


def make_amplicon_set(outfile, ref):
    assert len(ref) == 1000
    with open(outfile, "w") as out:
        print(
            "\t".join(
                [
                    "Amplicon_name",
                    "Primer_name",
                    "Left_or_right",
                    "Sequence",
                    "Position",
                ]
            ),
            file=out,
        )
        for i, (l_pos, r_pos) in enumerate(
            [(0, 300), (175, 475), (350, 650), (525, 825), (700, 1000)]
        ):
            primer_len = 20
            r_primer = revcomp(ref[r_pos, r_pos + primer_len])
            print(
                "\t".join(
                    [
                        f"amp{i+1}",
                        f"amp{i+1}_left",
                        "left",
                        ref[l_pos, l_pos + primer_len],
                        str(position),
                    ]
                ),
                file=out,
            )
            print(
                "\t".join(
                    [
                        f"amp{i+1}",
                        f"amp{i+1}_right",
                        "right",
                        r_primer,
                        str(position - primer_len),
                    ]
                ),
                file=out,
            )


def write_read(name, ref_seq, start, end, filehandle, revcomp=False):
    if revcomp:
        seq = pyfastaq.sequences.Fastq(name, ref_seq[start:end], "I" * (end - start))
        seq.revcomp()
        print(seq, file=filehandle)
    else:
        print(
            f"@{name}",
            ref_seq[start:end],
            "+",
            "I" * (end - start),
            file=filehandle,
            sep="\n",
        )


def make_test_unpaired_fastq(outfile, ref_seq):
    read_count = 0
    read_len = 100
    with open(outfile, "w") as f:
        for i in range(20):
            write_read(f"{read_count}", ref_seq, i, i + read_len, f)
            read_count += 1
            write_read(f"{read_count}", ref_seq, i, i + read_len, f, revcomp=True)
            read_count += 1


def make_test_paired_fastq(out1, out2, ref_seq):
    with open(out1, "w") as f1, open(out2, "w") as f2:
        read_len = 100
        pair_count = 0
        for i in range(20):
            write_read(f"{pair_count}/1", ref_seq, i, i + read_len, f1)
            write_read(
                f"{pair_count}/2",
                ref_seq,
                200 + i,
                200 + i + read_len,
                f2,
                revcomp=True,
            )
            pair_count += 1


def map_reads(ref, reads1, reads2, bam_out):
    command = (
        f"minimap2 -a -x sr {ref} {reads1} {reads2} | samtools sort -O BAM > {bam_out}"
    )
    subprocess.check_output(command, shell=True)
    subprocess.check_output(f"samtools index {bam_out}", shell=True)


@pytest.fixture(scope="session")
def test_data():
    random.seed(42)
    outdir = "tmp.sample_reads_test_data"
    subprocess.check_output(f"rm -rf {outdir}", shell=True)
    os.mkdir(outdir)
    data = {
        "dirname": outdir,
        "amplicons_tsv": os.path.join(outdir, "amplicons.tsv"),
        "ref_fasta": os.path.join(outdir, "ref.fasta"),
        "paired_fq1": os.path.join(outdir, "paired_reads_1.fq"),
        "paired_fq2": os.path.join(outdir, "paired_reads_2.fq"),
        "unpaired_fq": os.path.join(outdir, "unpaired_reads.fq"),
        "paired_bam": os.path.join(outdir, "paired_reads.bam"),
        "unpaired_bam": os.path.join(outdir, "unpaired_reads.bam"),
    }

    data["ref_seq"] = make_test_ref_genome(data["ref_fasta"])
    make_amplicons_tsv(data["amplicons_tsv"], data["ref_seq"])
    make_test_unpaired_fastq(data["unpaired_fq"], data["ref_seq"])
    make_test_paired_fastq(data["paired_fq1"], data["paired_fq2"], data["ref_seq"])
    map_reads(
        data["ref_fasta"], data["paired_fq1"], data["paired_fq2"], data["paired_bam"]
    )
    map_reads(data["ref_fasta"], data["unpaired_fq"], "", data["unpaired_bam"])
    yield data
    subprocess.check_output(f"rm -rf {outdir}", shell=True)


def test_sample_unpaired_reads(test_data):
    outprefix = "tmp.sample_unpaired_reads.out"
    subprocess.check_output(f"rm -f {outprefix}.*", shell=True)
    sampler = sample_reads.sample_reads(
        test_data["ref_fasta"],
        test_data["unpaired_bam"],
        outprefix,
        test_data["amplicons_bed"],
        5,
    )

    # TODO: check contents of files, and add more cases to the type of
    # reads in the BAM we're sampling from.
    assert sampler.number_of_amplicons() == 4
    assert os.path.exists(sampler.bam_out)
    assert os.path.exists(sampler.fq_out)
    assert not os.path.exists(sampler.fq_out1)
    assert not os.path.exists(sampler.fq_out2)
    subprocess.check_output(f"rm -f {outprefix}.*", shell=True)


def test_sample_paired_reads(test_data):
    outprefix = "tmp.sample_paired_reads.out"
    subprocess.check_output(f"rm -f {outprefix}.*", shell=True)
    amplicons = primers.AmpliconSet.from_tsv(test_data["amplicons_bed"])
    bam = readstore.Bam(test_data["paired_bam"])
    print(bam)
    print(amplicons)
    reads = readstore.Readstore(bam, amplicons)

    sampler = sample_reads.sample_reads(
        test_data["ref_fasta"],
        test_data["paired_bam"],
        outprefix,
        test_data["amplicons_bed"],
        5,
        min_sampled_depth_for_pass=5,
    )

    # TODO: check contents of files, and add more cases to the type of
    # reads in the BAM we're sampling from.
    assert sampler.number_of_amplicons() == 4
    assert os.path.exists(sampler.bam_out)
    assert not os.path.exists(sampler.fq_out)
    assert os.path.exists(sampler.fq_out1)
    assert os.path.exists(sampler.fq_out2)
    assert os.path.exists(sampler.failed_amps_file)
    with open(sampler.failed_amps_file) as f:
        got = [x.rstrip() for x in f]
    assert got == ["amp2", "amp3", "amp4"]
    assert sampler.failed_amplicons == 3
    subprocess.check_output(f"rm -f {outprefix}.*", shell=True)
