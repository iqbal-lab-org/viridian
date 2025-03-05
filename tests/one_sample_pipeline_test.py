import csv
import gzip
import itertools
import os
import pytest
import random

import pyfastaq

from viridian import one_sample_pipeline, utils

read_counter = 0


def write_amplicons_tsv(test_data, scheme_name, outfile):
    with open(outfile, "w") as f:
        print(
            "Amplicon_name",
            "Primer_name",
            "Left_or_right",
            "Sequence",
            "Position",
            sep="\t",
            file=f,
        )
        for name, d in test_data["amplicons"][scheme_name].items():
            test_data["reads_files"][name] = {}
            for i, (start, end) in enumerate(d["primers"]["left"]):
                print(
                    name,
                    f"{name}_L_{i}",
                    "left",
                    "A" * (end - start + 1),
                    start,
                    sep="\t",
                    file=f,
                )
            for i, (start, end) in enumerate(d["primers"]["right"]):
                print(
                    name,
                    f"{name}_L_{i}",
                    "right",
                    "A" * (end - start + 1),
                    start,
                    sep="\t",
                    file=f,
                )


def reads_for_one_amplicon(
    test_data,
    amp_name,
    left_primer,
    right_primer,
    outprefix,
    mutations={},
    coverage=100,
    scheme="scheme1",
):
    global read_counter
    amp = test_data["amplicons"][scheme][amp_name]
    amp_start = amp["primers"]["left"][left_primer][0]
    amp_end = amp["primers"]["right"][right_primer][1]
    amp_seq = test_data["ref_seq"][amp_start : amp_end + 1]

    if len(mutations) > 0:
        amp_seq = list(amp_seq)
        for ref_start, ref_end, alt in reversed(sorted(mutations)):
            assert amp_start < ref_start <= ref_end < amp_end
            amp_seq[ref_start - amp_start : 1 + ref_end - amp_start] = alt
        amp_seq = "".join(amp_seq)

    amp_seq_rev = pyfastaq.sequences.Fasta("x", amp_seq)
    amp_seq_rev.revcomp()

    with open(f"{outprefix}.fa", "w") as f:
        for i in range(coverage):
            print(f">{read_counter}.{amp_name}.{amp_start}.{amp_end}.f", file=f)
            print(amp_seq, file=f)
            print(f">{read_counter}.{amp_name}.{amp_start}.{amp_end}.r", file=f)
            print(amp_seq_rev.seq, file=f)
            read_counter += 2

    fwd_read = amp_seq[:200]
    rev_read = pyfastaq.sequences.Fasta("x", amp_seq[-200:])
    rev_read.revcomp()
    with open(f"{outprefix}.1.fa", "w") as f1, open(f"{outprefix}.2.fa", "w") as f2:
        for i in range(coverage):
            print(f">{read_counter}.{amp_name}.{amp_start}.{amp_end}/1", file=f1)
            print(fwd_read, file=f1)
            print(f">{read_counter}.{amp_name}.{amp_start}.{amp_end}/2", file=f2)
            print(rev_read.seq, file=f2)
            read_counter += 1


def fastas_of_all_perfect_reads(test_data, outprefix, amps_to_exclude=None):
    if amps_to_exclude is None:
        amps_to_exclude = set()
    all_ont_files = []
    all_ilm1_files = []
    all_ilm2_files = []
    for amp, d in test_data["reads_files"].items():
        if amp in amps_to_exclude:
            continue
        all_ont_files.extend(f"{x}.fa" for x in d.values())
        all_ilm1_files.extend(f"{x}.1.fa" for x in d.values())
        all_ilm2_files.extend(f"{x}.2.fa" for x in d.values())
    all_ont_files = " ".join(all_ont_files)
    all_ilm1_files = " ".join(all_ilm1_files)
    all_ilm2_files = " ".join(all_ilm2_files)
    utils.syscall(f"cat {all_ont_files} > {outprefix}.fa")
    utils.syscall(f"cat {all_ilm1_files} > {outprefix}.1.fa")
    utils.syscall(f"cat {all_ilm2_files} > {outprefix}.2.fa")


def load_qc_tsv(filename):
    results = {}
    with gzip.open(filename, "rt") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            key = (int(d["Ref_pos"]) - 1, int(d["Cons_pos"]) - 1)
            assert key not in results
            results[key] = d
    return results


@pytest.fixture(scope="session")
def test_data():
    random.seed(42)
    outdir = "tmp.one_sample_pipeline_test_data"
    utils.syscall(f"rm -rf {outdir}")
    os.mkdir(outdir)
    ref_list = ["A"] * 25 + random.choices(["A", "C", "G", "T"], k=1600) + ["A"] * 25
    ref_list[80] = "A"
    ref_list[100] = "T"
    ref_list[120] = "G"
    ref_list[140] = "G"
    ref_list[160] = "G"
    ref_list[180] = "A"
    ref_list[200] = "G"
    ref_list[220] = "G"
    ref_list[239] = "G"
    ref_list[240] = "A"
    ref_list[241] = "T"
    ref_list[259] = "G"
    ref_list[260] = "A"
    ref_list[261] = "T"
    ref_list[279] = "G"
    ref_list[280] = "A"
    ref_list[281] = "T"
    ref_seq = pyfastaq.sequences.Fasta("ref", "".join(ref_list))
    ref_seq_rev = pyfastaq.sequences.Fasta("ref_rev", "".join(ref_list))
    ref_seq_rev.revcomp()
    perfect_reads_prefix = os.path.join(outdir, "reads.perfect")
    data = {
        "ref_seq": ref_seq,
        "ref_fasta": os.path.join(outdir, "ref.fa"),
        "amps1_tsv": os.path.join(outdir, "amps.1.tsv"),
        "amps2_tsv": os.path.join(outdir, "amps.2.tsv"),
        "amp_schemes_tsv": os.path.join(outdir, "amp_schemes.tsv"),
        "reads_dir": os.path.join(outdir, "Reads"),
        "reads_files": {},
        "perfect_reads_ont": f"{perfect_reads_prefix}.fa",
        "perfect_reads_ilm1": f"{perfect_reads_prefix}.1.fa",
        "perfect_reads_ilm2": f"{perfect_reads_prefix}.2.fa",
    }
    with open(data["amp_schemes_tsv"], "w") as f:
        print("Name", "File", sep="\t", file=f)
        print("scheme1", data["amps1_tsv"], sep="\t", file=f)
        print("scheme2", data["amps2_tsv"], sep="\t", file=f)

    os.mkdir(data["reads_dir"])
    with open(data["ref_fasta"], "w") as f:
        print(ref_seq, file=f)

    amplicons1 = {
        "amp1": {
            "start": 50,
            "end": 350,
            "primers": {"left": [(50, 60)], "right": [(340, 350)]},
        },
        "amp2": {
            "start": 300,
            "end": 650,
            "primers": {
                "left": [(300, 310), (320, 330)],
                "right": [(620, 630), (640, 650)],
            },
        },
        "amp3": {
            "start": 600,
            "end": 950,
            "primers": {"left": [(590, 600)], "right": [(940, 950)]},
        },
        "amp4": {
            "start": 900,
            "end": 1250,
            "primers": {"left": [(900, 910)], "right": [(1240, 1250)]},
        },
        "amp5": {
            "start": 1200,
            "end": 1550,
            "primers": {"left": [(1200, 1210)], "right": [(1540, 1550)]},
        },
    }
    amplicons2 = {
        "amp1": {
            "start": 100,
            "end": 1000,
            "primers": {"left": [(90, 100)], "right": [(340, 350)]},
        },
        "amp2": {
            "start": 800,
            "end": 1500,
            "primers": {"left": [(790, 800)], "right": [(1490, 1500)]},
        },
    }
    data["amplicons"] = {"scheme1": amplicons1, "scheme2": amplicons2}
    write_amplicons_tsv(data, "scheme1", data["amps1_tsv"])
    write_amplicons_tsv(data, "scheme2", data["amps2_tsv"])
    for name, d in amplicons1.items():
        left_primers = d["primers"]["left"]
        right_primers = d["primers"]["right"]
        for t in itertools.product(range(len(left_primers)), range(len(right_primers))):
            out = os.path.join(data["reads_dir"], f"{name}.perfect.{t[0]}.{t[1]}")
            data["reads_files"][name][t] = out
            reads_for_one_amplicon(data, name, t[0], t[1], out)

    fastas_of_all_perfect_reads(data, perfect_reads_prefix)
    yield data
    utils.syscall(f"rm -r {outdir}")


def test_reads_file_not_found(test_data):
    reads_fq = "tmp.test_pipeline.reads_file_not_found.fq"
    outdir = "tmp.pipeline.reads_file_not_found.out"
    utils.syscall(f"rm -rf {outdir} {reads_fq}")
    with pytest.raises(FileNotFoundError):
        one_sample_pipeline.run_one_sample(
            "sars-cov-2",
            "ont",
            outdir,
            test_data["ref_fasta"],
            reads_file=reads_fq,
            gzip_files=False,
        )
    log = utils.load_json(os.path.join(outdir, "log.json"))
    assert log["run_summary"]["result"] == "Fail"
    assert log["run_summary"]["errors"][0] == "Error during stage: 1/10 Start pipeline"
    assert log["run_summary"]["last_stage_completed"] == "Finished"
    utils.syscall(f"rm -r {outdir}")


def test_empty_reads_file(test_data):
    reads_fq = "tmp.test_pipeline.reads.fq"
    outdir = "tmp.pipeline.empty_reads_file.out"
    with open(reads_fq, "w"):
        pass
    utils.syscall(f"rm -rf {outdir}")
    one_sample_pipeline.run_one_sample(
        "sars-cov-2",
        "ont",
        outdir,
        test_data["ref_fasta"],
        reads_file=reads_fq,
        gzip_files=False,
    )
    log = utils.load_json(os.path.join(outdir, "log.json"))
    assert log["run_summary"]["result"] == "Fail"
    assert log["run_summary"]["errors"][0].startswith("Not enough read coverage")
    assert log["run_summary"]["last_stage_completed"] == "Finished"
    utils.syscall(f"rm -r {outdir} {reads_fq}")


def test_scheme_not_match_reference(test_data):
    outdir = "tmp.pipeline.scheme_not_match_reference.out"
    utils.syscall(f"rm -rf {outdir}")
    one_sample_pipeline.run_one_sample(
        "sars-cov-2",
        "ont",
        outdir,
        test_data["ref_fasta"],
        reads_file=test_data["perfect_reads_ont"],
        gzip_files=False,
    )
    log = utils.load_json(os.path.join(outdir, "log.json"))
    assert log["run_summary"]["result"] == "Fail"
    assert log["run_summary"]["errors"][0].startswith(
        "Scheme does not match reference genome"
    )
    assert log["run_summary"]["last_stage_completed"] == "Finished"
    utils.syscall(f"rm -r {outdir}")


def test_all_perfect_reads_ont(test_data):
    outdir = "tmp.pipeline.all_perfect_reads_ont"
    utils.syscall(f"rm -rf {outdir}")
    one_sample_pipeline.run_one_sample(
        "sars-cov-2",
        "ont",
        outdir,
        test_data["ref_fasta"],
        reads_file=test_data["perfect_reads_ont"],
        tsv_of_amp_schemes=test_data["amp_schemes_tsv"],
        gzip_files=False,
    )
    log = utils.load_json(os.path.join(outdir, "log.json"))
    assert log["run_summary"]["result"] == "Success"
    assert log["run_summary"]["errors"] == []
    assert log["run_summary"]["last_stage_completed"] == "Finished"
    assert log["amplicon_scheme_name"] == "scheme1"
    cons = log["sequences"]["masked_consensus"]
    # Primers all have length 11.
    # Should have got first and last primer masked and no other Ns
    primer_len = 11
    assert cons.startswith("N" * primer_len)
    assert cons.endswith("N" * primer_len)
    assert "N" not in cons[primer_len:-primer_len]
    assert cons[primer_len:-primer_len] in test_data["ref_seq"].seq
    assert len(cons) >= 0.9 * len(test_data["ref_seq"])
    utils.syscall(f"rm -r {outdir}")


def test_all_perfect_reads_ilm(test_data):
    outdir = "tmp.pipeline.all_perfect_reads_ilm"
    utils.syscall(f"rm -rf {outdir}")
    one_sample_pipeline.run_one_sample(
        "sars-cov-2",
        "ont",
        outdir,
        test_data["ref_fasta"],
        reads_file1=test_data["perfect_reads_ilm1"],
        reads_file2=test_data["perfect_reads_ilm2"],
        tsv_of_amp_schemes=test_data["amp_schemes_tsv"],
        gzip_files=False,
    )
    log = utils.load_json(os.path.join(outdir, "log.json"))
    assert log["run_summary"]["result"] == "Success"
    assert log["run_summary"]["errors"] == []
    assert log["run_summary"]["last_stage_completed"] == "Finished"
    assert log["amplicon_scheme_name"] == "scheme1"
    cons = log["sequences"]["masked_consensus"]
    # Primers all have length 11.
    # Should have got first and last primer masked and no other Ns
    primer_len = 11
    assert cons.startswith("N" * primer_len)
    assert cons.endswith("N" * primer_len)
    assert "N" not in cons[primer_len:-primer_len]
    assert cons[primer_len:-primer_len] in test_data["ref_seq"].seq
    assert len(cons) >= 0.9 * len(test_data["ref_seq"])
    utils.syscall(f"rm -r {outdir}")


def test_all_perfect_reads_dropped_amp(test_data):
    outprefix = "tmp.pipeline.perfect_reads_dropped_amp"
    utils.syscall(f"rm -rf {outprefix}*")
    ref_length = len(test_data["ref_seq"])
    amplicons = test_data["amplicons"]["scheme1"]

    for amp in test_data["amplicons"]["scheme1"]:
        reads_prefix = f"{outprefix}.{amp}.reads"
        fastas_of_all_perfect_reads(test_data, reads_prefix, amps_to_exclude={amp})
        if amp in {"amp1", "amp5"}:
            expect_cons_len = ref_length - (
                1 + amplicons[amp]["end"] - amplicons[amp]["start"]
            )
        else:
            expect_cons_len = ref_length

        for tech in "ont", "illumina":
            out = f"{outprefix}.{amp}.{tech}"
            if tech == "ont":
                one_sample_pipeline.run_one_sample(
                    "sars-cov-2",
                    tech,
                    out,
                    test_data["ref_fasta"],
                    reads_file=f"{reads_prefix}.fa",
                    tsv_of_amp_schemes=test_data["amp_schemes_tsv"],
                    gzip_files=False,
                )
            else:
                one_sample_pipeline.run_one_sample(
                    "sars-cov-2",
                    tech,
                    out,
                    test_data["ref_fasta"],
                    reads_file1=f"{reads_prefix}.1.fa",
                    reads_file2=f"{reads_prefix}.2.fa",
                    tsv_of_amp_schemes=test_data["amp_schemes_tsv"],
                    gzip_files=False,
                )

            log = utils.load_json(os.path.join(out, "log.json"))
            assert log["run_summary"]["result"] == "Success"
            assert log["run_summary"]["errors"] == []
            assert log["run_summary"]["last_stage_completed"] == "Finished"
            assert log["amplicon_scheme_name"] == "scheme1"
            cons = log["sequences"]["masked_consensus"]
            assert len(cons) >= 0.9 * expect_cons_len
            for test_amp, test_amp_d in amplicons.items():
                amp_seq = test_data["ref_seq"][
                    test_amp_d["start"] : test_amp_d["end"] + 1
                ]
                end_trim = int(0.05 * len(amp_seq))
                amp_seq = amp_seq[end_trim:-end_trim]
                if test_amp == amp:
                    assert amp_seq not in cons
                else:
                    assert amp_seq in cons

            utils.syscall(f"rm -r {out}")

        utils.syscall(f"rm -r {reads_prefix}.*")


def test_failed_amps_or_too_many_Ns(test_data):
    outprefix = "tmp.pipeline.failed_amps_or_too_many_Ns"
    utils.syscall(f"rm -rf {outprefix}*")
    fastas_of_all_perfect_reads(
        test_data, outprefix, amps_to_exclude={"amp1", "amp2", "amp3"}
    )
    # Gradually change options (sometimes silly values needed to force erorrs!),
    # checking most of the points where the pipeline can stop do really happen
    options_to_change = [
        (None, "Not enough read coverage."),
        ({"coverage_min_pc": 30.0}, "Too many failed amplicons from assembler"),
        (
            {"max_percent_amps_fail": 0.0, "max_cons_n_percent": -1},
            "Too many Ns in assembler consensus:",
        ),
        (
            {
                "max_percent_amps_fail": 100.0,
                "max_cons_n_percent": 50.0,
                "masking_min_depth": 200,
            },
            "Too many Ns in final consensus",
        ),
    ]

    for tech in "ont", "illumina":
        options = {
            "tsv_of_amp_schemes": test_data["amp_schemes_tsv"],
            "coverage_min_pc": 50,
            "gzip_files": False,
        }
        out = f"{outprefix}.{tech}"
        if tech == "ont":
            options["reads_file"] = f"{outprefix}.fa"
            options["reads_file1"] = None
            options["reads_file2"] = None
        else:
            options["reads_file"] = None
            options["reads_file1"] = f"{outprefix}.1.fa"
            options["reads_file2"] = f"{outprefix}.2.fa"

        for opt_dict, error_message in options_to_change:
            if opt_dict is not None:
                options.update(opt_dict)
            one_sample_pipeline.run_one_sample(
                "sars-cov-2",
                tech,
                out,
                test_data["ref_fasta"],
                **options,
            )
            log = utils.load_json(os.path.join(out, "log.json"))
            assert log["run_summary"]["result"] == "Fail"
            assert log["run_summary"]["errors"][0].startswith(error_message)
            utils.syscall(f"rm -r {out}")

    utils.syscall(f"rm {outprefix}*")


def test_variants_and_force_consensus(test_data):
    outprefix = "tmp.pipeline.variants"
    utils.syscall(f"rm -rf {outprefix}*")
    reads_prefix = f"{outprefix}.reads.all"
    mut_reads_prefix = f"{outprefix}.reads.mut"
    fastas_of_all_perfect_reads(test_data, reads_prefix, amps_to_exclude={"amp1"})
    common_mutations = [
        (80, 80, "T"),
        (180, 180, "AG"),
        (240, 240, ""),
    ]
    mutations = [
        [
            (100, 100, "T"),
            (120, 120, "A"),
            (140, 140, "C"),
            (200, 200, "GA"),
            (260, 260, ""),
        ],
        [
            (100, 100, "G"),
            (120, 120, "C"),
            (160, 160, "C"),
            (220, 220, "GA"),
            (280, 280, ""),
        ],
    ]
    for i, muts in enumerate(mutations):
        out = f"{mut_reads_prefix}.{i}"
        reads_for_one_amplicon(
            test_data,
            "amp1",
            0,
            0,
            out,
            mutations=common_mutations + muts,
            coverage=25 + 50 * i,
            scheme="scheme1",
        )
        utils.syscall(f"cat {out}.1.fa >> {reads_prefix}.1.fa")
        utils.syscall(f"cat {out}.2.fa >> {reads_prefix}.2.fa")
        utils.syscall(f"cat {out}.fa >> {reads_prefix}.fa")
        utils.syscall(f"rm {out}.1.fa {out}.2.fa {out}.fa")

    for tech in "ont", "illumina":
        out = f"{outprefix}.run.{tech}"
        if tech == "ont":
            one_sample_pipeline.run_one_sample(
                "sars-cov-2",
                tech,
                out,
                test_data["ref_fasta"],
                reads_file=f"{reads_prefix}.fa",
                tsv_of_amp_schemes=test_data["amp_schemes_tsv"],
                fix_small_indels=False,
                gzip_files=False,
            )
        else:
            one_sample_pipeline.run_one_sample(
                "sars-cov-2",
                tech,
                out,
                test_data["ref_fasta"],
                reads_file1=f"{reads_prefix}.1.fa",
                reads_file2=f"{reads_prefix}.2.fa",
                tsv_of_amp_schemes=test_data["amp_schemes_tsv"],
                gzip_files=False,
            )

        log = utils.load_json(os.path.join(out, "log.json"))
        assert log["run_summary"]["result"] == "Success"
        assert log["run_summary"]["errors"] == []
        assert log["run_summary"]["last_stage_completed"] == "Finished"
        assert log["amplicon_scheme_name"] == "scheme1"
        qc = load_qc_tsv(os.path.join(out, "qc.tsv.gz"))
        cons = utils.load_single_seq_fasta(os.path.join(out, "consensus.fa"))

        assert qc[(80, 30)]["Ref_nt"] == "A"
        assert qc[(80, 30)]["Cons_nt"] == "T"
        assert cons[30] == "T"
        assert qc[(80, 30)]["Mask"] == "PASS"

        assert qc[(100, 50)]["Ref_nt"] == "T"
        assert qc[(100, 50)]["Cons_nt"] == "G"
        assert cons[50] == "K"
        assert qc[(100, 50)]["Mask"] == "HET"

        assert qc[(120, 70)]["Ref_nt"] == "G"
        assert qc[(120, 70)]["Cons_nt"] == "C"
        assert cons[70] == "M"
        assert qc[(120, 70)]["Mask"] == "HET"

        assert qc[(140, 90)]["Ref_nt"] == "G"
        assert qc[(140, 90)]["Cons_nt"] == "G"
        assert cons[90] == "S"
        assert qc[(140, 90)]["Mask"] == "HET"

        assert qc[(160, 110)]["Ref_nt"] == "G"
        assert qc[(160, 110)]["Cons_nt"] == "C"
        assert cons[110] == "S"
        assert qc[(160, 110)]["Mask"] == "HET"

        assert qc[(180, 130)]["Ref_nt"] == "A"
        assert qc[(180, 130)]["Cons_nt"] == "A"
        assert qc[(180, 131)]["Ref_nt"] == "-"
        assert qc[(180, 131)]["Cons_nt"] == "G"
        assert qc[(181, 132)]["Ref_nt"] == "C"
        assert qc[(181, 132)]["Cons_nt"] == "C"
        assert qc[(200, 151)]["Ref_nt"] == "G"
        assert qc[(200, 151)]["Cons_nt"] == "G"
        assert qc[(220, 171)]["Ref_nt"] == "G"
        assert qc[(220, 171)]["Cons_nt"] == "G"
        assert qc[(220, 172)]["Ref_nt"] == "-"
        assert qc[(220, 172)]["Cons_nt"] == "A"
        assert qc[(239, 191)]["Ref_nt"] == "G"
        assert qc[(239, 191)]["Cons_nt"] == "G"
        assert qc[(240, 191)]["Ref_nt"] == "A"
        assert qc[(240, 191)]["Cons_nt"] == "-"
        assert qc[(241, 192)]["Ref_nt"] == "T"
        assert qc[(241, 192)]["Cons_nt"] == "T"
        assert qc[(260, 211)]["Ref_nt"] == "A"
        assert qc[(260, 211)]["Cons_nt"] == "A"
        assert qc[(279, 230)]["Ref_nt"] == "G"
        assert qc[(279, 230)]["Cons_nt"] == "G"
        assert qc[(280, 230)]["Ref_nt"] == "A"
        assert qc[(280, 230)]["Cons_nt"] == "-"
        assert qc[(281, 231)]["Ref_nt"] == "T"
        assert qc[(281, 231)]["Cons_nt"] == "T"

        utils.syscall(f"rm -r {out}")

    # Now test force consensus option. We use the reference genome as the
    # forced consensus genome, but add in ambiguous bases.
    # Use the same reads as above. This means
    # that the reads will disagree with the consensus, so we can check
    # what happens at snps, indels, and ambiguous ref bases
    ref_seq = list(test_data["ref_seq"].seq)
    ref_seq[100] = "K"
    ref_seq[160] = "N"
    ref_seq = pyfastaq.sequences.Fasta("force", "".join(ref_seq))
    ref_fasta = f"{outprefix}.force_ref.fa"
    with open(ref_fasta, "w") as f:
        print(ref_seq, file=f)

    for tech in "ont", "illumina":
        out = f"{outprefix}.run.{tech}"
        if tech == "ont":
            one_sample_pipeline.run_one_sample(
                "sars-cov-2",
                tech,
                out,
                test_data["ref_fasta"],
                reads_file=f"{reads_prefix}.fa",
                tsv_of_amp_schemes=test_data["amp_schemes_tsv"],
                force_consensus=ref_fasta,
                gzip_files=False,
            )
        else:
            one_sample_pipeline.run_one_sample(
                "sars-cov-2",
                tech,
                out,
                test_data["ref_fasta"],
                reads_file1=f"{reads_prefix}.1.fa",
                reads_file2=f"{reads_prefix}.2.fa",
                tsv_of_amp_schemes=test_data["amp_schemes_tsv"],
                force_consensus=ref_fasta,
                gzip_files=False,
            )

        log = utils.load_json(os.path.join(out, "log.json"))
        assert log["run_summary"]["result"] == "Success"
        assert log["run_summary"]["errors"] == []
        assert log["run_summary"]["last_stage_completed"] == "Finished"
        assert log["amplicon_scheme_name"] == "scheme1"
        qc = load_qc_tsv(os.path.join(out, "qc.tsv.gz"))
        cons = utils.load_single_seq_fasta(os.path.join(out, "consensus.fa"))
        assert qc[(80, 80)]["Ref_nt"] == "A"
        assert qc[(80, 80)]["Cons_nt"] == "A"
        assert cons[80] == "N"
        assert qc[(80, 80)]["Mask"] == "FRS"

        assert qc[(100, 100)]["Ref_nt"] == "T"
        assert qc[(100, 100)]["Cons_nt"] == "K"
        assert cons[100] == "K"
        assert qc[(100, 100)]["Mask"] == "HET"

        assert qc[(120, 120)]["Ref_nt"] == "G"
        assert qc[(120, 120)]["Cons_nt"] == "G"
        assert cons[120] == "M"
        assert qc[(120, 120)]["Mask"] == "FRS;HET"

        assert qc[(140, 140)]["Ref_nt"] == "G"
        assert qc[(140, 140)]["Cons_nt"] == "G"
        assert cons[140] == "S"
        assert qc[(140, 140)]["Mask"] == "HET"

        assert qc[(160, 160)]["Ref_nt"] == "G"
        assert qc[(160, 160)]["Cons_nt"] == "N"
        assert cons[160] == "N"
        assert qc[(160, 160)]["Mask"] == "ASY;FRS;HET"

        # TO DO, check more when have decided how to report indels
        assert qc[(180, 180)]["Ref_nt"] == "A"
        assert qc[(180, 180)]["Cons_nt"] == "A"

        assert qc[(200, 200)]["Ref_nt"] == "G"
        assert qc[(200, 200)]["Cons_nt"] == "G"

        assert qc[(220, 220)]["Ref_nt"] == "G"
        assert qc[(220, 220)]["Cons_nt"] == "G"

        assert qc[(239, 239)]["Ref_nt"] == "G"
        assert qc[(239, 239)]["Cons_nt"] == "G"
        assert qc[(240, 240)]["Ref_nt"] == "A"
        assert qc[(240, 240)]["Cons_nt"] == "A"

        assert qc[(260, 260)]["Ref_nt"] == "A"
        assert qc[(260, 260)]["Cons_nt"] == "A"

        assert qc[(280, 280)]["Ref_nt"] == "A"
        assert qc[(280, 280)]["Cons_nt"] == "A"

        utils.syscall(f"rm -r {out}")

    utils.syscall(f"rm {reads_prefix}* {ref_fasta}")
