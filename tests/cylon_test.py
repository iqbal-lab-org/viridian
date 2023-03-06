import os
import pytest
import random

import pyfastaq

from viridian import cylon, utils


def perfect_reads_one_amp(ref, start, end, depth, outfile):
    read = ref[start : end + 1]
    read_rev = pyfastaq.sequences.Fasta("rev", ref[start : end + 1])
    read_rev.revcomp()
    with open(outfile, "w") as f:
        for i in range(int(depth / 2)):
            print(f">{i}", read, sep="\n", file=f)
            print(f">{i}.rev", read_rev.seq, sep="\n", file=f)


def test_run_cylon():
    root_dir = "tmp.cylon_test"
    utils.syscall(f"rm -rf {root_dir}")
    os.mkdir(root_dir)
    reads_dir = os.path.join(root_dir, "reads")
    os.mkdir(reads_dir)
    ref_bases = "".join(random.choices(["A", "C", "G", "T"], k=1000))
    ref_seq = pyfastaq.sequences.Fasta("ref", ref_bases)
    ref_fasta = os.path.join(root_dir, "ref.fa")
    with open(ref_fasta, "w") as f:
        print(ref_seq, file=f)

    amp1 = {"start": 50, "end": 400, "left_primer_end": 60, "right_primer_start": 390}
    amp2 = {"start": 300, "end": 700, "left_primer_end": 310, "right_primer_start": 690}
    amplicons = {
        "amplicons": {
            "amp1": amp1,
            "amp2": amp2,
        }
    }
    amp_json = os.path.join(root_dir, "amplicons.json")
    utils.write_json(amp_json, amplicons)

    manifest = {"amp1": "1.fa", "amp2": "2.fa"}
    perfect_reads_one_amp(ref_seq, 50, 400, 25, os.path.join(reads_dir, "1.fa"))
    perfect_reads_one_amp(ref_seq, 300, 700, 25, os.path.join(reads_dir, "2.fa"))
    utils.write_json(os.path.join(reads_dir, "manifest.json"), manifest)

    cylon_outdir = os.path.join(root_dir, "cylon")
    got_dict, got_error = cylon.run_cylon(
        reads_dir,
        "ont",
        ref_fasta,
        amp_json,
        cylon_outdir,
        debug=False,
    )

    # we're not fussed about the exact sequence here, just that cylon can run
    # ok and return something sane. The reads are perfect, so we should get a
    # subsequence of the original reference sequence
    assert got_error is None
    assert got_dict["run_summary"]["made_consensus"] is True
    cylon_cons_fa = os.path.join(cylon_outdir, "consensus.final_assembly.fa")
    got_seq = utils.load_single_seq_fasta(cylon_cons_fa)
    assert len(got_seq) > 500
    assert got_seq.seq in ref_seq.seq

    # replace the reads files with junk and check we get error message back
    with open(os.path.join(reads_dir, "1.fa"), "w") as f:
        print(">read", "A" * 300, sep="\n", file=f)
    with open(os.path.join(reads_dir, "2.fa"), "w") as f:
        print(">read", "A" * 300, sep="\n", file=f)

    utils.syscall(f"rm -r {cylon_outdir}")
    got_dict, got_error = cylon.run_cylon(
        reads_dir,
        "ont",
        ref_fasta,
        amp_json,
        cylon_outdir,
        debug=False,
    )
    assert got_error == "No consensus sequence from assembler"
    assert got_dict["run_summary"]["made_consensus"] is False

    utils.syscall(f"rm -r {cylon_outdir}")
    got_dict, got_error = cylon.run_cylon(
        reads_dir,
        "ont",
        "not_a_file",
        amp_json,
        cylon_outdir,
        debug=False,
    )
    assert got_error == "Error running cylon"
    assert got_dict == {}
    utils.syscall(f"rm -r {root_dir}")
