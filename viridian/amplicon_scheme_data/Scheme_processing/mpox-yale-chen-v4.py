#!/usr/bin/env python3

import csv
import json
import pyfastaq
import subprocess


ref = {}
pyfastaq.tasks.file_to_dict("../MT903345.fasta.gz", ref)
assert len(ref) == 1
ref = list(ref.values())[0]
amplicons = {}

outprefix = "mpox-yale-chen-v4"
bed_url = "https://content.protocols.io/files/h7pwbb8a7.bed"
bed_file = f"{outprefix}.bed"
subprocess.check_output(f"wget -O {bed_file} {bed_url}", shell=True)


with open(bed_file) as f_in, open("../mpox-yale-chen-v4.vwf.tsv", "w") as f_out:
    print("Amplicon_name",
        "Primer_name",
        "Left_or_right",
        "Sequence",
        "Position",
        sep="\t",
        file=f_out
    )

    for line in f_in:
        ref_name, start, end, primer_name, pool, strand  = line.rstrip().split("\t")
        start = int(start)
        end = int(end)
        assert start < end

        # primer names are like: MPXV_1_LEFT
        amplicon_name, l_or_r = primer_name.rsplit("_", maxsplit=1)
        assert l_or_r in ["LEFT", "RIGHT"]
        if primer_name.endswith("LEFT"):
            assert strand == "+"
        else:
            assert primer_name.endswith("RIGHT")
            assert strand == "-"
        seq = ref[start:end]

        print(
            amplicon_name,
            primer_name,
            l_or_r.lower(),
            seq,
            start,
            sep="\t",
            file=f_out
        )

