#!/usr/bin/env python3

import csv
import hashlib
import json
import pyfastaq
import subprocess


ref = {}
pyfastaq.tasks.file_to_dict("../MN908947.fasta", ref)
assert len(ref) == 1
ref = list(ref.values())[0]
amplicons = {}

outprefix = "covid-artic-v4.1"
bed_url = "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.primer.bed"
bed_file = f"{outprefix}.bed"
json_file = f"{outprefix}.json"
subprocess.check_output(f"wget -O {bed_file} {bed_url}", shell=True)
with open(bed_file, "rb") as f:
    bed_file_sha256 = hashlib.sha256(f.read()).hexdigest()


with open(bed_file) as f_in, open("../covid-artic-v4.1.vwf.tsv", "w") as f_out:
    print("Amplicon_name",
        "Primer_name",
        "Left_or_right",
        "Sequence",
        "Position",
        sep="\t",
        file=f_out
    )

    for line in f_in:
        fields = line.rstrip().split("\t")
        ref_name, start, end, primer_name, pool, strand, seq = fields
        start = int(start)
        end = int(end)

        # The coords of each primer are in the BED file, but we'll check they
        # match the reference sequence at the expected position anyway
        matches = ref.search(seq)
        assert len(matches) == 1
        match = matches[0]
        assert match[0] == start
        assert match[1] == strand

        print(line)
        if "_alt" in line:
            amplicon_name, number, l_or_r, alt = primer_name.split("_")
        else:
            amplicon_name, number, l_or_r = primer_name.split("_")
            alt = None

        assert l_or_r in ["LEFT", "RIGHT"]
        amplicon_name += "_" + number
        if amplicon_name not in amplicons:
            amplicons[amplicon_name] = {"left_primers": [], "right_primers": []}
        d = {
            "original_data": fields,
            "start": match[0],
            "end": match[0] + len(seq) - 1,
        }
        assert d["end"] + 1 == end
        amplicons[amplicon_name][l_or_r.lower() + "_primers"].append(d)
        print(
            amplicon_name,
            primer_name,
            l_or_r.lower(),
            seq,
            start,
            sep="\t",
            file=f_out
        )

for amplicon_name, d in amplicons.items():
    d["start"] = min([x["start"] for x in d["left_primers"]])
    d["end"] = max([x["end"] for x in d["right_primers"]])
    d["left_primer_end"] = max([x["end"] for x in d["left_primers"]])
    d["right_primer_start"] = min([x["start"] for x in d["right_primers"]])

json_data = {
    "name": "covid-artic-v4.1",
    "source_file": bed_url,
    "source_file_sha256": bed_file_sha256,
    "reference_accession": "MN908947.3",
    "amplicons": amplicons,
}

with open(json_file, "w") as f:
    json.dump(json_data, f, indent=2)
