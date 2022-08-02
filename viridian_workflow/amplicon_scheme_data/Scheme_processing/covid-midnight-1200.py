#!/usr/bin/env python3
# type: ignore

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

outprefix = "covid-midnight-1200"
tsv_url = "https://zenodo.org/record/3897530/files/SARS-CoV-2_primer_sets_RBK004_nanopore_sequencing.tab"
tsv_file = f"{outprefix}.tsv"
json_file = f"{outprefix}.json"
subprocess.check_output(f"wget -O {tsv_file} {tsv_url}", shell=True)
with open(tsv_file, "rb") as f:
    tsv_file_sha256 = hashlib.sha256(f.read()).hexdigest()


with open(tsv_file) as f_in, open("../covid-midnight-1200.vwf.tsv", "w") as f_out:
    reader = csv.DictReader(f_in, delimiter="\t")
    print(
        "Amplicon_name",
        "Primer_name",
        "Left_or_right",
        "Sequence",
        "Position",
        sep="\t",
        file=f_out,
    )

    for d in reader:
        matches = ref.search(d["Sequence"])
        assert len(matches) == 1
        match = matches[0]

        amplicon_name, number, l_or_r = d["Primer Name"].rsplit("_", maxsplit=2)
        assert l_or_r in ["LEFT", "RIGHT"]
        if l_or_r == "LEFT":
            assert match[0] == int(d["Start"])
        else:
            assert match[0] == int(d["End"])

        amplicon_name += "_" + number
        if amplicon_name not in amplicons:
            amplicons[amplicon_name] = {"left_primers": [], "right_primers": []}
        d = {
            "original_data": d,
            "start": match[0],
            "end": match[0] + len(d["Sequence"]) - 1,
        }
        amplicons[amplicon_name][l_or_r.lower() + "_primers"].append(d)
        print(
            amplicon_name,
            d["original_data"]["Primer Name"],
            l_or_r.lower(),
            d["original_data"]["Sequence"],
            d["start"],
            sep="\t",
            file=f_out,
        )

for amplicon_name, d in amplicons.items():
    d["start"] = min([x["start"] for x in d["left_primers"]])
    d["end"] = max([x["end"] for x in d["right_primers"]])
    d["left_primer_end"] = max([x["end"] for x in d["left_primers"]])
    d["right_primer_start"] = min([x["start"] for x in d["right_primers"]])

json_data = {
    "name": "covid-midnight-1200",
    "source_file": tsv_url,
    "source_file_sha256": tsv_file_sha256,
    "reference_accession": "MN908947.3",
    "amplicons": amplicons,
}

with open(json_file, "w") as f:
    json.dump(json_data, f, indent=2)
