#!/usr/bin/env python3

import csv
import hashlib
import json
import operator
import subprocess

import pyfastaq


def get_bed_and_md5(url, outfile):
    #subprocess.check_output(f"wget -O {outfile} {url}", shell=True)
    with open(outfile, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def file_to_set(infile):
    with open(infile) as f:
        return set(f)


ref = {}
pyfastaq.tasks.file_to_dict("../MN908947.fasta", ref)
assert len(ref) == 1
ref = list(ref.values())[0]
amplicons = {}

outprefix = "covid-varskip-v1a-2b"
bed_url_1a = "https://raw.githubusercontent.com/nebiolabs/VarSkip/main/schemes/NEB_VarSkip/V1a/NEB_VarSkip.scheme.bed"
bed_url_2b = "https://raw.githubusercontent.com/nebiolabs/VarSkip/main/schemes/NEB_VarSkip/V2b/NEB_VarSkip.scheme.bed"
bed_file_1a = f"{outprefix}.1a.bed"
bed_file_2b = f"{outprefix}.2b.bed"
json_file = f"{outprefix}.json"
bed_file_sha256_1a = get_bed_and_md5(bed_url_1a, bed_file_1a)
bed_file_sha256_2b = get_bed_and_md5(bed_url_2b, bed_file_2b)

bed_lines_1a = file_to_set(bed_file_1a)
bed_lines_2b = file_to_set(bed_file_2b)

bed_1a_only = []
for line in bed_lines_1a:
    if line in bed_lines_2b:
        continue

    assert "ALT" not in line
    fields = line.rstrip().split("\t")
    assert fields[3].endswith("_LEFT") or fields[3].endswith("_RIGHT")
    fields[3] += "_ALT-v1a"
    bed_1a_only.append(fields)


wanted_lines = bed_1a_only + list([x.rstrip().split("\t") for x in bed_lines_2b])
wanted_lines.sort(key=lambda x: int(x[1]))

with open("../covid-varskip-v1a-2b.vwf.tsv", "w") as f_out:
    print("Amplicon_name",
        "Primer_name",
        "Left_or_right",
        "Sequence",
        "Position",
        sep="\t",
        file=f_out
    )

    for fields in wanted_lines:
        ref_name, start, end, primer_name, pool, strand = fields
        start = int(start)
        end = int(end)

        if "_ALT" in primer_name:
            amplicon_name, number, l_or_r, alt = primer_name.split("_")
        else:
            amplicon_name, number, l_or_r = primer_name.split("_")
            alt = None

        assert l_or_r.split("_")[0] in ["LEFT", "RIGHT"]
        amplicon_name += "_" + number
        if amplicon_name not in amplicons:
            amplicons[amplicon_name] = {"left_primers": [], "right_primers": []}
        d = {
            "original_data": fields,
            "start": start,
            "end": end - 1,
        }
        amplicons[amplicon_name][l_or_r.lower() + "_primers"].append(d)
        print(
            amplicon_name,
            primer_name,
            l_or_r.lower(),
            ref[start:end],
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
    "name": "covid-varskip-v1a",
    "source_file_1a": bed_url_1a,
    "source_file_sha256_1a": bed_file_sha256_1a,
    "source_file_2b": bed_url_2b,
    "source_file_sha256_2b": bed_file_sha256_2b,
    "reference_accession": "MN908947.3",
    "amplicons": amplicons,
}

with open(json_file, "w") as f:
    json.dump(json_data, f, indent=2)

