#!/usr/bin/env python3

import csv
import hashlib
import json
import operator
import subprocess

import pyfastaq


def get_bed_and_md5(url, outfile):
    subprocess.check_output(f"wget -O {outfile} {url}", shell=True)
    with open(outfile, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def file_to_set(infile):
    with open(infile) as f:
        return set([x.rstrip() for x in f if len(x.strip()) > 0])


ref = {}
pyfastaq.tasks.file_to_dict("../MN908947.fasta", ref)
assert len(ref) == 1
ref = list(ref.values())[0]
amplicons = {}

outprefix = "covid-artic-v5.0-5.3.2_400"

scheme2bed_url = {
    "v5.0.0": "https://raw.githubusercontent.com/quick-lab/SARS-CoV-2/main/400/v5.0.0_400/SARS-CoV-2_v5.0.0_400.primer.bed",
    "v5.1.0": "https://raw.githubusercontent.com/quick-lab/SARS-CoV-2/main/400/v5.1.0_400/SARS-CoV-2_v5.1.0_400.primer.bed",
    "v5.2.0": "https://raw.githubusercontent.com/quick-lab/SARS-CoV-2/main/400/v5.2.0_400/SARS-CoV-2_v5.2.0_400.primer.bed",
    "v5.3.2": "https://raw.githubusercontent.com/quick-lab/SARS-CoV-2/main/400/v5.3.2_400/SARs-CoV-2_v5.3.2_400.primer.bed",
}

scheme2bed = {k: f"covid-artic-{k}_400.bed" for k in scheme2bed_url}
scheme2bed_sha = {k: get_bed_and_md5(scheme2bed_url[k], scheme2bed[k]) for k in scheme2bed}
scheme2bed_lines = {k: file_to_set(v) for k, v in scheme2bed.items()}


lines2scheme = {}
for scheme, lines in scheme2bed_lines.items():
    for line in lines:
        if line not in lines2scheme:
            lines2scheme[line] = []
        lines2scheme[line].append(scheme)


json_file = f"{outprefix}.json"
lines_out = []


for line, schemes in lines2scheme.items():
    fields = line.rstrip().split("\t")
    ref_name, start, end, primer_name, pool, strand, seq = fields
    start = int(start)
    end = int(end)

    # The coords of each primer are in the BED file, but we'll check they
    # match the reference sequence at the expected position anyway
    if primer_name == "SARS-CoV-2_400_84_RIGHT_2":
        assert seq == "TGTTCAACACCARTGTCTGTACTC"
        matches = ref.search(seq.replace("R", "G"))
    else:
        matches = ref.search(seq)
    assert len(matches) == 1
    match = matches[0]
    assert match[0] == start
    assert match[1] == strand


    assert "alt" not in primer_name.upper()
    assert primer_name.startswith("SARS-CoV-2_400_")
    amplicon_name, l_or_r, primer_number = primer_name.rsplit("_", maxsplit=2)

    assert l_or_r.split("_")[0] in ["LEFT", "RIGHT"]
    if amplicon_name not in amplicons:
        amplicons[amplicon_name] = {"left_primers": [], "right_primers": []}
    d = {
        "original_data": fields,
        "in_schemes": schemes,
        "start": start,
        "end": end - 1,
    }
    amplicons[amplicon_name][l_or_r.lower() + "_primers"].append(d)
    lines_out.append((
        amplicon_name,
        primer_name,
        l_or_r.lower(),
        seq,
        start,
    ))

lines_out.sort(key=lambda x: int(x[4]))

with open(f"../{outprefix}.vwf.tsv", "w") as f_out:
    print("Amplicon_name",
        "Primer_name",
        "Left_or_right",
        "Sequence",
        "Position",
        sep="\t",
        file=f_out
    )
    for line in lines_out:
        print(*line, sep="\t", file=f_out)


for amplicon_name, d in amplicons.items():
    d["start"] = min([x["start"] for x in d["left_primers"]])
    d["end"] = max([x["end"] for x in d["right_primers"]])
    d["left_primer_end"] = max([x["end"] for x in d["left_primers"]])
    d["right_primer_start"] = min([x["start"] for x in d["right_primers"]])

json_data = {
    "name": outprefix,
    "source_files": scheme2bed_url,
    "source_files_sha256": scheme2bed_sha,
    "reference_accession": "MN908947.3",
    "amplicons": amplicons,
}

with open(json_file, "w") as f:
    json.dump(json_data, f, indent=2)

