import csv
import json
import os

this_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(this_dir, "amplicon_scheme_data")


def get_built_in_schemes():
    json_index_file = os.path.join(DATA_DIR, "schemes.json")
    if not os.path.exists(json_index_file):
        raise Exception(
            f"Amplicons scheme index file not found. Something wrong with installation? Looked for: {json_index_file}"
        )
    with open(json_index_file) as f:
        schemes = json.load(f)

    for name in schemes:
        scheme_json = os.path.join(DATA_DIR, schemes[name])
        if not os.path.exists(scheme_json):
            raise Exception(
                f"Amplicons scheme file not found. Something wrong with installation? Looked for: {scheme_json}"
            )
        schemes[name] = scheme_json

    return schemes


def convert_tsv_to_viridian_json(tsv_in, json_out, scheme_name=None):
    amplicons = {}
    data = {
        "name": tsv_in if scheme_name is None else scheme_name,
        "source_file": tsv_in,
        "amplicons": amplicons,
    }
    with open(tsv_in) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for d in reader:
            if d["Amplicon_name"] not in amplicons:
                amplicons[d["Amplicon_name"]] = {
                    "left_primers": [],
                    "right_primers": [],
                }

            l_or_r = d["Left_or_right"].lower()
            amplicons[d["Amplicon_name"]][f"{l_or_r}_primers"].append(
                {
                    "original_data": d,
                    "start": int(d["Position"]),
                    "end": int(d["Position"]) + len(d["Sequence"]) - 1,
                }
            )

    for d in amplicons.values():
        d["start"] = min([x["start"] for x in d["left_primers"]])
        d["end"] = max([x["end"] for x in d["right_primers"]])
        d["left_primer_end"] = max([x["end"] for x in d["left_primers"]])
        d["right_primer_start"] = min([x["start"] for x in d["right_primers"]])

    with open(json_out, "w") as f:
        json.dump(data, f, indent=2)
