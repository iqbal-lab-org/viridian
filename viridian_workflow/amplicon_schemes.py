import csv
import json
import os

from viridian_workflow import primers

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


def load_list_of_amplicon_sets(built_in_names_to_use=None, tsv_others_to_use=None):
    assert built_in_names_to_use is not None or tsv_others_to_use is not None
    schemes = {}
    if built_in_names_to_use is not None:
        all_built_in_schemes = get_built_in_schemes()
        for name in built_in_names_to_use:
            if name not in all_built_in_schemes:
                names = ",".join(sorted(list(all_built_in_schemes.keys())))
                raise Exception(
                    f"No built-in amplicon scheme called {name}. Available names: {names}"
                )
            schemes[name] = all_built_in_schemes[name]

    if tsv_others_to_use is not None:
        with open(tsv_others_to_use) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for d in reader:
                if d["Name"] in schemes:
                    raise Exception(
                        f"Duplicate name '{d['Name']}' used. Cannot continue"
                    )
                if not os.path.exists(d["File"]):
                    raise FileNotFoundError(f"File not found: {d['File']}")
                schemes[d["Name"]] = d["File"]

    assert len(schemes) > 0
    return schemes, [
        primers.AmpliconSet(k, vwf_tsv_file=v) for k, v in sorted(schemes.items())
    ]
