import csv
import json
import os

from viridian_workflow import primers

this_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(this_dir, "amplicon_scheme_data")


def get_built_in_schemes():
    tsv_index_file = os.path.join(DATA_DIR, "schemes.tsv")
    if not os.path.exists(tsv_index_file):
        raise Exception(
            f"Amplicons scheme index file not found. Something wrong with installation? Looked for: {json_index_file}"
        )

    schemes = {}
    with open(tsv_index_file) as f:
        for line in f:
            name, tsv = line.strip().split()
            schemes[name] = tsv

    for name in schemes:
        scheme_json = os.path.join(DATA_DIR, schemes[name])
        if not os.path.exists(scheme_json):
            raise Exception(
                f"Amplicons scheme file not found. Something wrong with installation? Looked for: {scheme_json}"
            )
        schemes[name] = scheme_json

    return schemes


def load_list_of_amplicon_sets(built_in_names_to_use=None, tsv_others_to_use=None):
    assert built_in_names_to_use is not None or tsv_others_to_use is not None
    if isinstance(built_in_names_to_use, str):
        built_in_names_to_use = built_in_names_to_use.split(",")
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
    return (
        schemes,
        [primers.AmpliconSet.from_tsv(v, name=k) for k, v in sorted(schemes.items())],
    )
