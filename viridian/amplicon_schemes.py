import csv
import os

from viridian import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(this_dir, "amplicon_scheme_data")
REF_FASTA = os.path.join(DATA_DIR, "MN908947.fasta")
REF_FASTA_NO_POLY_A = os.path.join(DATA_DIR, "MN908947.no_polyA.fasta")


def get_built_in_schemes():
    json_index_file = os.path.join(DATA_DIR, "schemes.json")
    if not os.path.exists(json_index_file):
        raise Exception(
            f"Amplicons scheme index file not found. Something wrong with installation? Looked for: {json_index_file}"
        )
    schemes = utils.load_json(json_index_file)

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
    return schemes
