import csv
import copy
import os

from viridian import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(this_dir, "amplicon_scheme_data")
SCHEMES_JSON = os.path.join(DATA_DIR, "schemes.json")
try:
    SCHEMES = utils.load_json(SCHEMES_JSON)
except:
    raise Exception(
        f"Something wrong with installation! Schemes json file missing {SCHEMES_JSON}"
    )

REF_FASTAS = {k: os.path.join(DATA_DIR, v["ref_fasta"]) for k, v in SCHEMES.items()}
DECONTAM_REFS = copy.copy(REF_FASTAS)
DECONTAM_REFS["sars-cov-2"] = os.path.join(DATA_DIR, "MN908947.no_polyA.fasta")


def get_built_in_schemes(species):
    if species not in SCHEMES:
        allowed_species = ",".join(sorted(list(SCHEMES.keys())))
        raise Exception(
            f"Unknown species: {species}. Must be one of: {allowed_species}"
        )
    schemes = SCHEMES[species]["schemes"]

    for name in schemes:
        scheme_json = os.path.join(DATA_DIR, schemes[name])
        if not os.path.exists(scheme_json):
            raise Exception(
                f"Amplicons scheme file not found. Something wrong with installation? Looked for: {scheme_json}"
            )
        schemes[name] = scheme_json

    return schemes


def load_list_of_amplicon_sets(
    species, built_in_names_to_use=None, tsv_others_to_use=None
):
    assert built_in_names_to_use is not None or tsv_others_to_use is not None
    if isinstance(built_in_names_to_use, str):
        built_in_names_to_use = built_in_names_to_use.split(",")
    schemes = {}
    if built_in_names_to_use is not None:
        all_built_in_schemes = get_built_in_schemes(species)
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
