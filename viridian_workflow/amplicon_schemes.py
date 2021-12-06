import json
import os

this_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(this_dir, "amplicon_scheme_data")

def get_built_in_schemes():
    json_index_file = os.path.join(DATA_DIR, "schemes.json")
    if not os.path.exists(json_index_file):
        raise Exception(f"Amplicons scheme index file not found. Something wrong with installation? Looked for: {json_index_file}")
    with open(json_index_file) as f:
        schemes = json.load(f)

    for name in schemes:
        scheme_json = os.path.join(DATA_DIR, schemes[name])
        if not os.path.exists(scheme_json):
            raise Exception(f"Amplicons scheme file not found. Something wrong with installation? Looked for: {scheme_json}")
        schemes[name] = scheme_json

    return schemes
