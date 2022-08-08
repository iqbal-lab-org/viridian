"""
Collect the built-in amplicon scheme files
"""
from __future__ import annotations

import csv
from typing import Optional
from pathlib import Path

from viridian_workflow import primers

this_dir = Path(__file__).resolve()
DATA_DIR = this_dir / "amplicon_scheme_data"


def get_built_in_schemes() -> dict[str, Path]:
    """Read list of built-in schemes from schemes.tsv"""
    tsv_index_file = DATA_DIR / "schemes.tsv"
    if not tsv_index_file.exists():
        raise Exception(
            f"Amplicons scheme index file not found. Looked for: {tsv_index_file}"
        )

    schemes: dict[str, Path] = {}
    with open(tsv_index_file, encoding="utf-8") as tsv_fd:
        for line in tsv_fd:
            name, tsv = line.strip().split()
            schemes[name] = Path(tsv)

    for name, scheme in schemes.items():
        scheme_path = DATA_DIR / scheme
        if not scheme_path.exists():
            raise Exception(
                f"Amplicons scheme file not found. Looked for: {scheme_path}"
            )
        schemes[name] = scheme_path

    return schemes


def load_amplicon_index(
    index_tsv: Path, scheme_dir: Path, subset: dict[str, Path] = None
) -> dict[str, Path]:
    """Load either all or a subset of the built-in amplicon schemes"""
    index = {}
    with open(scheme_dir / index_tsv, encoding="utf-8") as index_fd:
        for record in index_fd:
            name, tsv = record.strip().split()
            if name == "Name" and tsv == "File":
                continue

            tsv_path = scheme_dir / tsv
            if not tsv_path.exists():
                raise Exception(f"Amplicon scheme {tsv_path} does not exist")
            index[name] = tsv_path

    if subset:
        for key in subset:
            if key not in index:
                raise Exception(
                    f"Selected subset of amplicon schemes ({','.join(subset)})\
                      are not in the builtin set: {','.join(index.keys())}"
                )
            del index[key]
    return index


def load_list_of_amplicon_sets(
    built_in_names_to_use: Optional[list[str]] = None,
    tsv_others_to_use: Optional[Path] = None,
) -> tuple[dict[str, Path], list[primers.AmpliconSet]]:
    """Load list of amplicon schemes from preferred source"""
    assert built_in_names_to_use is not None or tsv_others_to_use is not None
    if isinstance(built_in_names_to_use, str):
        built_in_names_to_use = built_in_names_to_use.split(",")
    schemes: dict[str, Path] = {}
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
        with open(tsv_others_to_use, encoding="utf-8") as tsv_fd:
            reader = csv.DictReader(tsv_fd, delimiter="\t")
            for row in reader:
                if row["Name"] in schemes:
                    raise Exception(
                        f"Duplicate name '{row['Name']}' used. Cannot continue"
                    )
                if not Path(row["File"]).exists():
                    raise FileNotFoundError(f"File not found: {row['File']}")
                schemes[row["Name"]] = Path(row["File"])

    assert len(schemes) > 0
    return (
        schemes,
        [primers.AmpliconSet.from_tsv(v, name=k) for k, v in sorted(schemes.items())],
    )
