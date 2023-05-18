from pkg_resources import get_distribution

try:
    __version__ = get_distribution("viridian").version
except:
    __version__ = "local"

__all__ = [
    "amplicon_schemes",
    "constants",
    "cylon",
    "ena",
    "maptools",
    "one_sample_pipeline",
    "qc",
    "read_it_and_keep",
    "reads",
    "scheme_id",
    "tasks",
    "utils",
    "varifier",
]

from viridian import *
