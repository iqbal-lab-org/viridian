from pkg_resources import get_distribution

try:
    __version__ = get_distribution("viridian_workflow").version
except:
    __version__ = "local"

__all__ = [
    "amplicon_schemes",
    "primers",
    "self_qc",
    "readstore",
    "tasks",
    "utils",
    "reads",
]

from viridian_workflow import *
