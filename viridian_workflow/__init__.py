from pkg_resources import get_distribution

try:
    __version__ = get_distribution("viridian_workflow").version
except:
    __version__ = "local"

__all__ = ["one_sample_pipeline", "tasks", "utils", "minimap", "qcovid", "config"]

from viridian_workflow import *
