"""
Viridian Workflow

Wrapper around cylon
"""
from pkg_resources import get_distribution
from viridian_workflow import (
    amplicon_schemes,
    primers,
    self_qc,
    readstore,
    tasks,
    utils,
    reads,
)

__version__ = get_distribution("viridian_workflow").version

__all__ = [
    "amplicon_schemes",
    "primers",
    "self_qc",
    "readstore",
    "tasks",
    "utils",
    "reads",
]
