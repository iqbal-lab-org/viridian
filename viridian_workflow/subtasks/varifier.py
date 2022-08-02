"""Varifier wrapper
"""
from __future__ import annotations

from typing import Optional
from pathlib import Path
from .task import Task
from viridian_workflow.utils import Index0, Index1


class Varifier(Task):
    def __init__(
        self,
        outdir: Path,
        ref: Path,
        consensus: Path,
        min_coord: Index0 = Index0(0),
        max_coord: Optional[Index0] = None,
        sanitise_gaps: bool = True,
        hp_min_fix_length: Optional[int] = None,
    ):
        vcf = outdir / "04.truth.vcf"
        msa = outdir / "04.msa"
        consensus_out = outdir / "04.qry_sanitised_gaps.fa"
        self.output: tuple[Path, Path, Path] = (vcf, msa, consensus_out)

        self.options: list[str] = ["--global_align"]

        self.options += ["--global_align_min_coord", str(min_coord + 1)]

        if max_coord is not None:
            self.options += ["--global_align_max_coord", str(max_coord + 1)]

        if hp_min_fix_length is not None:
            self.options += ["--hp_min_fix_length", str(hp_min_fix_length)]

        self.cmd = [
            "varifier",
            "make_truth_vcf",
        ]

        if sanitise_gaps:
            self.cmd.append("--sanitise_truth_gaps")

        self.cmd += [
            *self.options,
            str(consensus),
            str(ref),
            str(outdir),
        ]
