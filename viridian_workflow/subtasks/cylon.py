"""Cylon wrapper
"""
from __future__ import annotations

from typing import Any
import json
from pathlib import Path
from .task import Task


class Cylon(Task):
    def __init__(
        self,
        work_dir: Path,
        platform: str,
        ref: Path,
        amplicon_dir: Path,
        amplicon_manifest: dict[str, Any],
        amplicon_json: dict[str, Any],
    ):
        self.output: Path = (
            work_dir / "initial_assembly" / "consensus.final_assembly.fa"
        )
        self.work_dir: Path = work_dir

        with open(
            amplicon_dir / "manifest.json", "w", encoding="utf-8"
        ) as failed_amplicon_amps_fd:
            json.dump(amplicon_manifest, failed_amplicon_amps_fd, indent=2)

        with open(
            work_dir / "amplicons.json", "w", encoding="utf-8"
        ) as failed_amplicon_amps_fd:
            json.dump(amplicon_json, failed_amplicon_amps_fd, indent=2)

        self.cmd = [
            "cylon",
            "assemble",
            "--reads_per_amp_dir",
            str(amplicon_dir),
            platform,
            str(ref),
            str(work_dir / "amplicons.json"),
            str(work_dir / "initial_assembly"),
        ]
        super(Cylon, self).__init__(name="cylon")
