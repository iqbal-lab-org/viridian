"""Viridian wrapper
"""
import json
from .task import Task


class Viridian(Task):
    def __init__(
        self, work_dir, platform, ref, amplicon_dir, amplicon_manifest, amplicon_json
    ):
        self.output = work_dir / "viridian" / "consensus.final_assembly.fa"
        self.work_dir = work_dir

        with open(amplicon_dir / "manifest.json", "w") as failed_amplicon_amps_fd:
            json.dump(amplicon_manifest, failed_amplicon_amps_fd, indent=2)

        with open(work_dir / "amplicons.json", "w") as failed_amplicon_amps_fd:
            json.dump(amplicon_json, failed_amplicon_amps_fd, indent=2)

        self.cmd = [
            "viridian",
            "assemble",
            "--reads_per_amp_dir",
            amplicon_dir,
            platform,
            ref,
            work_dir / "amplicons.json",
            work_dir / "viridian",
        ]
