"""Viridian wrapper
"""
import os
from viridian_workflow.utils import run_process, check_file
from .task import Task


class Viridian(Task):
    def __init__(
        self, work_dir, platform, amplicon_dir, amplicon_manifest, amplicon_json
    ):
        output = work_dir / "viridian" / "consensus.final_assembly.fa"
        self.platform = platform
        self.work_dir = work_dir
        self.amplicon_dir = amplicon_dir

    def run(self):
        with open(self.amplicon_dir / "manifest.json", "w") as failed_amplicon_amps_fd:
            json.dump(self.amplicon_manifest, failed_amplicon_amps_fd, indent=2)

        with open(self.work_dir / "amplicons.json", "w") as failed_amplicon_amps_fd:
            json.dump(self.amplicon_json, failed_amplicon_amps_fd, indent=2)

        viridian_cmd = [
            "viridian",
            "assemble",
            "--reads_per_amp_dir",
            self.amplicon_dir,
            self.platform,
            self.ref,
            self.work_dir / "amplicons.json",
            self.work_dir / "viridian",
        ]
        utils.run_process(viridian_cmd)
        output = self.work_dir / "viridian" / "consensus.final_assembly.fa"
        utils.check_file(self.output)
        return self.output
