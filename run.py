import sys
from pathlib import Path
import tempfile
import json

from viridian_workflow import primers, readstore, minimap, utils


def run_viridian(amplicon_dir, amplicon_json):
    print(amplicon_json)
    with open(work_dir / "failed_amplicons.json", "w") as failed_amplicon_amps_fd:
        json.dump(amplicon_failures, failed_amplicon_amps_fd)

    viridian_cmd = [
        "viridian",
        "assemble",
        "--reads_per_amp_dir",
        amplicon_dir,
        "illumina",
        ref,
        work_dir / "failed_amplicons.json",
        work_dir / "viridian",
    ]
    utils.run_process(viridian_cmd)


# load amplicon sets
amplicon_sets = []
for name, tsv in [
    ("artic-v3", "viridian_workflow/amplicon_scheme_data/covid-artic-v3.vwf.tsv"),
    (
        "midnight-1200",
        "viridian_workflow/amplicon_scheme_data/covid-midnight-1200.vwf.tsv",
    ),
    (
        "covid-ampliseq-v1",
        "viridian_workflow/amplicon_scheme_data/covid-ampliseq-v1.vwf.tsv",
    ),
    ("artic-v4", "viridian_workflow/amplicon_scheme_data/covid-artic-v4.vwf.tsv"),
]:
    amplicon_sets.append(primers.AmpliconSet.from_tsv(tsv, name=name))

# input file args, working directory, and reference file path
fq1, fq2 = sys.argv[1], sys.argv[2]
ref = Path("../covid/MN908947.fasta")
with tempfile.TemporaryDirectory() as work_dir:
    work_dir = Path(work_dir)
    print(f"using {work_dir} {work_dir.exists()}")

    # generate name-sorted bam from fastqs
    unsorted_bam = minimap.run(work_dir / "testbam.bam", ref, fq1, fq2=fq2, sort=False)

    # pre-process input bam
    bam = readstore.Bam(unsorted_bam)

    # detect amplicon set
    amplicon_set = bam.detect_amplicon_set(amplicon_sets)

    # check it out: stats = bam.stats

    # construct readstore
    rs = readstore.ReadStore(amplicon_set, bam)

    # downsample to viridian assembly
    amp_dir = work_dir / "amplicons"
    amplicon_failures = rs.make_reads_dir_for_viridian(amp_dir, 1000)

    # run viridian
    run_viridian(amp_dir, amplicon_failures)
