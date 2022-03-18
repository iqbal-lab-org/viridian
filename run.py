import sys

from viridian_workflow import primers, readstore, minimap, utils

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

# generate name-sorted bam from fastqs
minimap.run(
    work_dir / "testbam.bam", sys.argv[1], sys.argv[1], fq2=sys.argv[2], sort=False
)

# pre-process input bam
bam = readstore.Bam(work_dir / "unsorted.bam")

# detect amplicon set
amplicon_set = bam.detect_amplicon_set(amplicon_sets)

stats = bam.stats

# construct readstore
rs = readstore.ReadStore(amplicon_set, bam)

# downsample to viridian assembly
failures = rs.make_reads_dir_for_viridian(work_dir / "amplicons", 1000)

for failure in failures:
    print(failure.name, file=failed_amplicon_amps_fd)

# run viridian
def run_viridian(self):
    logging.info("Making initial unmasked consensus using Viridian")
    viridian_cmd = " ".join(
        [
            "viridian",
            "assemble",
            "--bam",
            self.sampled_bam,
            "illumina",
            "--amplicons_to_fail_file",
            failed_amps_file,
            ref,
            amplicon_json,
            work_dir / "viridian",
        ]
    )
    utils.run_process(viridian_cmd)
