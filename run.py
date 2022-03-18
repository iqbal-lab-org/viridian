import sys

from viridian_workflow import primers
from viridian_workflow import readstore
from viridian_workflow import minimap

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
    f"{sys.argv[10]}/testbam.bam", sys.argv[5], sys.argv[7], fq2=sys.argv[8], sort=False
)

# pre-process input bam
bam = readstore.Bam(f"{sys.argv[10]}/testbam.bam")

# detect amplicon set
amplicon_set = bam.detect_amplicon_set(amplicon_sets)

stats = bam.stats

# construct readstore
rs = readstore.ReadStore(amplicon_set, bam)

# downsample to viridian assembly
failures = rs.make_reads_dir_for_viridian("/tmp/asdfg", 1000)

for failure in failures:
    print(failure.name)
