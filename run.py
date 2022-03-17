import sys

from viridian_workflow import primers
from viridian_workflow import detect_primers
from viridian_workflow import readstore
from viridian_workflow import minimap

print(f"{sys.argv}")


print("load amplicon sets")

amplicon_sets = {}

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
    amplicon_sets[name] = primers.AmpliconSet.from_tsv(tsv, name=name)


for name in amplicon_sets:
    print(len(amplicon_sets[name].tree))

bam = minimap.run(
    f"{sys.argv[10]}/testbam.bam", sys.argv[5], sys.argv[7], fq2=sys.argv[8], sort=False
)
stats = detect_primers.gather_stats_from_bam(bam, amplicon_sets.values())
print(stats)
# rs = readstore.ReadStore(amplicon_set, bam)
