import sys
import shutil
from pathlib import Path
import tempfile
import json

from viridian_workflow import primers, readstore, minimap, utils, self_qc, varifier


def run_viridian(work_dir, amplicon_dir, amplicon_manifest, amplicon_json):
    with open(amplicon_dir / "manifest.json", "w") as failed_amplicon_amps_fd:
        json.dump(amplicon_manifest, failed_amplicon_amps_fd, indent=2)

    with open(work_dir / "amplicons.json", "w") as failed_amplicon_amps_fd:
        json.dump(amplicon_json, failed_amplicon_amps_fd, indent=2)

    viridian_cmd = [
        "viridian",
        "assemble",
        "--reads_per_amp_dir",
        amplicon_dir,
        "illumina",
        ref,
        work_dir / "amplicons.json",
        work_dir / "viridian",
    ]
    utils.run_process(viridian_cmd)
    output = work_dir / "viridian" / "consensus.final_assembly.fa"
    utils.check_file(output)
    return output


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
# with tempfile.TemporaryDirectory() as work_dir:
work_dir = f"/tmp/vwf_test/"
if True:
    work_dir = Path(work_dir)
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir()
    print(f"using {work_dir} {work_dir.exists()}")

    # generate name-sorted bam from fastqs
    unsorted_bam = minimap.run(work_dir / "testbam.bam", ref, fq1, fq2=fq2, sort=False)

    # pre-process input bam
    bam = readstore.Bam(unsorted_bam)

    # detect amplicon set
    amplicon_set = bam.detect_amplicon_set(amplicon_sets)

    # check it out: stats = bam.stats

    # construct readstore
    # note: we could/should be subsampling at this point. we
    # could do this by passing in bam.stats
    rs = readstore.ReadStore(amplicon_set, bam)

    # downsample to viridian assembly
    amp_dir = work_dir / "amplicons"
    manifest_data = rs.make_reads_dir_for_viridian(amp_dir, 1000)

    # run viridian
    consensus = run_viridian(work_dir, amp_dir, manifest_data, rs.viridian_json)

    # varifier
    vcf = varifier.run(work_dir / "varifier", ref, consensus, min_coord=rs.start_pos, max_coord=rs.end_pos)

    # self qc
    position_stats = rs.remap(consensus)

    # annotate vcf
    annotated_vcf = self_qc.annotate_vcf(
        work_dir / "varifier" / "04.final.vcf",
        work_dir / "varifier" / "04.msa",
        position_stats
    )
    for rec in annotated_vcf:
        print(rec)
