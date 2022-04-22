import sys
import shutil
from pathlib import Path
import tempfile
import json

from viridian_workflow import primers, readstore, minimap, utils, self_qc, varifier
from viridian_workflow.tasks import minimap, varifier, viridian

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
    ("artic-v4.0", "viridian_workflow/amplicon_scheme_data/covid-artic-v4.vwf.tsv"),
]:
    amplicon_sets.append(primers.AmpliconSet.from_tsv(tsv, name=name))

# input file args, working directory, and reference file path
fq1, fq2 = sys.argv[1], sys.argv[2]
ref = Path("../covid/MN908947.fasta")
# with tempfile.TemporaryDirectory() as work_dir:
work_dir = sys.argv[3]

log = {"summary": {"version": "test-0.1", "status": "processing"}}

try:
    work_dir = Path(work_dir)
    if work_dir.exists():
        print(f"workdir {work_dir} exists, exiting")
        exit()
    #     shutil.rmtree(work_dir)
    work_dir.mkdir()
    print(f"using {work_dir} {work_dir.exists()}")

    # generate name-sorted bam from fastqs
    unsorted_bam = minimap.run(
        work_dir / "name_sorted.bam", ref, fq1, fq2=fq2, sort=False
    )

    # pre-process input bam
    bam = readstore.Bam(unsorted_bam)

    # detect amplicon set
    amplicon_set = bam.detect_amplicon_set(amplicon_sets)
    log["amplicons"] = bam.stats

    # check it out: stats = bam.stats

    # construct readstore
    # note: we could/should be subsampling at this point. we
    # could do this by passing in bam.stats
    rs = readstore.ReadStore(amplicon_set, bam)

    # save reads for viridian assembly
    amp_dir = work_dir / "amplicons"
    manifest_data = rs.make_reads_dir_for_viridian(amp_dir)

    # run viridian
    consensus = run_viridian(work_dir, amp_dir, manifest_data, rs.viridian_json)

    # varifier
    vcf, msa = Varifier(
        work_dir / "varifier",
        ref,
        consensus,
        min_coord=rs.start_pos,
        max_coord=rs.end_pos,
    ).run()

    log["varifier"] = varifier.log
    # self qc: remap reads to consensus
    pileup = rs.pileup(consensus, msa=msa)

    # mask output
    masked_fasta = pileup.mask()
    # log["self_qc"] = pileup.log
    log["summary"] = pileup.summary

    for i in pileup.summary:
        print(f"{i}\t{pileup.summary[i]}")

    with open(work_dir / "masked.fasta", "w") as fasta_out:
        print(masked_fasta, file=fasta_out)

    # annotate vcf
    print(f"msa: {msa}")
    annotated_vcf = pileup.annotate_vcf(vcf)

    with open(work_dir / "final.vcf", "w") as vcf_out:
        header, records = annotated_vcf
        for h in header:
            print(h, file=vcf_out)
        for rec in records:
            print("\t".join(map(str, rec)), file=vcf_out)

    print(log)
    with open(work_dir / "log.json", "w") as json_out:
        json.dump(log, json_out, indent=4)
except(e):
    log["summary"]["status"] = {"Failed", str(e)}

log["summary"]["status"] = "Success"
