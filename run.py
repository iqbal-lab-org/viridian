import sys
import shutil
from pathlib import Path
import tempfile
import json

from viridian_workflow import primers, readstore, minimap, utils, self_qc, varifier
from viridian_workflow.tasks import minimap, varifier, viridian


if __name__ == "__main__":
    # set up the pipeline

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

    # with tempfile.TemporaryDirectory() as work_dir:
    work_dir = "/tmp/vwf/"
    work_dir = Path(work_dir)
    if work_dir.exists():
        print(f"workdir {work_dir} exists, clobbering")
        shutil.rmtree(work_dir)
    work_dir.mkdir()
    print(f"using {work_dir}")

    log = {"summary": {"version": "test-0.1", "status": "processing"}}
    try:
        results = run_pipeline(work_dir, platform, fqs)
        log["results"] = results
    except (e):
        log["summary"]["status"] = {"Failed", str(e)}
        exit(0)  # ?

    log["summary"]["status"] = "Success"

    with open(work_dir / "log.json", "w") as json_out:
        json.dump(log, json_out, indent=4)


def run_pipeline(work_dir, platform, fqs):
    log = {}
    # generate name-sorted bam from fastqs
    ref = Path("../covid/MN908947.fasta")

    if platform == "illumina":
        fq1, fq2 = *fqs
        minimap = Minimap(work_dir, ref, fq1, fq2=fq2, sort=False)
    elif platform == "onp":
        fq = fqs[0]
        minimap = Minimap(work_dir, ref, fq, sort=False)
    unsorted_bam = minimap.run()
    # add minimap task log to result log
    # log["minimap"] = minimap.log

    # pre-process input bam
    bam = readstore.Bam(unsorted_bam)

    # detect amplicon set
    amplicon_set = bam.detect_amplicon_set(amplicon_sets)
    log["amplicons"] = bam.stats

    # construct readstore
    # this subsamples the reads
    rs = readstore.ReadStore(amplicon_set, bam)

    # save reads for viridian assembly
    amp_dir = work_dir / "amplicons"
    manifest_data = rs.make_reads_dir_for_viridian(amp_dir)

    # run viridian
    viridian = Viridian(work_dir, platform, amp_dir, manifest_data, rs.viridian_json)
    consensus = viridian.run()
    log["viridian"] = viridian.log

    # varifier
    varifier = Varifier(
        work_dir / "varifier",
        ref,
        consensus,
        min_coord=rs.start_pos,
        max_coord=rs.end_pos,
    )
    varifier.run()
    log["varifier"] = varifier.log

    # self qc: remap reads to consensus
    pileup = rs.pileup(consensus, msa=msa)

    # mask output
    masked_fasta = pileup.mask()
    # log["self_qc"] = pileup.log
    log["qc"] = pileup.summary

    with open(work_dir / "masked.fasta", "w") as fasta_out:
        print(masked_fasta, file=fasta_out)

    # annotate vcf
    annotated_vcf = pileup.annotate_vcf(vcf)

    with open(work_dir / "final.vcf", "w") as vcf_out:
        header, records = annotated_vcf
        for h in header:
            print(h, file=vcf_out)
        for rec in records:
            print("\t".join(map(str, rec)), file=vcf_out)

    return log
