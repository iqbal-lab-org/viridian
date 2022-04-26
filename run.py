import sys
import shutil
from pathlib import Path
import tempfile
import json

from viridian_workflow import primers, readstore, utils, self_qc
from viridian_workflow.subtasks import Minimap, Varifier, Viridian


def run_pipeline(work_dir, platform, fqs):
    log = {}
    # generate name-sorted bam from fastqs
    ref = Path("../covid/MN908947.fasta")

    if platform == "illumina":
        fq1, fq2 = fqs
        minimap = Minimap(work_dir / "name_sorted.bam", ref, fq1, fq2=fq2, sort=False)
    elif platform == "ont":
        fq = fqs[0]
        minimap = Minimap(work_dir / "name_sorted.bam", ref, fq, sort=False)
    elif platform == "iontorrent":
        raise NotImplementedError
    else:
        print(f"Platform {platform} is not supported.", file=sys.stderr)
        exit(1)

    unsorted_bam = minimap.run()
    # add minimap task log to result log
    # log["minimap"] = minimap.log

    # pre-process input bam
    bam = readstore.Bam(unsorted_bam)

    # detect amplicon set
    amplicon_set = bam.detect_amplicon_set(amplicon_sets)
    # log["amplicons"] = bam.stats

    # construct readstore
    # this subsamples the reads
    rs = readstore.ReadStore(amplicon_set, bam)
    log["amplcions"] = rs.summary

    # save reads for viridian assembly
    amp_dir = work_dir / "amplicons"
    manifest_data = rs.make_reads_dir_for_viridian(amp_dir)

    # run viridian
    viridian = Viridian(
        work_dir, platform, ref, amp_dir, manifest_data, rs.viridian_json
    )
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
    vcf, msa, varifier_consensus = varifier.run()
    log["varifier"] = varifier.log

    # self qc: remap reads to consensus
    pileup = rs.pileup(varifier_consensus, msa=msa)

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


if __name__ == "__main__":
    # set up the pipeline
    platform, *fqs = sys.argv[1:]

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

    work_dir = "/tmp/vwf/"
    work_dir = Path(work_dir)
    if work_dir.exists():
        print(f"work dir {work_dir} exists, clobbering", file=sys.stderr)
        shutil.rmtree(work_dir)
    work_dir.mkdir()

    log = {"summary": {"version": "test-0.1", "status": "Interrupted"}}
    # try:
    if True:
        results = run_pipeline(work_dir, platform, fqs)
        log["results"] = results
        log["summary"]["status"] = "Success"
    # except Exception as e:
    elif False:
        log["summary"]["status"] = {"Failure": str(e)}
        print(f"Pipeline failed with exception: {e}", file=sys.stderr)

    with open(work_dir / "log.json", "w") as json_out:
        json.dump(log, json_out, indent=4)
