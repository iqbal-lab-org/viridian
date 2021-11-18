import logging
import subprocess
from viridian_workflow import one_sample_pipeline


def run(options):
    if options.tech == "ont":
        assert options.reads1 is None and options.reads2 is None
        assert options.reads is not None
        fq1 = options.reads
        fq2 = None
    elif options.tech == "illumina":
        assert options.reads is None
        assert options.reads1 is not None and options.reads2 is not None
        fq1 = options.reads1
        fq2 = options.reads2

    if options.force:
        logging.info(f"--force option used, so deleting {options.outdir} if it exists")
        subprocess.check_output(f"rm -rf {options.outdir}", shell=True)

    one_sample_pipeline.run_one_sample(
        options.tech,
        options.outdir,
        options.ref_fasta,
        fq1,
        fq2,
        options.amplicon_json,
        keep_intermediate=options.debug,
        keep_bam=options.keep_bam,
        target_sample_depth=options.target_sample_depth,
        sample_name=options.sample_name,
        command_line_args=options,
    )
