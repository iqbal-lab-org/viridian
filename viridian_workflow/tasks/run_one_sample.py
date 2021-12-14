import logging
import subprocess
from viridian_workflow import one_sample_pipeline, utils


def run(options):
    fq1, fq2 = utils.check_tech_and_reads_opts_and_get_reads(options)

    if options.force:
        logging.info(f"--force option used, so deleting {options.outdir} if it exists")
        subprocess.check_output(f"rm -rf {options.outdir}", shell=True)

    one_sample_pipeline.run_one_sample(
        options.tech,
        options.outdir,
        options.ref_fasta,
        fq1,
        fq2=fq2,
        built_in_amp_schemes=options.built_in_amp_schemes,
        tsv_of_amp_schemes=options.amp_schemes_tsv,
        force_amp_scheme=options.force_amp_scheme,
        keep_intermediate=options.debug,
        keep_bam=options.keep_bam,
        target_sample_depth=options.target_sample_depth,
        sample_name=options.sample_name,
        command_line_args=options,
    )
