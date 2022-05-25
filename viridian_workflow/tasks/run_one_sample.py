import logging
import subprocess
from viridian_workflow import run, utils


def run(options):
    fq1, fq2 = utils.check_tech_and_reads_opts_and_get_reads(options)

    if options.force:
        logging.info(f"--force option used, so deleting {options.outdir} if it exists")
        subprocess.check_output(f"rm -rf {options.outdir}", shell=True)

    # TODO: old code function call was:
    #     one_sample_pipeline.run_one_sample(
    #       options.tech,
    #       options.outdir,
    #       options.ref_fasta,
    #       fq1,
    #       fq2=fq2,
    #       ... etc
    #     )
    # New function run.run_pipeline wants a list of fastq files
    fqs = "TODO"

    # TODO: reinstate functionality that handles the options:
    #  - options.built_in_amp_schemes
    #  - options.amp_schemes_tsv
    #  - options.force_amp_scheme
    # as per spec in the wiki:
    # https://github.com/iqbal-lab-org/viridian_workflow/wiki/Amplicon-schemes#using-any-combination-of-schemes
    #
    # Code that gathered the schemes to be considered was here:
    # https://github.com/iqbal-lab-org/viridian_workflow/blob/cb6c7a55cb145d74fc6352341a76a183b0c1b4e1/viridian_workflow/one_sample_pipeline.py#L134-L156
    # and forcing choice (if requested):
    # https://github.com/iqbal-lab-org/viridian_workflow/blob/cb6c7a55cb145d74fc6352341a76a183b0c1b4e1/viridian_workflow/one_sample_pipeline.py#L183-L191
    amplicon_sets = "TODO"

    # TODO: this needs to run and handle the keyword args
    run.run_pipeline(
        options.oudir,
        options.tech,
        fqs,
        amplicon_sets,
        ref=options.ref_fasta,
        built_in_amp_schemes=options.built_in_amp_schemes,
        tsv_of_amp_schemes=options.amp_schemes_tsv,
        force_amp_scheme=options.force_amp_scheme,
        keep_intermediate=options.debug,
        keep_bam=options.keep_bam,
        sample_name=options.sample_name,
        frs_threshold=options.frs_threshold,
        self_qc_depth=options.self_qc_depth,
        consensus_max_n_percent=options.max_cons_n_percent,
        max_percent_amps_fail=options.max_percent_amps_fail,
        command_line_args=options,
    )
