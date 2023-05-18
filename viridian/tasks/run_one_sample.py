import logging
from viridian import amplicon_schemes, constants, one_sample_pipeline, utils


def run(options):
    utils.check_tech_and_reads_options(options)

    if options.force:
        logging.info(f"--force option used, so deleting {options.outdir} if it exists")
        utils.syscall(f"rm -rf {options.outdir}")

    if options.masking_min_frs is None and options.tech is not None:
        options.masking_min_frs = constants.TECH2FRS[options.tech]

    if options.decontam == "COVID":
        options.decontam = amplicon_schemes.REF_FASTA_NO_POLY_A

    one_sample_pipeline.run_one_sample(
        options.tech,
        options.outdir,
        options.ref_fasta,
        ena_run=options.ena_run,
        keep_ena_reads=options.keep_ena_reads,
        reads_file1=options.reads1,
        reads_file2=options.reads2,
        reads_file=options.reads,
        reads_bam=options.reads_bam,
        decontam_ref_fa=options.decontam,
        built_in_amp_schemes=options.built_in_amp_schemes,
        tsv_of_amp_schemes=options.amp_schemes_tsv,
        force_amp_scheme=options.force_amp_scheme,
        detect_scheme_only=options.detect_scheme_only,
        force_consensus=options.force_consensus,
        debug=options.debug,
        msas_to_write=options.write_msa,
        keep_bam=options.keep_bam,
        qc_depth=options.qc_depth,
        min_scheme_score=options.min_scheme_score,
        max_scheme_ratio=options.max_scheme_ratio,
        sample_name=options.sample_name,
        max_percent_amps_fail=options.max_percent_amps_fail,
        max_cons_n_percent=options.max_cons_n_percent,
        het_min_pc=options.het_min_pc,
        assemble_depth=options.assemble_depth,
        coverage_min_x=options.coverage_min_x,
        coverage_min_pc=options.coverage_min_pc,
        masking_min_frs=options.masking_min_frs,
        masking_min_depth=options.masking_min_depth,
        command_line_args=options,
        temp_root=options.tmp_dir,
        gzip_files=not options.no_gzip,
    )
