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

    one_sample_pipeline.run_one_sample(
        options.tech,
        options.outdir,
        options.ref_fasta,
        options.amplicon_json,
        fq1,
        fq2,
        keep_intermediate=options.debug,
        keep_bam=options.keep_bam,
        target_sample_depth=options.target_sample_depth,
        sample_name=options.sample_name,
    )
