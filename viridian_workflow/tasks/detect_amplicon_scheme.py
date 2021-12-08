import logging
import os
from viridian_workflow import amplicon_schemes, detect_primers, minimap, utils


def run(options):
    fq1, fq2 = utils.check_tech_and_reads_opts_and_get_reads(options)

    logging.info("Gathering amplicons scheme files")
    if options.built_in_amp_schemes is None and options.amp_schemes_tsv is None:
        logging.info("No primer schemes provided. Using all built in schemes")
        options.built_in_amp_schemes = list(
            amplicon_schemes.get_built_in_schemes().keys()
        )
    (
        amplicon_scheme_name_to_tsv,
        amplicon_scheme_list,
    ) = amplicon_schemes.load_list_of_amplicon_sets(
        built_in_names_to_use=options.built_in_amp_schemes,
        tsv_others_to_use=options.amp_schemes_tsv,
    )

    temp_sam = f"{options.outprefix}.tmp.sam"
    logging.info("Mapping reads to reference")
    minimap.run(
        temp_sam,
        options.ref_fasta,
        fq1,
        fq2=fq2,
        sample_name=options.sample_name,
        sort=False,
    )

    logging.info("Detecting primers for each read")
    bam_out = f"{options.outprefix}.bam" if options.make_bam else None
    results = detect_primers.gather_stats_from_bam(
        temp_sam, bam_out, amplicon_scheme_list
    )
    results[
        "amplicon_scheme_set_matches"
    ] = detect_primers.amplicon_set_counts_to_json_friendly(
        results["amplicon_scheme_set_matches"]
    )

    json_out = f"{options.outprefix}.json"
    logging.info(f"Tidying files and writing final JSON {json_out}")
    if not options.debug:
        os.unlink(temp_sam)
    utils.write_json(json_out, results)
