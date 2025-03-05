import logging
from viridian import scheme_simulate, utils


def run(options):
    if options.force:
        logging.info(f"--force option used, so deleting {options.outdir} if it exists")
        utils.syscall(f"rm -rf {options.outdir}")

    scheme_simulate.simulate_all_schemes(
        options.species,
        options.ref_fasta,
        options.outdir,
        built_in_amp_schemes=options.built_in_amp_schemes,
        tsv_of_amp_schemes=options.amp_schemes_tsv,
        read_length=options.read_length,
    )
