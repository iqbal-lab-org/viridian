import logging
import os

from viridian import tasks, utils


def viridian_finished(indir):
    try:
        log = utils.load_json(os.path.join(indir, "log.json.gz"))
        last_stage = log["run_summary"]["last_stage_completed"]
    except:
        return False

    return last_stage == "Finished"


def run(options):
    with open(options.acc_file) as f:
        run_accessions = [x.rstrip() for x in f]
    logging.info(f"{len(run_accessions)} runs found in file {options.acc_file}")

    if not os.path.exists(options.outdir):
        logging.info(f"Making main output directory {options.outdir}")
        os.mkdir(options.outdir)
    else:
        logging.info(f"Existing output directory found, using it: {options.outdir}")

    options.keep_downloaded_reads = False
    options.reads1 = None
    options.reads2 = None
    options.reads = None
    options.reads_bam = None
    options.tech = None
    options.force = True

    os.chdir(options.outdir)

    for i, run in enumerate(run_accessions):
        message = f" STARTING {i+1}/{len(run_accessions)} {run} "
        logging.info(f"{message:#^60}")
        options.outdir = run
        options.run_accession = run
        options.sample_name = run
        if viridian_finished(options.outdir):
            logging.info(f"Already done, skipping: {run}")
        else:
            try:
                tasks.run_one_sample.run(options)
            except:
                logging.warn(f"Something went wrong: {run}")

        message = f" FINISHED {run} "
        logging.info(f"{message:#^60}")
