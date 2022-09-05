import logging
import sys
import os
import socket
import time
import subprocess
from pathlib import Path
from viridian_workflow import utils, primers, amplicon_schemes
from viridian_workflow.run import run_pipeline


def cuckoo(options):
    run(options, force_consensus=options.force_consensus)


def run(options, force_consensus=None):
    fq1, fq2 = utils.check_tech_and_reads_opts_and_get_reads(options)

    log: dict[str, Any] = {}

    log["Summary"]["command"] = " ".join(sys.argv)
    log["Summary"]["Version"] = ""
    log["Summary"]["Finished_running"] = False
    log["Summary"]["Success"] = False
    log["Summary"]["Progress"] = []
    log["Summary"]["cwd"] = os.getcwd()
    log["Summary"]["hostname"] = socket.gethostname()
    start_time = time.time()
    log["Summary"]["start_time"] = start_time

    log["Summary"]["options"] = {}
    for option, setting in options.__dict__.items():
        log["Summary"]["options"][str(option)] = setting

    if options.force:
        logging.info(f"--force option used, so deleting {options.outdir} if it exists")
        subprocess.check_output(f"rm -rf {options.outdir}", shell=True)

    # New function run.run_pipeline wants a list of fastq files
    fqs = [
        fq1,
    ]
    if fq2 is not None:
        fqs = [fq1, fq2]

    # Build the index of built-in schemes, possibly subsetted
    data_dir = Path(amplicon_schemes.__file__).resolve().parent / "amplicon_scheme_data"
    amplicon_index = amplicon_schemes.load_amplicon_index(
        Path("schemes.tsv"), data_dir, subset=options.built_in_amp_schemes
    )

    # If a set is forced, select it from the possibly subsetted built-ins
    chosen_amplicon_set = None
    if options.force_amp_scheme:
        # If they're forcing an amplicon scheme but have disabled all built-ins
        # this is an error. We may want to allow this to enable them to force
        # a custom scheme
        if options.amp_schemes_tsv and not options.built_in_schemes:
            raise Exception("Can only force amplicon scheme from built-in options")

        if options.force_amp_scheme in amplicon_index:
            chosen_amplicon_set = primers.AmpliconSet.from_tsv(
                amplicon_index[options.force_amp_scheme], name=options.force_amp_scheme
            )
        else:
            raise Exception(
                f"Chose to force amplicons scheme to be {options.force_amp_scheme}, but scheme not found. Found these: {','.join(amplicon_index.keys())}"
            )

    if options.amp_schemes_tsv:
        # if the user brings their own tsv index, ignore the built in set,
        # unless they also specified a subset from the built in set
        if options.built_in_amp_schemes:
            for name, scheme in load_amplicon_index(options.amp_schemes_tsv).items():
                amplicon_index[name] = scheme
        else:
            amplicon_index = load_amplicon_index(options.amp_schemes_tsv)

    amplicon_sets = [
        primers.AmpliconSet.from_tsv(tsv, name=name)
        for name, tsv in amplicon_index.items()
    ]

    try:
        pipeline_results = run_pipeline(
            options.outdir,
            options.tech,
            fqs,
            amplicon_sets,
            ref=options.ref_fasta,
            force_amp_scheme=chosen_amplicon_set,
            keep_intermediate=options.debug,
            keep_bam=options.keep_bam,
            dump_tsv=options.dump_tsv,
            sample_name=options.sample_name,
            frs_threshold=options.frs_threshold,
            self_qc_depth=options.self_qc_depth,
            consensus_max_n_percent=options.max_cons_n_percent,
            max_percent_amps_fail=options.max_percent_amps_fail,
            command_line_args=options,
            force_consensus=force_consensus,
            log=log,
        )
        log["Results"] = pipeline_log
        log["Summary"]["Success"] = True
    except Exception as e:
        log["Summary"]["Success"] = False
        # log["Summary"]["status"] = {"Failure": str(e)}
        print(f"Pipeline failed with exception: {e}", file=sys.stderr)

    with open(work_dir / "log.json", "w") as json_out:
        end_time = time.time()
        log["Summary"]["end_time"] = end_time
        log["Summary"]["run_time"] = end_time - start_time
        json.dump(log, json_out, indent=2)
