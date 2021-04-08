#!/usr/bin/env python3

import argparse
import logging
import viridian_workflow


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="viridian_workflow",
        usage="viridian_workflow <command> <options>",
        description="viridian_workflow: ...",  # FIXME
    )

    parser.add_argument("--version", action="version", version=viridian_workflow.__version__)
    parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ------------------------ run_one_sample ----------------------------
    subparser_run_one_sample = subparsers.add_parser(
        "run_one_sample",
        help="Help for run_one_sample",
        usage="viridian_workflow run_one_sample [options] <ref_fasta>",
        description="description of run_one_sample",
    )
    subparser_run_one_sample.add_argument(
        "ref_fasta",
        help="FASTA file of reference genome",
    )
    subparser_run_one_sample.add_argument(
        "amplicon_bed",
        help="Amplicon coordinate file",
    )
    subparser_run_one_sample.add_argument(
        "fastq1",
        help="fq1",
    )
    subparser_run_one_sample.add_argument(
        "fastq2",
        help="fq2",
    )
    subparser_run_one_sample.add_argument(
        "outdir",
        help="outdir",
    )

    subparser_run_one_sample.set_defaults(func=viridian_workflow.tasks.run_one_sample.run)

    args = parser.parse_args()

    logging.basicConfig(
        format="[%(asctime)s viridian_workflow %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    log = logging.getLogger()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
