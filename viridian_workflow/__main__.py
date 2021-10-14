#!/usr/bin/env python3

import argparse
import logging
import viridian_workflow


def check_reads_args(args):
    if args.tech == "ont":
        if args.reads1 is not None or args.reads2 is not None or args.reads is None:
            raise Exception(
                "When the --tech option is 'ont', please provide reads using the --reads option, and do not use --reads1 or --reads2."
            )
    elif args.tech == "illumina":
        if args.reads is not None or args.reads1 is None or args.reads2 is None:
            raise Exception(
                "When the --tech option is 'illumina', please provide reads using both options --reads1 and --reads2, and do not use the option --reads."
            )
    else:
        raise NotImplementedError(f"tech {args.tech} not implemented")


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="viridian_workflow",
        usage="viridian_workflow <command> <options>",
        description="viridian_workflow: amplicon-aware consensus maker",
    )
    parser.add_argument(
        "--version", action="version", version=viridian_workflow.__version__
    )
    parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ------------------------ run_one_sample ----------------------------
    tech_choices = ["illumina", "ont"]
    subparser_run_one_sample = subparsers.add_parser(
        "run_one_sample",
        help="Help for run_one_sample",
        usage=f"viridian_workflow run_one_sample --tech {'|'.join(tech_choices)} --ref_fasta ref.fasta --outdir out --amplicon_json <amplicon.json> <reads options (see help)>",
        description="run_one_sample: runs the pipeline on one sample",
        epilog="IMPORTANT: --tech, --ref_fasta, --outdir are REQUIRED. Reads files are required, and depend on the --tech option. Either use: 1) '--tech ont --reads reads.fq' or 2) '--tech illumina --reads1 reads1.fq --reads2 reads2.fq'.",
    )
    subparser_run_one_sample.add_argument(
        "--tech",
        choices=tech_choices,
        help=f"Sequencing technology, currently supported: {','.join(tech_choices)}",
        required=True,
    )
    subparser_run_one_sample.add_argument(
        "--ref_fasta",
        help="REQUIRED. FASTA file of reference genome",
        required=True,
        metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--amplicon_json",
        help="REQUIRED. JSON file of amplicons and primers",
        required=True,
        metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--reads1",
        help="Illumina reads file 1",
        metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--reads2",
        help="Illumina reads file 2",
        metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--reads",
        help="Unpaired reads (eg nanopore) file",
        metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--outdir",
        help="REQUIRED. Name of output directory (will be created). This must not exist already, unless the --force option is used to overwrite",
        required=True,
        metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output directory, if it already exists. Use with caution!",
    )
    subparser_run_one_sample.add_argument(
        "--keep_bam",
        action="store_true",
        help="Keep BAM file of reads mapped to reference genome (it is deleted by default)",
    )
    subparser_run_one_sample.add_argument(
        "--target_sample_depth",
        type=int,
        default=1000,
        help="Target coverage for amplicon depth normalisation [%(default)s]",
        metavar="INT",
    )

    subparser_run_one_sample.set_defaults(
        func=viridian_workflow.tasks.run_one_sample.run
    )

    args = parser.parse_args()
    check_reads_args(args)

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
