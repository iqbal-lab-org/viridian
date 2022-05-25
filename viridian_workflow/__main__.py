#!/usr/bin/env python3

import argparse
import logging
import sys
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

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ----------- general options common to all tasks ------------------------
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )
    tech_choices = ["illumina", "ont"]
    common_parser.add_argument(
        "--tech",
        choices=tech_choices,
        help=f"Sequencing technology, currently supported: {','.join(tech_choices)}",
        required=True,
    )

    # -------------- amplicons options ---------------------------------------
    amplicons_parser = argparse.ArgumentParser(add_help=False)
    scheme_names = ",".join(
        sorted(list(viridian_workflow.amplicon_schemes.get_built_in_schemes().keys()))
    )
    amplicons_parser.add_argument(
        "--built_in_amp_schemes",
        help=f"Comma-separated list of built in amplicon schemes to use [{scheme_names}]",
        metavar="scheme1,scheme2,...",
    )
    amplicons_parser.add_argument(
        "--amp_schemes_tsv",
        help="Tab-delimited file of amplicon schemes to use. Must have header line that includes column names 'Name' and 'File'",
        metavar="FILENAME",
    )

    # ----------------- reads and ref options --------------------------------
    reads_ref_parser = argparse.ArgumentParser(add_help=False)
    reads_ref_parser.add_argument(
        "--reads1", help="Illumina reads file 1", metavar="FILENAME",
    )
    reads_ref_parser.add_argument(
        "--reads2", help="Illumina reads file 2", metavar="FILENAME",
    )
    reads_ref_parser.add_argument(
        "--reads", help="Unpaired reads (eg nanopore) file", metavar="FILENAME",
    )
    reads_ref_parser.add_argument(
        "--ref_fasta",
        help="REQUIRED. FASTA file of reference genome",
        required=True,
        metavar="FILENAME",
    )
    reads_ref_parser.add_argument(
        "--sample_name",
        default="sample",
        help="Name of sample to put in header of final FASTA, VCF, and BAM files [%(default)s]",
        metavar="STRING",
    )
    reads_ref_epilog = "IMPORTANT: --tech, --ref_fasta, --outdir are REQUIRED. Reads files are required, and depend on the --tech option. Either use: 1) '--tech ont --reads reads.fq' or 2) '--tech illumina --reads1 reads1.fq --reads2 reads2.fq'."


    # ------------------------ run_one_sample ----------------------------
    subparser_run_one_sample = subparsers.add_parser(
        "run_one_sample",
        parents=[common_parser, amplicons_parser, reads_ref_parser],
        help="Run the complete pipeline on one sample",
        usage=f"viridian_workflow run_one_sample [options] --tech {'|'.join(tech_choices)} --ref_fasta ref.fasta --outdir out <reads options (see help)>",
        description="Run the complete pipeline on one sample",
        epilog=reads_ref_epilog,
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
        "--force_amp_scheme",
        help="Force choice of amplicon scheme. The value provided must exactly match a built-in name or a name in file given by --amp_schemes_tsv",
        metavar="STRING",
    )
    subparser_run_one_sample.add_argument(
        "--target_sample_depth",
        type=int,
        default=1000,
        help="Target coverage for amplicon depth normalisation [%(default)s]",
        metavar="INT",
    )
    subparser_run_one_sample.add_argument(
        "--frs_threshold",
        type=float,
        default=0.7,
        help="Masking threshold for consensus base support [%(default)s]",
        metavar="FLOAT",
    )
    subparser_run_one_sample.add_argument(
        "--self_qc_depth",
        type=int,
        default=10,
        help="Masking threshold for consensus base depth [%(default)s]",
        metavar="INT",
    )
    subparser_run_one_sample.add_argument(
        "--log_liftover", action="store_true", help=argparse.SUPPRESS,
    )
    subparser_run_one_sample.add_argument(
        "--test_amplicon_frs", action="store_true", help=argparse.SUPPRESS,
    )
    subparser_run_one_sample.add_argument(
        "--trim_5prime", action="store_true", help=argparse.SUPPRESS,
    )
    subparser_run_one_sample.add_argument(
        "--min_sample_depth",
        type=int,
        default=10,
        help="Minimum coverage required during amplicon depth normalisation. Any amplicon with depth below this is failed and it will have no consensus sequence generated [%(default)s]",
        metavar="INT",
    )
    subparser_run_one_sample.add_argument(
        "--max_percent_amps_fail",
        type=float,
        default=50.0,
        help="Maximum percent of amplicons allowed to fail during read sampling or making consensus for each amplicon. The pipeline is stopped as soon as too many failed amplicons are detected [%(default)s]",
        metavar="FLOAT",
    )
    subparser_run_one_sample.add_argument(
        "--max_cons_n_percent",
        type=float,
        default=50.0,
        help="Maximum allowed percentage of Ns in the consensus sequence. Pipeline is stopped if too many Ns [%(default)s]",
        metavar="FLOAT",
    )
    subparser_run_one_sample.set_defaults(
        func=viridian_workflow.tasks.run_one_sample.run
    )

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit()
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
