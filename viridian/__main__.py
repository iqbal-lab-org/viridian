#!/usr/bin/env python3

import argparse
import logging
import sys
import viridian


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="viridian",
        usage="viridian <command> <options>",
        description="viridian: amplicon-aware consensus maker",
    )
    parser.add_argument("--version", action="version", version=viridian.__version__)

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ----------- general options common to all tasks ------------------------
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )
    common_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output directory, if it already exists. Use with caution!",
    )

    # -------------- outdir --------------------------------------------------
    outdir_parser = argparse.ArgumentParser(add_help=False)
    outdir_parser.add_argument(
        "--outdir",
        help="REQUIRED. Name of output directory (will be created). This must not exist already, unless the --force option is used to overwrite",
        required=True,
        metavar="FILENAME",
    )

    # -------------- amplicons options ---------------------------------------
    amplicons_parser = argparse.ArgumentParser(add_help=False)
    scheme_names = ",".join(
        sorted(list(viridian.amplicon_schemes.get_built_in_schemes().keys()))
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

    # --------------------- reads options ------------------------------------
    reads_parser = argparse.ArgumentParser(add_help=False)
    reads_parser.add_argument(
        "--run_accession",
        help="INSDC run accession. Reads will be downloaded from the ENA and --tech determined from ENA metadata",
        metavar="RUN_ID",
    )
    reads_parser.add_argument(
        "--keep_downloaded_reads",
        help="Keep reads downloaded from the ENA. Only relevant if --run_accession is used, otherwise is ignored",
        action="store_true",
    )
    reads_parser.add_argument(
        "--reads1",
        help="Forwards reads file",
        metavar="FILENAME",
    )
    reads_parser.add_argument(
        "--reads2",
        help="Reverse reads file",
        metavar="FILENAME",
    )
    reads_parser.add_argument(
        "--reads",
        help="Unpaired reads file",
        metavar="FILENAME",
    )
    reads_parser.add_argument(
        "--reads_bam",
        help="Sorted by coordinate and indexed BAM file of reads mapped to the reference genome (reference must be same as that given by --ref_fasta). Incompatible with --reads* and --decontam options",
        metavar="FILENAME",
    )
    reads_parser.add_argument(
        "--sample_name",
        default="sample",
        help="Name of sample to put in various output files (eg final FASTA, VCF, and BAM) [%(default)s]",
        metavar="STRING",
    )
    tech_choices = sorted(list(viridian.constants.ALLOWED_TECH))
    reads_parser.add_argument(
        "--tech",
        choices=tech_choices,
        help=f"Sequencing technology, currently supported: {','.join(tech_choices)}",
    )
    reads_epilog = "IMPORTANT: --outdir is REQUIRED. Unless --run_accession is used, reads files are required, and depend on the --tech option. Use one of: 1) '--tech illumina|ont|iontorrent --reads reads.fq'; 2) '--tech illumina|iontorrent --reads1 reads1.fq --reads2 reads2.fq'; 3) (with any tech) --reads_bam reads.bam"

    # --------------------advanced options -----------------------------------
    advanced_parser = argparse.ArgumentParser(add_help=False)
    advanced_parser.add_argument(
        "--decontam",
        help="Decontaminate input reads at start of pipeline, using ReadItAndKeep with the provided reference genome FASTA file (use '--decontam COVID' to use the MN908947.3 SARS-CoV-2 reference genome with poly-A removed). Incompatible with --reads_bam",
        metavar="FILENAME",
    )
    advanced_parser.add_argument(
        "--keep_bam",
        action="store_true",
        help="Keep BAM file of original reads mapped to reference genome (it is deleted by default)",
    )
    msa_choices = sorted(list(viridian.constants.MSA_CHOICES))
    advanced_parser.add_argument(
        "--write_msa",
        action="append",
        choices=msa_choices,
        help=f"MSA(s) to write. Use this option once for each MSA type to write. Each written to separate FASTA file. Choose from: {';'.join(msa_choices)}",
        metavar="MSA_TYPE",
    )
    advanced_parser.add_argument(
        "--force_amp_scheme",
        help="Force choice of amplicon scheme. The value provided must exactly match a built-in name or a name in file given by --amp_schemes_tsv",
        metavar="STRING",
    )
    advanced_parser.add_argument(
        "--min_scheme_score",
        type=int,
        default=250,
        help="Minimum score required for a matching scheme. If all schemes are less than this cutoff, the pipeline is stopped unless --force_amp_scheme is used [%(default)s]",
        metavar="INT",
    )
    advanced_parser.add_argument(
        "--max_scheme_ratio",
        type=float,
        default=0.5,
        help="Maximum allowed value of (second best scheme score) / (best scheme score). See also --min_scheme_score [%(default)s]",
        metavar="FLOAT",
    )
    advanced_parser.add_argument(
        "--assemble_depth",
        type=int,
        help="Target read depth for each amplicon when assembling [%(default)s]",
        metavar="INT",
        default=100,
    )
    advanced_parser.add_argument(
        "--qc_depth",
        type=int,
        default=1000,
        help="Target coverage for each amplicon, when running QC and masking [%(default)s]",
        metavar="INT",
    )
    frs_default_str = "; ".join(
        [f"{k}:{v}" for k, v in viridian.constants.TECH2FRS.items()]
    )
    advanced_parser.add_argument(
        "--masking_min_frs",
        type=float,
        help=f"Consensus positions with a fraction of read support less than the given number are masked. Default depends on the --tech option [{frs_default_str}]",
        metavar="FLOAT",
    )
    advanced_parser.add_argument(
        "--het_min_pc",
        type=float,
        help="Cutoff for when to use ambiguous IUPAC codes. With --het_min_pc X, if two or more alleles each have support of more than X percent of reads, then an ambiguous IUPAC code is used and the position is given the 'HET' filter [%(default)s]",
        metavar="FLOAT",
        default=20.0,
    )
    advanced_parser.add_argument(
        "--masking_min_depth",
        type=int,
        default=20,
        help="Consensus positions with total read depth less than the given number are masked [%(default)s]",
        metavar="INT",
    )
    advanced_parser.add_argument(
        "--coverage_min_x",
        default=20,
        type=int,
        help="With options --coverage_min_x X --coverage_min_pc Y, after initial read mapping to the reference, if less than Y percent of the genome has at least X read coverage, then the pipeline is stopped. Default is to require at least 50 percent of the genome with at least 20X read depth",
        metavar="INT",
    )
    advanced_parser.add_argument(
        "--coverage_min_pc",
        default=50,
        type=float,
        help="Please see the help for --coverage_min_x",
        metavar="FLOAT",
    )
    advanced_parser.add_argument(
        "--max_percent_amps_fail",
        type=float,
        default=50.0,
        help="Maximum percent of amplicons allowed to fail during read sampling or making consensus for each amplicon. The pipeline is stopped as soon as too many failed amplicons are detected [%(default)s]",
        metavar="FLOAT",
    )
    advanced_parser.add_argument(
        "--max_cons_n_percent",
        type=float,
        default=50.0,
        help="Maximum allowed percentage of Ns in the consensus sequence. Pipeline is stopped as soon as too many Ns in current consensus [%(default)s]",
        metavar="FLOAT",
    )
    advanced_parser.add_argument(
        "--detect_scheme_only",
        action="store_true",
        help="Map the reads and detect the amplicon scheme, then stop the pipeline instead of making a consensus etc",
    )
    advanced_parser.add_argument(
        "--force_consensus",
        help="Use the provided FASTA for the consensus, instead of inferring it from the reads. This option is useful to evaluate the quality of an existing consensus sequence using the reads that made it",
        metavar="FILENAME",
    )
    advanced_parser.add_argument(
        "--tmp_dir",
        help="Directory in which to put a temporary directory for intermediate files that are not kept. The dir given by --tmp_dir must already exist. Default is platform dependent, and uses python's tempdir default locations (eg for linux, looks for /tmp first). If you use the --debug option, then the --tmp_dir option is ignored and a directory called 'Processing' is used in the output directory",
        metavar="DIRNAME",
    )
    advanced_parser.add_argument(
        "--no_gzip",
        help="Use this to not gzip the final output FASTA file and log.json file",
        action="store_true",
    )
    advanced_parser.add_argument(
        "--force_mafft",
        help="Force use of mafft for global align between consensus and ref, instead of nucmer-based method. By default, mafft only used when the reference genome is longer than 30kbp",
        action="store_true",
    )

    # --------------------------- ref options --------------------------------
    ref_parser = argparse.ArgumentParser(add_help=False)
    ref_parser.add_argument(
        "--ref_fasta",
        help="FASTA file of reference genome. Default is to use built-in MN908947.3 SARS-CoV-2 reference genome. Only use this if you know what you are doing, since the amplicon scheme(s) must match this reference [%(default)s]",
        metavar="FILENAME",
        default=viridian.amplicon_schemes.REF_FASTA,
    )

    # ------------------------- dl_and_run options -----------------------
    dl_and_run_parser = argparse.ArgumentParser(add_help=False)
    dl_and_run_parser.add_argument(
        "--acc_file",
        help="REQUIRED. File containing run accessions, one run per line",
        metavar="FILENAME",
        required=True,
    )
    dl_and_run_parser.add_argument(
        "--outdir",
        help="REQUIRED. Name of output directory. If it already exists, will assume it is from a previous run of `viridian dl_and_run`, keep existing runs of viridian and only run samples that did not finish",
        required=True,
        metavar="FILENAME",
    )

    # ------------------------ run_one_sample ----------------------------
    subparser_run_one_sample = subparsers.add_parser(
        "run_one_sample",
        parents=[
            common_parser,
            outdir_parser,
            amplicons_parser,
            reads_parser,
            ref_parser,
            advanced_parser,
        ],
        help="Run the complete pipeline on one sample",
        usage=f"viridian run_one_sample [options] --tech {'|'.join(tech_choices)} --outdir out <reads options (see help)>",
        description="Run the complete pipeline on one sample",
        epilog=reads_epilog,
    )
    subparser_run_one_sample.set_defaults(func=viridian.tasks.run_one_sample.run)

    # ------------------------ dl_and_run --------------------------------
    subparser_dl_and_run = subparsers.add_parser(
        "dl_and_run",
        parents=[
            dl_and_run_parser,
            common_parser,
            amplicons_parser,
            ref_parser,
            advanced_parser,
        ],
        help="Download reads from the ENA and run the pipeline, for multiple sequencing runs",
        usage="viridian dl_and_run [options] --acc_file foo.txt --outdir out",
        description="Download reads from the ENA and run the pipeline, for multiple sequencing runs",
    )
    subparser_dl_and_run.set_defaults(func=viridian.tasks.dl_and_run.run)

    # ------------------------ sim_schemes -------------------------------
    subparser_sim_schemes = subparsers.add_parser(
        "sim_schemes",
        parents=[amplicons_parser, common_parser, outdir_parser, ref_parser],
        help="For each scheme, simulate fragments and run scheme id code",
        usage="viridian sim_schemes [options] --outdir <out>",
        description="For each scheme, simulate fragments and run scheme id code",
    )
    subparser_sim_schemes.add_argument(
        "--read_length",
        help="Length of fragmented reads [%(default)s]",
        type=int,
        default=200,
        metavar="INT",
    )
    subparser_sim_schemes.set_defaults(func=viridian.tasks.sim_schemes.run)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit()

    logging.basicConfig(
        format="[%(asctime)s viridian %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
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
