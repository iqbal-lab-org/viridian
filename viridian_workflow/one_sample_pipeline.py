import argparse
import datetime
import json
import logging
import socket
import subprocess
import sys
import os

from viridian_workflow import (
    amplicon_schemes,
    detect_primers,
    minimap,
    qcovid,
    sample_reads,
    utils,
    varifier,
)
from viridian_workflow import __version__ as viridian_wf_version


def check_tech(tech):
    allowed_tech = {"ont", "illumina"}
    if tech not in allowed_tech:
        techs = ",".join(sorted(list(allowed_tech)))
        raise Exception(f"Tech '{tech}' not allowed, must be one of: {techs}")


def run_viridian(tech, outdir, ref_genome, amplicon_json, bam, bad_amplicons):
    check_tech(tech)
    viridian_cmd = " ".join(
        [
            "viridian",
            "assemble",
            "--bam",
            bam,
            tech,
            "--amplicons_to_fail_file",
            bad_amplicons,
            ref_genome,
            amplicon_json,
            outdir,
        ]
    )
    assembly = os.path.join(outdir, "consensus.final_assembly.fa")
    utils.run_process(viridian_cmd)
    utils.check_file(assembly)
    return assembly


def update_json_latest_stage(log_dict, latest_stage, outfile):
    log_dict["run_summary"]["last_stage_completed"] = latest_stage
    utils.write_json(outfile, log_dict)


def run_one_sample(
    tech,
    outdir,
    ref_genome,
    fq1,
    fq2=None,
    built_in_amp_schemes=None,
    tsv_of_amp_schemes=None,
    force_amp_scheme=None,
    keep_intermediate=False,
    keep_bam=False,
    target_sample_depth=1000,
    sample_name="sample",
    command_line_args=None,
):
    check_tech(tech)
    if tech == "ont":
        assert fq2 is None
        paired = False
    elif tech == "illumina":
        assert fq2 is not None
        paired = True
    else:
        raise NotImplementedError(f"tech not implemented: {tech}")

    logging.info(f"Start running viridian_workflow, output dir: {outdir}")
    os.mkdir(outdir)
    json_log = os.path.join(outdir, "log.json")
    # Make a dict of the command line options to go in the JSON output file.
    # The tests don't use argparse (they use Mock), which means convert to dict
    # doesn't work. Don't care about that case anyway in the final output, so
    # just set to None
    if isinstance(command_line_args, argparse.Namespace):
        options_dict = {k: v for k, v in vars(command_line_args).items() if k != "func"}
    else:
        options_dict = None
    start_time = datetime.datetime.now()
    log_info = {
        "run_summary": {
            "last_stage_completed": "start",
            "command": " ".join(sys.argv),
            "options": options_dict,
            "cwd": os.getcwd(),
            "version": viridian_wf_version,
            "finished_running": False,
            "start_time": start_time.replace(microsecond=0).isoformat(),
            "end_time": None,
            "hostname": socket.gethostname(),
        },
        "read_and_primer_stats": None,
        "read_sampling": None,
        "viridian": None,
    }
    utils.write_json(json_log, log_info)

    if built_in_amp_schemes is None and tsv_of_amp_schemes is None:
        logging.info("No primer schemes provided. Using all built in schemes")
        built_in_amp_schemes = list(amplicon_schemes.get_built_in_schemes().keys())
    (
        amplicon_scheme_name_to_tsv,
        amplicon_scheme_list,
    ) = amplicon_schemes.load_list_of_amplicon_sets(
        built_in_names_to_use=built_in_amp_schemes, tsv_others_to_use=tsv_of_amp_schemes
    )
    if force_amp_scheme is not None:
        if force_amp_scheme not in amplicon_scheme_name_to_tsv:
            names = ",".join(sorted(list(amplicon_scheme_name_to_tsv.keys())))
            raise Exception(
                f"Chose to force amplicons scheme to be {force_amp_scheme}, but scheme not found. Found these: {names}"
            )
    update_json_latest_stage(log_info, "Processed amplicon scheme files", json_log)

    processing_dir = os.path.join(outdir, "Processing")
    os.mkdir(processing_dir)
    logging.info("Mapping reads to reference")
    unsorted_bam = os.path.join(processing_dir, "map_reads.unsorted.bam")
    minimap.run(
        unsorted_bam, ref_genome, fq1, fq2=fq2, sample_name=sample_name, sort=False
    )
    update_json_latest_stage(log_info, "Initial map reads", json_log)

    logging.info("Detecting amplicon scheme and gathering read statistics")
    unsorted_read_tagged_bam = os.path.join(
        processing_dir, "map_reads.unsorted.read_tagged.bam"
    )
    log_info["read_and_primer_stats"] = detect_primers.gather_stats_from_bam(
        unsorted_bam, unsorted_read_tagged_bam, amplicon_scheme_list
    )
    log_info["read_and_primer_stats"][
        "amplicon_scheme_set_matches"
    ] = detect_primers.amplicon_set_counts_to_json_friendly(
        log_info["read_and_primer_stats"]["amplicon_scheme_set_matches"]
    )
    if force_amp_scheme is None:
        log_info["amplicon_scheme_name"] = log_info["read_and_primer_stats"][
            "chosen_amplicon_scheme"
        ]
    else:
        log_info["amplicon_scheme_name"] = force_amp_scheme
    amplicon_tsv = amplicon_scheme_name_to_tsv[log_info["amplicon_scheme_name"]]
    update_json_latest_stage(
        log_info, "Gather read stats and detect primer scheme", json_log
    )

    logging.info("Processing files for chosen amplicon scheme")
    amplicon_json = os.path.join(processing_dir, "amplicons.json")
    amplicon_schemes.convert_tsv_to_viridian_json(
        amplicon_tsv, amplicon_json, scheme_name=log_info["amplicon_scheme_name"]
    )
    amplicon_bed = os.path.join(processing_dir, "amplicons.bed")
    amplicons_start, amplicons_end = utils.amplicons_json_to_bed_and_range(
        amplicon_json, amplicon_bed
    )
    update_json_latest_stage(
        log_info, "Processed chosen amplicon scheme files", json_log
    )

    logging.info("Sorting and indexing BAM file of all reads mapped to reference")
    all_reads_bam = os.path.join(outdir, "reference_mapped.bam")
    utils.run_process(
        f"samtools sort -O BAM -o {all_reads_bam} {unsorted_read_tagged_bam}"
    )
    utils.run_process(f"samtools index {all_reads_bam}")
    update_json_latest_stage(
        log_info, "Sorted and indexed all reads BAM file", json_log
    )

    logging.info("Sampling reads")
    sample_outprefix = os.path.join(processing_dir, "sample_reads")
    sampler = sample_reads.sample_reads(
        ref_genome,
        all_reads_bam,
        sample_outprefix,
        amplicon_bed,
        target_depth=target_sample_depth,
    )
    log_info["read_sampling"] = utils.load_json(f"{sample_outprefix}.json")
    bam = sampler.bam_out
    update_json_latest_stage(log_info, "Sampled reads", json_log)

    logging.info(f"Running QC on all mapped reads in {bam}")
    if paired:
        bad_amplicons = qcovid.bin_amplicons(
            processing_dir, ref_genome, amplicon_bed, bam
        )
    else:
        bad_amplicons = qcovid.bin_amplicons_se(
            processing_dir, ref_genome, amplicon_bed, bam
        )
    update_json_latest_stage(log_info, "QC on all mapped reads", json_log)

    logging.info("Making initial unmasked consensus using Viridian")
    viridian_out = os.path.join(processing_dir, "viridian")
    assembly = run_viridian(
        tech, viridian_out, ref_genome, amplicon_json, bam, bad_amplicons
    )
    log_info["viridian"] = utils.load_json(os.path.join(viridian_out, "run_info.json"))
    update_json_latest_stage(log_info, "Viridian", json_log)

    logging.info("Mapping reads to consensus from Viridian")
    self_map_bam = os.path.join(processing_dir, "self_qc.bam")
    if paired:
        fq1 = sampler.fq_out1
        fq2 = sampler.fq_out2
    else:
        fq1 = sampler.fq_out
        fq2 = None
    minimap.run(
        self_map_bam,
        assembly,
        fq1,
        fq2=fq2,
        sample_name=sample_name,
        sort=True,
    )
    update_json_latest_stage(log_info, "Map reads to Viridian consensus", json_log)

    logging.info("Running QC on Viridian consensus to make masked FASTA")
    masked_fasta = qcovid.self_qc(processing_dir, assembly, self_map_bam)
    final_masked_fasta = os.path.join(outdir, "consensus.fa")
    utils.set_seq_name_in_fasta_file(masked_fasta, final_masked_fasta, sample_name)
    update_json_latest_stage(log_info, "Ran QC on reads mapped to consensus", json_log)

    logging.info("Making VCF file of variants")
    varifier_out = os.path.join(processing_dir, "varifier")
    varifier_vcf = varifier.run(
        varifier_out,
        ref_genome,
        final_masked_fasta,
        min_coord=amplicons_start,
        max_coord=amplicons_end,
    )
    utils.check_file(varifier_vcf)
    final_vcf = os.path.join(outdir, "variants.vcf")
    utils.set_sample_name_in_vcf_file(varifier_vcf, final_vcf, sample_name)
    utils.check_file(final_vcf)
    update_json_latest_stage(log_info, "Made VCF of variants", json_log)

    # clean up intermediate files
    if not keep_intermediate:
        logging.info("Deleting temporary files")
        if keep_bam:
            logging.info(
                f"Keeping BAM file {all_reads_bam} because --keep_bam option used"
            )
        else:
            utils.rm(all_reads_bam)
            utils.rm(all_reads_bam + ".bai")
        logging.info(f"Removing processing directory {processing_dir}")
        subprocess.check_output(f"rm -rf {processing_dir}", shell=True)
    else:
        logging.info("Debug mode: not deleting temporary files")

    logging.info(f"Writing JSON log file {json_log}")
    end_time = datetime.datetime.now()
    log_info["run_summary"]["end_time"] = end_time.replace(microsecond=0).isoformat()
    log_info["run_summary"]["run_time"] = str(end_time - start_time)
    log_info["run_summary"]["finished_running"] = True
    update_json_latest_stage(log_info, "Finished", json_log)
    logging.info(f"Finished running viridian_workflow, output dir: {outdir}")
