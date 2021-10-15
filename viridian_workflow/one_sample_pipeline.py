import argparse
import datetime
import json
import logging
import socket
import subprocess
import sys
import os

from viridian_workflow import minimap, qcovid, sample_reads, varifier
from viridian_workflow.utils import (
    amplicons_json_to_bed_and_range,
    check_file,
    load_json,
    run_process,
    rm,
    set_sample_name_in_vcf_file,
    set_seq_name_in_fasta_file,
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
    run_process(viridian_cmd)
    check_file(assembly)
    return assembly


def run_one_sample(
    tech,
    outdir,
    ref_genome,
    amplicon_json,
    fq1,
    fq2=None,
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
            "command": " ".join(sys.argv),
            "options": options_dict,
            "cwd": os.getcwd(),
            "version": viridian_wf_version,
            "finished_running": False,
            "start_time": start_time.replace(microsecond=0).isoformat(),
            "end_time": None,
            "hostname": socket.gethostname(),
        },
        "read_sampling": None,
        "viridian": None,
    }
    with open(json_log, "w") as f:
        json.dump(log_info, f, indent=2)

    processing_dir = os.path.join(outdir, "Processing")
    os.mkdir(processing_dir)
    amplicon_bed = os.path.join(processing_dir, "amplicons.bed")
    amplicons_start, amplicons_end = amplicons_json_to_bed_and_range(
        amplicon_json, amplicon_bed
    )
    if paired:
        all_reads_bam = minimap.run(
            outdir, ref_genome, fq1, fq2, sample_name=sample_name
        )
    else:
        all_reads_bam = minimap.run_se(outdir, ref_genome, fq1, sample_name=sample_name)
    logging.info("Mapping and sampling reads")
    sample_outprefix = os.path.join(processing_dir, "sample_reads")
    sampler = sample_reads.sample_reads(
        ref_genome,
        all_reads_bam,
        sample_outprefix,
        amplicon_bed,
        target_depth=target_sample_depth,
    )
    bam = sampler.bam_out
    logging.info(f"Running QC on all mapped reads in {bam}")
    if paired:
        bad_amplicons = qcovid.bin_amplicons(processing_dir, ref_genome, amplicon_bed, bam)
    else:
        bad_amplicons = qcovid.bin_amplicons_se(processing_dir, ref_genome, amplicon_bed, bam)

    logging.info("Making initial unmasked consensus using Viridian")
    viridian_out = os.path.join(processing_dir, "viridian")
    assembly = run_viridian(
        tech, viridian_out, ref_genome, amplicon_json, bam, bad_amplicons
    )

    logging.info("Mapping reads to consensus from Viridian")
    if paired:
        self_map = minimap.run(
            processing_dir,
            assembly,
            sampler.fq_out1,
            sampler.fq_out2,
            prefix="self_qc",
            sample_name=sample_name,
        )
    else:
        self_map = minimap.run_se(
            processing_dir, assembly, sampler.fq_out, prefix="self_qc", sample_name=sample_name
        )

    logging.info("Running QC on Viridian consensus to make masked FASTA")
    masked_fasta = qcovid.self_qc(processing_dir, assembly, self_map)
    final_masked_fasta = os.path.join(outdir, "consensus.fa")
    set_seq_name_in_fasta_file(masked_fasta, final_masked_fasta, sample_name)

    logging.info("Making VCF file of variants")
    varifier_out = os.path.join(processing_dir, "varifier")
    varifier_vcf = varifier.run(
        varifier_out,
        ref_genome,
        final_masked_fasta,
        min_coord=amplicons_start,
        max_coord=amplicons_end,
    )
    check_file(varifier_vcf)
    final_vcf = os.path.join(outdir, "variants.vcf")
    set_sample_name_in_vcf_file(varifier_vcf, final_vcf, sample_name)
    check_file(final_vcf)

    log_info["read_sampling"] = load_json(f"{sample_outprefix}.json")
    log_info["viridian"] = load_json(os.path.join(viridian_out, "run_info.json"))

    # clean up intermediate files
    if not keep_intermediate:
        if keep_bam:
            logging.info(f"Keeping BAM file {all_reads_bam} because --keep_bam option used")
        else:
            rm(all_reads_bam)
            rm(all_reads_bam + ".bai")
        logging.info(f"Removing processing directory {processing_dir}")
        subprocess.check_output(f"rm -rf {processing_dir}", shell=True)

    logging.info(f"Writing JSON log file {json_log}")
    end_time = datetime.datetime.now()
    log_info["run_summary"]["end_time"] = end_time.replace(microsecond=0).isoformat()
    log_info["run_summary"]["run_time"] = str(end_time - start_time)
    log_info["run_summary"]["finished_running"] = True
    with open(json_log, "w") as f:
        json.dump(log_info, f, indent=2)
    logging.info(f"Finished running viridian_workflow, output dir: {outdir}")
