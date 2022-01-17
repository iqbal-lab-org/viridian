import argparse
import datetime
import logging
import socket
import subprocess
import sys
import os

from viridian_workflow import (
    amplicon_schemes,
    detect_primers,
    primers,
    minimap,
    qcovid,
    self_qc,
    sample_reads,
    utils,
    varifier,
)
from viridian_workflow import __version__ as viridian_wf_version


class Pipeline:
    def __init__(
        self,
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
        min_sample_depth=10,
        max_percent_amps_fail=50.0,
        viridian_cons_max_n_percent=50.0,
        command_line_args=None,
    ):
        self.tech = tech
        self.outdir = outdir
        self.ref_genome = ref_genome
        self.fq1 = fq1
        self.fq2 = fq2
        self.built_in_amp_schemes = built_in_amp_schemes
        self.tsv_of_amp_schemes = tsv_of_amp_schemes
        self.force_amp_scheme = force_amp_scheme
        self.keep_intermediate = keep_intermediate
        self.keep_bam = keep_bam
        self.target_sample_depth = target_sample_depth
        self.sample_name = sample_name
        self.max_percent_amps_fail = max_percent_amps_fail
        self.min_sample_depth = min_sample_depth
        self.viridian_cons_max_n_percent = viridian_cons_max_n_percent
        self.command_line_args = command_line_args
        self.start_time = None
        self.amplicon_scheme_name_to_tsv = None
        self.amplicon_scheme_list = None
        self.set_command_line_dict()
        self.check_tech()
        self.set_paired_and_check_reads()
        self.json_log_file = os.path.join(self.outdir, "log.json")
        self.processing_dir = os.path.join(self.outdir, "Processing")
        self.unsorted_read_tagged_bam = os.path.join(
            self.processing_dir, "map_reads.unsorted.read_tagged.bam"
        )
        self.amplicon_json = os.path.join(self.processing_dir, "amplicons.json")
        self.amplicon_bed = os.path.join(self.processing_dir, "amplicons.bed")
        self.all_reads_bam = os.path.join(self.outdir, "reference_mapped.bam")
        self.viridian_outdir = os.path.join(self.processing_dir, "viridian")
        self.viridian_fasta = os.path.join(
            self.viridian_outdir, "consensus.final_assembly.fa"
        )
        self.bam_reads_v_viridian = os.path.join(self.processing_dir, "self_qc.bam")
        self.final_masked_fasta = os.path.join(self.outdir, "consensus.fa")
        self.varifier_vcf = None
        self.final_vcf = os.path.join(self.outdir, "variants.vcf")
        self.sampler = None
        self.amplicons_start = None
        self.amplicons_end = None
        self.sampled_bam = None
        self.amplicon_tsv = None
        self.amplicon_set = None
        self.log_dict = None

    def check_tech(self):
        allowed_tech = {"ont", "illumina"}
        if self.tech not in allowed_tech:
            techs = ",".join(sorted(list(allowed_tech)))
            raise Exception(f"Tech '{self.tech}' not allowed, must be one of: {techs}")

    def get_minimap_presets(self):
        if self.tech == "ont":
            return "map-ont"
        elif self.tech == "illumina":
            return "sr"

    def set_paired_and_check_reads(self):
        if self.tech == "ont":
            assert self.fq2 is None
            self.paired = False
        elif self.tech == "illumina":
            assert self.fq2 is not None
            self.paired = True
        else:
            raise NotImplementedError(f"tech not implemented: {self.tech}")

    def set_command_line_dict(self):
        # Make a dict of the command line options to go in the JSON output file.
        # The tests don't use argparse (they use Mock), which means convert to dict
        # doesn't work. Don't care about that case anyway in the final output, so
        # just set to None
        if isinstance(self.command_line_args, argparse.Namespace):
            self.command_line_dict = {
                k: v for k, v in vars(self.command_line_args).items() if k != "func"
            }
        else:
            self.command_line_dict = None

    def update_json_latest_stage(self, latest_stage):
        self.log_dict["run_summary"]["last_stage_completed"] = latest_stage
        utils.write_json(self.json_log_file, self.log_dict)

    def process_amplicon_schemes(self):
        if self.built_in_amp_schemes is None and self.tsv_of_amp_schemes is None:
            logging.info("No primer schemes provided. Using all built in schemes")
            self.built_in_amp_schemes = list(
                amplicon_schemes.get_built_in_schemes().keys()
            )
        (
            self.amplicon_scheme_name_to_tsv,
            self.amplicon_scheme_list,
        ) = amplicon_schemes.load_list_of_amplicon_sets(
            built_in_names_to_use=self.built_in_amp_schemes,
            tsv_others_to_use=self.tsv_of_amp_schemes,
        )
        if self.force_amp_scheme is not None:
            if self.force_amp_scheme not in self.amplicon_scheme_name_to_tsv:
                names = ",".join(sorted(list(self.amplicon_scheme_name_to_tsv.keys())))
                raise Exception(
                    f"Chose to force amplicons scheme to be {force_amp_scheme}, but scheme not found. Found these: {names}"
                )
        self.update_json_latest_stage("Processed amplicon scheme files")
        logging.info(
            f"Processed amplicon scheme files. Amplicon scheme names: {','.join(sorted(list(self.amplicon_scheme_name_to_tsv.keys())))}"
        )

    def initial_read_map_and_detect_amplicon_scheme(self):
        logging.info("Mapping reads to reference")
        unsorted_bam = os.path.join(self.processing_dir, "map_reads.unsorted.bam")
        minimap.run(
            unsorted_bam,
            self.ref_genome,
            self.fq1,
            fq2=self.fq2,
            sample_name=self.sample_name,
            sort=False,
        )
        self.update_json_latest_stage("Initial map reads")
        logging.info("Detecting amplicon scheme and gathering read statistics")

        primer_stats = detect_primers.gather_stats_from_bam(
            unsorted_bam, self.unsorted_read_tagged_bam, self.amplicon_scheme_list
        )

        self.log_dict["read_and_primer_stats"] = primer_stats

        self.log_dict["read_and_primer_stats"][
            "amplicon_scheme_set_matches"
        ] = detect_primers.amplicon_set_counts_to_json_friendly(
            self.log_dict["read_and_primer_stats"]["amplicon_scheme_set_matches"]
        )
        if self.force_amp_scheme is None:
            self.log_dict["amplicon_scheme_name"] = self.log_dict[
                "read_and_primer_stats"
            ]["chosen_amplicon_scheme"]
        else:
            self.log_dict["chosen_amplicon_scheme"] = self.force_amp_scheme

        chosen_scheme = primer_stats["chosen_amplicon_scheme"]
        self.amplicon_tsv = self.amplicon_scheme_name_to_tsv[chosen_scheme]

        self.amplicon_set = primers.AmpliconSet.from_tsv(
            self.amplicon_tsv, name=chosen_scheme
        )
        self.update_json_latest_stage("Gather read stats and detect primer scheme")

    def process_chosen_amplicon_scheme_files(self):
        scheme_name = self.log_dict["amplicon_scheme_name"]
        logging.info(f"Processing files for chosen amplicon scheme {scheme_name}")
        amplicon_schemes.convert_tsv_to_viridian_json(
            self.amplicon_tsv, self.amplicon_json, scheme_name=scheme_name
        )
        (
            self.amplicons_start,
            self.amplicons_end,
        ) = utils.amplicons_json_to_bed_and_range(self.amplicon_json, self.amplicon_bed)
        self.update_json_latest_stage("Processed chosen amplicon scheme files")

    def sample_the_reads(self):
        logging.info("Sampling reads")
        sample_outprefix = os.path.join(self.processing_dir, "sample_reads")
        self.sampler = sample_reads.sample_reads(
            self.ref_genome,
            self.all_reads_bam,
            sample_outprefix,
            self.amplicon_bed,
            target_depth=self.target_sample_depth,
            min_sampled_depth_for_pass=self.min_sample_depth,
        )
        self.log_dict["read_sampling"] = utils.load_json(f"{sample_outprefix}.json")
        self.sampled_bam = self.sampler.bam_out
        self.update_json_latest_stage("Sampled reads")
        if (
            100 * self.sampler.failed_amplicons / self.sampler.number_of_amplicons()
            >= self.max_percent_amps_fail
        ):
            message = "Too many amplicons are too low depth. STOPPING"
            logging.info(message)
            self.update_json_latest_stage(
                "Sampled reads, too many low coverage amplicons"
            )
            self.finalise_json_log([message])
            return False
        else:
            return True

    def check_viridian(self):
        try:
            self.log_dict["viridian"] = utils.load_json(
                os.path.join(self.viridian_outdir, "run_info.json")
            )
        except:
            return "Error getting viridian JSON file"

        try:
            consensus = self.log_dict["viridian"]["run_summary"]["consensus"]
            percent_n = 100.0 * consensus.count("N") / len(consensus)
            if percent_n > self.viridian_cons_max_n_percent:
                return f"Too many Ns in Viridian consensus: {round(percent_n, 2)}%"
        except:
            return "Error getting viridian consensus sequence and/or counting Ns"

        if not os.path.exists(self.viridian_fasta):
            return "No FASTA file made by Viridian"

        try:
            total_amps = self.log_dict["viridian"]["run_summary"]["total_amplicons"]
            bad_amps = (
                total_amps
                - self.log_dict["viridian"]["run_summary"]["successful_amplicons"]
            )
            if 100 * bad_amps / total_amps >= self.max_percent_amps_fail:
                return "Too many failed amplicons from Viridian"
        except:
            pass

    def run_viridian(self):
        logging.info("Making initial unmasked consensus using Viridian")
        viridian_cmd = " ".join(
            [
                "viridian",
                "assemble",
                "--bam",
                self.sampled_bam,
                self.tech,
                "--amplicons_to_fail_file",
                self.sampler.failed_amps_file,
                self.ref_genome,
                self.amplicon_json,
                self.viridian_outdir,
            ]
        )
        utils.run_process(viridian_cmd)
        message = self.check_viridian()

        if message is not None:
            self.update_json_latest_stage("Viridian. " + message)
            message += ". STOPPING"
            logging.info(message)
            self.finalise_json_log([message])
            return False

        self.update_json_latest_stage("Viridian")
        return True

    def map_reads_to_viridian_consensus(self):
        logging.info("Mapping reads to consensus from Viridian")
        if self.paired:
            fq1 = self.sampler.fq_out1
            fq2 = self.sampler.fq_out2
        else:
            fq1 = self.sampler.fq_out
            fq2 = None
        minimap.run(
            self.bam_reads_v_viridian,
            self.viridian_fasta,
            fq1,
            fq2=fq2,
            sample_name=self.sample_name,
            sort=True,
        )
        self.update_json_latest_stage("Map reads to Viridian consensus")

    def self_qc_and_make_masked_consensus(
        self, minimap_presets, amplicon_set, tagged_reads
    ):
        logging.info("Running QC on Viridian consensus to make masked FASTA")

        position_stats = self_qc.remap(
            self.viridian_fasta, self.get_minimap_presets(), amplicon_set, tagged_reads
        )
        masked_fasta, masking_log = self_qc.mask(
            self.viridian_fasta,
            position_stats,
            outpath=self.final_masked_fasta,
            name=self.sample_name,
        )

        self.log_dict["self_qc"] = masking_log
        self.update_json_latest_stage(
            "Ran QC on reads mapped to consensus and made final masked FASTA"
        )

    def make_vcf_wrt_reference(self):
        logging.info("Making VCF file of variants")
        varifier_out = os.path.join(self.processing_dir, "varifier")
        print(varifier_out)
        self.varifier_vcf = varifier.run(
            varifier_out,
            self.ref_genome,
            self.final_masked_fasta,
            min_coord=self.amplicons_start,
            max_coord=self.amplicons_end,
        )
        utils.check_file(self.varifier_vcf)
        utils.set_sample_name_in_vcf_file(
            self.varifier_vcf, self.final_vcf, self.sample_name
        )
        utils.check_file(self.final_vcf)
        self.update_json_latest_stage("Made VCF of variants")

    def clean_intermediate_files(self):
        if not self.keep_intermediate:
            logging.info("Deleting temporary files")
            if self.keep_bam:
                logging.info(
                    f"Keeping BAM file {self.all_reads_bam} because --keep_bam option used"
                )
            else:
                bai = f"{self.all_reads_bam}.bai"
                subprocess.check_output(f"rm -f {self.all_reads_bam} {bai}", shell=True)
            logging.info(f"Removing processing directory {self.processing_dir}")
            subprocess.check_output(f"rm -rf {self.processing_dir}", shell=True)
        else:
            logging.info("Debug mode: not deleting temporary files")

    def finalise_json_log(self, result):
        logging.info(f"Writing JSON log file {self.json_log_file}")
        end_time = datetime.datetime.now()
        self.log_dict["run_summary"]["end_time"] = end_time.replace(
            microsecond=0
        ).isoformat()
        self.log_dict["run_summary"]["run_time"] = str(end_time - self.start_time)
        self.log_dict["run_summary"]["finished_running"] = True
        self.log_dict["run_summary"]["result"] = result
        self.update_json_latest_stage("Finished")
        logging.info(f"Finished running viridian_workflow")

    def run(self):
        logging.info(f"Start running viridian_workflow, output dir: {self.outdir}")
        os.mkdir(self.outdir)
        self.start_time = datetime.datetime.now()
        self.log_dict = {
            "run_summary": {
                "last_stage_completed": "start",
                "command": " ".join(sys.argv),
                "options": self.command_line_dict,
                "cwd": os.getcwd(),
                "version": viridian_wf_version,
                "finished_running": False,
                "start_time": self.start_time.replace(microsecond=0).isoformat(),
                "end_time": None,
                "hostname": socket.gethostname(),
                "result": "Unknown",
            },
            "read_and_primer_stats": None,
            "read_sampling": None,
            "viridian": None,
            "self_qc": None,
        }
        self.update_json_latest_stage("start")
        self.process_amplicon_schemes()
        os.mkdir(self.processing_dir)
        self.initial_read_map_and_detect_amplicon_scheme()
        self.process_chosen_amplicon_scheme_files()

        logging.info("Sorting and indexing BAM file of all reads mapped to reference")
        utils.run_process(
            f"samtools sort -O BAM -o {self.all_reads_bam} {self.unsorted_read_tagged_bam}"
        )
        utils.run_process(f"samtools index {self.all_reads_bam}")
        self.update_json_latest_stage("Sorted and indexed all reads BAM file")
        if not self.sample_the_reads():
            self.clean_intermediate_files()
            return

        if not self.run_viridian():
            self.clean_intermediate_files()
            return

        self.self_qc_and_make_masked_consensus(
            self.viridian_fasta, self.amplicon_set, self.unsorted_read_tagged_bam
        )
        self.make_vcf_wrt_reference()
        self.clean_intermediate_files()
        self.finalise_json_log("Success")


def run_one_sample(
    tech, outdir, ref_genome, fq1, **kw,
):
    pipeline = Pipeline(tech, outdir, ref_genome, fq1, **kw)
    pipeline.run()
