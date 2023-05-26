import argparse
from collections import Counter
import datetime
import logging
import os
import shutil
import socket
import sys
import tempfile
import traceback

from viridian import (
    amplicon_schemes,
    cylon,
    constants,
    ena,
    maptools,
    qc,
    reads,
    read_it_and_keep,
    scheme_id,
    utils,
    varifier,
)
from viridian import __version__ as viridian_version


class Pipeline:
    def __init__(
        self,
        tech,
        outdir,
        ref_fasta,
        ena_run=None,
        keep_ena_reads=False,
        reads_file1=None,
        reads_file2=None,
        reads_file=None,
        reads_bam=None,
        decontam_ref_fa=None,
        built_in_amp_schemes=None,
        tsv_of_amp_schemes=None,
        force_amp_scheme=None,
        detect_scheme_only=False,
        force_consensus=None,
        debug=False,
        keep_bam=False,
        msas_to_write=None,
        qc_depth=1000,
        min_scheme_score=250,
        max_scheme_ratio=0.5,
        sample_name="sample",
        max_percent_amps_fail=50.0,
        max_cons_n_percent=50.0,
        coverage_min_x=20,
        coverage_min_pc=50.0,
        masking_min_frs=None,
        masking_min_depth=10,
        het_min_pc=10.0,
        assemble_depth=100,
        command_line_args=None,
        temp_root=None,
        fix_small_indels=True,
        gzip_files=True,
    ):
        self.tech = tech
        self.outdir = os.path.abspath(outdir)
        self.ena_run = ena_run
        self.ena_reads_dir = os.path.join(self.outdir, "ENA_download")
        self.keep_ena_reads = keep_ena_reads
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.all_reads_bam = os.path.join(self.outdir, "reference_mapped.bam")
        self.reads_bam = reads_bam
        self.reads_file = reads_file
        self.reads_file1 = reads_file1
        self.reads_file2 = reads_file2
        self.decontam_ref_fa = decontam_ref_fa
        self.built_in_amp_schemes = built_in_amp_schemes
        self.tsv_of_amp_schemes = tsv_of_amp_schemes
        self.force_amp_scheme = force_amp_scheme
        self.detect_scheme_only = detect_scheme_only
        self.debug = debug
        self.keep_bam = keep_bam
        self.msas_to_write = set() if msas_to_write is None else set(msas_to_write)
        self.qc_depth = qc_depth
        self.min_scheme_score = min_scheme_score
        self.max_scheme_ratio = max_scheme_ratio
        self.sample_name = sample_name
        self.max_percent_amps_fail = max_percent_amps_fail
        self.max_cons_n_percent = max_cons_n_percent
        self.coverage_min_x = coverage_min_x
        self.coverage_min_pc = coverage_min_pc
        self.masking_min_frs = masking_min_frs
        self.masking_min_depth = masking_min_depth
        self.het_min_pc = het_min_pc
        self.command_line_args = command_line_args
        self.set_command_line_dict()
        self.cylon_depth = assemble_depth
        self.primer_end_tolerance = constants.PRIMER_END_TOLERANCE

        self.json_log_file = os.path.join(self.outdir, "log.json")
        self.final_unmasked_fasta = os.path.join(self.outdir, "consensus.unmasked.fa")
        self.final_masked_fasta = os.path.join(self.outdir, "consensus.fa")
        self.final_vcf = os.path.join(self.outdir, "variants.vcf")
        self.qc_tsv = os.path.join(self.outdir, "qc.tsv.gz")
        self.log_dict = None
        self.force_consensus = force_consensus
        self.temp_root = temp_root
        self.fix_small_indels = fix_small_indels
        self.gzip_files = gzip_files
        if self.gzip_files:
            self.json_log_file += ".gz"
            self.final_masked_fasta += ".gz"
        self.minimap_x_opt = constants.TECH2MINIMAP_X[self.tech]


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

    def init_json_log(self):
        self.log_dict = {
            "run_summary": {
                "last_stage_completed": "start",
                "command": " ".join(sys.argv),
                "options": self.command_line_dict,
                "cwd": os.getcwd(),
                "version": viridian_version,
                "finished_running": False,
                "start_time": self.start_time.replace(microsecond=0).isoformat(),
                "end_time": None,
                "hostname": socket.gethostname(),
                "result": "Unknown",
                "errors": [],
                "temp_processing_dir": self.processing_dir,
            },
            "stages_completed": [],
            "amplicons": None,
            "scheme_choice": None,
            "self_qc": None,
            "debug": {},
        }

    def update_json_latest_stage(self, latest_stage):
        self.log_dict["stages_completed"].append(latest_stage)
        self.log_dict["run_summary"]["last_stage_completed"] = latest_stage
        utils.write_json(self.json_log_file, self.log_dict)

    def add_errors_to_log(self, errors):
        if self.log_dict is None:
            return
        if errors is not None:
            if isinstance(errors, list):
                self.log_dict["run_summary"]["errors"].extend(errors)
                logging.error("\n".join(errors))
            elif isinstance(errors, str):
                self.log_dict["run_summary"]["errors"].append(errors)
                logging.error(errors)
            else:
                raise NotImplementedError()

    def finalise_json_log(self, result, errors=None):
        if self.log_dict is None:
            return
        if "debug" in self.log_dict and not self.debug:
            del self.log_dict["debug"]
        logging.info(f"Writing JSON log file {self.json_log_file}")
        end_time = datetime.datetime.now(datetime.timezone.utc)
        self.log_dict["run_summary"]["end_time"] = end_time.replace(
            microsecond=0
        ).isoformat()
        self.log_dict["run_summary"]["run_time"] = str(end_time - self.start_time)
        self.log_dict["run_summary"]["finished_running"] = True
        self.log_dict["run_summary"]["result"] = result
        self.add_errors_to_log(errors)
        self.update_json_latest_stage("Finished")
        logging.info(f"Finished running viridian. Result: {result}")

    def setup_output_dirs(self):
        # tempfile tries the usual places for temp file (like /tmp), but
        # if all else fails it uses the current working directory. We Don't
        # want to use cwd, but instead make sure it's in the output directory
        if self.debug or tempfile.gettempdir() == os.getcwd():
            self.processing_dir = os.path.join(self.outdir, "Processing")
            os.mkdir(self.processing_dir)
        else:
            self.tempdir = tempfile.TemporaryDirectory(
                prefix="viridian.", dir=self.temp_root
            )
            self.processing_dir = self.tempdir.name

        logging.info(f"Putting temporary files in {self.processing_dir}")
        self.scheme_id_dir = os.path.join(self.processing_dir, "scheme_id")
        self.sample_reads_dir = os.path.join(self.processing_dir, "sample_reads")
        self.cylon_outdir = os.path.join(self.processing_dir, "cylon")
        if self.force_consensus is None:
            self.force_consensus = False
            self.cylon_fasta = os.path.join(
                self.cylon_outdir, "consensus.final_assembly.fa"
            )
        else:
            if not os.path.exists(self.force_consensus):
                raise FileNotFoundError(
                    f"Forced consensus file not found: {self.force_consensus}"
                )
            self.cylon_fasta = os.path.abspath(self.force_consensus)
            self.force_consensus = True
        self.qc_bams_dir = os.path.join(self.processing_dir, "qc_bams")
        self.varifier_dir = os.path.join(self.processing_dir, "varifier")

    def setup_reads_files(self):
        if self.ena_run is not None:
            got_reads, self.tech, ena_metadata = ena.download_run(
                self.ena_run, self.ena_reads_dir
            )
            logging.debug(f"ENA metadata: {ena_metadata}")
            self.log_dict["ena_metadata"] = ena_metadata
            if "reads1" in got_reads:
                self.reads_file1 = os.path.abspath(got_reads["reads1"])
                self.reads_file2 = os.path.abspath(got_reads["reads2"])
                self.paired = True
            else:
                self.reads_file1 = os.path.abspath(got_reads["unpaired_reads"])
                self.reads_file2 = None
                self.paired = False
        elif self.reads_bam is not None:
            self.reads_bam = os.path.abspath(self.reads_bam)
            self.all_reads_bam = self.reads_bam
            assert all(
                x is None for x in [self.reads_file, self.reads_file1, self.reads_file2]
            )
            self.paired = None
            assert self.decontam_ref_fa is None
        elif self.reads_file is not None:
            self.reads_file1 = os.path.abspath(self.reads_file)
            self.paired = False
            assert self.reads_file2 is None
        else:
            self.reads_file1 = os.path.abspath(self.reads_file1)
            self.reads_file2 = os.path.abspath(self.reads_file2)
            self.paired = True

    def start_pipeline(self):
        self.start_time = datetime.datetime.now(datetime.timezone.utc)
        logging.info(f"Start running viridian, output dir: {self.outdir}")
        try:
            os.mkdir(self.outdir)
        except:
            logging.error("Error making output directory. Cannot continue")
            raise Exception("Error making output directory. Cannot continue")
        self.setup_output_dirs()
        self.init_json_log()
        self.setup_reads_files()
        files_to_check = {
            "reference fasta": self.ref_fasta,
        }
        if self.decontam_ref_fa is not None:
            files_to_check["Decontamination reference fasta"] = self.decontam_ref_fa

        if self.reads_bam is not None:
            files_to_check["reads_bam"] = self.reads_bam
            files_to_check["reads_bam.bai"] = self.reads_bam + ".bai"
        elif self.paired:
            files_to_check["reads1"] = self.reads_file1
            files_to_check["reads2"] = self.reads_file2
        else:
            files_to_check["reads"] = self.reads_file1
        not_found = []
        for nice_name, filename in sorted(files_to_check.items()):
            found = os.path.exists(filename)
            if not found:
                not_found.append(f"{nice_name}: {filename}")
            logging.info(f"Check file exists {nice_name} {filename}: {found}")

        if len(not_found):
            not_found = "\n".join(not_found)
            raise FileNotFoundError(f"File(s) not found:\n{not_found}")

        return True

    def decontam_reads(self):
        logging.info("Decontaminate reads with readItAndKeep")
        outprefix = os.path.join(self.outdir, "decontaminate")
        self.log_dict["decontaminate"], files, error = read_it_and_keep.run_riak(
            outprefix, self.decontam_ref_fa, self.reads_file1, reads2=self.reads_file2
        )
        if error is not None:
            raise Exception(f"Error running readItAndKeep:\n{error}")
        self.reads_file1 = os.path.abspath(files["reads_1"])
        if self.paired:
            self.reads_file2 = os.path.abspath(files["reads_2"])
        return True

    def process_amplicon_schemes(self):
        if self.built_in_amp_schemes is None and self.tsv_of_amp_schemes is None:
            logging.info("No primer schemes provided. Using all built in schemes")
            self.built_in_amp_schemes = list(
                amplicon_schemes.get_built_in_schemes().keys()
            )
        self.amplicon_scheme_name_to_tsv = amplicon_schemes.load_list_of_amplicon_sets(
            built_in_names_to_use=self.built_in_amp_schemes,
            tsv_others_to_use=self.tsv_of_amp_schemes,
        )
        if self.force_amp_scheme is not None:
            if self.force_amp_scheme not in self.amplicon_scheme_name_to_tsv:
                names = ",".join(sorted(list(self.amplicon_scheme_name_to_tsv.keys())))
                raise Exception(
                    f"Chose to force amplicons scheme to be {self.force_amp_scheme}, but scheme not found. Found these: {names}"
                )
        logging.info(
            f"Processed amplicon scheme files. Amplicon scheme names: {','.join(sorted(list(self.amplicon_scheme_name_to_tsv.keys())))}"
        )
        return True

    def initial_read_map(self):
        logging.info("Mapping reads to reference")
        maptools.map_reads(
            self.all_reads_bam,
            self.ref_fasta,
            self.reads_file1,
            reads2=self.reads_file2,
            debug=self.debug,
            sample_name=self.sample_name,
            minimap_x_opt=self.minimap_x_opt,
        )
        return True

    def detect_amplicon_scheme(self):
        logging.info("Comparing mapped reads against amplicon schemes")
        id_results, error_message = scheme_id.analyse_bam(
            self.all_reads_bam,
            self.amplicon_scheme_name_to_tsv,
            self.scheme_id_dir,
            end_tolerance=self.primer_end_tolerance,
            sample_name=self.sample_name,
            debug=self.debug,
            min_depth_cutoff=self.coverage_min_x,
            min_percent_genome_cutoff=self.coverage_min_pc,
            max_primer_dist=constants.SCHEME_ID_PRIMER_WITHIN_END,
        )
        self.log_dict["reads"] = id_results["read_counts"]
        self.log_dict["amplicons"] = id_results["amplicons"]
        self.log_dict["scheme_choice"] = id_results["scheme_choice"]
        self.log_dict["read_depth"] = id_results["stats"]
        self.log_dict["debug"]["assembly_amplicons"] = id_results["cylon_amplicons"]
        try:
            self.log_dict["debug"]["read_depth"] = {
                "depth_per_position": id_results["stats"]["depth_per_position"],
                "depth_hist": id_results["stats"]["depth_hist"],
            }
            del self.log_dict["read_depth"]["depth_per_position"]
            del self.log_dict["read_depth"]["depth_hist"]
        except:
            pass
        if error_message is not None:
            self.add_errors_to_log(error_message)
            return False
        best_scheme = id_results["scheme_choice"]["best_scheme"]
        best_score = id_results["scheme_choice"]["best_score"]
        score_ratio = id_results["scheme_choice"]["score_ratio"]
        logging.info(
            f"Amplicon scheme that reads best match to: {best_scheme}, with score {best_score}, and (second best)/best = {score_ratio}"
        )

        if self.force_amp_scheme is None:
            score_ok = True
            if best_score < self.min_scheme_score:
                self.log_dict["amplicon_scheme_name"] = None
                self.add_errors_to_log(
                    f"Best scheme score is {best_score}, which is less than required minimum score {self.min_scheme_score}"
                )
                score_ok = False
            if score_ratio is None or score_ratio > self.max_scheme_ratio:
                self.log_dict["amplicon_scheme_name"] = None
                self.add_errors_to_log(
                    f"(second best scheme score) / (best score) is {score_ratio}, which is more than the required maximum {self.max_scheme_ratio}"
                )
                score_ok = False
            if not score_ok:
                return False
            self.log_dict["amplicon_scheme_name"] = best_scheme
            logging.info(f"Using amplicon scheme {best_scheme}")
        else:
            self.log_dict["amplicon_scheme_name"] = self.force_amp_scheme
            logging.info(
                f"Using amplicon scheme {self.force_amp_scheme} (forced by user, so may be different from the best matching scheme inferred from read mapping)"
            )

        files_to_copy = ["score_plot.pdf", "depth_across_genome.pdf"]
        for filename in files_to_copy:
            old = os.path.join(self.scheme_id_dir, filename)
            new = os.path.join(self.outdir, f"scheme_id.{filename}")
            shutil.copyfile(old, new)

        return True

    def sample_the_reads(self):
        logging.info("Sampling reads")
        self.read_sampler = reads.ReadSampler(
            self.log_dict["amplicons"],
            qc_depth=self.qc_depth,
            cylon_depth=self.cylon_depth,
            adapter_trim_tolerance=self.primer_end_tolerance,
        )
        self.read_sampler.sample_reads(self.all_reads_bam, self.sample_reads_dir)
        self.log_dict["reads"] = self.read_sampler.read_counts
        self.log_dict["debug"]["sample_reads"] = self.read_sampler.sample_stats
        logging.info("Finished sampling reads")
        return True

    def run_cylon(self):
        if self.force_consensus:
            message = f"Skip making initial consensus sequence because --force_consensus used with file {self.cylon_fasta}"
            self.log_dict["debug"]["assembly"] = message
            logging.info(message)
            return True

        logging.info("Making initial consensus sequence")
        reads_dir = os.path.join(self.sample_reads_dir, "cylon")
        amp_json = os.path.join(self.sample_reads_dir, "cylon.json")
        cylon_dict = {"amplicons": self.log_dict["debug"]["assembly_amplicons"]}
        utils.write_json(amp_json, cylon_dict)
        self.log_dict["debug"]["assembly"], error_message = cylon.run_cylon(
            reads_dir,
            constants.CYLON_TECH[self.tech],
            self.ref_fasta,
            amp_json,
            self.cylon_outdir,
            debug=self.debug,
            max_percent_n=self.max_cons_n_percent,
            max_percent_amps_fail=self.max_percent_amps_fail,
        )
        if error_message is not None:
            logging.warning(
                f"Error making initial consensus sequence. Cannot continue. Error: {error_message}"
            )
            self.add_errors_to_log(error_message)
            return False

        logging.info("Finished making initial consensus sequence")
        return True

    def run_varifier(self):
        logging.info("Making initial VCF file and multiple sequence alignment")
        cylon_amps = self.log_dict["debug"]["assembly_amplicons"]
        if (
            self.fix_small_indels
            and self.tech in constants.INDEL_FIX_LENGTH
            and not self.force_consensus
        ):
            indel_fix_length = constants.INDEL_FIX_LENGTH[self.tech]
        else:
            indel_fix_length = None

        (
            self.log_dict["sequences"],
            self.vcf_header,
            self.vcf_records,
            error_message,
        ) = varifier.run_varifier(
            self.cylon_fasta,
            self.ref_fasta,
            self.varifier_dir,
            min([x["start"] for x in cylon_amps.values()]) + 1,
            max([x["end"] for x in cylon_amps.values()]) + 1,
            sanitise_gaps=not self.force_consensus,
            indel_fix_length=indel_fix_length,
            debug=self.debug,
        )
        if error_message is not None:
            logging.warning(
                f"Error making MSA and/or initial VCF file. Cannot continue. Error: {error_message}"
            )
            self.add_errors_to_log(error_message)
            return False

        try:
            utils.write_fasta(
                f"{self.sample_name}.unmasked",
                self.log_dict["sequences"]["unmasked_consensus"],
                self.final_unmasked_fasta,
            )
        except:
            message = (
                f"Error making unmasked consensus sequence {self.final_unmasked_fasta}"
            )
            logging.warning(message)
            self.add_errors_to_log(message)
            return False

        logging.info("Finished initial VCF file and multiple sequence alignment")
        return True

    def self_qc(self):
        logging.info("Start QC using reads mapped to consensus")
        self.pileups = self.read_sampler.pileups(
            self.final_unmasked_fasta,
            self.qc_bams_dir,
            minimap_x_opt=self.minimap_x_opt,
            debug=self.debug,
        )
        if len(self.pileups) == 0:
            self.add_errors_to_log("No pileup data generated")
            return False

        if self.masking_min_frs is None:
            self.masking_min_frs = constants.TECH2FRS[self.tech]

        q = qc.Qc(
            self.log_dict["amplicons"],
            self.log_dict["sequences"]["msa_ref"],
            self.log_dict["sequences"]["msa_unmasked_consensus"],
            max_amp_n_percent=50,
            mask_min_depth=self.masking_min_depth,
            mask_min_frs=self.masking_min_frs,
            het_min_pc=self.het_min_pc,
        )

        for i, d in enumerate(self.log_dict["amplicons"]):
            d["dropped"] = i in q.dropped_amplicons
        total_amps = len(self.log_dict["amplicons"])
        succ_amps = total_amps - len(q.dropped_amplicons)
        self.log_dict["run_summary"]["total_amplicons"] = total_amps
        self.log_dict["run_summary"]["successful_amplicons"] = succ_amps
        if 100 * len(q.dropped_amplicons) / total_amps > self.max_percent_amps_fail:
            self.add_errors_to_log(
                f"Too many failed amplicons ({len(q.dropped_amplicons)}/{total_amps} failed)"
            )
            return False

        q.make_pileup(self.pileups)
        logging.info("Making per-position stats and masked consensus sequence")
        q.make_tsv_lines_and_masked_cons()
        q.write_qc_tsv_and_make_masked_cons_msa(self.qc_tsv)
        self.log_dict["self_qc"] = {
            "mask_counts": q.mask_counts,
            "masked_positions": q.masked_positions,
            "insertions": q.insertions,
        }
        logging.info("Making VCF file")
        q.annotated_vcf_file(
            self.vcf_header,
            self.vcf_records,
            self.final_vcf,
            sample_name=self.sample_name,
        )
        self.log_dict["sequences"]["masked_consensus"] = q.masked_cons
        self.log_dict["sequences"]["masked_consensus_msa"] = q.masked_cons_msa
        self.log_dict["sequences"][
            "masked_consensus_msa_indel_as_N"
        ] = q.masked_cons_msa_indel_as_N
        self.log_dict["sequences"][
            "masked_consensus_msa_indel_as_ref"
        ] = q.masked_cons_msa_indel_as_ref
        logging.info("Writing masked consensus FASTA")
        utils.write_fasta(
            f"{self.sample_name}.masked", q.masked_cons, self.final_masked_fasta
        )
        for msa_type in self.msas_to_write:
            logging.info(f"Writing {msa_type} MSA FASTA")
            utils.write_fasta(
                f"{self.sample_name}.masked_msa.{msa_type}",
                self.log_dict["sequences"][f"masked_consensus_msa_{msa_type}"],
                os.path.join(self.outdir, f"msa.{msa_type}.fa"),
            )

        logging.info("Finished QC using reads mapped to consensus")
        return True

    def final_qc_checks(self):
        masked_cons = self.log_dict["sequences"]["masked_consensus"]
        self.log_dict["run_summary"]["consensus_length"] = len(masked_cons)
        base_counts = Counter(masked_cons)
        N_count = base_counts["N"]
        acgt_count = sum(base_counts[x] for x in ["A", "C", "G", "T"])
        het_count = len(masked_cons) - N_count - acgt_count
        if len(masked_cons) == 0:
            percent_N = 0
            percent_acgt = 0
            percent_het = 0
        else:
            percent_N = round(100 * N_count / len(masked_cons), 2)
            percent_acgt = round(100 * acgt_count / len(masked_cons), 2)
            percent_het = round(100 * het_count / len(masked_cons), 2)
        self.log_dict["run_summary"]["consensus_N_count"] = N_count
        self.log_dict["run_summary"]["consensus_N_percent"] = percent_N
        self.log_dict["run_summary"]["consensus_ACGT_count"] = acgt_count
        self.log_dict["run_summary"]["consensus_ACGT_percent"] = percent_acgt
        self.log_dict["run_summary"]["consensus_het_count"] = het_count
        self.log_dict["run_summary"]["consensus_het_percent"] = percent_het

        logging.info(
            f"{percent_N}% ({N_count}/{len(masked_cons)}) of the consensus sequence is Ns"
        )
        if percent_N > self.max_cons_n_percent:
            self.add_errors_to_log(
                f"Too many Ns in final consensus: {percent_N}% ({N_count}/{len(masked_cons)})"
            )
            return False
        return True

    def final_tidy(self):
        logging.info("Tidying up intermediate files and contents of log")
        if "sample_reads" in self.log_dict["debug"]:
            sample_dict = self.log_dict["debug"]["sample_reads"]
            for amp in self.log_dict["amplicons"]:
                amp["reads"] = sample_dict.get(amp["name"], None)

        if not self.debug:
            utils.syscall(f"rm -rf {self.processing_dir} {self.final_unmasked_fasta}*")
            if not (self.keep_bam or self.reads_bam is not None):
                utils.syscall(f"rm {self.all_reads_bam}*")
            if self.ena_run is not None and not self.keep_ena_reads:
                utils.syscall(f"rm -r {self.ena_reads_dir}")
            if "debug" in self.log_dict:
                del self.log_dict["debug"]

        logging.info("Finished tidying")
        return True

    def run(self):
        stages = [(self.start_pipeline, "Start pipeline")]

        if self.decontam_ref_fa is not None:
            stages.append((self.decontam_reads, "Decontaminate reads"))

        stages.append((self.process_amplicon_schemes, "Process amplicon scheme files"))

        if self.reads_bam is None:
            stages.append((self.initial_read_map, "Map reads to reference"))

        stages.append((self.detect_amplicon_scheme, "Detect amplicon scheme"))

        if not self.detect_scheme_only:
            stages.extend(
                [
                    (self.sample_the_reads, "Sample reads"),
                    (self.run_cylon, "Initial consensus sequence"),
                    (self.run_varifier, "Initial VCF and MSA of consensus/reference"),
                    (self.self_qc, "QC using reads vs consensus sequence"),
                    (self.final_qc_checks, "Final QC checks"),
                ]
            )

        stages.append((self.final_tidy, "Tidy up final files and log"))

        for stage_number, (function, description) in enumerate(stages):
            stage_start_time = datetime.datetime.now()
            description = f"{stage_number + 1}/{len(stages)} " + description
            message = f" {description} "
            logging.info(f"{message:=^60}")
            function_ok = False
            traceback_lines = None
            error = None
            try:
                function_ok = function()
            except Exception as e:
                traceback_lines = traceback.format_exc().split("\n")
                error = e

            if traceback_lines is not None or not function_ok:
                logging.error("Stopping pipeline")
                self.add_errors_to_log(f"Error during stage: {description}")
                try:
                    self.final_tidy()
                except:
                    pass
                self.finalise_json_log("Fail", errors=traceback_lines)
                if traceback_lines is not None:
                    raise error
                return

            stage_end_time = datetime.datetime.now()
            stage_seconds = round(
                (stage_end_time - stage_start_time).total_seconds(), 1
            )
            self.update_json_latest_stage(f"{description} ({stage_seconds}s)")

        self.finalise_json_log("Success")


def run_one_sample(
    tech,
    outdir,
    ref_genome,
    **kw,
):
    pipeline = Pipeline(tech, outdir, ref_genome, **kw)
    pipeline.run()
