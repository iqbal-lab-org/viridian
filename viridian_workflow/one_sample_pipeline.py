import subprocess
import os

from viridian_workflow import minimap, qcovid, sample_reads, varifier
from viridian_workflow.utils import (
    amplicons_json_to_bed_and_range,
    check_file,
    run_process,
    rm,
    set_sample_name_in_vcf_file,
)


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

    os.mkdir(outdir)
    amplicon_bed = os.path.join(outdir, "amplicons.bed")
    amplicons_start, amplicons_end = amplicons_json_to_bed_and_range(
        amplicon_json, amplicon_bed
    )
    if paired:
        all_reads_bam = minimap.run(
            outdir, ref_genome, fq1, fq2, sample_name=sample_name
        )
    else:
        all_reads_bam = minimap.run_se(outdir, ref_genome, fq1, sample_name=sample_name)
    sample_outprefix = os.path.join(outdir, "sample")
    sampler = sample_reads.sample_reads(
        ref_genome,
        all_reads_bam,
        sample_outprefix,
        amplicon_bed,
        target_depth=target_sample_depth,
    )
    bam = sampler.bam_out
    if paired:
        bad_amplicons = qcovid.bin_amplicons(outdir, ref_genome, amplicon_bed, bam)
    else:
        bad_amplicons = qcovid.bin_amplicons_se(outdir, ref_genome, amplicon_bed, bam)

    viridian_out = os.path.join(outdir, "viridian")
    assembly = run_viridian(
        tech, viridian_out, ref_genome, amplicon_json, bam, bad_amplicons
    )

    varifier_out = os.path.join(outdir, "varifier")
    if paired:
        self_map = minimap.run(
            outdir,
            assembly,
            sampler.fq_out1,
            sampler.fq_out2,
            prefix="self_qc",
            sample_name=sample_name,
        )
    else:
        self_map = minimap.run_se(
            outdir, assembly, sampler.fq_out, prefix="self_qc", sample_name=sample_name
        )

    masked_fasta = qcovid.self_qc(outdir, assembly, self_map)

    varifier_vcf = varifier.run(
        varifier_out,
        ref_genome,
        masked_fasta,
        min_coord=amplicons_start,
        max_coord=amplicons_end,
    )
    check_file(varifier_vcf)
    final_vcf = os.path.join(outdir, "variants.vcf")
    set_sample_name_in_vcf_file(varifier_vcf, final_vcf, sample_name)
    check_file(final_vcf)

    # clean up intermediate files
    if not keep_intermediate:
        if not keep_bam:
            rm(all_reads_bam)
            rm(all_reads_bam + ".bai")
        rm(amplicon_bed)
        rm(self_map)
        rm(self_map + ".bai")
        rm(bam)
        rm(bam + ".bai")
        if paired:
            rm(sampler.fq_out1)
            rm(sampler.fq_out2)
        else:
            rm(sampler.fq_out)
