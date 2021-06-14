import subprocess
import os

from viridian_workflow import minimap, qcovid, sample_reads, varifier
from viridian_workflow.utils import check_file, run_process, rm


def run_viridian_illumina(outdir, ref_genome, amplicon_bed, bam, bad_amplicons):
    viridian_cmd = " ".join(
        [
            "viridian",
            "assemble",
            "--bam",
            bam,
            "illumina",
            "--amplicons_to_fail_file",
            bad_amplicons,
            ref_genome,
            amplicon_bed,
            outdir,
        ]
    )
    assembly = os.path.join(outdir, "consensus.final_assembly.fa")
    run_process(viridian_cmd)
    check_file(assembly)
    return assembly


def run_viridian_ont(outdir, ref_genome, amplicon_bed, bam, bad_amplicons):
    viridian_cmd = " ".join(
        [
            "viridian",
            "assemble",
            "--bam",
            bam,
            "ont",
            "--amplicons_to_fail_file",
            bad_amplicons,
            ref_genome,
            amplicon_bed,
            outdir,
        ]
    )
    assembly = os.path.join(outdir, "consensus.final_assembly.fa")
    run_process(viridian_cmd)
    check_file(assembly)
    return assembly


def run_one_sample_ont(outdir, ref_genome, amplicon_bed, fq, keep_intermediate=False, target_sample_depth=1000):
    os.mkdir(outdir)
    all_reads_bam = minimap.run_se(outdir, ref_genome, fq)
    sample_outprefix = os.path.join(outdir, "sample")
    sampler = sample_reads.sample_reads(ref_genome, all_reads_bam, sample_outprefix, amplicon_bed, target_depth=target_sample_depth)
    bam = sampler.bam_out
    bad_amplicons = qcovid.bin_amplicons(outdir, ref_genome, amplicon_bed, bam)

    viridian_out = os.path.join(outdir, "viridian")
    assembly = run_viridian_ont(
        viridian_out, ref_genome, amplicon_bed, bam, bad_amplicons
    )

    varifier_out = os.path.join(outdir, "varifier")
    self_map = minimap.run_se(outdir, assembly, sampler.fq_out, prefix="self_qc")
    masked_fasta = qcovid.self_qc(outdir, assembly, self_map)

    vcf = varifier.run(varifier_out, ref_genome, masked_fasta)
    check_file(vcf)

    # clean up intermediate files
    if not keep_intermediate:
        rm(all_reads_bam)
        rm(all_reads_bam + ".bai")
        rm(self_map)
        rm(self_map + ".bai")
        rm(bam)
        rm(bam + ".bai")
        rm(sampler.fq_out)


def run_one_sample(outdir, ref_genome, amplicon_bed, fq1, fq2, keep_intermediate=False, target_sample_depth=1000):
    os.mkdir(outdir)
    all_reads_bam = minimap.run(outdir, ref_genome, fq1, fq2)
    sample_outprefix = os.path.join(outdir, "sample")
    sampler = sample_reads.sample_reads(ref_genome, all_reads_bam, sample_outprefix, amplicon_bed, target_depth=target_sample_depth)
    bam = sampler.bam_out
    bad_amplicons = qcovid.bin_amplicons(outdir, ref_genome, amplicon_bed, bam)

    viridian_out = os.path.join(outdir, "viridian")
    assembly = run_viridian_illumina(
        viridian_out, ref_genome, amplicon_bed, bam, bad_amplicons
    )

    varifier_out = os.path.join(outdir, "varifier")
    self_map = minimap.run(outdir, assembly, sampler.fq_out1, sampler.fq_out2, prefix="self_qc")
    masked_fasta = qcovid.self_qc(outdir, assembly, self_map)

    vcf = varifier.run(varifier_out, ref_genome, masked_fasta)
    check_file(vcf)

    # clean up intermediate files
    if not keep_intermediate:
        rm(all_reads_bam)
        rm(all_reads_bam + ".bai")
        rm(self_map)
        rm(self_map + ".bai")
        rm(bam)
        rm(bam + ".bai")
        rm(sampler.fq_out1)
        rm(sampler.fq_out2)
