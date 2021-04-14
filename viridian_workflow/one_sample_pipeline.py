import subprocess
import os

from viridian_workflow import minimap, qcovid, varifier
from viridian_workflow.utils import check_file, run_process, rm


def run_viridian(outdir, ref_genome, amplicon_bed, bam, bad_amplicons):
    viridian_cmd = " ".join(
        [
            "viridian",
            "assemble",
            "--min_read_length",
            "50",
            "--bam",
            bam,
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


def run_one_sample(outdir, ref_genome, amplicon_bed, fq1, fq2, keep_intermediate=False):
    os.mkdir(outdir)
    bam = minimap.run(outdir, ref_genome, fq1, fq2)
    bad_amplicons = qcovid.bin_amplicons(outdir, ref_genome, amplicon_bed, bam)

    viridian_out = os.path.join(outdir, "viridian")
    assembly = run_viridian(viridian_out, ref_genome, amplicon_bed, bam, bad_amplicons)

    varifier_out = os.path.join(outdir, "varifier")
    self_map = minimap.run(outdir, assembly, fq1, fq2, prefix="self_qc")
    masked_fasta = qcovid.self_qc(outdir, assembly, self_map)

    vcf = varifier.run(varifier_out, ref_genome, masked_fasta)
    check_file(vcf)

    # clean up intermediate files
    if not keep_intermediate:
        rm(bam)
        rm(bam + ".bai")
        rm(self_map)
        rm(self_map + ".bai")
