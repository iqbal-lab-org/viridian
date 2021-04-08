import subprocess
import os

from viridian_workflow import minimap, qcovid

def run_viridian(outdir, ref_genome, amplicon_bed, bam, bad_amplicons):
    viridian_cmd = ['viridian', 'assemble', '--min_read_length', '50', '--bam', bam, '--amplicons_to_fail_file', bad_amplicons, ref_genome, amplicon_bed, outdir]
    assembly = os.path.join(outdir, 'assembly.fa')
    subprocess.Popen(viridian_cmd).wait()
    assert os.path.isfile(assembly)
    return assembly

def run_one_sample(outdir, ref_genome, amplicon_bed, fq1, fq2):
    bam = minimap.run(outdir, ref_genome, fq1, fq2)
    bad_amplicons = qcovid.bin_amplicons(outdir, ref_genome, amplicon_bed, bam)

    viridian_out = os.path.join(outdir, 'viridian')
    assembly = run_viridian(viridian_out, ref_genome, amplicon_bed, bam, bad_amplicons)

    self_map = minimap.run(outdir, assembly, fq1, fq2)
    #qcovid.self_qc(outdir, assembly, self_mapping)
