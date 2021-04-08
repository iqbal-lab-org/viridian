"""minimap2 and samtools utilities
"""
import os
# TODO: import logging
import subprocess

def run(outdir, ref_genome, fq1, fq2):
    minimap_cmd = ['minimap2', '-ax', 'sr', ref_genome, fq1, fq2]
    sort_cmd = ['samtools', 'sort', '-o', bam]
    bam = os.path.join(outdir, 'reference_mapped.bam')

    map_proc = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
    sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
    sort_proc.wait()
    assert os.path.isfile(bam)
    
    index_proc = subprocess.Popen(['samtools', 'index', bam])
    assert os.path.isfile(os.path.join(bam + '.bai'))
    return bam
