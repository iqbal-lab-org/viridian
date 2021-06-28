![Build Status](https://github.com/iqbal-lab-org/viridian_workflow/actions/workflows/build.yaml/badge.svg)
[![REUSE status](https://api.reuse.software/badge/github.com/iqbal-lab-org/viridian_workflow)](https://api.reuse.software/info/github.com/iqbal-lab-org/viridian_workflow)
 
# Viridian Workflow

## Installation

```
singularity build viridian_workflow.img Singularity.def
```

## Usage

```
singularity run viridian_workflow.img run_one_sample data/MN908947.fasta data/nCoV-artic-v3.bed sample_R1.fastq.gz sample_R2.fastq.gz sample_outdir/
```

Single-end read Nanopore:

```
singularity run viridian_workflow.img run_one_sample_ont data/MN908947.fasta data/nCoV-artic-v3.bed sample.fastq.gz sample_R2.fastq.gz sample_outdir/
```


