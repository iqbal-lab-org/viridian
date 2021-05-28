[![Build Status](https://www.travis-ci.com/iqbal-lab-org/viridian_workflow.svg?branch=master)](https://www.travis-ci.com/iqbal-lab-org/viridian_workflow)
# Viridian Workflow

## Installation

```
singularity build viridian_workflow.img Singularity.def
```

## Usage

```
singularity run viridian_workflow.img run_one_sample data/MN908947.fasta data/nCoV-artic-v3.bed sample_R1.fastq.gz sample_R2.fastq.gz sample_outdir/
```
