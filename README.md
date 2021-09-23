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

## Configuration

Quality control thresholds are configured in `viridian_workflow/config.ini`:

```INI
[minimap2]
threads = 1

[viridian]
illumina_read_lengths = 50

[qcovid]
min_depth = 50
min_template_coverage_75 = 0.80
freq_threshold = 0.80
```

Viridian workflow applies quality control metrics at the amplicon and base level, as follows:

* Reads are aligned to the reference (default `MN908947`). This removes reads that do not sufficiently match this sequence.
* Individual reads must be sufficiently long. For Illumina reads, `illumina_read_lengths` must be, by default, at least 50. 
* Templates are reconstructed from either paired or single end reads that have been aligned. The start and end positions are used to infer which amplicon they belong to (e.g. `nCoV-artic-v3.bed`). "off-target" templates are discarded. 
* Enough templates must cover an amplicon or else all calls inside the region covered solely by that amplicon are voided with `N`s: `min_depth`, default 50.
* Enough templates must span at least 75% of the amplicon they belong to, excluding internal indels, where template length is the `(alignment_end - alignment_start) / amplicon_length`. This threshold is controlled by `min_template_coverage_75`, default 80%. This filters out short fragments caused by PCR artefacts. If this threshold fails the entire amplicon is discarded.
* After the assembly is constructed the original reads are mapped to it. Individual bases are voided with `N` if there is less than `freq_threshold` agreement. Default is 80%.
