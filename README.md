![Build Status](https://github.com/iqbal-lab-org/viridian_workflow/actions/workflows/build.yaml/badge.svg)

# Viridian Workflow

Please see the [Viridian Workflow Wiki](https://github.com/iqbal-lab-org/viridian_workflow/wiki)
for full documentation.

## Installation

The recommended method is to use a pre-built Docker or Singularity container
(see the wiki for how to build your own).

Both the Docker and Singularity container have the main script
`viridian` installed.

### Docker
Get a Docker image of the latest release:
```
docker pull ghcr.io/iqbal-lab-org/cte:latest
```
All Docker images are listed in the
[packages page](https://github.com/iqbal-lab-org/viridian_workflow/pkgs/container/viridian_workflow).

### Singularity
[Releases](https://github.com/iqbal-lab-org/viridian_workflow/releases)
include a Singularity image to download.
Each release has a singularity image file called
`viridian_workflow_vX.Y.Z.img`, where `X.Y.Z` is the release version.


## Usage

These instructions assume that you are assembling SARS-CoV-2 data.

To run on paired Illumina reads:
```
viridian run_one_sample \
  --tech illumina \
  --reads1 reads_1.fastq.gz \
  --reads2 reads_2.fastq.gz \
  --outdir OUT
```
To run on unpaired nanopore reads:
```
viridian run_one_sample \
  --tech ont \
  --reads reads.fastq.gz \
   --outdir OUT
```

To run on paired or unpaired Ion Torrent reads, use either of the
above commands, but with the option `--tech iontorrent`.


## Output files

The default files in the output directory are:

* `consensus.fa.gz`: a gzipped FASTA file of the consensus sequence.
* `variants.vcf`: a VCF file of the identified variants between the consensus
  sequence and the reference genome.
* `log.json.gz`: a gzipped JSON file that contains logging information
  for the viridian workflow run.
* `qc.tsv.gz`: a gzipped tab-delimited file of per-base QC information
* `scheme_id.depth_across_genome.pdf`: a plot of the read depth across
  the genome, with amplicons coloured in the background.
* `scheme_id.score_plot.pdf`: a plot of the scoring for amplicon scheme
  identification.


If the option `--keep_bam` is used, then a sorted BAM file of the reads mapped
to the reference will also be present, called
`reference_mapped.bam` (and its index file `reference_mapped.bam.bai`).


## Useful options

* `--sample_name MY_NAME`: use this to change the sample name
  (default is "sample") that is put in the final FASTA file, BAM file, and
  VCF file.
* `--reads_bam MY_READS.bam`: instead of providing FASTQ (or FASTA) files of
  reads, you can provide a sorted by genome coordinate and indexed BAM file.
  The reference genome must be the same as that used by viridian
  (by default MN908947).
* `--keep_bam`: use this option to keep the BAM file of original input reads
  mapped to the reference genome.
* `--decontam COVID`: decontaminate the reads using ReadItAndKeep at the
  start of the pipeline (this is incompatible with `--reads_bam`)
* `--force`: use with caution - it will overwrite the output directory if
  it already exists.
* `--write_msa indel_as_ref`: this will write a FASTA file
  called `msa.indel_as_ref.fa` that can be
  used to build trees. It has the consensus sequence aligned to the
  reference genome, but with insertions in the consensus ignored and
  deletions replaced with the reference sequence.
* `--ena_run RUN_ID`: using this option will download the specified reads
  from the ENA, and infer the `--tech` option from the ENA metadata

