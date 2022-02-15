![Build Status](https://github.com/iqbal-lab-org/viridian_workflow/actions/workflows/build.yaml/badge.svg)

# Viridian Workflow

Please see the [Viridian Workflow Wiki](https://github.com/iqbal-lab-org/viridian_workflow/wiki)
for full documentation.

## Installation

Clone this repository and then from its root, build either a Singularity
or Docker container.

To build a singularity container:
```
singularity build viridian_workflow.img Singularity.def
```

To build a docker container:
```
docker build --network=host .
```
(without `--network=host` you will likely get `pip install` timing out and
the build failing).

Both the Docker and Singularity container will have the main script
`viridian_workflow` installed.

## Usage

To run on paired Illumina reads:
```
viridian_workflow run_one_sample \
  --tech illumina
  --ref_fasta data/MN908947.fasta \
  --reads1 reads_1.fastq.gz \
  --reads2 reads_2.fastq.gz \
  --outdir OUT
```
To run on unpaired nanopore reads:
```
viridian_workflow run_one_sample \
  --tech ont
  --ref_fasta data/MN908947.fasta \
  --reads reads.fastq.gz \
   --outdir OUT
```

The FASTA file in those commands can be found in the `viridian_workflow/amplicon_scheme_data/`
directory of this repository.

Other options:
* `--sample_name MY_NAME`: use this to change the sample name
  (default is "sample") that is put in the final FASTA file, BAM file, and
  VCF file.
* `--keep_bam`: use this option to keep the BAM file of original input reads
  mapped to the reference genome.
* `--force`: use with caution - it will overwrite the output directory if
  it already exists.


## Pipeline

```mermaid
flowchart TD
    A[Map] --> B[Identify Amplicons];
    B --> C[Downsample Reads];
    C --> D[Identify Primers];
    D --> E[Amplicon Assembly];
    E --> F[Consensus QC];
    D --> F;
    B -- reference_mapped.bam --> G;
    F -- consensus.fa --> G;
```

## Output files

The default files in the output directory are:

* `consensus.fa`: a FASTA file of the consensus sequence.
* `variants.vcf`: a VCF file of the identified variants between the consensus
  sequence and the reference genome.
* `log.json`: contains logging information for the viridian workflow run.

If the option `--keep_bam` is used, then a sorted BAM file of the reads mapped
to the reference will also be present, called
`reference_mapped.bam` (and its index file `reference_mapped.bam.bai`).

