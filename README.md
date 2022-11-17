![Build Status](https://github.com/iqbal-lab-org/minos/actions/workflows/build.yaml/badge.svg)

# minos
Variant call adjudication.

Minimal instructions are below. Please see the [minos wiki page](https://github.com/iqbal-lab-org/minos/wiki)
for more details.

## Installation

### Docker

Get a Docker image of the latest release:

```
docker pull ghcr.io/iqbal-lab-org/minos:latest
```

All Docker images are listed in the
[packages page](https://github.com/iqbal-lab-org/minos/pkgs/container/minos).

Alternatively, build your own Docker image:

```
sudo docker build --network=host .
```

### Singularity

[Releases](https://github.com/iqbal-lab-org/minos/releases)
include a Singularity image to download (from version 0.12.1 onwards).

Alternatively, build your own Singularity image:

```
singularity build minos.simg Singularity.def
```


### From source

Dependencies:

* Python 3 (tested on version 3.6.9)
* [gramtools](https://github.com/iqbal-lab-org/gramtools) commit
  04c4ba717399507b643fd4b77a61c048ef2ed83f
* [bcftools](https://samtools.github.io/bcftools/)
* [vt](https://github.com/atks/vt.git)
* [vcflib](https://github.com/vcflib/vcflib.git). Specifically,
  either `vcflib`, or all three of
  `vcfbreakmulti`, `vcfallelicprimitives`, and `vcfuniq` must be installed.
* Optionally, [nextflow](https://www.nextflow.io/) and [ivcfmerge](https://github.com/iqbal-lab-org/ivcfmerge) if you want to use the
  pipeline to regenotype a large number of samples.

Install by cloning this repository (or downloading the latest release), and
running:

```
pip3 install .
```


## Quick start

To run on one sample, you will need:
* A FASTA file of the reference genome.
* One or more VCF files of variant calls.
  The only requirement of these files is that they must contain the genotype field `GT`,
  and correspond to the reference FASTA file. All variants with a non-reference genotype
  call will be used (both alleles are considered for diploid calls)
* Illumina reads in FASTQ file(s).

For example, if you have two call sets in the files `calls1.vcf` and `calls2.vcf`,
then run:

```
minos adjudicate --reads reads1.fq --reads reads2.fq out ref.fasta calls1.vcf calls2.vcf
```

where `reads1.fq` and `reads2.fq` are FASTQ files of the reads and `ref.fasta`
is a FASTA of the reference corresponding to the two input VCF files.
The final call set will be `out/final.vcf`.


## Unit tests

Run `tox` to run all unit tests.
They require `nextflow`, `gramtools`, `vt`, `vcfbreakmulti`,
`vcfallelicprimitives`, `vcfuniq`  in your `$PATH`.

Run an individual test file with `tox tests/for_test.py::TestSpam::test_eggs`.

Run the main entry point with `python3 -m minos`.
