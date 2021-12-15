[![Build Status](https://travis-ci.org/iqbal-lab-org/minos.svg?branch=master)](https://travis-ci.org/iqbal-lab-org/minos)

# minos
Variant call adjudication.

Minimal instructions are below. Please see the [minos wiki page](https://github.com/iqbal-lab-org/minos/wiki)
for more details.

## Installation

Dependencies:

* Python 3 (tested on version 3.6.9)
* [gramtools](https://github.com/iqbal-lab-org/gramtools) commit
  8af53f6c8c0d72ef95223e89ab82119b717044f2
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

Alternatively, instead of running `pip3`, build a
[Singularity](https://sylabs.io/singularity/) container by running:

```
singularity build minos.simg Singularity.def
```

Build a Docker container by running:
```
sudo docker build --network=host .
```


## Quick start

Supposing you have variant calls for one sample, called by two different tools
in the files `calls1.vcf` and `calls2.vcf`. Run:

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
