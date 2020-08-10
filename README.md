[![Build Status](https://travis-ci.org/iqbal-lab-org/minos.svg?branch=master)](https://travis-ci.org/iqbal-lab-org/minos)

# minos
Variant call adjudication.


## Installation

Dependencies:

* Python 3 (tested on version 3.6.9)
* [vt](https://github.com/atks/vt.git)
* [vcflib](https://github.com/vcflib/vcflib.git). Specifically,
  `vcfbreakmulti`, `vcfallelicprimitives`, and `vcfuniq` must be installed.
* Optionally, [nextflow](https://www.nextflow.io/), if you want to use the
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
