# minos
Variant call adjudication.

## Dependencies

For python dependencies, see `setup.py`; 
can install them running `python3 setup.py install`.

On top of this, need:
* bgzip
* bwa
* gramtools
* mummer (for pymummer package)
* nextflow

Possibly the easiest way to get these is setting up a Vagrant machine using
[clockwork's vagrant file](https://github.com/iqbal-lab-org/clockwork/wiki/Information-for-developers)

## Unit tests

Run `python3 setup.py test`, or `python3 -m unittest discover -p "*test*"` to run all unit tests.
