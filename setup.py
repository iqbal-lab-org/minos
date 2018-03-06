import glob
from setuptools import setup, find_packages


setup(
    name='bio-minos',
    version='0.3.0',
    description='Variant call adjudication',
    packages = find_packages(),
    author='Martin Hunt',
    author_email='mhunt@ebi.ac.uk',
    url='https://github.com/iqbal-lab-org/minos',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
        'cluster_vcf_records >= 0.5.0',
        'gramtools',
        'pyfastaq >= 3.14.0',
        'pymummer >= 0.10.3',
        'pysam >= 0.12',
        'scipy >= 1.0.0',
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)

