from pkg_resources import get_distribution

try:
    __version__ = get_distribution('bio-minos').version
except:
    __version__ = 'local'


__all__ = [
    'adjudicator',
    'bam_read_extract',
    'dependencies',
    'genotyper',
    'genotype_confidence_simulator',
    'gramtools',
    'mapping_based_verifier',
    'multi_sample_pipeline',
    'plots',
    'tasks',
    'utils',
    'vcf_chunker',
    'vcf_file_split_deletions',
]

from minos import *
