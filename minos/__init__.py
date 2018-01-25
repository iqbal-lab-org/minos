from pkg_resources import get_distribution

try:
    __version__ = get_distribution('minos').version
except:
    __version__ = 'local'


__all__ = [
    'adjudicator',
    'dependencies',
    'genotyper',
    'gramtools',
    'mapping_based_verifier',
    'tasks',
    'utils',
]

from minos import *
