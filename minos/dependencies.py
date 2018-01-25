import json
import os
import shutil

import pyfastaq

from minos import utils, __version__

class Error (Exception): pass


def find_binary(program, allow_fail=False):
    binary = os.environ.get('MINOS_' + program.upper(), program)
    which_output = shutil.which(binary)
    if which_output is None and not allow_fail:
        raise Error('Error finding ' + program + ' in $PATH. Looked for "' + binary + '"')
    return which_output


def get_version_of_program(program, binary=None, allow_fail=False):
    if binary is None:
        binary = find_binary(program, allow_fail=allow_fail)

    if program == 'gramtools':
        gramtools_process = utils.syscall(binary + ' --version')
        gramtools_json = json.loads(gramtools_process.stdout)
        return gramtools_json.get('version_number', None)
    elif program == 'bwa':
        # To get version of BWA, need to run it with no options.
        # This returns an error code of 1, which we need to ignore
        bwa_process = utils.syscall(binary, allow_fail=True)
        for line in bwa_process.stderr.split('\n'):
            if line.strip().startswith('Version:'):
                try:
                    version = line.rstrip().split()[-1]
                except:
                    return None
                return version
        return None
    else:
        raise Error('Program name "' + program + '" not recognised. Cannot continue')


def find_python_packages():
    '''Returns a dict of data on required
    third-party python packeges, and minos itself:
    package_name => (version, path).
    Values will be None if package not found'''
    packages = {}
    for package in ['minos', 'pyfastaq', 'pysam']:
        try:
            exec('import ' + package)
            version = eval(package + '.__version__')
            path = eval(package + '.__file__')
        except:
            version = 'NOT_FOUND'
            path = 'NOT_FOUND'

        packages[package] = (version, path)

    return packages


def find_binaries_and_versions(programs=None):
    '''Returns a dict of data on required binaries:
    package_name => (version, path).
        Values will be None if package not found'''
    data = {}

    if programs is None:
        programs = ['bwa', 'gramtools']

    for program in programs:
        binary = find_binary(program, allow_fail=True)
        if binary is None:
            binary = 'NOT_FOUND'
            version = 'NOT_FOUND'
        else:
            version = get_version_of_program(program, binary, allow_fail=True)
            if version is None:
                version = 'NOT_FOUND'

        data[program] = (version, binary)

    return data


def dependencies_report(programs=None):
    '''Returns tuple:
    bool, list
    bool = True iff all dependencies are found
    list = a list of lines that are a
    report on the third party programs and Python
    packages required by minos'''
    all_ok = True
    lines = ['minos ' + __version__]

    if programs is None:
        programs = ['bwa', 'gramtools']
    programs_data = find_binaries_and_versions(programs=programs)

    for program in programs:
        lines.append(program + ' ' + ' '.join(programs_data[program]))
        if 'NOT_FOUND' in programs_data[program]:
            all_ok = False

    packages_data = find_python_packages()
    for package in sorted(packages_data):
        lines.append(package + ' ' + ' '.join(packages_data[package]))
        if 'NOT_FOUND' in packages_data[package]:
            all_ok = False

    return all_ok, lines


def check_and_report_dependencies(outfile, programs=None):
    '''Writes report of depndencies to file (could be stdout).
    Raises error if any depndency not found'''
    all_ok, report_lines = dependencies_report(programs=programs)
    f = pyfastaq.utils.open_file_write(outfile)
    print(*report_lines, sep='\n', file=f)
    print('All ok:', all_ok, file=f)
    pyfastaq.utils.close(f)
    if not all_ok:
        raise Error('At least one dependency not found')

