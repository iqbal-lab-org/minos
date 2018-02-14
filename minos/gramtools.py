import datetime
import json
import os

from cluster_vcf_records import vcf_file_read

from minos import dependencies, genotyper, utils
from minos import __version__ as minos_version

class Error (Exception): pass


def _get_quasimap_out_dir(gramtools_dir):
    '''When quasimap is run, it makes a new dir like gramtools.out/quasimap_outputs/1510738443_ksize15.
       This function returns that directory. Dies if finds more than one directory in quasimap_outputs/'''
    if not os.path.exists(gramtools_dir):
        raise Error('Gramtools directory not found "' + gramtools_dir + '"')
    quasimap_outputs = os.path.join(gramtools_dir, 'quasimap_outputs')
    if not os.path.exists(quasimap_outputs):
        raise Error('quasimap outputs directory not found "' + quasimap_outputs + '"')

    subdirs = [x for x in os.listdir(quasimap_outputs) if os.path.isdir(os.path.join(quasimap_outputs, x))]
    if len(subdirs) == 0:
        raise Error('No directories found in quasimap outputs directory "' + quasimap_outputs + '"')
    elif len(subdirs) > 1:
        raise Error(str(len(subdirs)) + ' directories found in quasimap outputs directory "' + quasimap_outputs + '"')

    return os.path.join(quasimap_outputs, subdirs[0])


def run_gramtools(output_dir, vcf_file, ref_file, reads, max_read_length):
    '''Runs gramtools build and quasimap. Returns quasimap output directory.
    "reads" can be one filename, or a list of filenames.'''
    if type(reads) is not list:
        assert type(reads) is str
        reads = [reads]

    gramtools_exe = dependencies.find_binary('gramtools')
    build_command = ' '.join([
        gramtools_exe,
        'build',
        '--gram-directory', output_dir,
        '--vcf', vcf_file,
        '--reference', ref_file,
        '--max-read-length', str(max_read_length),
    ])
    utils.syscall(build_command)

    quasimap_command = ' '.join([
        gramtools_exe,
        'quasimap',
        '--gram-directory', output_dir,
        ' '.join(['--reads ' + x for x in reads]),
    ])
    utils.syscall(quasimap_command)
    return _get_quasimap_out_dir(output_dir)


def load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir):
    '''Loads the perl_generated_vcf file and allele_coverage files.
    Sanity checks that they agree: 1) same number of lines (excluding header
    lines in vcf) and 2) number of alts agree on each line.
    Raises error at the first time somthing wrong is found.
    Returns a list of tuples: (VcfRecord, dict of allele -> coverage)'''
    allele_base_counts_file = os.path.join(quasimap_dir, 'allele_base_coverage.json')
    grouped_allele_counts_file = os.path.join(quasimap_dir, 'grouped_allele_counts_coverage.json')
    all_allele_coverage, allele_groups = load_allele_files(allele_base_counts_file, grouped_allele_counts_file)
    vcf_header, vcf_lines = vcf_file_read.vcf_file_to_list(vcf_file)
    total_coverage = 0

    if len(all_allele_coverage) != len(vcf_lines):
        raise Error('Number of records in VCF (' + str(len(vcf_lines)) + ') does not match number output from gramtools.(' + str(len(all_allele_coverage)) + '). Cannot continue')

    for i, (allele_combi_coverage, allele_per_base_coverage) in enumerate(all_allele_coverage):
        if len(allele_per_base_coverage) != 1 + len(vcf_lines[i].ALT):
            raise Error('Mismatch in number of alleles for this VCF record:\n' + str(vcf_lines[i]) + '\nLine number is ' + str(i+1))

        total_coverage += sum(allele_combi_coverage.values())

    return round(total_coverage/len(vcf_lines), 3), vcf_header, vcf_lines, all_allele_coverage, allele_groups


def update_vcf_record_using_gramtools_allele_depths(vcf_record, allele_combination_cov, allele_per_base_cov, allele_groups_dict, mean_depth, read_error_rate):
    '''allele_depths should be a dict of allele -> coverage.
    The REF allele must also be in the dict.
    So keys of dict must be equal to REF + ALTs sequences.
    This also changes all columns from QUAL onwards'''
    gtyper = genotyper.Genotyper(mean_depth, read_error_rate, allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    gtyper.run()
    genotype_indexes = set()

    if 0 in gtyper.genotype:
        genotype_indexes.add(0)
    for i in range(len(vcf_record.ALT)):
        if i + 1 in gtyper.genotype:
            genotype_indexes.add(i+1)

    if len(genotype_indexes) == 1:
        genotype_index = genotype_indexes.pop()
        genotype = str(genotype_index) + '/' + str(genotype_index)
    else:
        genotype = '/'.join([str(x) for x in sorted(list(genotype_indexes))])

    cov_string = ','.join([str(gtyper.singleton_alleles_cov.get(x, 0)) for x in range(1 + len(vcf_record.ALT))])
    vcf_record.QUAL = None
    vcf_record.FILTER = '.'
    vcf_record.INFO = {'DP': str(sum(allele_combination_cov.values()))}
    vcf_record.format_keys = ['GT', 'COV', 'GT_CONF']
    vcf_record.FORMAT = {
        'GT': genotype,
        'COV': cov_string,
        'GT_CONF': str(gtyper.genotype_confidence)
    }


def write_vcf_annotated_using_coverage_from_gramtools(mean_depth, vcf_records, all_allele_coverage, allele_groups, read_error_rate, outfile, sample_name='SAMPLE'):
    '''mean_depth, vcf_records, all_allele_coverage, allele_groups should be those
    returned by load_gramtools_vcf_and_allele_coverage_files().
    Writes a new VCF that has allele counts for all the ALTs'''
    assert len(vcf_records) == len(all_allele_coverage)

    with open(outfile, 'w') as f:
        print('##fileformat=VCFv4.2', file=f)
        print('##source=minos, version', minos_version, file=f)
        print('##fileDate=', datetime.date.today(), sep='', file=f)
        print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_name, sep='\t', file=f)

        for i in range(len(vcf_records)):
            update_vcf_record_using_gramtools_allele_depths(vcf_records[i], all_allele_coverage[i][0], all_allele_coverage[i][1], allele_groups, mean_depth, read_error_rate)
            print(vcf_records[i], file=f)


def load_allele_files(allele_base_counts_file, grouped_allele_counts_file):
    '''Loads the allele base counts and groupeed allele counts files
    made by gramtools qausimap.'''
    with open(allele_base_counts_file) as f:
        json_base_counts_data = json.load(f)
    with open(grouped_allele_counts_file) as f:
        json_allele_counts_data = json.load(f)

    try:
        allele_base_counts = json_base_counts_data['allele_base_counts']
    except:
        raise Error('Error in json file ' + allele_base_counts_file + '. allele_base_counts not found.')

    try:
        site_counts = json_allele_counts_data['grouped_allele_counts']['site_counts']
    except:
        raise Error('Error in json file ' + grouped_allele_counts_file + '. site_counts not found.')


    try:
        allele_groups = json_allele_counts_data['grouped_allele_counts']['allele_groups']
    except:
         raise Error('Error in json file ' + grouped_allele_counts_file + '. allele_groups not found.')

    if len(allele_base_counts) != len(site_counts):
        raise Error('Mismatch between number of records in json files ' + allele_base_counts_file + ' and ' + grouped_allele_counts_file)

    for key, value in allele_groups.items():
        allele_groups[key] = set(value)

    return list(zip(site_counts, allele_base_counts)), allele_groups
