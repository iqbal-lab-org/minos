import datetime
import json
import logging
import os

from cluster_vcf_records import vcf_file_read

from minos import dependencies, genotyper, utils
from minos import __version__ as minos_version

class Error (Exception): pass


def run_gramtools_build(outdir, vcf_file, ref_file, max_read_length, kmer_size=15):
    '''Runs gramtools build. Makes new directory called 'outdir' for
    the output'''
    gramtools_exe = dependencies.find_binary('gramtools')
    build_command = ' '.join([
        gramtools_exe,
        'build',
        '--gram-directory', outdir,
        '--vcf', vcf_file,
        '--reference', ref_file,
        '--max-read-length', str(max_read_length),
        '--kmer-size', str(kmer_size),
    ])
    logging.info('Running gramtools build: ' + build_command)
    utils.syscall(build_command)
    logging.info('Finished running gramtools build')


def run_gramtools(build_dir, quasimap_dir, vcf_file, ref_file, reads, max_read_length, kmer_size=15):
    '''If build_dir does not exist, runs runs gramtools build and quasimap.
    Otherwise, just runs quasimap. quasimap output is in new
    directory called quasimap_dir.
    "reads" can be one filename, or a list of filenames.
    Raises Error if either of the expected json coverage 
    files made by quasimap are not found.'''
    gramtools_exe = dependencies.find_binary('gramtools')
    if not os.path.exists(build_dir):
        run_gramtools_build(build_dir, vcf_file, ref_file, max_read_length, kmer_size=kmer_size)

    if type(reads) is not list:
        assert type(reads) is str
        reads = [reads]

    quasimap_command = ' '.join([
        gramtools_exe,
        'quasimap',
        '--gram-directory', build_dir,
        '--run-directory', quasimap_dir,
        ' '.join(['--reads ' + x for x in reads]),
    ])
    logging.info('Running gramtools quasimap: ' + quasimap_command)
    utils.syscall(quasimap_command)
    logging.info('Finished running gramtools quasimap')

    build_report = os.path.join(build_dir, 'build_report.json')
    quasimap_report = os.path.join(quasimap_dir, 'report.json')
    allele_base_counts_file = os.path.join(quasimap_dir, 'allele_base_coverage.json')
    grouped_allele_counts_file = os.path.join(quasimap_dir, 'grouped_allele_counts_coverage.json')
    files_ok = True
    for filename in build_report, quasimap_report, allele_base_counts_file, grouped_allele_counts_file:
        if not(os.path.exists(filename)):
            files_ok = False
            logging.error('gramtools file not found: ' + filename)

    if not files_ok:
        error_message = 'Looks like something went wrong duing gramtools run. At least one output file not present. Cannot continue.'
        logging.error(error_message)
        raise Error(error_message)

    with open(build_report) as f:
        json_build_report = json.load(f)
    with open(quasimap_report) as f:
        json_quasimap_report = json.load(f)

    return json_build_report, json_quasimap_report


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


def update_vcf_record_using_gramtools_allele_depths(vcf_record, allele_combination_cov, allele_per_base_cov, allele_groups_dict, mean_depth, read_error_rate, kmer_size):
    '''allele_depths should be a dict of allele -> coverage.
    The REF allele must also be in the dict.
    So keys of dict must be equal to REF + ALTs sequences.
    This also changes all columns from QUAL onwards'''
    gtyper = genotyper.Genotyper(mean_depth, read_error_rate, allele_combination_cov, allele_per_base_cov, allele_groups_dict)
    gtyper.run()
    genotype_indexes = set()

    if '.' in gtyper.genotype:
        genotype = './.'
    else:
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
    vcf_record.INFO = {'KMER': str(kmer_size)}
    vcf_record.format_keys = ['DP', 'GT', 'COV', 'GT_CONF']
    vcf_record.FORMAT = {
        'DP': str(sum(allele_combination_cov.values())),
        'GT': genotype,
        'COV': cov_string,
        'GT_CONF': str(gtyper.genotype_confidence)
    }


def write_vcf_annotated_using_coverage_from_gramtools(mean_depth, vcf_records, all_allele_coverage, allele_groups, read_error_rate, outfile, kmer_size, sample_name='SAMPLE', max_read_length=None):
    '''mean_depth, vcf_records, all_allele_coverage, allele_groups should be those
    returned by load_gramtools_vcf_and_allele_coverage_files().
    Writes a new VCF that has allele counts for all the ALTs'''
    assert len(vcf_records) == len(all_allele_coverage)

    with open(outfile, 'w') as f:
        print('##fileformat=VCFv4.2', file=f)
        print('##source=minos, version', minos_version, file=f)
        print('##fileDate=', datetime.date.today(), sep='', file=f)
        print('##FORMAT=<ID=COV,Number=R,Type=Integer,Description="Number of reads on ref and alt alleles">', file=f)
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=f)
        print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="total kmer depth from gramtools",Source="minos">', file=f)
        print('##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence. Difference in log likelihood of most likely and next most likely genotype">', file=f)
        print('##INFO=<ID=KMER,Number=1,Type=Integer,Description="Kmer size at which variant was discovered (kmer-size used by gramtools build)">', file=f)

        if max_read_length is not None:
            print('##minos_max_read_length=' + str(max_read_length), file=f)

        print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_name, sep='\t', file=f)

        for i in range(len(vcf_records)):
            logging.debug('Genotyping: ' + str(vcf_records[i]))
            update_vcf_record_using_gramtools_allele_depths(vcf_records[i], all_allele_coverage[i][0], all_allele_coverage[i][1], allele_groups, mean_depth, read_error_rate, kmer_size)
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
