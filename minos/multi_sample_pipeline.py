import logging
import os
import shutil

from cluster_vcf_records import vcf_file_read

from minos import dependencies, utils, vcf_file_split_deletions

class Error (Exception): pass

class MultiSamplePipeline:
    def __init__(self,
        ref_fasta,
        input_data_tsv,
        output_dir,
        max_alleles_per_cluster=5000,
        min_large_ref_length=50,
        gramtools_max_read_length=0,
        gramtools_kmer_size=5,
        gramtools_build_threads=1,
        nextflow_config_file=None,
        nextflow_work_dir=None,
        force=False,
        no_run=False,
        clean=True,
        testing=False,
        variants_per_split=None,
        alleles_per_split=None,
        total_splits=None,
        nf_ram_cluster_small_vars=2,
        nf_ram_gramtools_build_small=12,
        nf_ram_minos_small_vars=5,
        nf_ram_merge_small_vars=4,
        use_unmapped_reads=False,
    ):
        self.ref_fasta = os.path.abspath(ref_fasta)
        if not os.path.exists(self.ref_fasta):
            raise Error('Reference FASTA file not found: ' + ref_fasta)

        self.input_data_tsv = os.path.abspath(input_data_tsv)
        if not os.path.exists(self.input_data_tsv):
            raise Error('Data TSV file not found: ' + input_data_tsv)

        self.output_dir = os.path.abspath(output_dir)
        self.nextflow_config_file = None if nextflow_config_file is None else os.path.abspath(nextflow_config_file)
        self.max_alleles_per_cluster = max_alleles_per_cluster
        self.min_large_ref_length = min_large_ref_length
        self.gramtools_max_read_length = gramtools_max_read_length
        self.gramtools_kmer_size = gramtools_kmer_size
        self.gramtools_build_threads = gramtools_build_threads

        if nextflow_work_dir is None:
            self.nextflow_work_dir = os.path.join(self.output_dir, 'nextflow.work')
        else:
            self.nextflow_work_dir = os.path.abspath(nextflow_work_dir)

        self.force = force
        self.nextflow_input_tsv = os.path.join(self.output_dir, 'nextflow.input.tsv')
        self.log_file = os.path.join(self.output_dir, 'log.txt')
        self.no_run = no_run
        self.clean = clean
        self.testing = testing
        self.variants_per_split = variants_per_split
        self.alleles_per_split = alleles_per_split
        self.total_splits = total_splits
        self.nf_ram_cluster_small_vars =  nf_ram_cluster_small_vars
        self.nf_ram_gramtools_build_small = nf_ram_gramtools_build_small
        self.nf_ram_minos_small_vars = nf_ram_minos_small_vars
        self.nf_ram_merge_small_vars = nf_ram_merge_small_vars
        self.use_unmapped_reads = use_unmapped_reads



    @classmethod
    def _load_input_data_tsv(cls, infile):
        logging.info('Start reading file ' + infile)
        data = []
        with open(infile) as f:
            for line in f:
                try:
                    vcf_file, *reads_files = line.rstrip().split('\t')
                except:
                    raise Error('Bad line in input TSV file: ' + line.rstrip())

                if not(os.path.exists(vcf_file)):
                    raise Error('VCF file not found: ' + vcf_file)
                for reads_file in reads_files:
                    if not(os.path.exists(reads_file)):
                        raise Error('Reads file not found: ' + reads_file)

                data.append((os.path.abspath(vcf_file), [os.path.abspath(x) for x in reads_files]))

        logging.info('Finish reading file ' + infile + '. Loaded ' + str(len(data)) + ' samples')
        return data


    @classmethod
    def _merge_vcf_files(cls, infiles_list, outfile):
        '''Reimplementation of bcftools merge. Load all files into
        memory, then write output. bcftools opens all files at the same
        time, which doesn't work for lots of files'''
        vcf_file_data_list_per_sample = []
        sample_names = []
        header_lines = []
        common_first_column_data = []
        first_file = True

        for filename in infiles_list:
            new_data = []

            with open(filename) as f_vcf:
                for vcf_line in f_vcf:
                    if vcf_line.startswith('#'):
                        if first_file and vcf_line.startswith('##'):
                            header_lines.append(vcf_line.rstrip())
                        elif vcf_line.startswith('#CHROM'):
                            fields = vcf_line.rstrip().split('\t')
                            assert len(fields) == 10
                            sample_names.append(fields[-1])
                    else:
                        first_columns, last_column = vcf_line.rstrip().rsplit('\t', maxsplit=1)
                        new_data.append(last_column)
                        if first_file:
                            common_first_column_data.append(first_columns)

            vcf_file_data_list_per_sample.append(new_data)
            first_file = False

        with open(outfile, 'w') as f:
            print(*header_lines, sep='\n', file=f)
            print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *sample_names, sep='\t', file=f)
            for i, common_first_data in enumerate(common_first_column_data):
                sample_data_cols = [vcf_file_data_list_per_sample[j][i] for j in range(len(vcf_file_data_list_per_sample))]
                print(common_first_data, *sample_data_cols, sep='\t', file=f)


    @classmethod
    def _filter_input_file_for_clustering(cls, infile, outfile):
        header_lines, vcf_records = vcf_file_read.vcf_file_to_dict(infile, sort=True, homozygous_only=False, remove_asterisk_alts=True, remove_useless_start_nucleotides=True)
        with open(outfile, 'w') as f:
            print(*header_lines, sep='\n', file=f)
            for ref_name in vcf_records:
                for vcf_record in vcf_records[ref_name]:
                    if 'MISMAPPED_UNPLACEABLE' in vcf_record.FILTER:
                        continue
                    if vcf_record.FORMAT is None or 'GT' not in vcf_record.FORMAT:
                        logging.warning('No GT in vcf record:' + str(vcf_record))
                        continue

                    genotype = vcf_record.FORMAT['GT']
                    genotypes = genotype.split('/')

                    called_alleles = set(genotypes)
                    if called_alleles == {'0'} or '.' in called_alleles:
                        continue

                    genotypes = sorted([int(x) for x in genotypes])

                    if len(called_alleles) == 1:
                        assert 0 not in genotypes
                        vcf_record.set_format_key_value('GT', '1/1')
                        vcf_record.ALT = [vcf_record.ALT[int(genotypes[0]) - 1]]
                    else:
                        assert len(called_alleles) == 2
                        vcf_record.set_format_key_value('GT', '0/1')
                        if 0 in genotypes:
                            vcf_record.set_format_key_value('GT', '0/1')
                            vcf_record.ALT = [vcf_record.ALT[genotypes[1] - 1]]
                        else:
                            vcf_record.set_format_key_value('GT', '1/2')
                            vcf_record.ALT = [vcf_record.ALT[genotypes[0] - 1], vcf_record.ALT[genotypes[1] - 1]]

                    print(vcf_record, file=f)


    @classmethod
    def _nextflow_helper_process_input_vcf_file(cls, infile, out_small_vars, out_big_vars, out_sample_name, min_large_ref_length):
        splitter = vcf_file_split_deletions.VcfFileSplitDeletions(infile, out_small_vars, out_big_vars, min_large_ref_length=min_large_ref_length)
        splitter.run()
        header_lines = vcf_file_read.get_header_lines_from_vcf_file(infile)
        sample_name = vcf_file_read.get_sample_name_from_vcf_header_lines(header_lines)
        assert sample_name is not None
        max_read_length = None
        for line in header_lines:
            if line.startswith('##minos_max_read_length='):
                max_read_length = int(line.rstrip().split('=')[1])

        with open(out_sample_name, "w") as f:
            sample_name = vcf_file_read.get_sample_name_from_vcf_file(infile)
            assert sample_name is not None
            print(sample_name, file=f)

        return max_read_length


    @classmethod
    def _write_nextflow_data_tsv(cls, data, outfile):
        with open(outfile, 'w') as f:
            print('sample_id', 'vcf_file', 'reads_files', sep='\t', file=f)
            for i, (vcf_file, reads_files) in enumerate(data):
                print(i, vcf_file, ' '.join(reads_files), sep='\t', file=f)


    @classmethod
    def _write_nextflow_script(cls, outfile):
        with open(outfile, 'w') as f:
            print(r'''params.data_in_tsv = ""
params.ref_fasta = ""
params.min_large_ref_length = 0
params.gramtools_max_read_length = 0
params.gramtools_kmer_size = 0
params.gramtools_build_threads = 1
params.final_outdir = ""
params.testing = false
params.cluster_small_vars_ram = 2
params.gramtools_build_small_vars_ram = 12
params.minos_small_vars_ram = 5
params.pre_cluster_small_vars_merge_ram = 8
params.pre_cluster_small_vars_merge_threads = 10
params.merge_small_vars_ram = 4
params.variants_per_split = 0
params.alleles_per_split = 0
params.total_splits = 0
params.max_alleles_per_cluster = 5000
params.use_unmapped_reads = false


if (params.testing) {
    params.pre_cluster_small_vars_merge_threads = 2
}

if (params.use_unmapped_reads) {
    use_unmapped_reads = "--use_unmapped_reads"
}
else {
    use_unmapped_reads = ""
}

data_in_tsv = file(params.data_in_tsv).toAbsolutePath()
ref_fasta = file(params.ref_fasta).toAbsolutePath()
final_outdir = file(params.final_outdir).toAbsolutePath()

if (!data_in_tsv.exists()) {
    exit 1, "Input data TSV file not found: ${params.data_in_tsv} -- aborting"
}

if (!ref_fasta.exists()) {
    exit 1, "Reference FASTA file not found: ${params.ref_fasta} -- aborting"
}

if (params.min_large_ref_length < 1) {
    exit 1, "Must use option --min_large_ref_length -- aborting"
}

if (params.gramtools_kmer_size < 1) {
    exit 1, "Must use option --gramtools_kmer_size -- aborting"
}

if (!final_outdir.exists()) {
    exit 1, "Output directory not found: ${params.final_outdir} -- aborting"
}

split_tsv = Channel.from(data_in_tsv).splitCsv(header: true, sep:'\t')


process process_input_vcf_file {
    memory '0.5 GB'
    input:
    val tsv_fields from split_tsv

    output:
    file("small_vars.${tsv_fields['sample_id']}.vcf") into process_input_vcf_file_out_small
    set(val(tsv_fields), file("big_vars.${tsv_fields['sample_id']}.vcf")) into merge_small_and_large_vars_in
    set(val(tsv_fields), file("sample_name.${tsv_fields['sample_id']}")) into minos_all_small_vars_tsv_in
    stdout into max_read_lengths

    """
    #!/usr/bin/env python3
    from minos import multi_sample_pipeline
    multi_sample_pipeline.MultiSamplePipeline._filter_input_file_for_clustering(
        "${tsv_fields.vcf_file}",
        'filtered.vcf'
    )
    max_read_length = multi_sample_pipeline.MultiSamplePipeline._nextflow_helper_process_input_vcf_file(
        'filtered.vcf',
        "small_vars.${tsv_fields['sample_id']}.vcf",
        "big_vars.${tsv_fields['sample_id']}.vcf",
        "sample_name.${tsv_fields['sample_id']}",
        ${params.min_large_ref_length}
    )
    if max_read_length is None:
        print(0, end='')
    else:
        print(max_read_length, end='')
    """
}


process pre_cluster_small_vars_merge {
    errorStrategy {task.attempt < 3 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 1.GB * params.pre_cluster_small_vars_merge_ram * task.attempt}
    maxRetries 3
    cpus {params.testing? 2 : params.pre_cluster_small_vars_merge_threads}

    input:
    val(file_list) from process_input_vcf_file_out_small.collect()

    output:
    file('pre_cluster_small_vars_merge.vcf') into pre_cluster_small_vars_merge_out

    """
    #!/usr/bin/env python3
    from cluster_vcf_records import vcf_merge
    import pyfastaq
    ref_seqs = dict()
    pyfastaq.tasks.file_to_dict("${ref_fasta}", ref_seqs)

    file_list = ["${file_list.join('", "')}"]

    vcf_merge.merge_vcf_files(file_list, ref_seqs, "pre_cluster_small_vars_merge.vcf", threads=${params.pre_cluster_small_vars_merge_threads})
    """

}

process cluster_small_vars_vcf {
    errorStrategy {task.attempt < 3 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 1.GB * params.cluster_small_vars_ram * task.attempt}
    maxRetries 3

    input:
    file('pre_cluster_small_vars_merge.vcf') from pre_cluster_small_vars_merge_out

    output:
    file('small_vars_clustered.vcf') into cluster_small_vars_vcf_out

    """
    #!/usr/bin/env python3
    from cluster_vcf_records import vcf_clusterer
    clusterer = vcf_clusterer.VcfClusterer(["pre_cluster_small_vars_merge.vcf"], "${ref_fasta}", "small_vars_clustered.vcf", max_alleles_per_cluster=${params.max_alleles_per_cluster})
    clusterer.run()
    """
}


process gramtools_build_small_vars {
    errorStrategy {task.attempt < 3 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 1.GB * params.gramtools_build_small_vars_ram * task.attempt}
    maxRetries 3
    cpus params.gramtools_build_threads

    input:
    file('small_vars_clustered.vcf') from cluster_small_vars_vcf_out
    val(max_read_length) from max_read_lengths.max()

    output:
    set(file('small_vars_clustered.vcf'), file('small_vars_clustered.gramtools.build')) into gramtools_build_small_vars_out

    """
    #!/usr/bin/env python3
    import sys
    from minos import gramtools, vcf_chunker
    if ${params.gramtools_max_read_length} == 0:
        max_read_length = ${max_read_length}
    else:
        max_read_length = ${params.gramtools_max_read_length}

    if max_read_length == 0:
        print('Error! max read length could not be inferred from input VCF files. Must use option --gramtools_max_read_length')
        sys.exit(1)

    total_splits = ${params.total_splits} if ${params.total_splits} > 0 else None
    variants_per_split = ${params.variants_per_split} if ${params.variants_per_split} > 0 else None
    alleles_per_split = ${params.alleles_per_split} if ${params.alleles_per_split} > 0 else None
    print("total_splits", total_splits)
    print("variants_per_split", variants_per_split)
    print("alleles_per_split", alleles_per_split)

    if total_splits is None and variants_per_split is None and alleles_per_split is None:
        gramtools.run_gramtools_build(
            "small_vars_clustered.gramtools.build",
            "small_vars_clustered.vcf",
            "${params.ref_fasta}",
            max_read_length,
            kmer_size=${params.gramtools_kmer_size},
        )
    else:
        chunker = vcf_chunker.VcfChunker(
            "small_vars_clustered.gramtools.build",
            vcf_infile="small_vars_clustered.vcf",
            ref_fasta="${params.ref_fasta}",
            variants_per_split=variants_per_split,
            alleles_per_split=alleles_per_split,
            max_read_length=max_read_length,
            total_splits=total_splits,
            flank_length=max_read_length,
            gramtools_kmer_size=${params.gramtools_kmer_size},
            threads=${params.gramtools_build_threads},
        )
        chunker.make_split_files()
    """
}



process minos_all_small_vars {
    errorStrategy {task.attempt < 3 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 1.GB * params.minos_small_vars_ram * task.attempt}
    maxRetries 3

    input:
    set(file('small_vars_clustered.vcf'), file('small_vars_clustered.gramtools.build')) from gramtools_build_small_vars_out
    set(val(tsv_fields), file("sample_name.${tsv_fields['sample_id']}")) from minos_all_small_vars_tsv_in

    output:
    file("small_vars.minos.${tsv_fields['sample_id']}") into minos_all_small_vars_out

    """
    sample_name=\$(cat sample_name.${tsv_fields['sample_id']})
    minos_outdir=small_vars.minos.${tsv_fields['sample_id']}
    minos adjudicate ${use_unmapped_reads} --sample_name \$sample_name --gramtools_build_dir "small_vars_clustered.gramtools.build" --reads ${tsv_fields['reads_files'].replaceAll(/ /, " --reads ")} \$minos_outdir ${ref_fasta} "small_vars_clustered.vcf"
    """
}


// This takes the list of files output by minos_all_small_vars
// and merges them.
process merge_small_vars_vcfs {
    memory '1 GB'
    publishDir path: final_outdir, mode: 'move', overwrite: true

    input:
    val(minos_dir_list) from minos_all_small_vars_out.collect()

    output:
    file('combined_calls.vcf')

    """
    #!/usr/bin/env python3
    # Files end with .N (N=0,1,2,3,...) Sort numerically on this N
    import os
    from minos import multi_sample_pipeline

    minos_dir_list = ["${minos_dir_list.join('", "')}"]
    tuple_list = []
    for filename in minos_dir_list:
        fields = filename.rsplit('.', maxsplit=1)
        tuple_list.append((int(fields[1]), filename))
    tuple_list.sort()
    filenames = [os.path.join(x[1], 'debug.calls_with_zero_cov_alleles.vcf')  for x in tuple_list]
    multi_sample_pipeline.MultiSamplePipeline._merge_vcf_files(filenames, 'combined_calls.vcf')
    """
}

''', file=f)


    def _make_output_dir(self):
        if os.path.exists(self.output_dir):
            if self.force:
                shutil.rmtree(self.output_dir)
            else:
                raise Error('Error! Output directory already exists. ' + self.output_dir)
        os.mkdir(self.output_dir)


    def _prepare_nextflow_input_files(self):
        input_data = MultiSamplePipeline._load_input_data_tsv(self.input_data_tsv)
        MultiSamplePipeline._write_nextflow_data_tsv(input_data, self.nextflow_input_tsv)


    def run(self):
        self._make_output_dir()
        fh = logging.FileHandler(self.log_file, mode='w')
        log = logging.getLogger()
        formatter = logging.Formatter('[minos %(asctime)s %(levelname)s] %(message)s', datefmt='%d-%m-%Y %H:%M:%S')
        fh.setFormatter(formatter)
        log.addHandler(fh)
        dependencies.check_and_report_dependencies(programs=['nextflow'])

        self._prepare_nextflow_input_files()
        original_dir = os.getcwd()
        os.chdir(self.output_dir)
        nextflow_script = 'nextflow.run.nf'
        MultiSamplePipeline._write_nextflow_script(nextflow_script)
        logging.info('Prepared nextflow files. cd ' + self.output_dir)

        nextflow = dependencies.find_binary('nextflow')
        nextflow_command = [
            nextflow, 'run',
            '-work-dir', self.nextflow_work_dir,
            '-with-dag', 'nextflow.out.dag.pdf',
            '-with-trace', 'nextflow.out.trace.txt',
        ]

        if self.nextflow_config_file is not None:
            nextflow_command.extend(['-c', self.nextflow_config_file])

        nextflow_command += [
            nextflow_script,
            '--ref_fasta', self.ref_fasta,
            '--data_in_tsv', self.nextflow_input_tsv,
            '--max_alleles_per_cluster', str(self.max_alleles_per_cluster),
            '--min_large_ref_length', str(self.min_large_ref_length),
            '--final_outdir', self.output_dir,
            '--gramtools_max_read_length', str(self.gramtools_max_read_length),
            '--cluster_small_vars_ram', str(self.nf_ram_cluster_small_vars),
            '--gramtools_build_small_vars_ram', str(self.nf_ram_gramtools_build_small),
            '--gramtools_kmer_size', str(self.gramtools_kmer_size),
            '--gramtools_build_threads', str(self.gramtools_build_threads),
            '--minos_small_vars_ram', str(self.nf_ram_minos_small_vars),
            '--merge_small_vars_ram', str(self.nf_ram_merge_small_vars),
        ]

        if self.testing:
            nextflow_command.append('--testing')

        if self.use_unmapped_reads:
            nextflow_command.append('--use_unmapped_reads')

        if self.variants_per_split is not None:
            nextflow_command.append('--variants_per_split ' + str(self.variants_per_split))
        if self.alleles_per_split is not None:
            nextflow_command.append('--alleles_per_split ' + str(self.alleles_per_split))
        elif self.total_splits is not None:
            nextflow_command.append('--total_splits ' + str(self.total_splits))

        nextflow_command = ' '.join(nextflow_command)

        if self.no_run:
            print('Prepared nextflow pipeline. --no_run used, so not running. The nextflow command to run is:')
            print(nextflow_command)
            return
        else:
            logging.info('Start running nextflow: ' + nextflow_command)
            syscall_process = utils.syscall(nextflow_command)
            logging.info('Finish running nextflow. Writing nextflow stdout/stderr to files')
            with open('nextflow.stdout', 'w') as f:
                print(syscall_process.stdout.rstrip(), file=f)
            with open('nextflow.stderr', 'w') as f:
                print(syscall_process.stderr.rstrip(), file=f)

            logging.info('cd ' + original_dir)

        if self.clean:
            logging.info('Delete nextflow work directory ' + self.nextflow_work_dir)
            shutil.rmtree(self.nextflow_work_dir)
            logging.info('Delete .nextflow directory')
            shutil.rmtree('.nextflow')

        logging.info('Rename .nextflow.log -> nextflow.log')
        os.rename('.nextflow.log', 'nextflow.log')
        os.chdir(original_dir)
