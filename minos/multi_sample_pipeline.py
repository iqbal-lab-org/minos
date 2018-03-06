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
        min_large_ref_length=50,
        gramtools_max_read_length=0,
        nextflow_config_file=None,
        nextflow_work_dir=None,
        force=False,
        no_run=False,
        clean=True,
        testing=False,
        nf_ram_cluster_small_vars=2,
        nf_ram_gramtools_build_small=12,
        nf_ram_minos_small_vars=5,
        nf_ram_bcftools_merge=2,
    ):
        self.ref_fasta = os.path.abspath(ref_fasta)
        if not os.path.exists(self.ref_fasta):
            raise Error('Reference FASTA file not found: ' + ref_fasta)

        self.input_data_tsv = os.path.abspath(input_data_tsv)
        if not os.path.exists(self.input_data_tsv):
            raise Error('Data TSV file not found: ' + input_data_tsv)

        self.output_dir = os.path.abspath(output_dir)
        self.nextflow_config_file = None if nextflow_config_file is None else os.path.abspath(nextflow_config_file)
        self.min_large_ref_length = min_large_ref_length
        self.gramtools_max_read_length = gramtools_max_read_length

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
        self.nf_ram_cluster_small_vars =  nf_ram_cluster_small_vars
        self.nf_ram_gramtools_build_small = nf_ram_gramtools_build_small
        self.nf_ram_minos_small_vars = nf_ram_minos_small_vars
        self.nf_ram_bcftools_merge = nf_ram_bcftools_merge


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
params.final_outdir = ""
params.bcftools = "bcftools"
params.testing = false
params.cluster_small_vars_ram = 2
params.gramtools_build_small_vars_ram = 12
params.minos_small_vars_ram = 5
params.bcftools_merge_ram = 2


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
    max_read_length = multi_sample_pipeline.MultiSamplePipeline._nextflow_helper_process_input_vcf_file(
        "${tsv_fields.vcf_file}",
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


process cluster_small_vars_vcf {
    errorStrategy {task.attempt < 3 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 1.GB * params.cluster_small_vars_ram * task.attempt}
    maxRetries 3

    input:
    val(file_list) from process_input_vcf_file_out_small.collect()

    output:
    file('small_vars_clustered.vcf') into cluster_small_vars_vcf_out

    """
    #!/usr/bin/env python3
    from cluster_vcf_records import vcf_clusterer
    file_list = ["${file_list.join('", "')}"]
    clusterer = vcf_clusterer.VcfClusterer(file_list, "${ref_fasta}", "small_vars_clustered.vcf")
    clusterer.run()
    """
}


process gramtools_build_small_vars {
    errorStrategy {task.attempt < 3 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 1.GB * params.gramtools_build_small_vars_ram * task.attempt}
    maxRetries 3

    input:
    file('small_vars_clustered.vcf') from cluster_small_vars_vcf_out
    val(max_read_length) from max_read_lengths.max()

    output:
    set(file('small_vars_clustered.vcf'), file('small_vars_clustered.gramtools.build')) into gramtools_build_small_vars_out

    """
    #!/usr/bin/env python3
    import sys
    from minos import gramtools
    if ${params.gramtools_max_read_length} == 0:
        max_read_length = ${max_read_length}
    else:
        max_read_length = ${params.gramtools_max_read_length}

    if max_read_length == 0:
        print('Error! max read length could not be inferred from input VCF files. Must use option --gramtools_max_read_length')
        sys.exit(1)

    gramtools.run_gramtools_build("small_vars_clustered.gramtools.build", "small_vars_clustered.vcf", "${params.ref_fasta}", max_read_length)
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
    minos adjudicate --sample_name \$sample_name --gramtools_build_dir "small_vars_clustered.gramtools.build" --reads ${tsv_fields['reads_files'].replaceAll(/ /, " --reads ")} \$minos_outdir ${ref_fasta} "small_vars_clustered.vcf"
    bgzip \$minos_outdir/final.vcf
    tabix -p vcf \$minos_outdir/final.vcf.gz
    """
}


// This just takes the list of files output by minos_all_small_vars
// and writes a new file, with them sorted. Have this as a separate
// process from the bash script that runs bcftools because the number
// of files could go over bash's character/length limits
process make_bcftools_small_vars_fofn {
    memory '1 GB'

    input:
    val(minos_dir_list) from minos_all_small_vars_out.collect()

    output:
    file('vcf_file_list_for_bcftools.txt') into make_bcftools_small_vars_fofn_out

    """
    #!/usr/bin/env python3
    # Files end with .N (N=0,1,2,3,...) Sort numerically on this N
    import os
    minos_dir_list = ["${minos_dir_list.join('", "')}"]
    tuple_list = []
    for filename in minos_dir_list:
        fields = filename.rsplit('.', maxsplit=1)
        tuple_list.append((int(fields[1]), fields[0]))
    tuple_list.sort()

    with open('vcf_file_list_for_bcftools.txt', 'w') as f:
        for t in tuple_list:
            vcf = os.path.join(t[1] + '.' + str(t[0]), 'final.vcf.gz')
            print(vcf, file=f)
    """
}


process bcftools_merge {
    errorStrategy {task.attempt < 3 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 1.GB * params.bcftools_merge_ram * task.attempt}
    maxRetries 3

    publishDir path: final_outdir, mode: 'move', overwrite: true
    input:
    file('vcf_file_list_for_bcftools.txt') from make_bcftools_small_vars_fofn_out

    output:
    file('combined_calls.vcf')

    """
    ${params.bcftools} merge -l vcf_file_list_for_bcftools.txt -o combined_calls.vcf
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
        dependencies.check_and_report_dependencies(programs=['bcftools', 'nextflow'])

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

        bcftools = dependencies.find_binary('bcftools')

        nextflow_command += [
            nextflow_script,
            '--bcftools', bcftools,
            '--ref_fasta', self.ref_fasta,
            '--data_in_tsv', self.nextflow_input_tsv,
            '--min_large_ref_length', str(self.min_large_ref_length),
            '--final_outdir', self.output_dir,
            '--gramtools_max_read_length', str(self.gramtools_max_read_length),
            '--cluster_small_vars_ram', str(self.nf_ram_cluster_small_vars),
            '--gramtools_build_small_vars_ram', str(self.nf_ram_gramtools_build_small),
            '--minos_small_vars_ram', str(self.nf_ram_minos_small_vars),
            '--bcftools_merge_ram', str(self.nf_ram_bcftools_merge),
        ]

        if self.testing:
            nextflow_command.append('--testing')

        nextflow_command = ' '.join(nextflow_command)

        if self.no_run:
            print('Prepared nextflow pipeline. --no_run used, so not running. The nextflow command to run is:')
            print(nextflow_command)
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
