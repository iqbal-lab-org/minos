from minos import multi_sample_pipeline

def run(options):
    pipeline = multi_sample_pipeline.MultiSamplePipeline(
        options.ref_fasta,
        options.data_tsv,
        options.outdir,
        min_large_ref_length=options.min_large_ref_length,
        gramtools_max_read_length=options.max_read_length,
        nextflow_config_file=options.nextflow_config_file,
        nextflow_work_dir=options.nextflow_work_dir,
        force=options.force,
        no_run=options.no_run,
    )
    pipeline.run()
