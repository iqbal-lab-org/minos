from minos import multi_sample_pipeline

def run(options):
    pipeline = multi_sample_pipeline.MultiSamplePipeline(
        options.ref_fasta,
        options.data_tsv,
        options.outdir,
        max_alleles_per_cluster=options.max_alleles_per_cluster,
        min_large_ref_length=options.min_large_ref_length,
        gramtools_max_read_length=options.max_read_length,
        gramtools_kmer_size=options.gramtools_kmer_size,
        gramtools_build_threads=options.gramtools_build_threads,
        nextflow_config_file=options.nextflow_config_file,
        nextflow_work_dir=options.nextflow_work_dir,
        force=options.force,
        no_run=options.no_run,
        clean=not options.no_clean,
        variants_per_split=options.variants_per_split,
        alleles_per_split=options.alleles_per_split,
        total_splits=options.total_splits,
        nf_ram_cluster_small_vars=options.nf_ram_cluster_small_vars,
        nf_ram_gramtools_build_small=options.nf_ram_gramtools_build_small,
        nf_ram_minos_small_vars=options.nf_ram_minos_small_vars,
        nf_ram_merge_small_vars=options.nf_ram_merge_small_vars,
        testing=options.testing,
        use_unmapped_reads=options.use_unmapped_reads,
    )
    pipeline.run()
