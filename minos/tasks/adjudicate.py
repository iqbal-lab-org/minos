from minos import adjudicator


def run(options):
    adj = adjudicator.Adjudicator(
        options.outdir,
        options.ref_fasta,
        options.reads,
        options.vcf_files,
        read_error_rate=options.read_error_rate,
        overwrite_outdir=options.force,
        max_alleles_per_cluster=options.max_alleles_per_cluster,
        gramtools_build_dir=options.gramtools_build_dir,
        sample_name=options.sample_name,
        variants_per_split=options.variants_per_split,
        alleles_per_split=options.alleles_per_split,
        total_splits=options.total_splits,
        clean=not options.debug,
        gramtools_kmer_size=options.gramtools_kmer_size,
        use_unmapped_reads=options.use_unmapped_reads,
        filter_min_dp=options.filter_min_dp,
        filter_min_gcp=options.filter_min_gcp,
        filter_max_dp=options.filter_max_dp,
        filter_min_frs=options.filter_min_frs,
        call_hets=options.include_het_calls,
        debug=options.debug,
    )
    adj.run()
