from minos import adjudicator

def run(options):
    adj = adjudicator.Adjudicator(
        options.outdir,
        options.ref_fasta,
        options.reads,
        options.vcf_files,
        max_read_length=options.max_read_length,
        read_error_rate=options.read_error_rate,
        overwrite_outdir=options.force,
        max_alleles_per_cluster=options.max_alleles_per_cluster,
        gramtools_build_dir=options.gramtools_build_dir,
        sample_name=options.sample_name,
    )
    adj.run()

