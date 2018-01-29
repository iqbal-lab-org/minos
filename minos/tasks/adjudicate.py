from minos import adjudicator

def run(options):
    adj = adjudicator.Adjudicator(
        options.outdir,
        options.ref_fasta,
        options.reads_file,
        options.vcf_files,
        max_read_length=options.max_read_length,
        read_error_rate=options.read_error_rate,
        overwrite_outdir=options.force,
        max_snps_per_cluster=options.max_snps_per_cluster,
        max_ref_len=options.max_ref_len,
    )
    adj.run()

