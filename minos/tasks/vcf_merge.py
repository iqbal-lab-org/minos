from cluster_vcf_records import variant_tracking


def run(options):
    with open(options.vcf_fofn) as f:
        vcf_files = [line.rstrip() for line in f]
    tracker = variant_tracking.VariantTracker(options.outdir, options.ref_fasta)
    tracker.merge_vcf_files(
        vcf_files,
        temp_dir=options.temp_dir,
        cpus=options.cpus,
        mem_limit=options.mem_limit,
        force=options.force,
        sample_limit=options.sample_limit,
    )
