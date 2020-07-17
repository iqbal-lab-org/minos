from cluster_vcf_records import variant_tracking


def run(options):
    tracker = variant_tracking.VariantTracker(options.merge_dir, options.ref_fasta)
    tracker.cluster(
        options.outprefix,
        options.max_ref_len,
        max_alleles=options.max_alleles,
        cpus=options.cpus,
    )
