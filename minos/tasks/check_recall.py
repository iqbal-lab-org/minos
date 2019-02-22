from minos import evaluate_recall

def run(options):
    verifier = evaluate_recall.EvaluateRecall(
        options.truth_vcf_file,
        options.truth_vcf_ref,
        options.query_vcf_file,
        options.query_vcf_ref,
        options.outprefix,
        flank_length=options.flank_length,
        merge_length=options.variant_merge_length,
        filter_and_cluster_vcf=not options.no_filter_cluster,
        discard_ref_calls=not options.include_ref_calls,
        allow_flank_mismatches=options.allow_flank_mismatches,
        exclude_regions_bed_file=options.exclude_bed,
        max_soft_clipped=options.max_soft_clipped,
    )
    verifier.run()

