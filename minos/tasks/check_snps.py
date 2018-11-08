from minos import dnadiff_mapping_based_verifier

def run(options):
    verifier = dnadiff_mapping_based_verifier.DnadiffMappingBasedVerifier(
        options.dnadiff_snps_file,
        options.dnadiff_ref,
        options.dnadiff_query,
        options.vcf_file1,
        options.vcf_file2,
        options.vcf_ref,
        options.outprefix,
        flank_length=options.flank_length,
        merge_length=options.variant_merge_length,
        filter_and_cluster_vcf=not options.no_filter_cluster,
        discard_ref_calls=not options.include_ref_calls,
        allow_flank_mismatches=options.allow_flank_mismatches,
    )
    verifier.run()

