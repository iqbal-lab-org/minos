from minos import mapping_based_verifier

def run(options):
    verifier = mapping_based_verifier.MappingBasedVerifier(
        options.vcf_file,
        options.vcf_ref,
        options.truth_ref,
        options.outprefix,
        flank_length=options.flank_length,
        expected_variants_vcf=options.expected_variants_vcf,
        run_dnadiff=options.run_dnadiff,
        filter_and_cluster_vcf=not options.no_filter_cluster,
    )
    verifier.run()

