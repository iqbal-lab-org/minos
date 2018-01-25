from minos import mapping_based_verifier

def run(options):
    verifier = mapping_based_verifier.MappingBasedVerifier(
        options.vcf_file,
        options.vcf_ref,
        options.truth_ref,
        options.outprefix,
        flank_length=options.flank_length
    )
    verifier.run()

