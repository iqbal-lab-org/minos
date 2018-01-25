from cluster_vcf_records import vcf_clusterer

def run(options):
    clusterer = vcf_clusterer.VcfClusterer(
        options.vcf_files,
        options.ref_fasta,
        options.outfile,
        max_distance_between_variants=options.max_var_dist,
        max_REF_len=options.max_ref_len,
        max_snps_per_cluster=options.max_snps_per_cluster,
    )
    clusterer.run()

