from cluster_vcf_records import variant_tracking
from minos import utils


def run(options):
    tracker = variant_tracking.VariantTracker(options.merge_dir, options.ref_fasta)
    tracker.cluster(
        options.outprefix,
        options.max_ref_len,
        max_alleles=options.max_alleles,
        cpus=options.cpus,
    )
    clustered_vcf = f"{options.outprefix}.vcf"
    utils.remove_vars_from_vcf_at_contig_ends(
        clustered_vcf, clustered_vcf, ref_fasta=options.ref_fasta
    )
