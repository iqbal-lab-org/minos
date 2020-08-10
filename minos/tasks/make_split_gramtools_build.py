from minos import vcf_chunker


def run(options):
    chunker = vcf_chunker.VcfChunker(
        options.outdir,
        vcf_infile=options.vcf_file,
        ref_fasta=options.ref_fasta,
        variants_per_split=options.variants_per_split,
        alleles_per_split=options.alleles_per_split,
        total_splits=options.total_splits,
        flank_length=200,
        gramtools_kmer_size=options.gramtools_kmer_size,
        threads=options.threads,
    )
    chunker.make_split_files()
