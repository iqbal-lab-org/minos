nextflow.enable.dsl=2

def check_parameter(params, parameter_name){
    if (!params[parameter_name]){
      error "You must specifiy a " + parameter_name
    } else {
      variable = params[parameter_name]
      return variable
    }
}


def build_vcf_from_gramtools_build_dir(build_dir) {
    build_vcf = file("${build_dir}/build.vcf")
    if (!build_vcf.exists()) {
        exit 1, "Error finding VCF in gramtools build dir ${build_dir} -- aborting"
    }
    return build_vcf
}


process parse_manifest {
    memory "1 GB"
    cpus params.parse_manifest_cpus

    input:
      file manifest_tsv
      file ref_fasta


    output:
      file "vcfs.fofn"
      file "sample_data.tsv"

    script:
    """
    #!/usr/bin/env python3
    from minos import regeno_helper
    regeno_helper.parse_manifest_file(
        "${manifest_tsv}",
        "vcfs.fofn",
        "sample_data.tsv",
        "${ref_fasta}",
        cpus=${params.parse_manifest_cpus},
        max_number_of_records=${params.max_variants_per_sample},
        max_ref_proportion=${params.max_genome_proportion_per_sample},
    )
    """
}


process vcf_merge {
    errorStrategy {task.attempt < 3 ? "retry" : "terminate"}
    maxRetries 3
    cpus params.vcf_merge_cpus

    input:
      file vcfs_fofn
      file ref_fasta


    output:
      path "outdir"

    script:
      """
      # Call faidx via pysam to avoid needing samtools installed
      python3 -c 'import pysam; pysam.faidx("${ref_fasta}")'
      minos vcf_merge --cpus ${params.vcf_merge_cpus} --mem_limit 2 ${vcfs_fofn} $ref_fasta outdir
      """
}


process vcf_cluster {
    cpus params.vcf_cluster_cpus
    errorStrategy {task.attempt < 3 ? "retry" : "terminate"}
    maxRetries 3

    input:
      file ref_fasta
      file merge_dir
      val max_ref_len
      val max_alleles

    output:
      file "cluster.vcf"

    script:
      """
      minos vcf_cluster --cpus ${params.vcf_cluster_cpus} --max_ref_len $max_ref_len --max_alleles $max_alleles $ref_fasta $merge_dir cluster
      """
}


process gramtools_build {
    errorStrategy {task.attempt < 3 ? "retry" : "terminate"}
    maxRetries 3
    cpus params.gramtools_build_cpus

    input:
      file vcf_in
      file ref_fasta
      val kmer
      val total_splits

    output:
      file "gramtools_build"

    script:
      """
      if [[ ${total_splits} -eq 1 ]]
      then
        gramtools build --vcf $vcf_in --gram-dir gramtools_build --kmer-size $kmer --reference $ref_fasta
      else
        minos make_split_gramtools_build --total_splits $total_splits --gramtools_kmer_size $kmer --threads ${params.gramtools_build_cpus} gramtools_build $vcf_in $ref_fasta
      fi
      """
}


process minos {
    errorStrategy {task.attempt < 3 ? 'retry' : "ignore"}
    maxRetries 3

    input:
      val sample_dict
      file build_dir
      file ref_fasta
      file vcf

    output:
      file "data_for_per_sample_dir.tsv"
      file "for_ivcfmerge.fofn"


    script:
      """
      minos adjudicate --sample_name ${sample_dict.name} --gramtools_build_dir $build_dir ${sample_dict.reads} minos_out $ref_fasta $vcf
      echo "${sample_dict.name}	\$PWD/minos_out" > data_for_per_sample_dir.tsv
      echo \$PWD/minos_out/debug.calls_with_zero_cov_alleles.vcf > for_ivcfmerge.fofn
      """
}


process make_per_sample_vcfs_dir {
    memory "1 GB"
    errorStrategy {task.attempt < 3 ? 'retry' : "terminate"}
    maxRetries 3

    input:
      file samples_tsv
      file manifest_tsv
      file outdir

    script:
    """
    #!/usr/bin/env python3
    from minos import regeno_helper
    regeno_helper.make_per_sample_vcfs_dir(
      "${samples_tsv}",
      "${outdir}",
      original_manifest="${manifest_tsv}",
      samples_per_dir=1000,
      cpus=1
    )
    """
}


process ivcf_merge_chunks {
    memory {1.GB * task.attempt}
    errorStrategy {task.attempt < 3 ? 'retry' : "terminate"}
    maxRetries 3

    input:
      val vcf_fofn

    output:
      file "merged.vcf"

    script:
    """
    ivcfmerge ${vcf_fofn} merged.vcf
    """
}


process ivcf_final_merge {
    memory {1.GB * task.attempt}
    errorStrategy {task.attempt < 3 ? 'retry' : "terminate"}
    maxRetries 3

    input:
      val file_list
      file outdir

    output:
      file "done_file"

    script:
    """
    cat <<"EOF" > files.txt
${file_list.join('\n')}
EOF
    ivcfmerge files.txt ${outdir}/merged.vcf
    touch done_file
    """
}


process distance_matrix {
    memory {5.GB * task.attempt}
    errorStrategy {task.attempt < 3 ? 'retry' : "terminate"}
    maxRetries 3

    input:
      file outdir
      file final_merge_done
      file mask_bed

    script:
    def mask_opt = mask_bed.name == "NO_FILE" ? "mask_bed_file=None" : "mask_bed_file='${mask_bed}'"
    """
    #!/usr/bin/env python3
    import os
    from minos import dist_matrix
    dist_matrix.distance_matrix_from_vcf_file(
        os.path.join("${outdir}", "merged.vcf"),
        os.path.join("${outdir}", "distance_matrix.txt"),
        ${mask_opt},
    )
    """
}


workflow make_vcf_for_gramtools {
    take:
        vcfs_fofn
        ref_fasta
        max_ref_len
        max_alleles

    main:
        vcf_merge(vcfs_fofn, ref_fasta)
        vcf_cluster(ref_fasta, vcf_merge.out, max_ref_len, max_alleles)

    emit:
        vcf_cluster.out
}


workflow {
    mask_bed = file(params.mask_bed_file)
    if (params.mask_bed_file != "NO_FILE" && !mask_bed.exists()) {
        exit 1, "Mask BED file not found: ${params.mask_bed_file} -- aborting"
    }

    parse_manifest(manifest, ref_fasta)
    if (params.vcf == "" && params.gramtools_build_dir == "") {
        gramtools_vcf = make_vcf_for_gramtools(
            parse_manifest.out[0],
            ref_fasta,
            params.max_ref_allele_len,
            params.max_alleles_per_site
        )
    }
    else if (params.vcf != "") {
        gramtools_vcf = file(params.vcf)
        if (!gramtools_vcf.exists()) {
            exit 1, "VCF file not found: ${params.vcf} -- aborting"
        }
    }

    if (params.gramtools_build_dir == "") {
        gramtools_build_dir = gramtools_build(
            gramtools_vcf,
            ref_fasta,
            params.gramtools_kmer,
            params.number_of_ref_chunks
        )
    }
    else {
        gramtools_build_dir = file(params.gramtools_build_dir)
        gramtools_vcf = build_vcf_from_gramtools_build_dir(params.gramtools_build_dir)
    }

    samples = parse_manifest.out[1].splitCsv(header: true, sep:'\t')
    minos(samples, gramtools_build.out, ref_fasta, gramtools_vcf)
    per_sample_tsv = minos.out[0].collectFile(name:"per_sample_data_for_final_outdir.tsv", newLine: false)
    make_per_sample_vcfs_dir(per_sample_tsv, manifest, outdir)
    ivcf_chunks = minos.out[1].collectFile().splitText(by: params.vcf_combine_batch_size, file: "ivcf_merge_chunks")
    ivcf_merge_chunks(ivcf_chunks)
    ivcf_final_merge(ivcf_merge_chunks.out.collect(), outdir)
    if (params.make_distance_matrix) {
        distance_matrix(outdir, ivcf_final_merge.out, mask_bed)
    }
}


params.help = false
params.ref_fasta = ""
params.root_out = ""
params.vcf = ""
params.gramtools_build_dir = ""
params.mask_bed_file = "NO_FILE"
params.gramtools_kmer = 7
params.parse_manifest_cpus = 5
params.max_variants_per_sample = 5000
params.max_genome_proportion_per_sample = 0.05
params.vcf_merge_cpus = 5
params.vcf_cluster_cpus = 5
params.gramtools_build_cpus = 1
params.max_ref_allele_len = 50
params.max_alleles_per_site = 500
params.number_of_ref_chunks = 10
params.vcf_combine_batch_size = 100
params.make_distance_matrix = false

if (params.help){
    log.info"""
        minos regenotyping pipeline: genotypes all samples at the same
        positions.

        Usage: nextflow run regenotype.nf <arguments>

        Required arguments:
          --ref_fasta FILENAME  FASTA file of reference genome
          --manifest FILENAME   TSV file containing sample information
          --outdir DIRNAME      Output directory (must not already exist)

        Please see template config file for other options.
    """.stripIndent()

    exit 0
}



if (params.vcf != "" && params.gramtools_build_dir != "") {
    exit 1, "The options --vcf and --gramtools_build_dir are mutually exclusive -- aborting"
}

ref_fasta = file(check_parameter(params, "ref_fasta"))
if (!ref_fasta.exists()) {
    exit 1, "Reference FASTA file not found: ${params.ref_fasta} -- aborting"
}

manifest = file(check_parameter(params, "manifest"))
if (!manifest.exists()) {
    exit 1, "Manifest TSV file not found: ${params.manifest} -- aborting"
}

outdir = file(check_parameter(params, "outdir"))
if (!outdir.mkdirs()) {
    error "Error making output directory " + outdir
}

