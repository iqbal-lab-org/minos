import filecmp
import shutil
import os
import pytest

from minos import adjudicator, genotype_confidence_simulator

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "adjudicator")


def _load_vcf_record_lines(infile):
    with open(infile) as f:
        return [x.rstrip() for x in f if not x.startswith("#")]


def test_get_gramtools_kmer_size():
    """test _get_gramtools_kmer_size"""
    build_dir = os.path.join(data_dir, "get_gramtools_kmer_size.build")
    assert adjudicator.Adjudicator._get_gramtools_kmer_size(build_dir, None) == 42
    assert adjudicator.Adjudicator._get_gramtools_kmer_size(build_dir, 20) == 42
    assert adjudicator.Adjudicator._get_gramtools_kmer_size(None, 20) == 20
    assert adjudicator.Adjudicator._get_gramtools_kmer_size(None, None) == 10


def test_run_clean_is_false():
    """test run when not cleaning up files afterwards"""
    # We're just testing that it doesn't crash.
    # Check the output files exist, but not their contents (except VCF)
    # First run using splitting of VCF file.
    # Then run without splitting.
    outdir = "tmp.adjudicator.noclean.out"
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    ref_fasta = os.path.join(data_dir, "run.ref.fa")
    reads_file = os.path.join(data_dir, "run.bwa.bam")
    vcf_files = [
        os.path.join(data_dir, x) for x in ["run.calls.1.vcf", "run.calls.2.vcf"]
    ]
    adj = adjudicator.Adjudicator(
        outdir,
        ref_fasta,
        [reads_file],
        vcf_files,
        variants_per_split=3,
        clean=False,
        gramtools_kmer_size=5,
        genotype_simulation_iterations=1000,
    )
    adj.run()
    assert os.path.exists(outdir)
    assert os.path.exists(adj.log_file)
    assert os.path.exists(adj.final_vcf)
    assert os.path.exists(adj.clustered_vcf)
    expect_vcf = os.path.join(data_dir, "run.expect_final.vcf")
    expect_lines = _load_vcf_record_lines(expect_vcf)
    got_lines = _load_vcf_record_lines(adj.final_vcf)
    assert got_lines == expect_lines

    # Clean up and then run without splitting
    shutil.rmtree(outdir)
    adj = adjudicator.Adjudicator(
        outdir,
        ref_fasta,
        [reads_file],
        vcf_files,
        clean=False,
        gramtools_kmer_size=5,
        genotype_simulation_iterations=1000,
    )
    adj.run()
    assert os.path.exists(outdir)
    assert os.path.exists(adj.log_file)
    assert os.path.exists(adj.final_vcf)
    assert os.path.exists(adj.gramtools_build_dir)
    assert os.path.exists(adj.gramtools_quasimap_dir)
    assert os.path.exists(adj.clustered_vcf)

    # Now we've run the adjudicator, we have a gramtools
    # build directory. Rerun, but this time use the build
    # directory, so we test the gramtools_build_dir option
    outdir2 = "tmp.adjudicator.out.2"
    gramtools_build_dir = adj.gramtools_build_dir
    if os.path.exists(outdir2):
        shutil.rmtree(outdir2)
    ref_fasta = os.path.join(data_dir, "run.ref.fa")
    reads_file = os.path.join(data_dir, "run.bwa.bam")
    # When gramtools build dir supplied, the Adjudicator assumes
    # one clsutered VCF file that matches the gramtools build run.
    # This is the clustered VCF made by the Adjudicator, so we
    # use that instead of the list of original VCF files
    vcf_files = [adj.clustered_vcf]
    adj = adjudicator.Adjudicator(
        outdir2,
        ref_fasta,
        [reads_file],
        vcf_files,
        gramtools_build_dir=gramtools_build_dir,
        clean=False,
        gramtools_kmer_size=5,
        genotype_simulation_iterations=1000,
    )
    adj.run()
    assert os.path.exists(outdir2)
    assert os.path.exists(adj.log_file)
    assert os.path.exists(adj.final_vcf)
    assert os.path.exists(adj.gramtools_build_dir)
    assert os.path.exists(adj.gramtools_quasimap_dir)
    assert os.path.exists(adj.clustered_vcf)
    assert not os.path.exists(os.path.join(outdir2, "gramtools.build"))
    shutil.rmtree(outdir)
    shutil.rmtree(outdir2)


def test_run_clean_is_true():
    """test run when we do clean files afterwards"""
    # We're just testing that it doesn't crash.
    # Check the output files exist, but not their contents.
    # First run using splitting of VCF file.
    # Then run without splitting.
    outdir = "tmp.adjudicator.clean.out"
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    ref_fasta = os.path.join(data_dir, "run.ref.fa")
    reads_file = os.path.join(data_dir, "run.bwa.bam")
    vcf_files = [
        os.path.join(data_dir, x) for x in ["run.calls.1.vcf", "run.calls.2.vcf"]
    ]
    adj = adjudicator.Adjudicator(
        outdir,
        ref_fasta,
        [reads_file],
        vcf_files,
        clean=True,
        gramtools_kmer_size=5,
        genotype_simulation_iterations=1000,
    )
    adj.run()
    assert os.path.exists(outdir)
    assert os.path.exists(os.path.join(outdir, "final.vcf"))
    shutil.rmtree(outdir)


def test_run_empty_vcf_input_files():
    """test run when input files have no variants"""
    outdir = "tmp.adjudicator.out"
    if os.path.exists(outdir):
        shutil.rmtree(outdir)

    ref_fasta = os.path.join(data_dir, "run.ref.fa")
    reads_file = os.path.join(data_dir, "run.bwa.bam")
    vcf_files = [
        os.path.join(data_dir, x)
        for x in ["run.calls.empty.1.vcf", "run.calls.empty.2.vcf"]
    ]
    adj = adjudicator.Adjudicator(
        outdir,
        ref_fasta,
        [reads_file],
        vcf_files,
        clean=False,
        gramtools_kmer_size=5,
        genotype_simulation_iterations=1000,
    )
    with pytest.raises(Exception):
        adj.run()
    assert os.path.exists(outdir)
    assert os.path.exists(adj.log_file)
    assert not os.path.exists(adj.final_vcf)
    assert not os.path.exists(adj.gramtools_build_dir)
    assert not os.path.exists(adj.gramtools_quasimap_dir)
    assert os.path.exists(adj.clustered_vcf)
    shutil.rmtree(outdir)


def test_add_gt_conf_percentile_and_filters_to_vcf_file():
    """test _add_gt_conf_percentile_and_filters_to_vcf_file"""
    original_file = os.path.join(data_dir, "add_gt_conf_percentile_to_vcf_file.in.vcf")
    tmp_file = "tmp.adjudicator.add_gt_conf_percentile_to_vcf_file.vcf"
    expect_file = os.path.join(
        data_dir, "add_gt_conf_percentile_to_vcf_file.expect.vcf"
    )
    shutil.copyfile(original_file, tmp_file)
    error_rate = 0.00026045894282438386
    simulations = genotype_confidence_simulator.GenotypeConfidenceSimulator(
        60, 100, error_rate, iterations=1000, call_hets=False,
    )
    simulations.run_simulations()
    adjudicator.Adjudicator._add_gt_conf_percentile_and_filters_to_vcf_file(
        tmp_file, simulations, min_dp=2, min_gcp=2.5, min_frs=0.9
    )
    assert filecmp.cmp(tmp_file, expect_file, shallow=False)
    os.unlink(tmp_file)


def test_0MeanDepth_stillRuns():
    """
    When mean depth is 0, we can get math errors: math.log(0) in genotype likelihood computation,
    and division by 0 in genotype confidence simulation.

    Note the former only actually occurs if there is a variant site with non-zero coverage;
    in this case, mean depth can get set to 0 due to rounding imprecision. This is tested in genotyper unit tests.
    """

    outdir = "tmp.adjudicator.out"
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    ref_fasta = os.path.join(data_dir, "run.ref.fa")
    reads_file = os.path.join(data_dir, "no_map_reads.fastq")
    vcf_files = [os.path.join(data_dir, "run.calls.1.vcf")]

    adj = adjudicator.Adjudicator(
        outdir, ref_fasta, [reads_file], vcf_files, clean=False, gramtools_kmer_size=5,
    )
    adj.run()
    # Make sure the coverage is 0
    assert adj.mean_depth == 0
    # And also the test passes if it raises no math related errors.
    shutil.rmtree(outdir)
