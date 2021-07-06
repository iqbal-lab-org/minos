import datetime
import shutil
import os
import pytest

from cluster_vcf_records import vcf_file_read, vcf_record

from minos import genotyper, gramtools
from minos import __version__ as minos_version

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "gramtools")


def test_build_json_file_is_good():
    """test _build_json_file_is_good"""
    build_file = os.path.join(data_dir, "build_json_file_is_good.good.json")
    assert gramtools._build_json_file_is_good(build_file)
    build_file = os.path.join(data_dir, "build_json_file_is_good.bad.1.json")
    assert not gramtools._build_json_file_is_good(build_file)
    build_file = os.path.join(data_dir, "build_json_file_is_good.bad.2.json")
    assert not gramtools._build_json_file_is_good(build_file)
    build_file = os.path.join(
        data_dir, "build_json_file_is_good.bad.2.json_thisfiledoesnotexist"
    )
    assert not gramtools._build_json_file_is_good(build_file)


def test_run_gramtools_build():
    """test run_gramtools_build"""
    tmp_out_build = "tmp.run_gramtools.out.build"
    if os.path.exists(tmp_out_build):
        shutil.rmtree(tmp_out_build)
    vcf_file = os.path.join(data_dir, "run_gramtools.calls.vcf")
    ref_file = os.path.join(data_dir, "run_gramtools.ref.fa")
    gramtools.run_gramtools_build(tmp_out_build, vcf_file, ref_file, kmer_size=5)
    assert os.path.exists(tmp_out_build)
    shutil.rmtree(tmp_out_build)


def test_run_gramtools():
    """test run_gramtools"""
    tmp_out_build = "tmp.run_gramtools.out.build"
    tmp_out_quasimap = "tmp.run_gramtools.out.quasimap"
    if os.path.exists(tmp_out_build):
        shutil.rmtree(tmp_out_build)
    if os.path.exists(tmp_out_quasimap):
        shutil.rmtree(tmp_out_quasimap)
    vcf_file = os.path.join(data_dir, "run_gramtools.calls.vcf")
    ref_file = os.path.join(data_dir, "run_gramtools.ref.fa")
    reads_file = os.path.join(data_dir, "run_gramtools.reads.fq")
    build_report, quasimap_report = gramtools.run_gramtools(
        tmp_out_build, tmp_out_quasimap, vcf_file, ref_file, reads_file, kmer_size=5,
    )
    # We're trusing gramtools output for this test. The point here is to check
    # that gramtools can run. Parsing its output is checked elsewhere.
    assert os.path.exists(tmp_out_build)
    assert os.path.exists(tmp_out_quasimap)
    assert os.path.exists(
        os.path.join(tmp_out_quasimap, "quasimap_outputs", "allele_base_coverage.json")
    )
    assert os.path.exists(
        os.path.join(
            tmp_out_quasimap, "quasimap_outputs", "grouped_allele_counts_coverage.json"
        )
    )
    assert "gramtools_build" in build_report
    assert "gramtools_cpp_quasimap" in quasimap_report
    shutil.rmtree(tmp_out_build)
    shutil.rmtree(tmp_out_quasimap)


def test_run_gramtools_fails():
    """test run_gramtools when fails"""
    # Don't trust error code. Instead, we check
    # that gramtools wrote the quesimap files we expected, as this
    # is a good proxy for success. One way to stop these files
    # from being written is to have no variants in the input VCF,
    # so that's what we do here
    tmp_out_build = "tmp.run_gramtools.fail.out.build"
    tmp_out_quasimap = "tmp.run_gramtools.fail.out.quasimap"
    if os.path.exists(tmp_out_build):
        shutil.rmtree(tmp_out_build)
    if os.path.exists(tmp_out_quasimap):
        shutil.rmtree(tmp_out_quasimap)
    vcf_file = os.path.join(data_dir, "run_gramtools.empty.vcf")
    ref_file = os.path.join(data_dir, "run_gramtools.ref.fa")
    reads_file = os.path.join(data_dir, "run_gramtools.reads.fq")
    with pytest.raises(Exception):
        gramtools.run_gramtools(
            tmp_out_build,
            tmp_out_quasimap,
            vcf_file,
            ref_file,
            reads_file,
            kmer_size=5,
        )
    shutil.rmtree(tmp_out_build)


def test_run_gramtools_two_reads_files():
    """test run_gramtools"""
    tmp_out_build = "tmp.run_gramtools.2files.out.build"
    tmp_out_quasimap = "tmp.run_gramtools.2files.out.quasimap"
    if os.path.exists(tmp_out_build):
        shutil.rmtree(tmp_out_build)
    if os.path.exists(tmp_out_quasimap):
        shutil.rmtree(tmp_out_quasimap)
    vcf_file = os.path.join(data_dir, "run_gramtools.calls.vcf")
    ref_file = os.path.join(data_dir, "run_gramtools.ref.fa")
    reads_file1 = os.path.join(data_dir, "run_gramtools.reads_1.fq")
    reads_file2 = os.path.join(data_dir, "run_gramtools.reads_2.fq")
    gramtools.run_gramtools(
        tmp_out_build,
        tmp_out_quasimap,
        vcf_file,
        ref_file,
        [reads_file1, reads_file2],
        kmer_size=5,
    )
    # We're trusing gramtools output for this test. The point here is to check
    # that gramtools can run. Parsing its output is checked elsewhere.
    assert os.path.exists(tmp_out_build)
    assert os.path.exists(tmp_out_quasimap)
    assert os.path.exists(
        os.path.join(tmp_out_quasimap, "quasimap_outputs", "allele_base_coverage.json")
    )
    assert os.path.exists(
        os.path.join(
            tmp_out_quasimap, "quasimap_outputs", "grouped_allele_counts_coverage.json"
        )
    )
    shutil.rmtree(tmp_out_build)
    shutil.rmtree(tmp_out_quasimap)


def test_grouped_allele_counts_coverage_json_to_cov_list():
    infile = os.path.join(
        data_dir, "grouped_allele_counts_coverage_json_to_cov_list.json"
    )
    gramtools.grouped_allele_counts_coverage_json_to_cov_list(infile) == [11, 10]


def test_coverage_list_from_allele_coverage():
    """test _coverage_list_from_allele_coverage"""
    all_allele_coverage = [
        ({"0": 1, "1": 15}, [[1], [15]]),
        ({"0": 15, "1": 13}, [[14, 2], [13]]),
    ]
    expect = [16, None]
    got = gramtools._coverage_list_from_allele_coverage(all_allele_coverage)
    assert expect == got
    expect = [15, 14]
    got = gramtools._coverage_list_from_allele_coverage(
        all_allele_coverage, use_indels=True
    )
    assert expect == got


def test_load_gramtools_vcf_and_allele_coverage_files():
    """test load_gramtools_vcf_and_allele_coverage_files"""
    vcf_file = os.path.join(data_dir, "load_gramtools_vcf_and_allele_coverage.vcf")
    quasimap_dir = os.path.join(
        data_dir, "load_gramtools_vcf_and_allele_coverage_files.quasimap"
    )
    (
        got_mean_depth,
        got_depth_variance,
        got_vcf_header,
        got_vcf_records,
        got_allele_coverage,
        got_allele_groups,
    ) = gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir)

    expected_header, expected_vcf_records = vcf_file_read.vcf_file_to_list(vcf_file)
    assert got_vcf_header == expected_header
    assert got_vcf_records == expected_vcf_records
    assert got_mean_depth == 6.5
    assert got_depth_variance == 24.5

    # now test bad files cause error to be raised
    vcf_file = os.path.join(
        data_dir, "load_gramtools_vcf_and_allele_coverage.short.vcf"
    )
    with pytest.raises(Exception):
        gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir)

    vcf_file = os.path.join(data_dir, "load_gramtools_vcf_and_allele_coverage.long.vcf")
    with pytest.raises(Exception):
        gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir)

    vcf_file = os.path.join(
        data_dir, "load_gramtools_vcf_and_allele_coverage.bad_allele_count.vcf"
    )
    with pytest.raises(Exception):
        gramtools.load_gramtools_vcf_and_allele_coverage_files(vcf_file, quasimap_dir)


def test_update_vcf_record_using_gramtools_allele_depths_homozygous():
    """test update_using_gramtools_allele_depths homozygous"""
    record = vcf_record.VcfRecord(
        "ref\t4\t.\tT\tTC,G\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0"
    )
    allele_combination_cov = {"1": 1, "2": 80}
    allele_groups_dict = {"1": {0}, "2": {2}, "3": {1, 2}}
    allele_per_base_cov = [[1], [0, 0], [80]]
    expected = vcf_record.VcfRecord(
        "ref\t4\t.\tT\tTC,G\t.\t.\t.\tGT:DP:ALLELE_DP:FRS:COV_TOTAL:COV:GT_CONF\t2/2:80:1,0,80:0.9877:81:1,0,80:646.06"
    )
    mean_depth = 85
    depth_variance = 200
    error_rate = 0.001
    gtyper = genotyper.Genotyper(
        mean_depth, depth_variance, error_rate, call_hets=False,
    )
    got_filtered = gramtools.update_vcf_record_using_gramtools_allele_depths(
        record, gtyper, allele_combination_cov, allele_per_base_cov, allele_groups_dict,
    )
    assert record == expected
    expected_filtered = vcf_record.VcfRecord(
        "ref\t4\t.\tT\tG\t.\t.\t.\tGT:DP:ALLELE_DP:FRS:COV_TOTAL:COV:GT_CONF\t1/1:80:1,80:0.9877:81:1,80:646.06"
    )
    assert got_filtered == expected_filtered


def test_update_vcf_record_using_gramtools_allele_depths_not_callable():
    """test update_using_gramtools_allele_depths not callable"""
    record = vcf_record.VcfRecord(
        "ref\t4\t.\tT\tG,TC\t228\t.\tINDEL;IDV=54;IMF=0.885246;DP=61;VDB=7.33028e-19;SGB=-0.693147;MQSB=0.9725;MQ0F=0;AC=2;AN=2;DP4=0,0,23,31;MQ=57\tGT:PL\t1/1:255,163,0"
    )
    allele_combination_cov = {"1": 0, "2": 0}
    allele_groups_dict = {"1": {0}, "2": {0}}
    allele_per_base_cov = [[0], [0], [0, 0]]
    expected = vcf_record.VcfRecord(
        "ref\t4\t.\tT\tG,TC\t.\t.\t.\tGT:DP:ALLELE_DP:FRS:COV_TOTAL:COV:GT_CONF\t./.:.:0,0,0:.:0:0,0,0:0.0"
    )
    mean_depth = 85
    error_rate = 0.001
    depth_variance = 200
    gtyper = genotyper.Genotyper(
        mean_depth, depth_variance, error_rate, call_hets=False,
    )
    got_filtered = gramtools.update_vcf_record_using_gramtools_allele_depths(
        record, gtyper, allele_combination_cov, allele_per_base_cov, allele_groups_dict,
    )
    assert record == expected
    assert got_filtered == expected


def test_write_vcf_annotated_using_coverage_from_gramtools():
    """test write_vcf_annotated_using_coverage_from_gramtools"""
    vcf_file_in = os.path.join(
        data_dir, "write_vcf_annotated_using_coverage_from_gramtools.in.vcf"
    )
    quasimap_dir = os.path.join(
        data_dir, "write_vcf_annotated_using_coverage_from_gramtools.quasimap"
    )
    (
        mean_depth,
        depth_variance,
        vcf_header,
        vcf_records,
        allele_coverage,
        allele_groups,
    ) = gramtools.load_gramtools_vcf_and_allele_coverage_files(
        vcf_file_in, quasimap_dir
    )
    tmp_outfile = "tmp.gramtools.write_vcf_annotated_using_coverage_from_gramtools.vcf"
    tmp_outfile_filtered = tmp_outfile + ".filter.vcf"
    error_rate = 0.001
    gramtools.write_vcf_annotated_using_coverage_from_gramtools(
        mean_depth,
        depth_variance,
        vcf_records,
        allele_coverage,
        allele_groups,
        error_rate,
        tmp_outfile,
        sample_name="sample_42",
        filtered_outfile=tmp_outfile_filtered,
        ref_seq_lengths={"ref": "10000"},
        call_hets=False,
    )
    expected_vcf = os.path.join(
        data_dir, "write_vcf_annotated_using_coverage_from_gramtools.out.vcf"
    )
    expected_vcf_filtered = os.path.join(
        data_dir,
        "write_vcf_annotated_using_coverage_from_gramtools.out.vcf.filter.vcf",
    )
    # Today's date and the verison of minos get added to the header.
    # We'll have to take account
    # of those by fixing what we get from the expected file
    def check_vcfs(expected_vcf, got_vcf):
        expected_header, expected_vcf_records = vcf_file_read.vcf_file_to_list(
            expected_vcf
        )
        got_header, got_vcf_records = vcf_file_read.vcf_file_to_list(got_vcf)
        for i in range(len(expected_header)):
            if expected_header[i].startswith("##fileDate="):
                expected_header[i] = "##fileDate=" + str(datetime.date.today())
            elif expected_header[i].startswith("##source=minos"):
                expected_header[i] = "##source=minos, version " + minos_version

        assert got_header == expected_header
        assert got_vcf_records == expected_vcf_records

    check_vcfs(expected_vcf, tmp_outfile)
    check_vcfs(expected_vcf_filtered, tmp_outfile_filtered)
    os.unlink(tmp_outfile)
    os.unlink(tmp_outfile_filtered)


def test_load_allele_files():
    """test load_allele_files"""
    expected_counts_list = [
        ({"0": 10, "1": 3, "2": 9}, [[1, 1, 1], [1, 1, 0]]),
        ({"0": 30, "3": 2, "4": 1}, [[0, 0, 1], [2, 2, 0], [2, 2, 0, 1, 3]]),
    ]

    expected_groups_dict = {"0": {0}, "1": {1}, "2": {0, 1}, "3": {2}, "4": {1, 2}}

    allele_base_counts_file = os.path.join(
        data_dir, "load_allele_files.allele_base_counts.json"
    )
    grouped_allele_counts_file = os.path.join(
        data_dir, "load_allele_files.grouped_allele_counts.json"
    )
    got_counts_list, got_groups_dict = gramtools.load_allele_files(
        allele_base_counts_file, grouped_allele_counts_file
    )
    assert got_counts_list == expected_counts_list
    assert got_groups_dict == expected_groups_dict
