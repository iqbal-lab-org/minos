import filecmp
import shutil
import os
import unittest

import cluster_vcf_records

from minos import vcf_chunker

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "vcf_chunker")


class TestVcfChunker(unittest.TestCase):
    def test_total_variants_and_alleles_in_vcf_dict(self):
        """test _total_variants_and_alleles_in_vcf_dict"""

        class FakeVcf:
            def __init__(self, alt):
                self.ALT = alt

        test_dict = {
            "chrom1": [FakeVcf("123"), FakeVcf("1"), FakeVcf("123456789")],
            "chrom2": [FakeVcf("12"), FakeVcf("1234")],
        }
        expect_variants = 5
        expect_alleles = 24
        (
            got_variants,
            got_alleles,
        ) = vcf_chunker.VcfChunker._total_variants_and_alleles_in_vcf_dict(test_dict)
        self.assertEqual(expect_variants, got_variants)
        self.assertEqual(expect_alleles, got_alleles)

    def test_chunk_end_indexes_from_vcf_record_list(self):
        """test _chunk_end_indexes_from_vcf_record_list"""
        record_list = [
            cluster_vcf_records.vcf_record.VcfRecord("ref\t1\t.\tA\tG\t.\t.\t.\t."),
            cluster_vcf_records.vcf_record.VcfRecord(
                "ref\t2\t.\tC\tT,A,G,TA\t.\t.\t.\t."
            ),
            cluster_vcf_records.vcf_record.VcfRecord("ref\t3\t.\tT\tA,C\t.\t.\t.\t."),
            cluster_vcf_records.vcf_record.VcfRecord(
                "ref\t5\t.\tAGAGTCACGTA\tG\t.\t.\t.\t."
            ),
            cluster_vcf_records.vcf_record.VcfRecord("ref\t18\t.\tA\tG\t.\t.\t.\t."),
            cluster_vcf_records.vcf_record.VcfRecord("ref\t21\t.\tG\tT\t.\t.\t.\t."),
        ]

        self.assertEqual(
            (0, 0, 1),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=1
            ),
        )
        self.assertEqual(
            (0, 0, 1),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=2
            ),
        )
        self.assertEqual(
            (0, 0, 1),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=3
            ),
        )
        self.assertEqual(
            (0, 0, 1),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=4
            ),
        )
        self.assertEqual(
            (0, 0, 1),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=5
            ),
        )
        self.assertEqual(
            (0, 0, 1),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=6
            ),
        )
        self.assertEqual(
            (0, 1, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=7
            ),
        )
        self.assertEqual(
            (0, 1, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=8
            ),
        )
        self.assertEqual(
            (0, 1, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=9
            ),
        )
        self.assertEqual(
            (0, 2, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=10
            ),
        )
        self.assertEqual(
            (0, 2, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=11
            ),
        )
        self.assertEqual(
            (0, 3, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_alleles=12
            ),
        )

        self.assertEqual(
            (0, 0, 1),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 1, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_sites=2
            ),
        )
        self.assertEqual(
            (0, 2, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_sites=3
            ),
        )
        self.assertEqual(
            (0, 3, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_sites=4
            ),
        )
        self.assertEqual(
            (0, 4, 4),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_sites=5
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_sites=6
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_sites=7
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 1, total_sites=8
            ),
        )

        self.assertEqual(
            (0, 0, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 2, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 1, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 2, total_sites=2
            ),
        )
        self.assertEqual(
            (0, 2, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 2, total_sites=3
            ),
        )
        self.assertEqual(
            (0, 3, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 2, total_sites=4
            ),
        )
        self.assertEqual(
            (0, 4, 4),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 2, total_sites=5
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 2, total_sites=6
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 2, total_sites=7
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 2, total_sites=8
            ),
        )

        self.assertEqual(
            (0, 0, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 3, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 1, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 3, total_sites=2
            ),
        )
        self.assertEqual(
            (0, 2, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 3, total_sites=3
            ),
        )
        self.assertEqual(
            (0, 3, 4),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 3, total_sites=4
            ),
        )
        self.assertEqual(
            (0, 4, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 3, total_sites=5
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 3, total_sites=6
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 3, total_sites=7
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 3, total_sites=8
            ),
        )

        self.assertEqual(
            (0, 0, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 4, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 1, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 4, total_sites=2
            ),
        )
        self.assertEqual(
            (0, 2, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 4, total_sites=3
            ),
        )
        self.assertEqual(
            (0, 3, 4),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 4, total_sites=4
            ),
        )
        self.assertEqual(
            (0, 4, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 4, total_sites=5
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 4, total_sites=6
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 4, total_sites=7
            ),
        )

        self.assertEqual(
            (0, 1, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 1, 1, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 1, 2),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 1, 2, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 1, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 1, 3, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 1, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 1, 15, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 1, 4),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 1, 16, total_sites=1
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 1, 1, total_sites=6
            ),
        )

        self.assertEqual(
            (4, 4, 4),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 4, 1, total_sites=1
            ),
        )
        self.assertEqual(
            (4, 4, 4),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 4, 2, total_sites=1
            ),
        )
        self.assertEqual(
            (3, 4, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 4, 3, total_sites=1
            ),
        )
        self.assertEqual(
            (4, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 4, 1, total_sites=2
            ),
        )

        self.assertEqual(
            (5, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 1, total_sites=1
            ),
        )
        self.assertEqual(
            (5, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 1, total_sites=2
            ),
        )
        self.assertEqual(
            (5, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 2, total_sites=2
            ),
        )
        self.assertEqual(
            (4, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 3, total_sites=2
            ),
        )
        self.assertEqual(
            (4, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 4, total_sites=2
            ),
        )
        self.assertEqual(
            (4, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 5, total_sites=2
            ),
        )
        self.assertEqual(
            (3, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 6, total_sites=2
            ),
        )
        self.assertEqual(
            (3, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 7, total_sites=2
            ),
        )
        self.assertEqual(
            (3, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 17, total_sites=2
            ),
        )
        self.assertEqual(
            (2, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 18, total_sites=2
            ),
        )
        self.assertEqual(
            (1, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 19, total_sites=2
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 20, total_sites=2
            ),
        )
        self.assertEqual(
            (0, 5, 5),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 5, 21, total_sites=2
            ),
        )

        # These records caused minos error because variant at 800
        # was included in the last split file, but the use_end_index was at
        # position of the variant at 610. So the one at 800 was not getting used.
        record_list = [
            cluster_vcf_records.vcf_record.VcfRecord("ref\t75\t.\tA\tG\t.\t.\t.\t."),
            cluster_vcf_records.vcf_record.VcfRecord("ref\t150\t.\tG\tA,T\t.\t.\t.\t."),
            cluster_vcf_records.vcf_record.VcfRecord("ref\t450\t.\tT\tC\t.\t.\t.\t."),
            cluster_vcf_records.vcf_record.VcfRecord("ref\t610\t.\tA\tG\t.\t.\t.\t."),
            cluster_vcf_records.vcf_record.VcfRecord("ref\t800\t.\tC\tCA\t.\t.\t.\t."),
        ]

        self.assertEqual(
            (0, 1, 1),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 0, 100, total_sites=2
            ),
        )
        self.assertEqual(
            (2, 3, 3),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 2, 100, total_sites=2
            ),
        )
        self.assertEqual(
            (4, 4, 4),
            vcf_chunker.VcfChunker._chunk_end_indexes_from_vcf_record_list(
                record_list, 4, 100, total_sites=2
            ),
        )

    def test_make_split_files(self):
        """test make_split_files"""
        infile = os.path.join(data_dir, "make_split_files.in.vcf")
        tmp_out = "tmp.vcf_chunker.make_split_files"
        ref_fa = os.path.join(data_dir, "make_split_files.in.ref.fa")
        if os.path.exists(tmp_out):
            shutil.rmtree(tmp_out)

        vcf1 = cluster_vcf_records.vcf_record.VcfRecord(
            "ref1\t1\t.\tG\tT\t.\tPASS\t.\t.\t."
        )
        vcf2 = cluster_vcf_records.vcf_record.VcfRecord(
            "ref1\t2\t.\tC\tT\t.\tPASS\t.\t.\t."
        )
        vcf3 = cluster_vcf_records.vcf_record.VcfRecord(
            "ref1\t3\t.\tT\tA\t.\tPASS\t.\t.\t."
        )
        vcf4 = cluster_vcf_records.vcf_record.VcfRecord(
            "ref1\t5\t.\tAGAGTCACGTA\tG\t.\tPASS\t.\t.\t."
        )
        vcf5 = cluster_vcf_records.vcf_record.VcfRecord(
            "ref1\t18\t.\tA\tG\t.\tPASS\t.\t.\t."
        )
        vcf6 = cluster_vcf_records.vcf_record.VcfRecord(
            "ref1\t21\t.\tG\tT\t.\tPASS\t.\t.\t."
        )
        vcf7 = cluster_vcf_records.vcf_record.VcfRecord(
            "ref2\t42\t.\tC\tG\t.\tPASS\t.\t.\t."
        )
        header_lines = [
            "##header1",
            "##header2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_name",
        ]

        chunker = vcf_chunker.VcfChunker(
            tmp_out,
            vcf_infile=infile,
            ref_fasta=ref_fa,
            variants_per_split=2,
            flank_length=1,
            gramtools_kmer_size=5,
        )
        chunker.make_split_files()
        self.assertTrue(os.path.exists(chunker.metadata_pickle))

        got_header, got_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            os.path.join(tmp_out, "split.0.in.vcf")
        )
        self.assertEqual(header_lines, got_header)
        self.assertEqual([vcf1, vcf2, vcf3], got_records)

        got_header, got_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            os.path.join(tmp_out, "split.1.in.vcf")
        )
        self.assertEqual(header_lines, got_header)
        self.assertEqual([vcf2, vcf3, vcf4], got_records)

        got_header, got_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            os.path.join(tmp_out, "split.2.in.vcf")
        )
        self.assertEqual(header_lines, got_header)
        self.assertEqual([vcf5, vcf6], got_records)

        got_header, got_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            os.path.join(tmp_out, "split.3.in.vcf")
        )
        self.assertEqual(header_lines, got_header)
        self.assertEqual([vcf7], got_records)

        self.assertFalse(os.path.exists(os.path.join(tmp_out, "split.4.in.vcf")))
        shutil.rmtree(tmp_out)

        chunker = vcf_chunker.VcfChunker(
            tmp_out,
            vcf_infile=infile,
            ref_fasta=ref_fa,
            variants_per_split=4,
            flank_length=3,
            gramtools_kmer_size=5,
        )
        chunker.make_split_files()
        self.assertTrue(os.path.exists(chunker.metadata_pickle))

        got_header, got_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            os.path.join(tmp_out, "split.0.in.vcf")
        )
        self.assertEqual(header_lines, got_header)
        self.assertEqual([vcf1, vcf2, vcf3, vcf4, vcf5], got_records)

        got_header, got_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            os.path.join(tmp_out, "split.1.in.vcf")
        )
        self.assertEqual(header_lines, got_header)
        self.assertEqual([vcf4, vcf5, vcf6], got_records)

        got_header, got_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            os.path.join(tmp_out, "split.2.in.vcf")
        )
        self.assertEqual(header_lines, got_header)
        self.assertEqual([vcf7], got_records)

        self.assertFalse(os.path.exists(os.path.join(tmp_out, "split.3.in.vcf")))

        chunker2 = vcf_chunker.VcfChunker(tmp_out, gramtools_kmer_size=5)
        self.assertEqual(chunker.vcf_infile, chunker2.vcf_infile)
        self.assertEqual(chunker.ref_fasta, chunker2.ref_fasta)
        self.assertEqual(chunker.variants_per_split, chunker2.variants_per_split)
        self.assertEqual(chunker.total_splits, chunker2.total_splits)
        self.assertEqual(chunker.flank_length, chunker2.flank_length)
        self.assertEqual(chunker.gramtools_kmer_size, chunker2.gramtools_kmer_size)
        self.assertEqual(chunker.total_split_files, chunker2.total_split_files)
        self.assertEqual(chunker.vcf_split_files, chunker2.vcf_split_files)
        shutil.rmtree(tmp_out)

    def test_make_split_files_2(self):
        """test make_split_files with different input from previous test"""
        # These records cause a minos bug. Last record was not being used
        # when merging because the index was wrong.
        # They are test data from multi_sample_pipeline tests
        infile = os.path.join(data_dir, "make_split_files2.in.vcf")
        tmp_out = "tmp.vcf_chunker.make_split_files2"
        ref_fa = os.path.join(data_dir, "make_split_files2.in.ref.fa")
        if os.path.exists(tmp_out):
            shutil.rmtree(tmp_out)

        chunker = vcf_chunker.VcfChunker(
            tmp_out,
            vcf_infile=infile,
            ref_fasta=ref_fa,
            variants_per_split=2,
            flank_length=200,
            gramtools_kmer_size=5,
        )
        chunker.make_split_files()
        self.assertTrue(os.path.exists(chunker.metadata_pickle))
        chunker2 = vcf_chunker.VcfChunker(tmp_out, gramtools_kmer_size=5)
        self.assertEqual(1, len(chunker2.vcf_split_files))
        self.assertEqual(3, len(chunker2.vcf_split_files["ref.0"]))
        self.assertEqual(4, chunker2.vcf_split_files["ref.0"][-1].use_end_index)
        shutil.rmtree(tmp_out)

        # Test with two threads
        chunker = vcf_chunker.VcfChunker(
            tmp_out,
            vcf_infile=infile,
            ref_fasta=ref_fa,
            variants_per_split=2,
            flank_length=200,
            threads=2,
            gramtools_kmer_size=5,
        )
        chunker.make_split_files()
        self.assertTrue(os.path.exists(chunker.metadata_pickle))
        chunker2 = vcf_chunker.VcfChunker(tmp_out, gramtools_kmer_size=5)
        self.assertEqual(1, len(chunker2.vcf_split_files))
        self.assertEqual(3, len(chunker2.vcf_split_files["ref.0"]))
        self.assertEqual(4, chunker2.vcf_split_files["ref.0"][-1].use_end_index)
        shutil.rmtree(tmp_out)

    def test_merge_files(self):
        """test merge_files"""
        vcf_to_split = os.path.join(data_dir, "merge_files.in.vcf")
        ref_fasta = os.path.join(data_dir, "merge_files.in.ref.fa")
        tmp_outdir = "tmp.vcf_chunker.merge_files"
        chunker = vcf_chunker.VcfChunker(
            tmp_outdir,
            vcf_infile=vcf_to_split,
            ref_fasta=ref_fasta,
            variants_per_split=4,
            flank_length=3,
            gramtools_kmer_size=5,
        )
        chunker.make_split_files()
        to_merge = {}
        for ref, split_list in chunker.vcf_split_files.items():
            to_merge[ref] = [x.filename for x in split_list]
        tmp_vcf_out = "tmp.vcf_chunker.merge_files.out.vcf"
        chunker.merge_files(to_merge, tmp_vcf_out)
        self.assertTrue(filecmp.cmp(vcf_to_split, tmp_vcf_out, shallow=False))
        os.unlink(tmp_vcf_out)
        shutil.rmtree(tmp_outdir)
