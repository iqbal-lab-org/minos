import pyfastaq
import pymummer

from minos import utils

class Error (Exception): pass


class Dnadiff:
    def __init__(self, ref_fasta, query_fasta, outdir):
        self.ref_fasta = ref_fasta
        self.qry_fasta = qry_fasta
        self.outfile = outfile


    @classmethod
    def run_dnadiff(cls, ref_fasta, query_fasta, outprefix):
        command = ' '.join(['dnadiff -p', outprefix, ref_fasta, query_fasta])
        utils.syscall(command)


    @classmethod
    def load_qdiff_file(cls, qdiff_file):
        '''Loads the out.qdiff file made by dnadiff.
        Returns dict of seq name -> list of sorted non-overlapping
        regions that contain variants compared to the ref'''
        diff_regions = {}

        with open(qdiff_file) as f:
            for line in f:
                fields = line.rstrip().split()
                start, end = sorted([int(fields[2]) - 1, int(fields[3]) - 1])
                interval = pyfastaq.intervals.Interval(start, end)
                if fields[0] not in diff_regions:
                    diff_regions[fields[0]] = []
                diff_regions[fields[0]].append(interval)

        for interval_list in diff_regions.values():
            pyfastaq.intervals.merge_overlapping_in_list(interval_list)

        return diff_regions

