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


