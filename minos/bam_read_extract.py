import pysam


def get_read_names(infile):
    '''Returns set of read names from input bam file'''
    samfile = pysam.AlignmentFile(infile, 'rb')
    reads = [x.query_name for x in samfile.fetch(until_eof=True)]
    samfile.close()
    return reads


def get_unmapped_reads(infile, outfile):
    '''Writes BAM file of unmapped reads from infile'''
    pysam.view('-b', '-f', '0x4', '-o', outfile, infile, catch_stdout=False)


def get_region(infile, ref_name, start, end, outfile):
    '''Writes BAM file of the given region'''
    region = ref_name + ':' + str(start + 1) + '-' + str(end + 1)
    pysam.view('-b', '-F', '0x4', '-o', outfile, infile, region, catch_stdout=False)

