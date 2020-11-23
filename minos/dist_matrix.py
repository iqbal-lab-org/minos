import itertools
import json
import os

import numpy as np

from cluster_vcf_records import vcf_file_read, vcf_record


def bed_mask_file_to_dict(bed_file):
    mask = {}
    with open(bed_file) as f:
        for line in f:
            chrom, start, end = line.rstrip().split("\t")
            if chrom not in mask:
                mask[chrom] = set()
            mask[chrom].update(range(int(start), int(end)))

    return mask


def chrom_pos_in_mask(chrom, pos, ref_string, mask):
    if chrom in mask:
        for i in range(pos, pos + len(ref_string)):
            if i in mask[chrom]:
                return True
    return False


def gt_str_to_number(gt):
    gt = set(gt.split("/"))
    if len(gt) != 1 or "." in gt:
        return 0
    else:
        try:
            gt = int(gt.pop())
            return gt + 1
        except:
            return 0


def load_genotypes_from_vcf_file(filename, mask=None):
    mask = {} if mask is None else mask

    with open(filename) as f:
        genotypes = None
        samples = None

        for line in f:
            if line.startswith("#CHROM\t"):
                samples = line.rstrip().split("\t")[9:]
                genotypes = {sample: [] for sample in samples}
            elif not line.startswith("#"):
                fields = line.rstrip().split("\t")
                assert len(fields) == len(samples) + 9
                if chrom_pos_in_mask(fields[0], int(fields[1]) - 1, fields[3], mask):
                    continue
                else:
                    info_keys = fields[8].split(":")
                    for i, sample in enumerate(samples):
                        info = dict(zip(info_keys, fields[i + 9].split(":")))
                        if info.get("FT", "") != "PASS":
                            genotypes[sample].append(0)
                        else:
                            genotypes[sample].append(
                                gt_str_to_number(info.get("GT", "."))
                            )

    for sample, gt_list in genotypes.items():
        genotypes[sample] = np.array(gt_list, dtype=np.uint16)
    return genotypes


def distance_matrix_from_vcf_file(vcf_file, outfile, mask_bed_file=None):
    if mask_bed_file is None or mask_bed_file == "":
        mask = {}
    else:
        mask = bed_mask_file_to_dict(mask_bed_file)
    genotypes = load_genotypes_from_vcf_file(vcf_file, mask=mask)
    samples = sorted(list(genotypes.keys()))
    distances = {}

    for s1, s2 in itertools.combinations(samples, 2):
        key = tuple(sorted([s1, s2]))
        distances[key] = np.sum(
            genotypes[s1] * genotypes[s2] * (genotypes[s1] - genotypes[s2]) != 0
        )

    for s in samples:
        distances[(s, s)] = 0

    with open(outfile, "w") as f:
        print(len(genotypes), file=f)
        for sample in samples:
            out = [distances[tuple(sorted([sample, x]))] for x in samples]
            print(sample, *out, sep="\t", file=f)
