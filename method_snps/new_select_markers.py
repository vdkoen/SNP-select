# -*- coding: utf-8 -*-
# !/bin/env python3.4

"""
Script reads a (filtered) vcf file and output the locations and SNP alleles into a data table (csv)
Additionally snps are sorted on information content (entropy), most informative to least
informative.
"""

from optparse import OptionParser
from operator import itemgetter
import vcf
import numpy as np

__author__ = "Koen van Diemen"
__date__ = "18/10/2017"

def parse_options():
    """Option Parser"""
    parser = OptionParser()

    parser.add_option("-v", "--vcf", dest="vcf_file",
                      help="input vcf file", metavar="VCF", type="string")

    parser.add_option("-o", "--out", dest="output",
                      help="Output name", metavar="OUT", type="string")
    parser.add_option("-p", "--pairs", dest="pair_file",
                      help="File specifying pairs every line is one pair is comma seperated",
                      metavar="PAIR_FILe", type="string", default=None)

    opts, args = parser.parse_args()
    return opts, args


"""only main function, reads the data and sorts the records on Minor allel frequency
After collecting informative records output it in csv format"""


# make list as input loop through the list should be fine
def calc_entropy(record):
    """Calculates Shannon entropy (information theory)"""
    genotypes = {}

    for snp in record:
        if snp != "-":
            if snp in genotypes:
                genotypes[snp] += 1
            else:
                genotypes[snp] = 1
    #
    probs = [f / sum(genotypes.values()) for f in genotypes.values()]
    return np.sum(-p * np.log2(p) for p in probs)


def get_genotypes(record):
    """Get the genotype calls for every loci, when nog pairs are specified"""
    snp_list = []

    for sample in range(0, len(record.samples)):
        if record.samples[sample].called:
            snp_list.append(record.samples[sample].gt_bases)
        else:
            snp_list.append("-")

    return snp_list


def get_pair_genotypes(record, pairs, sample_names):
    """Get the genoype calls if pairs are specified"""
    snp_dict = {}
    snp_list = []

    for pair in pairs:
        pair_hits = []
        for individual in pair:
            pair_hits.append(record.samples[record._sample_indexes[individual]].gt_bases)
        hits = list(set(pair_hits))
        if len(hits) > 1 or hits[0] is None:
            hit = "-"
        else:
            hit = hits[0]

        snp_dict[individual[:-1]] = hit

    for name in sample_names[2:]:
        snp_list.append(snp_dict[name])

    return snp_list


def get_pairs(pair_file):
    """Read pairs from input file and put them in a matrix (nested list)"""
    pairs = []

    with open(pair_file, "r") as file:
        for line in file:
            clean_line = line.rstrip()
            splitted_line = clean_line.split(",")
            pairs.append(splitted_line)

    return pairs


def main(opts):
    """Main function, read, parse and output data"""

    write_file = []
    sample_names = ["Location", "Entropy"]

    vcf_reader = vcf.Reader(open(opts.vcf_file, "rb"))

    if opts.pair_file is not None:
        pairs = get_pairs(opts.pair_file)

        for id in vcf_reader.samples:
            if id[:-1] not in sample_names:
                sample_names.append(id[:-1])
    else:
        for id in vcf_reader.samples:
            sample_names.append(id)

    # print(sample_names)
    for record in vcf_reader:
        if opts.pair_file is not None:
            snps = get_pair_genotypes(record, pairs, sample_names)
        else:
            snps = get_genotypes(record)

        entropy = calc_entropy(snps)
        if entropy != 0:
            loci = record.CHROM + ":" + str(record.POS)
            snps.insert(0, entropy)
            snps.insert(0, loci)

            write_file.append(snps)

    sorted_write_file = (sorted(write_file, key=itemgetter(1), reverse=True))
    sorted_write_file.insert(0, sample_names)

    with open(opts.output, 'w') as f:
        for row in zip(*sorted_write_file):
            f.write('\t'.join(list(map(str, row))) + '\n')


if __name__ == '__main__':
    opts, args = parse_options()
    main(opts)
