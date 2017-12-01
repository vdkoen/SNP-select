# -*- coding: utf-8 -*-

"""
Script reads a (filtered) vcf file and output the locations and SNP alleles into a data table (csv)
This data table is easy to use and to play with in R.
"""

import vcf
from optparse import OptionParser
from operator import itemgetter
import numpy as np


def parse_options():
    parser = OptionParser()

    parser.add_option("-v", "--vcf", dest="vcf_file",
                      help="input vcf file", metavar="VCF", type="string")

    parser.add_option("-o", "--out", dest="output",
                      help="Output name", metavar="OUT", type="string")
    parser.add_option("-p", "--pairs", dest="pair_file",
                      help="File specifying pairs every line is one pair is komma seperated",
                      metavar="PAIR_FILe", type="string", default=None)

    opts, args = parser.parse_args()
    return opts, args

"""only main function, reads the data and sorts the records on Minor allel frequency
After collecting informative records output it in csv format"""


def calc_entropy(record):
    genotypes = {}

    for sample in range(0, len(record.samples)):
        if record.samples[sample].called:
            # print(record.samples[sample].gt_bases)

            if record.samples[sample].gt_bases in genotypes:
                genotypes[record.samples[sample].gt_bases] += 1
            else:
                genotypes[record.samples[sample].gt_bases] = 1

    probs = [f/sum(genotypes.values()) for f in genotypes.values()]
    return(np.sum(-p * np.log2(p) for p in probs))



def main(opts):

    snp_list = []
    sampleid = []
    write_file = []

    vcf_reader = vcf.Reader(filename=opts.vcf_file)
    for record in vcf_reader:
        entropy = calc_entropy(record)

        if entropy != 0:
            temp_list = [record.CHROM, record.POS, round(entropy, 3)]
            snp_list.append(temp_list)

    sorted_snp_list = (sorted(snp_list, key=itemgetter(2), reverse=True))

    sampleid.append("location")
    sampleid.append("Entropy")

    for id in vcf_reader.samples:
        sampleid.append(id)

    write_file.append(sampleid)


    # print(sorted_snp_list)
    for snp in sorted_snp_list:
        info_snp = vcf_reader.fetch(snp[0], snp[1] - 1)

        write_var = []
        write_var.append(str(snp[0]) + ":" + str(snp[1]))
        write_var.append(str(snp[2]))
        variant = next(info_snp)
        for sample in variant:
            if sample.gt_bases == None:
                write_var.append("-")
            else:
                write_var.append(sample.gt_bases)
        write_file.append(write_var)

    with open(opts.output, 'w') as f:
        for row in zip(*write_file):
            f.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    opts, args = parse_options()
    main(opts)