"""
Script will select markers calculate the allel frequency between the groups,
then it will calculate the entropy over those frequency to determine which have information to separate the groups.
Then output it in csv format. for the algorithm to procces.

input: pairs definition and the filtered vcf file.
output: csv, row's are groups, with columns as loci, data are the allel frequencies with reference. sorted on
"""

import numpy as np
import vcf
from optparse import OptionParser
from operator import itemgetter


def parse_options():
    parser = OptionParser()

    parser.add_option("-v", "--vcf", dest="vcf_file",
                      help="Directory containing all GBS data", metavar="VCF", type="string")

    parser.add_option("-p", "--pairs", dest="pair_file",
                      help="File specifying pairs every line is one pair tab separated",
                      metavar="PAIR_FILe", type="string", default=None)

    parser.add_option("-o", "--out", dest="output",
                      help="Output name", metavar="OUT", type="string")

    opts, args = parser.parse_args()
    return opts, args


def get_pair_freq(record, pairs):
    all_freq = []

    for pair in pairs:

        allel_counts = {}

        for individual in pair:
            # total_reads += mapped_reads[individual][-1]
            if record.samples[record._sample_indexes[individual]].gt_bases is not None:
                bases = record.samples[record._sample_indexes[individual]].gt_bases.split("/")
                amount = record.samples[record._sample_indexes[individual]]["AD"]

                for base in range(0, len(bases)):
                    if bases[base] in allel_counts:
                        allel_counts[bases[base]] += amount[base]
                    else:
                        allel_counts[bases[base]] = amount[base]

        if sum(allel_counts.values()) == 0:
            freq = "-"
        else:
            try:
                freq = allel_counts[record.REF] / sum(allel_counts.values())
            except KeyError:
                freq = float(0)

        all_freq.append(freq)

    # for pair in pairs:
    #     pair_hits = []
    #     for individual in pair:
    #         if record.samples[record._sample_indexes[individual]].gt_bases is not None:
    #             pair_hits.append(record.samples[record._sample_indexes[individual]].gt_bases)
    #
    #     allels = [allels for genotypes in pair_hits for allels in genotypes.split("/")]
    #
    #     try:
    #         ref_freq = allels.count(record.REF) / len(allels)
    #         all_freq.append(ref_freq)
    #
    #
    #     # you cannot divide 6/0 meaning len(allels) is empty so no frequency information is possible.
    #     except ZeroDivisionError:
    #         all_freq.append("-")

    return all_freq


def get_freq(record):
    all_freq = []

    for sample in range(0, len(record.samples)):
        if record.samples[sample].called:
            bases = record.samples[sample].gt_bases.split("/")
            amount = record.samples[sample]["AD"]

            try:
                pos = bases.index(record.REF)
                try:
                    freq = amount[pos] / sum(amount)
                except ZeroDivisionError:
                    freq = "-"
            except ValueError:
                freq = float(0)
        else:
            freq = "-"

        all_freq.append(freq)

    return all_freq

    # if record.samples[record._sample_indexes[individual]].gt_bases is not None:
    #     bases = record.samples[record._sample_indexes[individual]].gt_bases.split("/")
    #     amount = record.samples[record._sample_indexes[individual]]["AD"]


def calc_entropy(record):
    ranges = {"0.25": 0, "0.50": 0, "0.75": 0, "1.0": 0}

    for freq in record:
        try:
            if 0.0 <= freq <= 0.25:
                ranges["0.25"] += 1
            elif 0.25 <= freq <= 0.50:
                ranges["0.50"] += 1
            elif 0.50 <= freq <= 0.75:
                ranges["0.75"] += 1
            elif 0.75 <= freq <= 1.0:
                ranges["1.0"] += 1
        except TypeError:
            pass

    probs = [f / sum(ranges.values()) for f in ranges.values()]
    probs = [x for x in probs if x != 0.0]

    return (np.sum(-p * np.log2(p) for p in probs))


def get_pairs(pair_file):
    pairs = []

    with open(opts.pair_file, "r") as file:
        for line in file:
            clean_line = line.rstrip()
            splitted_line = clean_line.split(",")
            pairs.append(splitted_line)

    return pairs


def main(opts):
    counter = 0

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

    for record in vcf_reader:
        snp_data = []

        if opts.pair_file is not None:
            frequencies = get_pair_freq(record, pairs)
        else:
            frequencies = get_freq(record)

        entropy = calc_entropy(frequencies)
        if entropy != 0.0:
            loci = record.CHROM + ":" + str(record.POS)

            frequencies.insert(0, loci)
            frequencies.insert(1, entropy)
            write_file.append(frequencies)


            # write all info

            # print(frequencies)

    sorted_write_file = (sorted(write_file, key=itemgetter(1), reverse=True))
    sorted_write_file.insert(0, sample_names)

    with open(opts.output, 'w') as f:
        for row in zip(*sorted_write_file):
            f.write('\t'.join(list(map(str, row))) + '\n')


if __name__ == '__main__':
    opts, args = parse_options()
    main(opts)
