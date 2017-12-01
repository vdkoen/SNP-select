"""
Script will calculate the allel frequency which is equal to the reference genome
over all the groups. the idea is to get a representation of allel frequencies instead of actual snp variants.
since you can report frequencies for every allel we will be reporting the frequencies for all the reference allels.

it calculates allele frequency in multiple ways from allele frequency and from the individual calls this allows
for a comparison which method is the best.
"""

import vcf
import pysam
from optparse import OptionParser
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt

def parse_options():
    parser = OptionParser()

    parser.add_option("-v", "--vcf", dest="vcf_file",
                      help="Directory containing all GBS data", metavar="VCF", type="string")

    parser.add_option("-p", "--pairs", dest="pair_file",
                      help="File specifying pairs every line is one pair tab separated",
                      metavar="PAIR_FILe", type="string", default=None)

    opts, args = parser.parse_args()
    return opts, args


def check_pairs(record, pairs):

    allow_record = True

    all_freq = []

    reference = record.REF
    for pair in pairs:
        pair_hits = []
        for individual in pair:
            if record.samples[record._sample_indexes[individual]].gt_bases is not None:
                pair_hits.append(record.samples[record._sample_indexes[individual]].gt_bases)

        # print(pair)

        allels = [allels for genotypes in pair_hits for allels in genotypes.split("/")]
        # print(allels)
        try:
            ref_freq = allels.count(record.REF) / len(allels)
            all_freq.append(ref_freq)
        except ZeroDivisionError:
            all_freq.append("-")


        # print(reference, ref_freq, pair[0][:-1], allels)

    return(record.CHROM, record.POS, all_freq, pair[0][:-1], allels)

        # for genotype in pair_hits:
        #     try:
        #         allel = genotype.split("/")
        #         all_allels.append()print(allel)
        #     except AttributeError:
        #         pass



                ##############################
    #     if allow_record:
    #         if len(set(pair_hits)) > 1:
    #             for sample in pair:
    #                 record.samples[record._sample_indexes[sample]] = edit_sample(record.samples[record._sample_indexes[sample]], vcf_type)
    #                 record.samples[record._sample_indexes[sample]].called = False
    #
    # return record
def calc_pool_freq(snp, pooled_vcf):

    try:
        info_snp = pooled_vcf.fetch(snp[0], int(snp[1]) - 1)
        variant = next(info_snp)
    except (ValueError, StopIteration):
        return False

    all_pool = []

    for sample in range(0, len(variant.samples)):
        individual = variant.samples[sample]
        if individual.called:
            if individual.is_het:
                freq = individual["AD"][0] / sum(individual["AD"])
            else:
                base = individual.gt_bases.split("/")
                if base[0] == variant.REF:
                    freq = float(1)
                else:
                    freq = float(0)
        else:
            freq = "-"

        all_pool.append(freq)

    return all_pool

####CONTINU HERE#### calculate weighted, weight times the allel freq
def weight_freq(record, pairs, total_weights):

    weighted_freq = []
    ids = []
    for pair in pairs:

        pair_hits = []
        allel_counts = {}

        for individual in pair:
            # total_reads += mapped_reads[individual][-1]
            if record.samples[record._sample_indexes[individual]].gt_bases is not None:
                bases = record.samples[record._sample_indexes[individual]].gt_bases.split("/")
                amount = record.samples[record._sample_indexes[individual]]["AD"]

                # print(bases)
                # print(amount)

                for base in range(0, len(bases)):
                    if bases[base] in allel_counts:
                        # append the new number to the existing array at this slot
                        allel_counts[bases[base]] += amount[base]
                    else:
                        allel_counts[bases[base]] = amount[base]

        try:
            freq = allel_counts[record.REF] / sum(allel_counts.values())
        except KeyError:
            freq = float(0)

        weighted_freq.append(freq)
        ids.append(pair[0][:-1])
    return (weighted_freq, ids)


def get_weights(samfile, pairs):

    indv_reads = {}
    counter = 0
    for read in samfile:

        counter += 1
        if counter > 100000:
            break

        if not read.is_unmapped:
            RG = read.tags[-1][-1]
            RG_splitted = RG.split("_")
            sample = RG_splitted[-1]

            if sample in indv_reads:
                # append the new number to the existing array at this slot
                indv_reads[sample] += 1
            else:
                # create a new array in this slot
                indv_reads[sample] = 1
    # print(indv_reads)
    total_reads = {}

    for pair in pairs:
        for indivudual in pair:
            if indivudual[:-1] in total_reads:
                # append the new number to the existing array at this slot
                total_reads[indivudual[:-1]] += indv_reads[indivudual]
            else:
                # create a new array in this slot
                total_reads[indivudual[:-1]] = indv_reads[indivudual]

    final_weigth = {}

    for pair in pairs:
        for sample in pair:
            weight = indv_reads[sample] / total_reads[sample[:-1]]
            fraction = 1/len(pair)
            final_weigth[sample] = fraction / weight

    return final_weigth



def main(opts):
    counter = 0

    #mapped reads = 275096963

    pooled_vcf = vcf.Reader(filename="/nfs/BigData01/Big_Data/Lolium/results/grassen_pool/grassen_pool_samtools.vcf.gz")
    vcf_reader = vcf.Reader(open(opts.vcf_file, "rb"))
    samfile = pysam.Samfile("/nfs/BigData01/Big_Data/Lolium/results/grassen_named/sorted_grassen_name.bam", 'r')

    write_file = []

    # only if a pair file is specified
    if opts.pair_file is not None:
        # print("its not none")
        pairs = []
        with open(opts.pair_file, "r") as file:
            for line in file:
                clean_line = line.rstrip()
                splitted_line = clean_line.split(",")
                pairs.append(splitted_line)

    final_weight = get_weights(samfile, pairs)

    # print(final_weight)
    # print(len(final_weight))
    # f.write("pooled" + "," + "weighted" + "," + "raw_freq" + "," + "sample" + "," + "loci" + "\n")


    all_names = ["names"]

    for pair in pairs:
        print(pair[0][:-1])
        print(pair[0][:-1] + ".P")
        all_names.append(pair[0][:-1])
    for pair in pairs:
        all_names.append(pair[0][:-1] + ".P")

    write_file.append(all_names)
    print(all_names)

    for record in vcf_reader:
        frequency = check_pairs(record, pairs)
        w_frequency = weight_freq(record, pairs, final_weight)
        pooled_results = calc_pool_freq(frequency, pooled_vcf)

        snp = []


        if pooled_results is False:
            pass
        else:
            # if w_frequency[0].count(1.0) > 5:
            #     pass
            # elif w_frequency[0].count(0.0) > 0:
            #     pass
            # elif pooled_results.count(1.0) > 5:
            #     pass
            # elif pooled_results.count(0.0) > 0:
            #     pass
            # else:

                # print("indv", frequency[2])
                # print("pool", pooled_results)

            snp.append(record.CHROM + ":" + str(record.POS))
            for call in range(0,len(w_frequency[0])):
                snp.append(str(w_frequency[0][call]))
            for pool in range(0,len(pooled_results)):
                snp.append(str(pooled_results[pool]))

            write_file.append(snp)

                # checker = [pooled_results[call], w_frequency[0][call], frequency[2][call]]
                # if 0 in checker or 1 in checker:
                #     pass
                # else:
                #     f.write(str(pooled_results[call]) + "," + str(w_frequency[0][call]) + "," + str(frequency[2][call]) + "," + w_frequency[1][call] + "," + record.CHROM + ":" + str(record.POS) + "\n")

        counter += 1
        # print(counter)
        # if counter > 700:
        #     break



    with open("mega_test.csv", 'w') as f:
        for row in zip(*write_file):
            f.write('\t'.join(row) + '\n')






if __name__ == '__main__':
    opts, args = parse_options()
    main(opts)