"""
Script filters VCF file to only obtain informative sides.
"""

import vcf
from optparse import OptionParser


def parse_options():
    parser = OptionParser()

    parser.add_option("-v", "--vcf", dest="vcf_file",
                      help="Directory containing all GBS data", metavar="VCF", type="string")

    parser.add_option("-o", "--out", dest="output",
                      help="Output name", metavar="OUT", type="string")

    parser.add_option("-q", "--qual", dest="record_qual",
                      help="Minimal quality a record needs to have to pass the filtering", metavar="QUAL", type="int",
                      default="20")

    parser.add_option("-d", "--depth", dest="record_depth",
                      help="Minimal amount of reads mapped before calling a variant at that location", metavar="DEPTH",
                      type="int", default="30")

    parser.add_option("-s", "--sample_depth", dest="sample_depth",
                      help="Minimal depth per sample", metavar="SAMPLE_DEPTH", type="int", default="3")

    parser.add_option("-r", "--call", dest="call_rate",
                      help="Amount of genotypes which need to be present, before including in downstream analysis",
                      metavar="CALL_RATE", type="float", default="0.8")

    parser.add_option("-f", "--fraction_het", dest="fraction_het",
                      help="Fraction of heterozygous in a record before filtering it, the higher it is the larger the "
                           "chance of not removing faulty heterozygous records. But if its to low it could remove "
                           "a record only based a single bad heterozygous call.",
                      metavar="HET", type="float", default="0.2")

    opts, args = parser.parse_args()
    return opts, args


# Here we want to chance the GT tag to "./." for samples we don't want to pass our filter step.
# Since pyvcf makes a Calldata object which is not assignable, we make a new one with the correct information.

def edit_sample(sample):
    header = sample.site.FORMAT
    header_list = header.split(':')

    call_data = vcf.model.make_calldata_tuple(header_list)
    values = ["./.", sample["DP"], sample["AD"], sample["RO"],
              sample["QR"], sample["AO"], sample["QA"], sample["GL"]]

    model = vcf.model._Call(sample.site,
                            sample.sample,
                            call_data(*values))

    return model


# some filter
def filter_sample_depth(record, sample_depth):
    for sample in range(0, len(record.samples)):
        try:
            if record.samples[sample]["DP"] < sample_depth:
                record.samples[sample] = edit_sample(record.samples[sample])
                record.samples[sample].called = False

            # heterozygoot calls which are 1/2 are not trustworthy so filter out.
            if record.samples[sample].is_het:
                if sum(record.samples[sample]["AD"]) < 4:
                    record.samples[sample] = edit_sample(record.samples[sample])
                    record.samples[sample].called = False
        except TypeError:

            record.samples[sample] = edit_sample(record.samples[sample])
            record.samples[sample].called = False

    return record


"""
This function checks if a record is trust worthy based on the heterozygous calls, you can expect those to be around
0.5/0.5 for diploide samples. If this ratio is off for all the calls combined, it might indicate a faulty side 
and therefore those sides are filtered out and excluded. 
"""


def filter_site(record, fraction_het):
    het = 0
    hom = 0
    some_value = 0

    for sample in range(0, len(record.samples)):
        # if not called dont check it
        if record.samples[sample].called == False:
            pass

        else:
            # if it is a heterozygote with more then 2 different bases, we don't look at it yet. let the record pass
            # can also let it fail... not sure what to do.
            if len(record.samples[sample]["AD"]) > 2:
                return False
            else:
                # if heterozygote only has two allels. count het
                if record.samples[sample].is_het:
                    het += 1
                # if not hetero count it as a homozygote.
                elif record.samples[sample].is_het == False:
                    hom += 1

    try:
        # fraction of the called heterozygotes must be atleast 20 procent before checking the ratio.
        if het / (het + hom) > fraction_het:
            for sample in range(0, len(record.samples)):
                # if sample called and heterozygote count the fraction. All it to some_value
                if record.samples[sample].called:
                    if record.samples[sample].is_het:
                        some_value += max(record.samples[sample]["AD"]) / sum(record.samples[sample]["AD"])

            # Check if some_value has a logical ratio, if not remove the entire record.
            if not 0.3 <= some_value / het <= 0.7:
                return False
            else:
                return True

        # if ratio het /het + hom is not 20 percent we just let it pass.
        else:
            return True

    except ZeroDivisionError:
        # return false since something/0 is not possible. meaning that something is called but not a het or hom.
        # can not trust these sites.
        return False


# Record filter on quality
def filter_qual(record, minimal_qual):
    if record.QUAL >= minimal_qual:
        return True
    else:
        return False


# Record filter on depth
def filter_depth(record, minimal_depth):
    if record.INFO["DP"] >= minimal_depth:
        return True
    else:
        return False


def sample_het(record):
    for sample in range(0, len(record.samples)):
        if record.samples[sample].called:
            if record.samples[sample].is_het:
                try:
                    if not 0.3 <= max(record.samples[sample]["AD"]) / sum(record.samples[sample]["AD"]) <= 0.7:
                        record.samples[sample] = edit_sample(record.samples[sample])
                        record.samples[sample].called = False
                except ZeroDivisionError:
                    pass
    return record


"""Keep only records with atleast information about a user specified fraction"""


def genotype_depth(record, call_rate):
    if record.call_rate >= call_rate:
        return True
    else:
        return False


def main(opts):
    counter = 0

    vcf_reader = vcf.Reader(open(opts.vcf_file, "rb"))
    writer_willem = vcf.Writer(open(opts.output, 'w'), vcf_reader, lineterminator='\n')

    for record in vcf_reader:
        counter += 1

        genotype_info_first = genotype_depth(record, opts.call_rate)
        if genotype_info_first:
            qual_pass = filter_qual(record, opts.record_qual)
            depth_pass = filter_depth(record, opts.record_depth)

            if qual_pass and depth_pass:
                # print("pass qual and dept")
                new_record = filter_sample_depth(record, opts.sample_depth)
                some_site = filter_site(new_record, opts.fraction_het)

                if some_site:
                    # print("passed heterozygot filter")
                    new_record = sample_het(new_record)
                    genotype_info_second = genotype_depth(new_record, opts.call_rate)
                    if genotype_info_second:
                        # print("passed all writing")
                        writer_willem.write_record(new_record)

        # print(counter)
        # if counter > 10000:
        #     break

    writer_willem.close()


if __name__ == '__main__':
    opts, args = parse_options()
    main(opts)



    # for filtering the heterozygous ration on sample depth only.

    # total_allel_depth = sum(record.samples[sample]["AD"])
    #
    # if len(record.samples[sample]["AD"]) == 2:
    #     if not 0.3 <= record.samples[sample]["AD"][0] / total_allel_depth <= 0.7:  # 10/0
    #         record.samples[sample] = edit_sample(record.samples[sample])
    #         record.samples[sample].called = False
    # if record.samples[sample].is_het:
