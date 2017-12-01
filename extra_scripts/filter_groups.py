import vcf


def edit_sample(sample, vcf_type):

    header = sample.site.FORMAT
    header_list = header.split(':')
    call_data = vcf.model.make_calldata_tuple(header_list)

    if vcf_type == "samtools":
        values = ["./.", sample["PL"], sample["DP"], sample["SP"],
                  sample["ADF"], sample["ADR"], sample["AD"]]
    elif vcf_type == "freeBayes":
        values = ["./.", sample["DP"], sample["AD"], sample["RO"],
                  sample["QR"], sample["AO"], sample["QA"], sample["GL"]]

    model = vcf.model._Call(sample.site,
                            sample.sample,
                            call_data(*values))

    return model


def check_pairs(record, pairs):

    check = []

    for pair in pairs:
        if len(pair) == 2:
            sample1 = record.samples[record._sample_indexes[pair[0]]]
            sample2 = record.samples[record._sample_indexes[pair[1]]]

            if sample1.gt_bases == sample2.gt_bases:
                check.append("0")
            else:
                check.append("1")

        else:
            sample1 = record.samples[record._sample_indexes[pair[0]]]
            sample2 = record.samples[record._sample_indexes[pair[1]]]
            sample3 = record.samples[record._sample_indexes[pair[2]]]

            if sample1.gt_bases == sample2.gt_bases == sample3.gt_bases:
                check.append("0")
            else:
                check.append("1")

    if "1" in check:
        return False
    else:
        return True

def main():
    counter = 0

    pairs = []

    vcf_reader = vcf.Reader(open("framboos_only_freebayes_raw.vcf.gz", "rb"))
    writer_willem = vcf.Writer(open("pairs_framboos_only_freebayes_raw.vcf", 'w'), vcf_reader, lineterminator='\n')

    with open("pairs.txt", "r") as file:
        for line in file:
            splitted_line = line.split()
            pairs.append(splitted_line)

    for record in vcf_reader:

        wow = check_pairs(record, pairs)
        if wow:
            writer_willem.write_record(record)

    writer_willem.close()




if __name__ == '__main__':
    main()