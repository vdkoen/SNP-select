import vcf
import csv

"""Select samples which are not used in the making of the snp set and get the information of those
now they can be clustered to see if the snp set works."""

def main():
    # pass
    # Kwanza3,GlenDee3,Meeker3,Versailles3,CascadeDelight3,Imara3,Paris3,Tulameen3,Nr3,Lagorai3,Willamette3,Heritage3

    new_samples = ["Ace4","Ace5","Adeline4","Alison4","Antonella4","Bardrina4","Bareuro4","Barlicum4","Bizet4","Brio4","Calgary4","Candice4","Chardin4","Greenway4","Hayley4","LP4","Marcellina4","Montreux4","Nagano4","Nautica4","Platinum4","Prunella4","Sakura4","Titus4","Zenia4","Adeline5","Alison5","Antonella5","Bardrina5","Bareuro5","Barlicum5","Bizet5","Brio5","Calgary5","Candice5","Chardin5","Greenway5","Hayley5","LP5","Marcellina5","Montreux5","Nagano5","Nautica5","Platinum5","Prunella5","Sakura5","Titus5","Zenia5"]

    write_file = []
    write_row_name = [""]

    with open('/nfs/BigData01/Big_Data/Lolium/results/grassen_subset/Results/snps_2_filtered_grassen_subset.csv', 'r') as f:
        first_line = f.readline()

    print(first_line)

    vcf_reader = vcf.Reader(filename="/nfs/BigData01/Big_Data/Lolium/results/grassen_named/grassen_name_samtools.vcf.gz")

    snp_index = first_line.split()

    for snp in snp_index[1:]:
        # coordinate = snp.split(".") # for framboos
        coordinate = snp.rsplit(".", 1)

        write_var = [snp]

        info_snp = vcf_reader.fetch(coordinate[0][1:], int(coordinate[1][:-1]) -1 )
        variant = next(info_snp)

        for sample in range(0, len(variant.samples)):
            if variant.samples[sample].sample in new_samples:
                # print(variant.samples[sample].sample, variant.samples[sample].gt_bases)

                if variant.samples[sample].sample not in write_row_name:
                    write_row_name.append(variant.samples[sample].sample)

                if variant.samples[sample].gt_bases == None:
                    write_var.append("-")
                else:
                    write_var.append(variant.samples[sample].gt_bases)
        write_file.append(write_var)


    write_file.insert(0, write_row_name)
    write_file[0][0] = "-"
    print(write_file)

    with open("/nfs/BigData01/Big_Data/L15-217_framboos/results/koen/framboos_subset/Results/classify_subset _2_grassen.csv", 'w') as f:
        for row in zip(*write_file):
            f.write('\t'.join(row) + '\n')


# /nfs/BigData01/Big_Data/L15-217_framboos/data/raw_sequences/merged_framboos_RGtag/framboos_freebayes_raw.vcf.gz
# B4, C5, C9, C12, D7, E12, F2, F5, F8, F11, G4, G11


if __name__ == '__main__':
    main()