#!usr/bin/env python3

"""
Script merges multiple GBS fastq.gz files into one large file .gz file, it does this while adding
its origin to the to the fastq header to maintain sample information.

used input(fastq header): @HWI-D00148:384:H55W3ADXX:1:1101:1384:2194 1:N:0:ACCTT
file name: A10_ACCTT_merged_R1_001.fastq.gz
"""

import os
import gzip
from optparse import OptionParser


def parse_options():
    """Option parser"""

    parser = OptionParser()
    parser.add_option("-s", "--sep", dest="separator",
                      help="the separator sign in the file name, example: "
                           "C2_CAGAC_merged_R1_001.fastq.gz therefore separator is an underscore",
                      metavar="SEP", type="string", default="_")

    parser.add_option("-d", "--dir", dest="directory",
                      help="Directory containing all GBS data", metavar="DIR", type="string")

    parser.add_option("-o", "--out", dest="output",
                      help="Output name of the file in gz format", metavar="OUT", type="string")

    parser.add_option("-l", "--loc", dest="location",
                      help="Output location of the file, use without / on the end", metavar="LOC",
                      type="string")

    opts, args = parser.parse_args()
    return opts, args


def merge_files(path_to_file, files, new_gz_file):
    """
    Merges the file together while adding the RG:Z: identifier to the fastq header,
    and writes the correct merged header. It also filters out the uniq read group tag identifiers
    and returns those.
    """

    RG_tag_list = []
    # merger = 1
    for single_file in files:

        # merger = 1
        counter = 0

        file_name = get_file_name(single_file, opts.separator)

        with gzip.open(path_to_file + "/" + single_file, 'rt') as file:
            for line in file:

                if counter == 0:
                    clean_line = line.strip()
                    splitted_line = clean_line.split(":")

                    clean_line = clean_line.replace(" ", "_", 1)

                    new_gz_file.write(
                        clean_line + "\t" + "RG:Z:" + splitted_line[2] + "_" + splitted_line[
                            3] + "_" + file_name + "\n")
                    counter += 1

                    if "RG:Z:" + splitted_line[2] + "_" + splitted_line[3] + \
                            "_" + file_name not in RG_tag_list:
                        RG_tag_list.append(
                            "RG:Z:" + splitted_line[2] + "_" + splitted_line[3] + "_" + file_name)

                elif counter == 3:
                    # print(line.strip())
                    new_gz_file.write(line.strip() + "\n")
                    counter = 0

                else:
                    # print(line.strip())
                    new_gz_file.write(line.strip() + "\n")
                    counter += 1

                    # merger += 1
                    # if merger == 51010:
                    #     break

    new_gz_file.close()
    return RG_tag_list


def get_file_name(fasta_file, separator):
    """Filters and isolates the file name, returns the file name"""

    splitted = fasta_file.split(separator)
    name = splitted[0]
    return name


def main(opts):
    """Main function, output the file with all unique read group identifiers."""

    if opts.location is None:
        new_gz_file = gzip.open(opts.output, 'wt')
    else:
        new_gz_file = os.path.join(opts.location, opts.output)

    all_files = os.listdir(opts.directory)

    # Check if the files have right extension, these files are selected for merging in the folder
    files = [f for f in all_files if f.endswith('fq.gz') or f.endswith('.fastq.gz')]

    RG_tags = merge_files(opts.directory, files, new_gz_file)

    output_dir = os.path.dirname(opts.output)
    filename = os.path.basename(opts.output)
    basename = filename.split(".")

    writer_jack = open(output_dir + "/RG_tag_list_" + basename[0] + ".txt", "w")

    for tag in RG_tags:
        writer_jack.write(tag + "\n")
    writer_jack.close()

if __name__ == '__main__':
    opts, args = parse_options()
    main(opts)
