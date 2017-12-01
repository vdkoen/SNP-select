#!/usr/bin/python3.4

"""
Script reads two files a sam file and a file containing the unique read group identifiers those
unique identifiers are added to the sam header. Script also adds flags to bad regions to exclude
for variant calling.

The created bam file is also sorted and outputted.
"""

import pysam
import os
from optparse import OptionParser


def parse_options():
    """Option parser"""
    parser = OptionParser()

    parser.add_option("-s", "--sam", dest="sam_file",
                      help="Sam file you want to process", metavar="SAM", type="string")

    parser.add_option("-r", "--rgt", dest="read_group_tags",
                      help="Output name of the file in gz format", metavar="RGT", type="string")

    parser.add_option("-o", "--out", dest="output",
                      help="Output name", metavar="OUT", type="string")

    opts, args = parser.parse_args()
    return opts, args


def complete_header(read_group_tags, samfile, new_header):
    """Make correct object to parse to the header"""

    all_RG = []

    # Collect all RG headers and put them in a list, add this list to the original header (dict).
    with open(read_group_tags, "r") as file:
        for line in file:
            splitted_1 = line.strip().split(":")
            splitted_2 = line.strip().split("_")
            ID = splitted_1[-1]
            SM = splitted_2[-1]
            new_RG = add_RG_header(ID, SM, samfile)
            all_RG.append(new_RG)

    new_header['RG'] = all_RG
    return new_header


def add_RG_header(ID, SM, samfile):
    """RG header template (dict) for adding the correct format in the sam header"""

    RG_template = {'ID': '',
                   # Read group identifier. e.g., Illumina flowcell + lane name and number
                   'SM': '',  # Sample. Use pool name where a pool is being sequenced.
                   'PL': 'ILLUMINA'}  # Platform/technology used to produce the reads.

    # fill the template, with correct info
    RG_template = RG_template.copy()
    RG_template['ID'] = ID
    RG_template['SM'] = SM

    return RG_template


def main(opts):
    """
    Main function open samfile, collect the correct header and write output as bam file.
    A sorted Bamfile is also outputted.
    Function also flags the XA marked regions in the samfile,
    which can then be excluded in the variant calling.
    """

    samfile = pysam.Samfile(opts.sam_file, 'r')
    new_header = samfile.header.copy()

    final_header = complete_header(opts.read_group_tags, samfile, new_header)

    # "wb" here means write as bam file.
    if not opts.output.endswith(".bam"):
        opts.output = opts.output + ".bam"
    writer_piet = pysam.AlignmentFile(opts.output, "wb", header=final_header)

    for read in samfile:
        tags = dict(read.tags)
        if "XA" in tags:
            read.is_qcfail = True
        writer_piet.write(read)

    writer_piet.close()
    samfile.close()

    print("Bam file created, continuing to sort bam file....")

    path_and_name = os.path.split(opts.output)
    path_and_name = list(path_and_name)

    sorted_bam_out = os.path.join(path_and_name[0], "sorted_" + path_and_name[1])
    pysam.sort("-o", sorted_bam_out, opts.output)
    print("Done, Sorted bam file created")


if __name__ == '__main__':
    opts, args = parse_options()
    main(opts)
