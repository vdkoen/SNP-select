"""
Script takes the selected snp positions and gets the flanking sequences from the reference genome,
during this process it will find all other variants in the regions and denote those as uipac
nucleotides in the sequence.

It also outputs the final vcf file containing the selected variants.
"""

from optparse import OptionParser
import vcf
import pyfasta

__author__ = "Koen van Diemen"
__date__ = "18/10/2017"


def parse_options():
    """
    Option parser
    :return: returns selected parameters/settings.
    """

    parser = OptionParser()

    parser.add_option("-v", "--vcf", dest="vcf_in",
                      help="input raw vcf file in .gz format also needs to be tabix indexed",
                      metavar="VCF",
                      type="string")

    parser.add_option("-r", "--ref", dest="reference",
                      help="reference genome used to make the alignment, for fetching the "
                           "sequences",
                      metavar="REF", type="string")

    parser.add_option("-s", "--snp", dest="snp_set",
                      help="snp set produced by configure_snp_set.R", metavar="SNP", type="string")

    parser.add_option("-f", "--fasta", dest="fasta_out",
                      help="flanking sequences fasta output with snp locations marked in it to "
                           "design primers on",
                      metavar="FASTA", type="string")

    parser.add_option("-c", "--vcf_out", dest="vcf_out",
                      help="selected snp set output in vcf format fetched from the raw vcf file.",
                      metavar="VCF_OUT",
                      type="string")

    parser.add_option("-l", "--length", dest="length",
                      help="Length of the flanking sequences", metavar="LENGTH", type="int",
                      default="100")

    parser.add_option("-m", "--maf", dest="maf",
                      help="Minor allele frequency of which a variant atleast needs to have before "
                           "being reported as iupac nucleotide", metavar="MAF", type="float",
                      default="0.2")

    parser.add_option("-p", "--call_rate", dest="call_rate",
                      help="minimal amount of genotypes which need to participate before reporting "
                           "the variant as iupac nucleotide", metavar="MAF", type="float",
                      default="0.5")

    opts, args = parser.parse_args()
    return opts, args


def parse_sequence(coordinate, snp_locs, snp_seq, snp_call, main_snp, start):
    """
    Function get the snps and converts those in to iupac nucleotides returns
    the sequence with iupac nucleotides.
    :param coordinate: Position of the snp
    :param snp_locs: Locations of the snps in the flanking sequences.
    :param snp_seq: Full sequence around the main snp.
    :param snp_call: Nucleotides of the snps at the snp_locs.
    :param main_snp: selected snp.
    :param start: Start position of the flanking sequence.
    :return: String if a sequence containg all the flanking sequences in iupac nucleotides.
    """

    lib = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M", "CGT": "B", "AGT": "D",
           "ACT": "H", "ACG": "V", "ACGT": "N"}

    seq_lijst = list(snp_seq)
    snp_pos = int(coordinate) - (start + 1)
    seq_lijst[snp_pos] = "[" + str(main_snp[0]) + "/" + str(main_snp[1]) + "]"

    for loc in range(0, len(snp_locs)):
        snps = [str(i) for i in snp_call[loc]]
        uni = sorted(snps)

        iupac = "".join(uni)
        iupac_code = lib[iupac]

        pos = snp_locs[loc] - (start + 1)
        seq_lijst[pos] = iupac_code

    return "".join(seq_lijst)


def get_all_snps(vcf_reader, coordinate, start, stop, writer_willem, maf_value, call_rate):
    """
    function collects all the variants around the main snp, and returns those.
    :param vcf_reader: Raw vcf file object.
    :param coordinate: Position of the snp.
    :param start: start position of the flanking sequence.
    :param stop: end position of the flanking sequence.
    :param writer_willem: vcf out writer.
    :return: all snps positions in flanking sequences, all variants at the the positions, and the
    actural snp.
    """

    snp_locs = []
    snp_call = []

    for record in vcf_reader.fetch(coordinate[0], start, stop):
        allels = {}
        if record.POS == int(coordinate[1]):
            writer_willem.write_record(record)
            main_snp = record.alleles

        # Small filter whether to report a snp, filter on maf value and call rate
        if record.call_rate > call_rate and record.POS != int(coordinate[1]):
            for sample in range(0, len(record.samples)):
                if record.samples[sample].called:
                    bases = record.samples[sample].gt_bases.split("/")
                    for nucleotide in bases:
                        if nucleotide in allels:
                            allels[nucleotide] += 1
                        else:
                            allels[nucleotide] = 1

            # filter on maf value, if this is small variant is rare so its not that interesting.
            maf = min(allels.values()) / sum(allels.values())
            if maf > maf_value:  # maf value make input variable?
                # ignore indels
                if len(max(record.alleles, key=len)) < 2:
                    snp_locs.append(record.POS)
                    snp_call.append(record.alleles)

    return (snp_locs, snp_call, main_snp)


def main(opts):
    """
    main function
    :param opts: input parameters
    :return: file containing flanking sequences with mutations as iupac nucleotides, and vcf file
    with the selected sequences.
    """
    # input files
    vcf_reader = vcf.Reader(filename=opts.vcf_in)  # read in raw vcf
    reference = pyfasta.Fasta(opts.reference)  # read in reference genome
    with open(opts.snp_set, 'r') as f:  # open selected markers but only first line
        first_line = f.readline()

    # output files
    primer_seq = open(opts.fasta_out, "w")  # output iupack nucleotide file

    # output vcf with only the selected sequences from raw vcf file. (contains all the info)
    writer_willem = vcf.Writer(open(opts.vcf_out, 'w'), vcf_reader, lineterminator='\n')

    # For each of the selected snps
    snp_index = first_line.split()
    for snp in snp_index:

        coordinate = snp.rsplit(".", 1)
        scaffold_len = len(str(reference[coordinate[0]]))
        start = int(coordinate[1]) - opts.length
        stop = int(coordinate[1]) + opts.length

        # if reference sequence is not long enough ajust the lengths of start and stop.
        if start < 0:
            start = 0
        if stop > scaffold_len:
            stop = scaffold_len

        snp_seq = reference[coordinate[0]][start:stop]
        snp_locs, snp_call, main_snp = get_all_snps(vcf_reader, coordinate, start, stop,
                                                    writer_willem, opts.maf, opts.call_rate)
        new_seq = parse_sequence(coordinate[1], snp_locs, snp_seq, snp_call, main_snp, start)

        primer_seq.write(coordinate[0] + ":" + str(start) + "-" + str(
            stop) + ":" + coordinate[1] + "\t")
        primer_seq.write(new_seq + "\n")

    primer_seq.close()
    writer_willem.close()


if __name__ == '__main__':
    opts, args = parse_options()
    main(opts)
