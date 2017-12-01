==============
Output summary
==============

The pipeline will output several results and some examples and how to interpret those are shown below.

**Pairwise distance table**

Below you see the pairwise distance table value's. For the snps method (shown below) the values indicate
the sum of different genotype calls over all selected snps. For the frequency method the values are
the absolute difference in allele frequency.

+-----+-----+-----+-----+----+----+----+----+
|     | A10 | A11 | A12 | A1 | A2 | A3 | A4 |
+-----+-----+-----+-----+----+----+----+----+
| A10 | 0   |     |     |    |    |    |    |
+-----+-----+-----+-----+----+----+----+----+
| A11 | 7   | 0   |     |    |    |    |    |
+-----+-----+-----+-----+----+----+----+----+
| A12 | 13  | 10  | 0   |    |    |    |    |
+-----+-----+-----+-----+----+----+----+----+
| A1  | 14  | 15  | 14  | 0  |    |    |    |
+-----+-----+-----+-----+----+----+----+----+
| A2  | 10  | 10  | 10  | 6  | 0  |    |    |
+-----+-----+-----+-----+----+----+----+----+
| A3  | 16  | 18  | 14  | 14 | 10 | 0  |    |
+-----+-----+-----+-----+----+----+----+----+
| A4  | 14  | 14  | 11  | 17 | 11 | 8  | 0  |
+-----+-----+-----+-----+----+----+----+----+

**Selected snps**

This is an overview of the selected snps. *Snp_1* will always have the format contig + position of the snp.
``Entropy`` is the information gain value, the higher this value the better the snp can distinguish between the varieties.
The rest of the table contains the genotype calls for that position.

+---------+-------+-------+-------+-------+
|         | Snp.1 | Snp.2 | Snp.3 | Snp.4 |
+---------+-------+-------+-------+-------+
| Entropy | 1.54  | 1.25  | 0.98  | 1.32  |
+---------+-------+-------+-------+-------+
| A10     |       | G/G   | C/G   | T/T   |
+---------+-------+-------+-------+-------+
| A11     | A/A   | G/G   | C/G   | T/T   |
+---------+-------+-------+-------+-------+
| A12     | A/G   | C/G   | G/G   | C/C   |
+---------+-------+-------+-------+-------+
| A1      | A/A   | G/G   | G/G   | T/T   |
+---------+-------+-------+-------+-------+
| A2      | A/A   | G/G   | G/G   | T/T   |
+---------+-------+-------+-------+-------+
| A3      | A/G   | C/G   | G/G   | T/C   |
+---------+-------+-------+-------+-------+
| A4      | A/G   | C/G   | G/G   | C/C   |
+---------+-------+-------+-------+-------+

**VCF file**

A vcf file containing the selected snps, with all their information.

**Flanking sequences**

In the example output you see two columns, the first, left column contains information over the location of the snp.
``Contig1:267879-268079:267979:``.

    * ``Contig1`` is the name of the chromosome/scaffold the snp is located.
    * ``267879-268079`` is the range of the fetched sequence(in the right column).
    * ``267979`` is the exact position of the selected snp.

The second column is the flanking sequences around the main snp, indicated between the brackets ``[C/T]``. The other snps in the flanking sequences
contain iupac nucleotides.
More information about the iupac nucleotides can be found `here <https://www.bioinformatics.org/sms/iupac.html>`_ and on `wikipedia <https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation>`_.

Example output::

    Contig1:267879-268079:267979    AACAAGCGGAGCCTGTTT[C/T]RAWYCWGKTCAAGAAATT
    Contig2:99417-99617:99517       GAGAGCTCGGACGGACTC[C/G]GAGTTGCTCTGTTTCTTC
    Contig5:324044-324244:324144    GTTCCTTCTGGAGCAGAW[G/T]CTCCAGAAGCATATGCCA
    Contig2:1055761-1055961:1055861 GGCTACCATAATGAGAAC[A/T]CCAACCACTGCGCCAATG


http://www.tablesgenerator.com/text_tables
