========
Pipeline
========

Introduction
------------

Snp-select is a pipeline for finding minimal amount of snps for accurate identification of plant variety. Snp select is
designed with GBS data but also works with WGS data. The end results is a snp panel in ``vcf format`` and ``flanking sequences``
with iupac nucleotides for accurate primer design.

Their are two version of the pipeline available:

    #. A full pipeline including: alignment, variant calling, variant filtering, and snp selection.
    #. A pipeline from a raw vcf file including: variant filtering and snp selection.

The pipeline has two methods: ``snps`` or ``frequency``, which one you choose is depended on how your species of interest reproduces:

    #. ``frequency``, recommended to be used with cross-pollination species. (Grass, Corn, Spinach, etc.)
    #. ``snps``, recommended to be used with self-pollination/hybrid or cuttings reproduced species. (Rassperry, Tomato, Wheat, Rice, etc.)

For a overview of the pipeline(s) see :numref:`mylabel`

.. _mylabel:
.. figure:: contents/full_workflow_3.png
   :scale: 80 %
   :figclass: align-center

   Snakemake workflow, two version exist the left side is the full pipeline. On the right side a smaller version which takes a raw vcf file incase you have your own.


Usage
-----

Pipeline is designed with `snakemake <http://snakemake.readthedocs.io/en/stable/>`_

Some important ``snakemake`` parameters:

    * ``-p`` Printout shell commands
    * ``-j`` Amount of threats you want to use, default is 8
    * ``-n`` Do dry run

**Full pipeline:**

To run snakemake on the full pipeline edit the **"snakemake_config.yaml"** config file to your liking, or you can also make one yourself see: (:ref:`config_file`) and run::

   snakemake -s Snakefile --configfile snakemake_config.yaml -p -n

This will perform a dry run to see if all the paths are correctly specified. If no errors occur run::

   snakemake -s Snakefile --configfile snakemake_config.yaml -p

**VCF pipeline:**

To run snakemake on the vcf pipeline edit the **"snakemake_config_vcf.yaml"** config file to your liking, or you can also make one yourself see: (:ref:`config_file`) and run::

   snakemake -s Snakefile_vcf --configfile snakemake_config_vcf.yaml -p -n

This will perform a dry run to see if all the paths are correctly specified. If no errors occur run::

   snakemake -s Snakefile_vcf --configfile snakemake_config_vcf.yaml -p

Snakemake has not a very easy way of logging used commands and settings. Here ``snakemake.log`` is a specified file by the user::

   snakemake -s Snakefile_vcf --configfile snakemake_config_vcf.yaml -p | tee snakemake.log

.. note::

   Using ``tee`` will remove the colour output from Snakemake, if you want don't want this you can use ``unbuffer`` to install unbuffer: ``sudo apt-get install expect-dev``


Configuration
-------------

The configuration file for the pipeline can be in ``YAML`` or ``JSON`` format, recommended format is ``YAML``
see (:ref:`config_file`) for an example.

Configuration file can contain the following fields:

``reference_genome:*``
    * Full path to the reference genome, this much be in .gz format!
``output_dir:*``
    * Output directory where you want your results.
``basename:*``
    * A general name for you analysis, this name will be shown in the names of the generated files, a short name is recommended.
``input_dir:*``
    * Input directory containing the input files.
``vcf_file``
    * input raw vcf file (only use this with short pipeline version, from vcf to snp set)
``pair_file:``
    * File specifying if you have multiple individuals from the same variety see (:ref:`pair_file`) for an example and more information.
``variant_caller:*``
    * (String) Variant caller you would like to use. You have two options: **samtools** or **freebayes**.
``method:*``
    * (String) Method you would like to use. You have two options: **snps** and **frequency**.

        #. **snps** will tell the pipeline to use the genotype calls to make a snp panel.
        #. **frequency** will tell the pipeline to calculate allele frequencies using the reference allele, and make a snp panel.
``vcf_filter:*``
    * ``-q:`` (Numeric) Minimal quality a record needs to have to pass the filtering
    * ``-d:`` (Numeric) Minimal amount of reads that need to be mapped in a record
    * ``-s:`` (Numeric) Minimal depth per sample.
    * ``-r:`` (Numeric) The minimal call rate in percentage (between 0 - 1). So the amount of genotypes that need to be present in a record.
    * ``-f:`` (Numeric) Percentage of heterozygotes in a record (between 0 - 1). it checks if a record is trust worthy based on the heterozygous calls, you can expect those to be around 0.5/0.5 for diploide samples. If this ratio is off for all the calls combined it might indicate a faulty side.
``configure_snp_set:*``
    * ``genetic_distance:*``
        * (Numeric) The minimal genetic distance you want to build the snp panel with.
``flanking_sequences:*``
    * ``-l:`` (Numeric) Length of the flanking sequences for primer design.
    * ``-m:`` (Numeric) Minimal minor allele frequency before reporting this variant with iupac nucleotides.
    * ``-p:`` (Numeric) Minimal percentage of the genotypes which have to participate in a snp, before reporting the variants.


.. note::

   Every field with a * behind it is a compulsory field! This much be specified in the config file.

.. _pair_file:                                 

Pair file
^^^^^^^^^
The pair file is not compulsory. Specifying a pair file does have a great influence on how the pipeline functions.
It also depends on what method you are using **snps** or **frequency**

**snps:**
 So when a pair is specified using the snps method the pipeline will search for variants which are stable within the groups specified in the pair file.
 Therefore will only select those variants which have the same genotype call within the specified groups. Since all the
 members of the group now have the same variants they are merged into one sample. If a pair file is simply not specified
 in the configuration this filtering step is not applied.

**frequency:**
 For the frequency method the pair file has different function. If a pair file is specified the allele frequency is
 calculated over all the group members instead if for every group member. This means that all the group members will be merged
 into one single sample having one allele frequency. If a pair file is not specified allele frequency are calculated over
 all the samples instead and nothing will be merged.  


Here you see an example of how the pair file should look, every line represents a group. The individual members should be comma separated.
you can specify as many groups as you want. The names in this file represent the RG identifiers / sample names in the VCF file.


.. include:: contents/pair_example.txt
    :literal:

.. _config_file:

Configuration example full pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example of the configuration file in ``YAML``

.. code-block:: yaml

    reference_genome:
      /nfs/BigData01/Big_Data/Genomes/Rubus_Occidentalis/BWA_index/Rubus_occidentalis_v1.0.a1.scaffolds.fasta.gz
    output_dir:
      /nfs/BigData01/Big_Data/L15-217_framboos/results/koen/framboos_meh_test_2/
    basename:
      framboos_name # how your data is called, plz do not use a extension here.
    input_dir:
      /nfs/BigData01/Big_Data/L15-217_framboos/data/raw_sequences_reduced
    pair_file:
      /home/koenvd/SNP_selection/framboos_pairs_names.txt
    variant_caller:
      samtools
    #  freebayes
    method:
    #  frequency
      snps
    vcf_filter:
      -q: 20
      -d: 30
      -s: 4
      -r: 0.8
    configure_snp_set:
      genetic_distance: 1
    flanking_sequences:
      -l: 100
      -m: 0.2
      -p: 0.5

.. _config_file_vcf:

Configuration example vcf pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    reference_genome:
      /nfs/BigData01/Big_Data/Genomes/Lolium_perenne/clean_genome_lolium.fna.gz
    output_dir:
      /nfs/BigData01/Big_Data/Lolium/results/grassen_pipeline_2/
    vcf_file:
      /nfs/BigData01/Big_Data/Lolium/results/grassen_named/grassen_name_samtools.vcf.gz
    pair_file:
      /home/koenvd/SNP_selection/grassen_pairs.txt
    basename:
      grassen_pair_test
    method:
      frequency
    #  snps
    vcf_filter:
      -q: 20
      -d: 30
      -s: 30
      -r: 0.8
    configure_snp_set:
      genetic_distance: 5
    flanking_sequences:
      -l: 100
      -m: 0.2
      -p: 0.5


Tips and tricks
===============

Reference genome must have a .gz variant and one without. kinda strange but pyfasta cannot work with .gz files
where other tools need to have a .gz variant

for fetch the reference sequences the header must an exact match with the contig or scaffold name.
example::

    >MEHO01000001.1
    CAAGTATCGTTCTGTCATCTAAACTAGATAGTTCTAATATCTCACATGCATAGCAGAAGCATTTTCATTATCTTCGACTT
    ATCTTATATTTTGCTTCTTCAGGAATTTGTGCAAGATGACAAACGTGGTTCAGAGTTTCCAGGCAAACCATCAGGAGTTT

Something like this will not work for now::

    >MEHO01000001.1 Lolium perenne isolate 4540-9 Ryegrass_Norm_contig_4785, whole genome shotgun sequence
    CAAGTATCGTTCTGTCATCTAAACTAGATAGTTCTAATATCTCACATGCATAGCAGAAGCATTTTCATTATCTTCGACTT
    ATCTTATATTTTGCTTCTTCAGGAATTTGTGCAAGATGACAAACGTGGTTCAGAGTTTCCAGGCAAACCATCAGGAGTTT


More comming soon

