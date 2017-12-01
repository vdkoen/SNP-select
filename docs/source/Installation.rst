Requirements
============

Clone the pipeline from one of the following locations::

    git clone git@gitlab.naktuinbouw.net:bioinformatics/GBS-SNP-selection.git
    git clone https://github.com/vdkoen/SNP-select.git

Its recommended to use at least the versions on which the pipeline is designed::

    Python: 3.4.3
    R version: 3.4.2
    Snakemake: 4.0.0
    Bwa: 0.7.15-r1144-dirty
    Freebayes: v1.1.0-44-gd784cf8
    Samtools: 1.5
    Htslib: 1.5

To check if you have missing dependencies run::

    bash check_dependencies.sh

Snakemake
---------
Installation instructions can be found `here: <http://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_

Burrow-Wheeler Aligner
----------------------
Download here: https://sourceforge.net/projects/bio-bwa/

To install::

    tar -xvf bwa-x.x.x.tar.bz2
    cd bwa-x.x.x

    ./configure
    make
    make install

Freebayes
---------
To install::

    git clone --recursive git://github.com/ekg/freebayes.git
    cd
    make
    make install

Samtools and Htslib
-------------------
Installation instructions can be found `here <http://www.htslib.org/download/>`_

Python & Python packages
------------------------
To install python3::

    sudo apt-get install python3

For the python packages you can use the requirements.txt::

    pip3 install -r requirements.txt

Or you can install everything manually

numpy::

    pip3 install numpy

pyvcf::

    pip3 install pyvcf

pyfasta::

    pip3 install pyfasta

pysam::

    pip3 install pysam

R & R Libraries
---------------

Install R::

    sudo apt-get install r-base

Ape and Reshape::

    R
    install.packages("ape")
    install.packages("reshape")



Sometimes installing "ape" will give the following error::

    /usr/bin/ld: cannot find -llapack
    /usr/bin/ld: cannot find -lblas
    collect2: error: ld returned 1 exit status

To solve this you will need the following packages::

    sudo apt-get install liblapack-dev
    sudo apt-get install liblapack3
    sudo apt-get install libopenblas-base
    sudo apt-get install libopenblas-dev


