Requirements
============

For the software to work some python and R packages are required

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

Python packages
---------------
numpy::

    pip install numpy

pyvcf::

    pip install pyvcf

pyfasta::

    pip install pyfasta

pysam::

    pip install pysam

R Libraries
-----------
Ape and Reshape::

    install.packages("ape", "reshape")
