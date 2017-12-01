#!/usr/bin/env bash

command -v freebayes >/dev/null 2>&1 || { echo "freebayes not founds, please install before usage.  " >&2; }
command -v snakemake >/dev/null 2>&1 || { echo "snakemake not founds, please install before usage.  " >&2; }
command -v bwa >/dev/null 2>&1 || { echo "bwa not founds, please install before usage.  " >&2; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not founds, please install before usage.  " >&2; }
command -v bgzip >/dev/null 2>&1 || { echo "bgzip not founds, please install htslib before usage.  " >&2; }
command -v tabix >/dev/null 2>&1 || { echo "tabix not founds, please install htslib before usage.  " >&2; }
command -v python3 >/dev/null 2>&1 || { echo "python 3 not founds, please install before usage.  " >&2; }
command -v R >/dev/null 2>&1 || { echo "R not founds, please install before usage.  " >&2; }

echo "Checking python packages..."
pip freeze | egrep "numpy">/dev/null 2>&1 || { echo "numpy might not be installed.  " >&2; }
pip freeze | egrep "pyvcf">/dev/null 2>&1 || { echo "pyvcf might not be installed.  " >&2; }
pip freeze | egrep "pyfasta">/dev/null 2>&1 || { echo "pyfasta might not be installed.  " >&2; }
pip freeze | egrep "pysam">/dev/null 2>&1 || { echo "pysam might not be installed.  " >&2; }

echo "Checking R packages..."
R -e "installed.packages()" | grep "ape">/dev/null 2>&1 || { echo "ape might not be installed.  " >&2; }
R -e "installed.packages()" | grep "reshape">/dev/null 2>&1 || { echo "reshape might not be installed.  " >&2; }
echo "DONE"