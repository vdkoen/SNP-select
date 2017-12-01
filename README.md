[![Documentation Status](https://readthedocs.org/projects/snp-select/badge/?version=latest)](http://snp-select.readthedocs.io/en/latest/?badge=latest)
[Full documentation can be found here](http://snp-select.readthedocs.io/en/latest/)

#### Root folder

- This folder / repository root contains the pipeline files Snakefile and Snakefile_vcf including
their configuration files and the scripts which are used with both methods **frequency** and **snps**.
 
#### docs

- Documentation made using Sphinx

#### methods_frequency

- Scripts which are unique for the frequency method. Also includes the script "select_pool_frequency.py"
which is able to fetch snp locations and calculate pooled frequencies. 

#### method_snps 

- Scripts which are unique to the snps method, also includes the script "get_selected_markers.py"
which allows you to fetch the genotype calls from specific snp locations for new samples.  

#### extra scripts

- Some extra scripts designed in the procces of building this pipeline, which might be useful.