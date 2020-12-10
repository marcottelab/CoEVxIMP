The directories and scripts included were used to calculate evolutionary couplings between subunits of protein complexes and used the intermolecular coupling results to build integrative models.

################# EVCouplings ##########################

A basic configuration file is provided. This can be used
in combination with the batchComplexCoupling.py script to
compute pairwise evolutionary couplings between all
proteins provided. For use:

python batchComplexCoupling.py -h

usage: batchComplexCoupling.py [-h] proteinlist [proteinlist ...]

Calculates pairwaise coevolutionary scores for protein pairs in an interaction
group.

positional arguments:
  proteinlist  list of uniprot ids of proteins in a complex groups.

optional arguments:
  -h, --help   show this help message and exit

All protein scripts used are provided in their prospective
directories.

The CoEv_analysis.ipynb notebook can be used for selection of
appropriate intermolecular couplings for integrative modeling.

############## Integrative Modeling ####################

All input files and modeling scripts are provided in their
prospective directories. Analysis scripts and exhaustive sampling
scripts are provided as well.
