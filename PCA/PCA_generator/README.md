Generating PCA plots based on genotype data.

Two published databases used: 1000GP and HGDP

Correction for these two platfroms is optional in code. 

Two versions are available: 1.unphased version for unphased data
                            2.phased version for phased data


-----------Usage-----------

python2.7 script.py -f inputfile [options]


-----------options-----------

-e:

If included, evenized all samples for all sub-populations to increase population sensitivity of PCA

-c:

if included, correction discrepancy between two platforms using common sub-populations (TSI) as the reference

-p:

if included, generate PCs first then project target sample on PCs. if not included, generate PCs with target sample. 

