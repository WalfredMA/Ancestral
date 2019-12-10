
Use Conditional Binomial Probability (CBP) to determine local ethnicity

-----------Usage-----------

python2.7 call.py -f inputfile -t target_sample_name [options]

-----------options-----------

-p:

If included, enter phased mode, used for phased loci. Default mode is for unphase loci.

-r [string]:

the target region for cbp, for example '0_1000000' means determine first 1Mb 

-m:

the max score each SNP can have.

-i [int]:

To split the target region into inputted pieces.

-o [string]:

Output path

-c:

If included, correction discrepancy between two platforms using common population (TSI) as the reference.

