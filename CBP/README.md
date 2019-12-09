
Use Conditional Binomial Probability to determine local ethnicity

-----------Usage-----------

python2.7 call.py -f inputfile -t target_sample_name [options]

-----------options-----------

-p:

If included, enter phased model

-c:

if included, correction discrepancy between two platforms using common sub-populations (TSI) as the reference

-p:

if included, generate PCs first then project target sample on PCs. if not included, generate PCs with target sample.



test='Sample'
phased=0
mp='2.0'
rangei='0_10000000000'
piece='1'
allele=''
out=''
correct=0
opts,args=getopt.getopt(sys.argv[1:],"o:f:t:r:m:i:cpa")
for op, value in opts:
	if op=='-f':
	    inputpath=value
        if op=='-t':
            test=value
        if op=='-p':
            phased=1
        if op=='-r':
            rangei=value
	if op=='-m':
            mp==float(value)
        if op=='-i':
            piece=value
        if op=='-a':
            allele=' -a'
        if op=='-o':
            out=value
        if op=='-c':
            correct=1
