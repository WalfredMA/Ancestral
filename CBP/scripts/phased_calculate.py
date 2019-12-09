#!/usr/bin/python

import pandas as pd
import collections as cl
import os
import Queue
import sys
import getopt
import os

opts,args=getopt.getopt(sys.argv[1:],"f:")
for op, value in opts:
	if op=='-f':
		inputfile=value

try:
    os.mkdir('phased_data_cal')
except:
    pass

t=pd.read_csv(inputfile,sep='\t')

allsamples=t.columns.tolist()[3:]

genotype=list(t['Genotype'])
martrix=[]

for x in allsamples:
	
	
	martrix.append([])


headers=[x for x in allsamples]
selected=[]

num0=len(allsamples)
for k in xrange(len(t)):
	
	allfreqs=list(t.iloc[k][allsamples])
   
        g=genotype[k]

        p0=(1.0+sum(allfreqs))/(num0+2)
        q0=1-p0

        if g==2:
            allsum=p0*p0
        elif g==1:
            allsum=p0
        else:
            allsum=q0

	for i,x in enumerate(allsamples):
	    #allsum=1	
            if (allfreqs[i]*allfreqs[i])<=0:
                martrix[i].append(1.)
            elif g==2:
	        martrix[i].append(allfreqs[i]*allfreqs[i])
	    elif g==1:
                martrix[i].append(allfreqs[i]/allsum)
            else:
                martrix[i].append((1-allfreqs[i])/allsum)


for i,x in enumerate(allsamples):
	
    t[headers[i]]=martrix[i]

	
		
t.to_csv(inputfile.replace('_freq','_cal'),sep='\t',index=False, mode='w')



