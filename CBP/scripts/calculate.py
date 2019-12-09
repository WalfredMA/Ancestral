#!/usr/bin/python

import pandas as pd
import collections as cl
import os
import Queue
import sys
import getopt
import os

out=''
opts,args=getopt.getopt(sys.argv[1:],"f:o:")
for op, value in opts:
	if op=='-f':
		inputfile=value
        if op=='-o':
                out=value
prefix=inputfile.split('/')[-1].split('_')[-2]

try:
    os.mkdir(out+prefix+'_cal')
except:
    pass

t=pd.read_csv(inputfile,sep='\t')

allsamples=t.columns.tolist()[3:]

genotype=list(t['Genotype'])
martrix=[]

for x in allsamples:
	
	
	martrix.append([[] for x in allsamples])


headers=[[x+'_'+y for y in allsamples] for x in allsamples]
selected=[]

num0=len(allsamples)
for k in xrange(len(t)):
	
	allfreqs=list(t.iloc[k][allsamples])
   
        g=genotype[k]

        p0=(sum(allfreqs)+1.00)/(num0+2.00)
        q0=1-p0

        if g==2:
            allsum=p0*p0
        elif g==1:
            allsum=2*p0*q0
        else:
            allsum=q0*q0

	for i,x in enumerate(allsamples):
		
            for j,y in enumerate(allsamples[i:]):
	            j=j+i
                    if allfreqs[i]*allfreqs[j]<=0:
                        martrix[i][j].append(1.0)  
                    elif g==2:
			martrix[i][j].append(allfreqs[i]*allfreqs[j]/allsum)
	            elif g==1:
                        martrix[i][j].append(((1-allfreqs[i])*allfreqs[j]+(1-allfreqs[j])*allfreqs[i])/allsum)
                    else:
                        martrix[i][j].append((1-allfreqs[i])*(1-allfreqs[j])/allsum)
	
for i,x in enumerate(allsamples):
	
    for j,y in enumerate(allsamples[i:]):
                j=j+i
	        head=headers[i][j]	
		t[head]=martrix[i][j]

	
for x in allsamples:
	
	del t[x]
	
		
t.to_csv(out+prefix+'_cal/'+inputfile.split('/')[-1],sep='\t',index=False, mode='w')



