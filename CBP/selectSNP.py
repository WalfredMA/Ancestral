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
    os.mkdir('data_selected')

except:
    pass
t=pd.read_csv(inputfile,sep='\t')

allsamples=t.columns.tolist()[3:]

selected=[]
for i in xrange(len(t)):
	
	allfreqs=list(t.iloc[i][allsamples])
	
    	if min(allfreqs)>=0 and max(allfreqs)-min(allfreqs)>0:
		
		selected.append(i)
		
t.iloc[selected].to_csv(inputfile,sep='\t',index=False, mode='w')
