#!/usr/bin/python

import pandas as pd
import numpy
import collections as cl
import Queue
import sys
import getopt
import os

test='Sample'
allele=0
correct=0
opts,args=getopt.getopt(sys.argv[1:],"f:t:a")
for op, value in opts:
	if op=='-f':
		inputfile=value

        if op=='-t':
                test=value

        if op=='-a':
                allele=-1
        
        if op=='-c':
                correct=1


try:		
    os.mkdir('phased_data_freq')
except:
    pass

def calculation(genotype, allgenotypes, pop_index,missinginfor):

	
	change=[]
	for i,pop0 in enumerate(pop_index.keys()):
                
                if (len(pop_index[pop0])-sum([missinginfor[i] for i in pop_index[pop0]])*0.5)==0:
                    maf0=-1
                else:
                    num0=2*len(pop_index[pop0])-sum([missinginfor[i] for i in pop_index[pop0]])
                    maf0=(sum([allgenotypes[i] for i in pop_index[pop0]])+1.00)/(num0+2.00)
		
			
		change.append(maf0)
	
	return [change[i] for i in xrange(len(change))]
		

with open('allEU.txt', mode='r') as f:
    read=[x for x in f.readlines() if 'Russia' in x]
f.close()

russians=[x.split('\t')[0] for x in read]

with open('phase_index.txt',mode='r') as f:
    read=[int(x) for x in f.read().split(';')]
f.close()

index00=[x for x in read]

t0=pd.read_csv(inputfile,sep='\t',dtype='str')

#t0=t0.loc[(t0['QUAL']!='-')]

posi=[int(x) for x in list(t0['POS'])]
index01=[i for i,x in enumerate(posi) if x >= min(index00) and x <=max(index00)]
sample=t0.columns.tolist()[10:]

t0.loc[index01][['Sample']+russians].to_csv('log.txt',mode='w',sep='\t')

sample_geno=list(t0['Sample'])
index02=[i for i,x in enumerate(sample_geno) if '/' not in x]
t0=t0.iloc[index02]

t1=pd.read_csv('allpop.txt',sep='\t',header=None)

sample_name=list(t1[0])
pop=list(t1[1])
superpop=list(t1[2])

sample=[x for x in sample if x in sample_name and superpop[sample_name.index(x)] in ['EU','EUR']]

t0=t0[['ID','Sample','POS']+sample]

sample_pop=[pop[sample_name.index(x)] for x in sample]

pop_index={}

for i,x in enumerate(sample_pop):
	
        if x not in pop_index.keys():
            pop_index[x]=[i]
        else:
	    pop_index[x].append(i)

patience_genotype=[x[allele] for x in list(t0[test])]
pop_chance=[]
for i in xrange(len(t0)):
    
	allgenotypes=[x.count('1') for x in list(t0.iloc[i][sample])]
        missing=[x.count('-') for x in list(t0.iloc[i][sample])]
	pop_chance.append(calculation(patience_genotype[i], allgenotypes,pop_index,missing))

pos=list(t0['POS'])
SNPs=list(t0['ID'])
out=[[pos[i],SNPs[i],patience_genotype[i]]+pop_chance[i] for i in xrange(len(t0))]

pd.DataFrame.from_records(out).to_csv('phased_data_freq/'+inputfile.split('/')[-1]+'_'+test+"_"+str(allele+1), sep='\t', index=False, header=['POS','ID','Genotype']+pop_index.keys(),mode='w')


