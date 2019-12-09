#!/usr/bin/python

import numpy as np
import pandas as pd
import collections as cl
import os
import Queue
import sys
import getopt
from scipy import linalg

even=0
correct=0
projecton=0
source=2
opts,args=getopt.getopt(sys.argv[1:],"f:ecps:")
for op, value in opts:
	if op=='-f':
	    inputfile=value
        if op=='-e':
            even=1
        if op=='-c':
            correct=1
        if op=='-p':
            projecton=1
try:
    os.mkdir('ref_phased')
except:
    pass




def checkerror(diff_index,eigen_snps):


    #print 100*sum([s.tolist()[i] for i in [x for x in diff_index if x<len(s.tolist())]])/sum(s.tolist())


    v0=np.absolute(eigen_snps[:,0]).tolist()
    v1=np.absolute(eigen_snps[:,1]).tolist()

    snps=np.absolute(eigen_snps[diff_index,:]).sum()

    print 100*snps/np.absolute(eigen_snps).sum()
    print 100*sum([v0[i] for i in diff_index])/sum([x for x in v0]), 100*sum([v1[i] for i in diff_index])/sum([x for x in v1])


def runpca(allsamples, group ,even,sample,pop,superpop):

    allsuperpop=[superpop[sample.index(x)] for x in allsamples]

    allsamples=[x for i,x in enumerate(allsamples) if allsuperpop[i] in group]
    allsuperpop=[superpop[sample.index(x)] for x in allsamples]
    allpop=[pop[sample.index(x)] for x in allsamples]


    pcamatrix=pd.DataFrame.as_matrix(t[allsamples]).tolist()
    pcamatrix=[[float(x.count('1')) for x in vec] for vec in pcamatrix]

    if correct==1:

        TSI=[i for i,x in enumerate(allpop) if x=='TSI']
        tuscan=[i for i,x in enumerate(allpop) if x=='Tuscan']

        m=np.matrix(pcamatrix)

        TSI_samples=m[:,TSI]
        tuscan_samples=m[:,tuscan]
        combine_samples=m[:,TSI+tuscan]

        TSI_mean=TSI_samples.mean(1)
        tuscan_mean=tuscan_samples.mean(1)

        combine_mean=combine_samples.mean(1)

        TSI_adj=combine_mean/TSI_mean
        tuscan_adj=combine_mean/tuscan_mean

        tuscan_adj[np.isinf(TSI_adj)] = 0
        tuscan_adj[np.isnan(TSI_adj)] = 0
        TSI_adj[np.isinf(tuscan_adj)]=0
        TSI_adj[np.isnan(tuscan_adj)]=0

        TSI_adj[np.isinf(TSI_adj)]=0
        TSI_adj[np.isnan(TSI_adj)]=0
        tuscan_adj[np.isinf(tuscan_adj)] = 0
        tuscan_adj[np.isnan(tuscan_adj)] = 0

        eur_index=[i for i,x in enumerate(allsamples) if allsuperpop[i] =='EUR']
        eu_index=[i for i,x in enumerate(allsamples) if allsuperpop[i] =='EU']

        m[:,eur_index]=np.array(m[:,eur_index])*np.array(TSI_adj)
        m[:,eu_index]=np.array(m[:,eu_index])*np.array(tuscan_adj)

        pcamatrix=m.tolist()

    pcamatrix=pd.DataFrame.as_matrix(t[allsamples]).tolist()
    pcamatrix=[[x.count('1') for x in vec] for vec in pcamatrix]
    means=[(sum(vec)+0.500000)/(len(vec)+1) for vec in pcamatrix]
    stds=[np.std(vec) for vec in pcamatrix]

    matrix_std=[[(x-means[i]) for x in vec] for i,vec in enumerate(pcamatrix)]

    if even==1:

	pop_list=list(set(allpop))
	mean_pop={x:[np.mean([y for i,y in enumerate(vec) if allpop[i]==x]) for vec in matrix_std] for x in pop_list}

	matrix_std=np.matrix(matrix_std)

	for i in xrange(matrix_std.shape[1]):
		matrix_std[:,i]=np.transpose(np.matrix(mean_pop[allpop[i]])).tolist()
    else:
        matrix_std=np.matrix(matrix_std)

    print 'svd'


    eigen_snps,s,v=np.linalg.svd(matrix_std.tolist())

    #checkerror(diff_index,eigen_snps)

    print 'diagsvd'
    eigen_val=np.matrix(linalg.diagsvd(s, len(eigen_snps), len(v)))

    eigen_vec=np.matrix.transpose(np.matrix(v))

    projection=np.matrix(eigen_snps)*np.transpose(linalg.diagsvd([1.0/x for x in s.tolist()], len(v), len(eigen_snps)))

    return eigen_snps, s,projection,eigen_vec

def N50(j1):
    sum0=sum(j1)
    sumhalf=sum(j1)/2
    j1_sort=sorted(j1)
    num0=len(j1)

    max0=num0-1
    min0=0
    if max0-min0>1:
        nextindex=(max0+min0)/2

        if sum(j1[nextindex:])>sumhalf:
            max0=nextindex
        else:
            min0=nextindex

    return j1[max0]

        

def findprojection(j1,j2):


    j1N50=[N50([abs(x) for x in vec]) for vec in j1]
    j2N50=[N50([abs(x) for x in vec]) for vec in j2]

    j1=[[0.0 if abs(x) <j1N50[i] else x for x in vec] for i,vec in enumerate(j1)]
    j2=[[0.0 if abs(x) <j2N50[i] else x for x in vec] for i,vec in enumerate(j2)] 

    matrix=np.matrix([[max(j1[i][j], j2[i][j]) for j in xrange(len(j1[i]))] for i in xrange(len(j1))])

    return np.transpose(matrix)



t=pd.read_csv(inputfile,sep='\t',dtype='str')

table2=pd.read_csv('allpop.txt',sep='\t',header=None)

with open('phase_index.txt', mode='r') as f:
    read=[int(x) for x in f.read().split(';')]
f.close()



pos=list(t['POS'])

snpindex=[i for i,x in enumerate(pos) if int(x) >= min(read) and int(x) <= max(read)]

if projecton ==0:
    t=t.iloc[snpindex]

print len(t)
sample=list(table2[0])+['Sample']
pop=list(table2[1])+['XX']
superpop=list(table2[2])+['XX']

if source==2:
    source=['EU','EUR']
elif source==1:
    source=['EUR']
else:
    source=['EU']


eu=[sample0 for sample0 in t.columns.tolist()[10:] if superpop[sample.index(sample0)] in source or pop[sample.index(sample0)]=='Russian']

t=t[t.columns.tolist()[:10]+eu]

allsamples=t.columns.tolist()[10:]


eigen_snps, s,j1,v1=runpca(allsamples,source ,even,sample,pop,superpop)

projection=j1

if projecton==1:
    for x in [i for i in xrange(len(t)) if i not in snpindex]:
        projection[x,:]=0


pcamatrix=pd.DataFrame.as_matrix(t[allsamples]).tolist()
pcamatrix=[[x.count('1') for x in vec] for vec in pcamatrix]
means=[(sum(vec)+0.500000)/(len(vec)+1) for vec in pcamatrix]
matrix_std=[[(x-means[i]) for x in vec] for i,vec in enumerate(pcamatrix)]

eigen_vec=np.transpose(np.matrix(matrix_std))*projection
sample_vec1=[(2*(x[0].count('1'))-means[i]) for i,x in enumerate(list(t['Sample']))]
sample_vec2=[(2*(x[-1].count('1'))-means[i]) for i,x in enumerate(list(t['Sample']))]
sample_project1=np.matrix(sample_vec1)*projection
sample_project2=np.matrix(sample_vec2)*projection

out=pd.DataFrame.from_records(sample_project1.tolist()+sample_project2.tolist()).append(pd.DataFrame.from_records(eigen_vec.tolist()),ignore_index=True)

out=pd.DataFrame(out,dtype='str')


allsuperpop=[superpop[sample.index(x)] for x in allsamples]
#allsamples=[x for i,x in enumerate(allsamples) if allsuperpop[i] in ['Russian']]
namelist=['allele0','allele1']+allsamples

out.insert(0,'name',namelist)
out.to_csv('ref_phased/'+inputfile.split('/')[-1]+'_phasevec',mode='w',sep='\t',header=None, index=False)

out=pd.DataFrame.from_records(eigen_snps.tolist())
snps=list(t['ID'])
out.insert(0,'id',snps)

out.to_csv('ref_phased/'+inputfile.split('/')[-1]+'_phasedir',mode='w',sep='\t',header=None, index=False)


out=pd.DataFrame.from_records([[x*x] for x in s.tolist()])
out.to_csv('ref_phased/'+inputfile.split('/')[-1]+'_phaseval',mode='w',sep='\t',header=None, index=False)

os.system('python plotpca.py -f {:s}'.format('ref_phased/'+inputfile.split('/')[-1]+'_phasevec'))
