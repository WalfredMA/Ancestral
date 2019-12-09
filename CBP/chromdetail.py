#!/usr/bin/python

import pandas as pd
import collections as cl
import os
import Queue
import sys
import getopt
import math
import time
import numpy as np

piece=1
rangei=['0','1000000000']
mp=3.0
out='./'
opts,args=getopt.getopt(sys.argv[1:],"f:i:r:m:o:")
for op, value in opts:
	if op=='-f':
		inputfile=value
        if op=='-i':
                piece=int(value)
        if op=='-r':
                rangei=value.split('_')
        if op=='-m':
                mp=float(value)
        if op=='-o':
                out=value


prefix=inputfile.split('/')[-1].split('_')[-2]

try:
    os.mkdir(out+prefix+'_fullchrom')
except:
    pass

range0=(int(rangei[0]),int(rangei[1]))
t=pd.read_csv(inputfile,sep='\t')

ids=list(t['ID'])
pos=[int(x) for x in list(t['POS'])]
posi=[i for i,x in enumerate(pos) if x>range0[0] and x<range0[1]]

sortedindex=sorted(posi, key=lambda x: pos[x])

allsamples=t.columns.tolist()[3:]

t=t.iloc[sortedindex]
t.reset_index(drop=True, inplace=True)

def findethinic(t0,allsamples,mp,f,b):
    allsamples=t0.columns.tolist()[3:]

    if f==1:
        results=[x/2.0 for x in list(t0.iloc[0][3:])]
    else:
        results=[0 for x in allsamples]

    for i in xrange(f,len(t0)-b):
	
        allchances=list(t0.iloc[i][3:])

        if max(allchances)==0:
            max0=mp+1.0
        else:
            max0=math.log10(float(max(allchances)))
        for j,k in enumerate(allchances): 
            if k==0:
                k=max0-mp
            else:
                k=max(max0-mp,math.log10(float(k)))
	    results[j]=results[j]+k


    l0=len(t0)
    

    back=[x/2.0 for x in list(t0.iloc[l0-1])[3:]]


    if b==1:

        results=[x+back[i] for i,x in enumerate(results)]

    min0=min(results)
    d=[]

    results_fix=[x-max(results) for x in results]

    #sortindex=sorted(range(len(results_fix)), key=lambda x:results_fix[x],reverse=True)[:3]
    return {x:results_fix[i] for i,x in enumerate(allsamples)}


def binarycut(l10, p):

    l2=bin(l10)[2:]

    piece_size=int(l2[:p],2)
    remainder=l2[p:]

    remainder_mark=[i for i,x in enumerate(remainder) if x=='1']

    return remainder,remainder_mark,piece_size

def checkmark(mark, i,p):

    i2=len((bin(i)[2:]).split('1')[-1])
    
    if i2 in mark:
        return 1
    else:
        return 0

def dictadd(d1,d2):

    sumd={}
    allkeys=list(set(d1.keys()+d2.keys()))

    for key in allkeys:

        sumd[key]=0
        if key in d1.keys() and key in d2.keys():
            sumd[key]=d1[key]+d2[key]
    
    return sumd

def cleandict(d,cut=2.0):

    d={x:d[x]-max(d.values()) for x in d.keys()}

    return [x+':'+str(d[x])[:5] for x in sorted([x for x in d.keys() if d[x]>-cut], key=lambda x: d[x], reverse=True)]

t.sort_values(by=['POS'])
pos=list(t['POS'])
pops=list(set([x for y in allsamples for x in y.split('_')]))
inter=len(t)/piece

remainder,mark,piece_size=binarycut(len(t),6)

allsamples=t.columns.tolist()[3:]

piece_num=2**len(remainder)

start=0
frontborder=0

cuttime=len(remainder)

records=[[] for x in xrange(len(bin(piece_num)[2:]))]

cordi=[]

gbr=0

sumborder=0
print 'num',piece_num
for i in xrange(piece_num):

    border=checkmark(mark,i+1,cuttime)

    length=piece_size+border

    result=findethinic(t.iloc[(start-frontborder):(start+length)],allsamples ,mp, frontborder,border)

    pos0=pos[start-1]

    start=start+length

    cordi.append([pos0, pos[start-1]])

    frontborder=border

    records[0].append(result)

cordi[0][0]=min([int(x) for x in pos])

lastrecord=records[0]

sumall={x:0 for x in records[0][0].keys()}

for x in records[0]:

    for y in x.keys():

        sumall[y]=sumall[y]+x[y]


for record in records[1:]:

    iter0=0


    for x in lastrecord:

        if iter0==0:

            lastvalue=x

        else:
            
            add=dictadd(x,lastvalue)
            record.append(add)
            lastvalue=0

        iter0=1-iter0


    lastrecord=record

results=[]
for record in records:
    
    results.append([cleandict(x) for x in record])


allpiecelength=[int(x[1])-int(x[0]) for x in cordi]
percentages=[]
for i,result in enumerate(results):


    percent={}
    sum0=0


    for j,x in enumerate(result):

        percent0={}

        size={}

    	if 0 >= 0:

            for y in x:

                len0=sum(allpiecelength[(j*len(allpiecelength)/len(result)):((j+1)*len(allpiecelength)/len(result))])
                races=y.split(':')[0].split('_')
                value=float(y.split(':')[1])
                if races[0] not in percent0.keys():
                    percent0[races[0]]=0
                    size[races[0]]=0

                percent0[races[0]]=percent0[races[0]]+(10**value)
                size[races[0]]=len0*(10**value)+size[races[0]]

                if races[1] not in percent0.keys():
                    percent0[races[1]]=0
                    size[races[1]]=0

                percent0[races[1]]=percent0[races[1]]+(10**value)
                size[races[1]]=len0*(10**value)+size[races[1]]


        allkeys=sorted(percent0.keys(), key=lambda x: percent0[x], reverse=True)[:3]

        allkeys2=[x for x in percent0.keys() if percent0[x]>=0.4995*sum(percent0.values())]
       

        if 1.0*sum([percent0[x] for x in allkeys])/sum(percent0.values())>=0.999 or len(allkeys2)>0:
                
            for y in list(set(allkeys+allkeys2)):


                if y not in percent.keys():

                    percent[y]=0

                percent[y]=percent[y]+size[y]

                sum0=sum0+size[y]

    


    percent=sorted([(x,100*percent[x]/sum0) for j,x in enumerate(percent.keys())],key=lambda y:y[1],reverse=True)
    percentages.append(percent)

    result=[y for x in result  for y in [x for a in xrange(2**i)]]
    results[i]=result

    
print percentages

for sumall in percentages:

    fullsize=cordi[-1][-1]-cordi[0][0]

    #outfile=';'.join(['{:s}:{:f}'.format(sumall[i][0], sumall[i][1]) for i in xrange(len(sumall))])+'\n'

    outfile=','.join([x[0] for x in sumall])+"\n"+','.join([str(x[1]) for x in sumall])+'\n'

    with open('currentsum.txt', mode='a') as f:
        f.write(outfile)
    f.close()


    d=[cordi]+results

    d=[[x[i] for x in d] for i in xrange(len(d[0]))]


    #pd.DataFrame.from_records(d).to_csv(out+prefix+'_fullchrom/'+inputfile.split('/')[-1]+'_sum.csv',sep=',',index=False, mode='w')



