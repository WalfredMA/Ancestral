#!/usr/bin/python

import pandas as pd
import numpy
import collections as cl
import Queue
import sys
import getopt
import os




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

inputfile=inputpath.split('/')[-1]
prefix=inputfile.split('_')[-2]

if out=='':
    out='/'.join(inputpath.split('/')[:-2])+'/'


if phased==0:
    os.system('python freq.py -f {:s} -t {:s} -o {:s}'.format(inputpath,test,out))
    os.system('python phased_correct.py -f {:s}{:s}_freq/{:s} -c {:d}'.format(out,prefix,inputfile, correct))
    os.system('python selectSNP.py -f {:s}{:s}_freq/{:s}'.format(out,prefix,inputfile))
    os.system('python calculate.py -f {:s}{:s}_freq/{:s} -o {:s}'.format(out,prefix,inputfile,out))
    os.system('python chromdetail.py -f {:s}{:s}_cal/{:s} -r {:s} -m {:s}'.format(out,prefix,inputfile, rangei, mp))

else:
    os.system('python phased_freq.py -f {:s} -t {:s}'.format(inputpath,test))
    os.system('python phased_correct.py -f phased_data_freq/{:s}_{:s}_{:d} -c {:d}'.format(inputfile, test,1, correct))
    os.system('python phased_calculate.py -f phased_data_freq/{:s}_{:s}_{:d}'.format(inputfile, test,1))
    os.system('python phased_fullchrom.py -f phased_data_cal/{:s}_{:s}_{:d} -m {:s} -n {:s}_0'.format(inputfile,test,1, mp, test))


    os.system('python phased_freq.py -f {:s} -t {:s} -a'.format(inputpath,test))
    os.system('python phased_correct.py -f phased_data_freq/{:s}_{:s}_{:d} -c {:d}'.format(inputfile, test,0, correct))
    os.system('python phased_calculate.py -f phased_data_freq/{:s}_{:s}_{:d}'.format(inputfile,test, 0))
    os.system('python phased_fullchrom.py -f phased_data_cal/{:s}_{:s}_{:d} -m {:s} -n {:s}_1'.format(inputfile,test,0, mp, test))
