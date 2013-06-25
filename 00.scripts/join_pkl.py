#!/usr/bin/env python
import sys,os,pickle,shutil,fnmatch,itertools
import os.path,datetime,time
from glob import glob
from sys import argv
import numpy as np
from random import *

my_dir = os.path.abspath(os.path.dirname(__file__))

num = sys.argv[1]
data = []
collected_work = {}
def load_pkl(st):
    collected_work[st]={}
    acc = []
    for path in glob(os.path.join(my_dir,'*-sfwf.pkl*')):
        print path
        wrk_pkl={}
        wrk_pkl= pickle.load(open(path,'rb'))
        print len(wrk_pkl[st])
        seeds = wrk_pkl[st].keys()
        count = 0
        for s in seeds:
            count += 1
            sample_i = wrk_pkl[st][s][1]
            acc.append(sample_i)
            collected_work[st][s] = {}
            collected_work[st][s][1] = sample_i
    data = np.array(acc)
    print data.shape
    print len(collected_work[st])
    pickle.dump(collected_work,open('%s-sfwf.pkl' % st,'w'))
if '00' == num:
    print num,'no data acquired'
    sys.exit()
load_pkl(num)


'''
collect_deltaf = []
def main_call(st,w_c,d_cp):
    phase = int(st)-1
    acc = []
    wrk_pkl={}
    wrk_pkl= pickle.load(open('%s-sfwf.pkl' % st,'rb'))
    pcorr = pickle.load(open('%s-pmfc.pkl' % st,'rb'))
    p_arr = np.array(pcorr.values()) 
    w_c   = np.cumsum(p_arr)
    print w_c
    print len(wrk_pkl[st])
    seeds = wrk_pkl[st].keys()
    for s in seeds:
        sample_i = wrk_pkl[st][s][1]
        acc.append(sample_i)
    data = np.array(acc)
    print data.shape
    quota.append(data.shape[0])

w_c={}
w_c[0]=0
d_cp={}

# main call
dirs = []
for i in range(1,int(num)+1):
    dirs.append(str(i).zfill(2))
[main_call(st,w_c,d_cp) for st in sorted(dirs)]
'''
