#!/usr/bin/env python
import sys,os,fnmatch,itertools,pickle
import os.path,datetime,time
from glob import glob
import numpy as np
from random import *

#removed matplotlib here for keeneland, but it is not
#needed in general...
#import matplotlib
#import matplotlib.pyplot as plt

my_dir = os.path.abspath(os.path.dirname(__file__))
num=sys.argv[1]

def pack_pkl(stage):
    ''' Combine bonding pickles into 1 file.
        hb_protein-protein.pkl  =>  01-sd_hb.pkl
    '''
    with open('%s-sd_hb.pkl' % stage,'w') as hpk:
        for path in glob(os.path.join(my_dir,'%s/*/*-hb_pr*pr*.pkl.*' % stage)):
            dct_s_hb={}
            print path
            seed = path.split('.')[-1]
            sample_i = pickle.load(open(path,'rb'))
            if len(sample_i)==xxlenbpklxx:
                dct_s_hb[seed]=[sample_i]
                os.remove(path)
            if len(dct_s_hb)>0:
                pickle.dump(dct_s_hb,hpk) 

    ''' Combine bonding pickles into 1 file.
        protein-water           =>  01-sd_wp.pkl
    '''
    if my_dir.split('/')[-2].split('.')[1]=='exp':
        with open('%s-sd_wp.pkl' % stage,'w') as wpk:
            for path in glob(os.path.join(my_dir,'%s/*/*-hb_pr*segid*.pkl.*' % stage)):
                dct_s_whb={}
                print path
                seed = path.split('.')[-1]
                sample_i = pickle.load(open(path,'rb'))
                if len(sample_i)==xxlenbpklxx:
                    dct_s_whb[seed]=[sample_i]
                    os.remove(path)
                if len(dct_s_hb)>0:
                    pickle.dump(dct_s_whb,wpk)

# main call
pack_pkl(num)
