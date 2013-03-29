#!/usr/bin/env python
import sys,os,fnmatch,itertools,pickle
import os.path
import datetime,time
from glob import glob
import numpy as np
from random import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import LSQUnivariateSpline

my_dir = os.path.abspath(os.path.dirname(__file__))
num=sys.argv[1]

class mdict(dict):
    def __setitem__(self,key,value):
        self.setdefault(key,[]).append(value)
def print_dict(dct):
    for key,val in dct.items():
        print key,val
        print ''
    return key

config = pickle.load(open('config.pkl','rb'))
key = print_dict(config)

vel  = key
dist = config[vel][0][0]
ts   = config[vel][0][1]
path_seg   = config[vel][0][2]
path_svel  = config[vel][0][3]
path_vel   = config[vel][0][4]
path_steps = config[vel][0][5]
dct        = config[vel][0][6]       # 'freq'  50*ts/1000
dt         = dct['freq']*ts/1000
path_v_aps = path_vel/ts*1000
domain     = np.cumsum(((path_steps*ts)/1000)*path_v_aps)

spos=xxsposxx
beta=-0.6
quota=xxquotaxx*xxhowmanyxx

def pack_pkl(stage):
    dct_s_hb={}
    dct_s_whb={}
    phase = int(stage)-1
    # protein - protein
    for path in glob(os.path.join(my_dir,'%s/*/*-hb_pr*pr*.pkl.*' % stage)):
        print path
        seed = path.split('.')[-1]
        sample_i = pickle.load(open(path,'rb'))
        if len(sample_i)==100:
            dct_s_hb[seed]=[sample_i]
            os.remove(path)
    if len(dct_s_hb)>0:
        pickle.dump(dct_s_hb,open('%s-sd_hb.pkl' % stage,'w'))
    if my_dir.split('/')[-2].split('.')[1]=='exp':
        for path in glob(os.path.join(my_dir,'%s/*/*-hb_pr*segid*.pkl.*' % stage)):
            print path
            seed = path.split('.')[-1]
            sample_i = pickle.load(open(path,'rb'))
            if len(sample_i)==100:
                dct_s_whb[seed]=[sample_i]
                os.remove(path)
        if len(dct_s_hb)>0:
            pickle.dump(dct_s_whb,open('%s-sd_wp.pkl' % stage,'w'))

# main call
pack_pkl(num)
