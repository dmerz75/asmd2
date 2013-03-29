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
num =str(len(path_steps)).zfill(2)

def plot_hb(data,stage):
    phase = int(stage)-1
    if stage =='01':
        d = np.linspace(spos,spos+domain[phase],data.shape[0])
        plt.plot(d,data,'k-',label="hydrogen bonds",linewidth=1)
    elif stage !='01':
        d = np.linspace(spos+domain[phase-1],spos+domain[phase],data.shape[0])
        plt.plot(d,data,'k-',linewidth=1)

def plot_pkl(stage,sel):
    dct_sd_hb = pickle.load(open('%s-sd_%s.pkl' % (stage,sel),'rb'))
    #print '%s-sd_%s.pkl' % (stage,sel),dct_sd_hb.items()[0]
    seeds = dct_sd_hb.keys()
    acclens=[]
    for s in seeds:
        acc=[]
        sample_i = dct_sd_hb[s]
        for c in range(len(sample_i[0])):
            hbc=len(sample_i[0][c])
            acc.append(hbc)
        acclens.append(acc)
    idata= np.array(acclens)
    data = idata.mean(axis=0)
    plot_hb(data,stage)

def main_bond(sel):
    # PLOT - pmf
    fig=plt.figure()
    plt.clf()

    # main call
    dirs = []
    for i in range(1,int(num)+1):
        dirs.append(str(i).zfill(2))
    #[pack_pkl(st) for st in sorted(dirs)]
    [plot_pkl(st,sel) for st in sorted(dirs)]

    # matplotlib
    plt.xlabel('end-to-end distance (A)')
    plt.ylabel('average H-bond count')
    plt.title('xxmoleculexx - xxngnxx - ASMD \n \
               xxenvironxx xxvelxx A/ns')
    plt.legend()
    plt.gca().set_ylim(ymin=-0.2)
    plt.draw()

    texdir = os.path.join(('/'.join(my_dir.split('/')[0:-2])), \
                       'tex_%s/fig_bond' % my_dir.split('/')[-3])
    if not os.path.exists(texdir): os.makedirs(texdir)
    plotname = 'danamdimpvv100_%s' % sel
    plt.savefig('%s/%s.png' % (texdir,plotname))
    plt.savefig('%s/%s.eps' % (texdir,plotname))

main_bond('hb')
solv_test = my_dir.split('/')[-2].split('.')[1]
print 'solvent test',solv_test
if solv_test == 'exp':
    main_bond('wp')
