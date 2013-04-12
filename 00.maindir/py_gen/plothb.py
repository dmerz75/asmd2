#!/usr/bin/env python
import sys,os,fnmatch,itertools,pickle,re
import os.path,datetime,time
from glob import glob
import numpy as np
from random import *
import matplotlib
import matplotlib.pyplot as plt

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

spos=13.0
beta=-0.6
num =str(len(path_steps)).zfill(2)
#_____________________________________________________________________________

def plot_hb(data,stage,sel):
    phase = int(stage)-1
    if sel != 'ihb':
        if stage =='01':
            d = np.linspace(spos,spos+domain[phase],data.shape[0])
            print d[0:-1:len(d)]
            plt.plot(d,data,'k-',label="hydrogen bonds",linewidth=1)
        elif stage !='01':
            d = np.linspace(spos+domain[phase-1],spos+domain[phase],data.shape[0])
            plt.plot(d,data,'k-',linewidth=1)
            print d[0:-1:len(d)]
    if sel == 'ihb':
        if stage=='01':
            d = np.linspace(spos,spos+domain[phase],len(data.mean(axis=1)[0]))
            plt.plot(d,data.mean(axis=1)[0],'r-',label="i->i+3",linewidth=1.5)
            d = np.linspace(spos,spos+domain[phase],len(data.mean(axis=1)[1]))
            plt.plot(d,data.mean(axis=1)[1],'k-',label="i->i+4",linewidth=1.5)
            d = np.linspace(spos,spos+domain[phase],len(data.mean(axis=1)[2]))
            plt.plot(d,data.mean(axis=1)[2],'g-',label="i->i+5",linewidth=1.5)
        if stage!='01':
            d = np.linspace(spos+domain[phase-1],spos+domain[phase], \
                                             len(data.mean(axis=1)[0]))
            plt.plot(d,data.mean(axis=1)[0],'r-',linewidth=1.5)
            d = np.linspace(spos+domain[phase-1],spos+domain[phase], \
                                             len(data.mean(axis=1)[1]))
            plt.plot(d,data.mean(axis=1)[1],'k-',linewidth=1.5)
            d = np.linspace(spos+domain[phase-1],spos+domain[phase], \
                                             len(data.mean(axis=1)[2]))
            plt.plot(d,data.mean(axis=1)[2],'g-',linewidth=1.5)
def w_plot_hb(avgB,st,color,lw):
    phase = int(st)-1
    if (st=='01'):
        d = np.linspace(spos,spos+domain[phase],avgB.shape[0])
        plt.plot(d,avgB,color,label="hydrogen bonds",linewidth=lw)
    elif st !='01':
        d = np.linspace(spos+domain[phase-1],spos+domain[phase],avgB.shape[0])
        plt.plot(d,avgB,color,linewidth=lw)
def w_plot_hb_bluedot(avgB,st,color,lw):
    phase = int(st)-1
    if (st=='01'):
        d = np.linspace(spos,spos+domain[phase],avgB.shape[0])
        plt.plot(d,avgB,color,linewidth=lw)
    elif st !='01':
        d = np.linspace(spos+domain[phase-1],spos+domain[phase],avgB.shape[0])
        plt.plot(d,avgB,color,linewidth=lw)

#___________________________________________________________________________
def pack(stage):
    seed_bond={}
    wght_bond={}
    wrk_pkl={}
    wrk_pkl=pickle.load(open('%s-sfwf.pkl' % stage,'rb'))
    for path in glob(os.path.join(my_dir,'%s/*/*-hb_pr*pr*.pkl.*' % stage)):
        seed = path.split('.')[-1]
        sample_i = pickle.load(open(path,'rb'))
        bond_clist=[]
        for i in range(len(sample_i)):
            cnt = len(sample_i[i])
            bond_clist.append(cnt)
        seed_bond[seed]=np.array(bond_clist)
    seeds = wrk_pkl[stage].keys()
    for s in seeds:
        sample_w = np.exp(wrk_pkl[stage][s][1][::,3]*beta).astype(float)
        sample_b = (seed_bond[s]).astype(float)
        lenf_w = len(sample_w)/100
        lenf_b = len(sample_b)/100
        print lenf_w,lenf_b
        B_list=[]
        W_list=[]
        for b in range(len(sample_b)):
            wv = int(((b+1)/lenf_b)*lenf_w)
            sum_B=(sample_b[b]*np.exp(beta*sample_w[wv]))
            sum_W=(np.exp(beta*sample_w[wv]))
            print sum_B,sum_W
            B_list.append(sum_B)
            W_list.append(sum_W)
        avg_B=np.array(B_list).cumsum()/np.array(W_list).cumsum()
        wght_bond[s]=avg_B
        plot_hb_bluedot(avg_B[::2],stage,'b.',0.1)
    wb_vecs = np.array(wght_bond.values()).mean(axis=0)
    plot_hb(wb_vecs,stage,'k-',2)
    print type(wb_vecs)
    print len(wb_vecs)
#____________________________________________________________________________
def plot_pkl(stage,sel):
    def residue_index(label):
        return int(re.sub("[^0-9]","",label))
    def charac_bond2(trajectory,distance_target):
        acc_count_frames = []
        for frame in trajectory:
            acc_count = 0
            for bond in frame:
                distance = residue_index(bond[2])-residue_index(bond[3])
                if distance == distance_target:
                    acc_count += 1
            acc_count_frames.append(acc_count)
        return acc_count_frames

    if sel != 'ihb':
        dct_sd_hb=pickle.load(open('%s-sd_%s.pkl' % (stage,sel),'rb'))
        print '%s-sd_%s.pkl' %(stage,sel)
        seeds = dct_sd_hb.keys()
        traj_s= dct_sd_hb.values()
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
        plot_hb(data,stage,sel)
    else:
        dct_sd_hb=pickle.load(open('%s-sd_%s.pkl' % (stage,sel[1:3]),'rb'))
        print '%s-sd_%s.pkl' %(stage,sel[1:3])
        seeds = dct_sd_hb.keys()
        b_data = np.array([[charac_bond2(dct_sd_hb[s][0],n) for s in seeds] \
                     for n in [3,4,5]])
        # len 3 list, with len ~100,200 length 100 lists inside, 3x100,100_long
        plot_hb(b_data,stage,sel)

def main_bond(sel):
    # matplotlib
    fig=plt.figure()
    plt.clf()
    # sub calls - list dirs 01,02, ... 10; call plot_pkl
    dirs = []
    for i in range(1,int(num)+1):
        dirs.append(str(i).zfill(2))
    [plot_pkl(st,sel) for st in sorted(dirs)]
    # matplotlib - label axes on plots
    plt.xlabel('end-to-end distance (A)')
    plt.ylabel('average H-bond count')
    plt.title('xxmoleculexx - xxngnxx - ASMD \n xxenvironxx xxvelxx A/ns')
    plt.legend()
    plt.gca().set_ylim(ymin=-0.2)
    plt.draw()
    texdir = os.path.join(('/'.join(my_dir.split('/')[0:-2])), \
                       'tex_%s/fig_bond' % my_dir.split('/')[-3])
    if not os.path.exists(texdir): os.makedirs(texdir)
    plotname = 'xxplotnamexx_%s' % sel
    plt.savefig('%s/%s.png' % (texdir,plotname))
    plt.savefig('%s/%s.eps' % (texdir,plotname))
# main call  >>>  'hb','wp','ihb' ____________________________________
main_bond('hb')
if my_dir.split('/')[-2].split('.')[1] == 'exp':
    main_bond('wp')
main_bond('ihb')
