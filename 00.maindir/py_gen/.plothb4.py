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
        #plot_hb_bluedot(avg_B[::2],stage,'b.',0.1)
    wb_vecs = np.array(wght_bond.values()).mean(axis=0)
    plot_hb(wb_vecs,stage,'k-',2)
    print type(wb_vecs)
    print len(wb_vecs)
#____________________________________________________________________________
def plot_pkl(stage,sel,acc_d,acc_b,index=0,color='k-',b_label='hydrogen bonds'):
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

    if sel != 'ihb':     # sel ==  'wp', 'hb'
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
        #plot_hb(data,stage,sel)      # formerly plotted by stage
        phase=int(stage)-1            # now, append, line up domains, plot
        if stage=='01':
            d = np.linspace(spos,spos+domain[phase],data.shape[0])
        elif stage !='01':
            d = np.linspace(spos+domain[phase-1],spos+domain[phase],data.shape[0])
        acc_d.append(d[2:-2])
        acc_b.append(data[2:-2])
        if stage == '10':
            dd = np.array(acc_d)
            print dd.shape
            d3 = np.reshape(dd,dd.shape[0]*dd.shape[1])
            print d3.shape
            bb = np.array(acc_b)
            b3 = np.reshape(bb,bb.shape[0]*bb.shape[1])
            print b3.shape
            plt.plot(d3,b3,'k-',linewidth=1.5,label=b_label)
    else:                # sel == 'ihb'
        dct_sd_hb=pickle.load(open('%s-sd_%s.pkl' % (stage,sel[1:3]),'rb'))
        print '%s-sd_%s.pkl' %(stage,sel[1:3])
        seeds = dct_sd_hb.keys()
        b_data = np.array([[charac_bond2(dct_sd_hb[s][0],n) for s in seeds] \
                     for n in [3,4,5]])
        print b_data.shape
        data  = b_data[index].mean(axis=0)
        # len 3 list, with len ~100,200 length 100 lists inside, 3x100,100_long
        #plot_hb(b_data,stage,sel)
        print data.shape
        phase=int(stage)-1            # now, append, line up domains, plot
        if stage=='01':
            d = np.linspace(spos,spos+domain[phase],data.shape[0])
        elif stage !='01':
            d = np.linspace(spos+domain[phase-1],spos+domain[phase],data.shape[0])
        acc_d.append(d[2:-2])
        print len(acc_d)
        acc_b.append(data[2:-2])
        if stage == '10':
            dd = np.array(acc_d)
            d3 = np.reshape(dd,dd.shape[0]*dd.shape[1])
            print d3.shape
            bb = np.array(acc_b)
            print bb
            b3 = np.reshape(bb,bb.shape[0]*bb.shape[1])
            plt.plot(d3,b3,color,linewidth=1.5,label=b_label)

def main_bond(sel,indx_clr=(0,'k-')):
    acc_domain=[]
    acc_bond  =[]
    # matplotlib
    fig=plt.figure()
    plt.clf()
    # TITLE
    plt.title('xxmoleculexx - xxngnxx - ASMD \n xxenvironxx xxvelxx $\AA$/ns')
    # FONTSIZE
    fpropxxl=matplotlib.font_manager.FontProperties(size='xx-large')
    fpropxl=matplotlib.font_manager.FontProperties(size='x-large')
    # AXES labels
    plt.xlabel('end-to-end distance ($\AA$)',fontproperties=fpropxxl)
    plt.ylabel('Average H-bond count',fontproperties=fpropxxl)
    # TICKS
    #list = [0,10,20,30]
    #plt.yticks(list,fontproperties=fpropxl)
    #plt.yticks((0,10,20,30),fontproperties=fpropxl)
    #plt.xticks((15,20,25,30),fontproperties=fpropxl)
    # sub calls - list dirs 01,02, ... 10; call plot_pkl
    dirs = []
    for i in range(1,int(num)+1):
        dirs.append(str(i).zfill(2))
    if sel == 'ihb':
        for i in range(0,len(indx_clr)):
            acc_domain=[]
            acc_bond  =[]
            [plot_pkl(st,sel,acc_domain,acc_bond,indx_clr[i][0],indx_clr[i][1],\
                      indx_clr[i][2]) for st in sorted(dirs)]
    else:
        [plot_pkl(st,sel,acc_domain,acc_bond,index,color) for st in sorted(dirs)]
    # matplotlib - label axes on plots                      # MAIN CALL _ (2)
    # LEGEND
    plt.legend(loc='upper right')
    plt.gca().set_ylim(ymin=-0.2)
    leg = plt.gca().get_legend()
    leg.draw_frame(False)
    # DRAW
    plt.draw()

    texdir = os.path.join(('/'.join(my_dir.split('/')[0:-2])), \
                       'tex_%s/fig_bond' % my_dir.split('/')[-3])
    if not os.path.exists(texdir): os.makedirs(texdir)
    plotname = 'xxplotnamexx_%s' % sel
    plt.savefig('%s/%s.png' % (texdir,plotname))
    plt.savefig('%s/%s.eps' % (texdir,plotname))

# main call  >>>  'hb','wp','ihb' ____________________________________
main_bond('hb')                                  # MAIN CALL _ (1)
if my_dir.split('/')[-2].split('.')[1] == 'exp':
    main_bond('wp')
main_bond('ihb',[(0,'r-',r"i$\rightarrow$i+3"),(1,'k-',r"i$\rightarrow$i+4"), \
                                     (2,'g-',r"i$\rightarrow$i+5")])
