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

spos=xxsposxx
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
    phase=int(stage)-1
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
    #_________________________________________________________________________
    if sel != 'ihb':     # sel ==  'wp', 'hb'
        dct_sd_hb=pickle.load(open('%s-sd_%s.pkl' % (stage,sel),'rb'))
        print '%s-sd_%s.pkl' %(stage,sel)     # pkl: sd_hb or sd_wp
        seeds = dct_sd_hb.keys()
        acclens=[]
        for s in seeds:
            acc=[]
            sample_i = dct_sd_hb[s]  # trajectory,dcd-length list with
            for c in range(len(sample_i[0])):     # width of bonds per frame
                hbc=len(sample_i[0][c]) # hbc-hydrogen-bond-count, 1 frame
                acc.append(hbc)      # acc: counts over full trajectory
            acclens.append(acc)      # acclens: all trajectories
        data = np.array(acclens).mean(axis=0)
        if stage=='01':
            d = np.linspace(spos,spos+domain[phase],data.shape[0])
        elif stage !='01':
            d = np.linspace(spos+domain[phase-1],spos+domain[phase],data.shape[0])
        # establish domain, by linspacing - vector same length as data(frames)
        acc_d.append(d)              # acc_d.append(d[2:-2])
        acc_b.append(data)           # acc_b.append(data[2:-2])
    else: # sel == 'ihb'
        dct_sd_hb=pickle.load(open('%s-sd_%s.pkl' % (stage,sel[1:3]),'rb'))
        print '%s-sd_%s.pkl' %(stage,sel[1:3])
        seeds = dct_sd_hb.keys()
        b_data = np.array([[charac_bond2(dct_sd_hb[s][0],n) for s in seeds] \
                     for n in [3,4,5]])
        '''
        if stage=='01':
            d = np.linspace(spos,spos+domain[phase],b_data.shape[2])
        elif stage !='01':
            d = np.linspace(spos+domain[phase-1],spos+domain[phase], \
                                 b_data.shape[2])
        acc_d.append(d)
        '''
        acc_b.append(b_data)
#_____________________________________________________________________________
def main_bond(sel,indx_clr=[(0,'k-','')]):
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
    # list = [0,10,20,30]
    # plt.yticks(list,fontproperties=fpropxl)
    # plt.yticks((0,10,20,30),fontproperties=fpropxl)
    plt.yticks((0,1,2,3,4,5,6),fontproperties=fpropxl)
    plt.xticks((10,15,20,25,30,35),fontproperties=fpropxl)
    # sub calls - list dirs 01,02, ... 10; call plot_pkl
    dirs = []
    for i in range(1,int(num)+1):
        dirs.append(str(i).zfill(2))
    acc_domain=[]
    acc_bond  =[]
    if sel == 'ihb':
        [plot_pkl(st,sel,acc_domain,acc_bond,indx_clr[0][0],indx_clr[0][1],\
                  indx_clr[0][2]) for st in sorted(dirs)]
        print acc_bond[0].shape
        acc_count_traj = []
        for ab in acc_bond:
            print ab.shape
            print ab.shape[1]
            acc_count_traj.append(ab.shape[1])
        min_val = min(acc_count_traj)
        print type(min_val),'min_val',min_val
        acc_count_trim_traj = []
        for ab in acc_bond:
            b = ab[:,:min_val,:]
            print b.shape
            acc_count_trim_traj.append(b)
        bnd = np.concatenate(acc_count_trim_traj,axis=2)
        print bnd.shape
        ext = np.linspace(spos,spos+dist,bnd.shape[2])
        plt.plot(ext,bnd[0,::,::].mean(axis=0),'r-',linewidth=1.0, \
                                              label=r"i$\rightarrow$i+3")
        plt.plot(ext,bnd[1,::,::].mean(axis=0),'k-',linewidth=1.0, \
                                              label=r"i$\rightarrow$i+4")
        plt.plot(ext,bnd[2,::,::].mean(axis=0),'g-',linewidth=1.0, \
                                              label=r"i$\rightarrow$i+5")
    else: # 'hb','wp'
        [plot_pkl(st,sel,acc_domain,acc_bond,indx_clr[0][0],indx_clr[0][1],\
                   ) for st in sorted(dirs)]
        ext = np.concatenate(acc_domain,axis=0)
        bnd = np.concatenate(acc_bond,axis=0)
        plt.plot(ext,bnd,'k-',linewidth=1.5,label='hydrogen bonds')
    # LEGEND
    plt.legend(loc='upper right')
    plt.gca().set_ylim(ymin=-0.2)
    leg = plt.gca().get_legend()
    leg.draw_frame(False)
    # DRAW
    plt.draw()
    # SAVE:  ../../tex_workdir/fig_bond/plotname.png|.eps
    texdir = os.path.join(('/'.join(my_dir.split('/')[0:-2])), \
                       'tex_%s/fig_bond' % my_dir.split('/')[-3])
    if not os.path.exists(texdir): os.makedirs(texdir)
    plotname = 'xxplotnamexx_%s' % sel
    plt.savefig('%s/%s.png' % (texdir,plotname))
    plt.savefig('%s/%s.eps' % (texdir,plotname))
    # matplotlib - end

#___main_call_'hb','wp','ihb'_________________________________________________
main_bond('hb')
if my_dir.split('/')[-2].split('.')[1]=='exp':
    main_bond('wp')
main_bond('ihb')
