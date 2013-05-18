#!/usr/bin/env python
import sys,os,fnmatch,itertools,pickle,re
import os.path,datetime,time
from glob import glob
import numpy as np
from random import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from pylab import *

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
#num =str(len(path_steps)).zfill(2)
count = 0
for path in glob(os.path.join(my_dir,'*-sd_hb.pkl*')):
    count+=1
num = str(count).zfill(2)
quota = []
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
def iter_pickle(filename):
    with open(filename) as fp:
        while True:
            try:
                entry = pickle.load(fp)
            except EOFError:
                break
            yield entry

def find_seed(seq,seed):
    for dct in seq:
        if dct['seed']==seed:
            return dct
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
        dct_sd_hb={}
        pfile = '%s-sd_%s.pkl' %(stage,sel)     # pkl: sd_hb or sd_wp
        for hb_seq in iter_pickle(pfile):
            dct_sd_hb.update(hb_seq)
        seeds = dct_sd_hb.keys()
        print seeds
        quota.append(len(seeds))
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
        pfile = '%s-sd_%s.pkl' %(stage,sel[1:3])
        dct_sd_hb={}
        for hb_seq in iter_pickle(pfile):
            dct_sd_hb.update(hb_seq)
        seeds = dct_sd_hb.keys()
        print seeds
        quota.append(len(seeds))
        b_data = np.array([[charac_bond2(dct_sd_hb[s][0],n) for s in seeds] \
                     for n in [3,4,5]])
        acc_b.append(b_data)
#_____________________________________________________________________________
def main_bond(sel,indx_clr=[(0,'k-','')]):
    # matplotlib
    fig,(ax1)=plt.subplots(1)
    plt.clf()
    subplot(111)
    # FONTSIZE
    fpropxxl=matplotlib.font_manager.FontProperties(size='xx-large')
    fpropxl=matplotlib.font_manager.FontProperties(size='x-large')
    fpropl=matplotlib.font_manager.FontProperties(size='large')
    fpropm=matplotlib.font_manager.FontProperties(size='medium')
    # AXES labels
    plt.xlabel('end-to-end distance ($\AA$)',fontproperties=fpropxl)
    plt.ylabel('Average H-bond count',fontproperties=fpropxl)
    # TICKS
    # list = [0,10,20,30]
    # plt.yticks(list,fontproperties=fpropxl)
    # plt.yticks((0,10,20,30),fontproperties=fpropxl)
    plt.yticks((0,1,2,3,4,5,6),fontproperties=fpropl)
    plt.xticks((15,20,25,30),fontproperties=fpropl)
    plt.xlim([spos,spos+dist])
    #plt.ylim([-.1,6.5])
    # sub calls - list dirs 01,02, ... 10; call plot_pkl
    dirs = []
    for i in range(1,int(num)+1):
        dirs.append(str(i).zfill(2))
    acc_domain=[]
    acc_bond  =[]
    sel = 'hb'
    #else: # 'hb','wp'
    if sel == 'hb':
        [plot_pkl(st,sel,acc_domain,acc_bond,indx_clr[0][0],indx_clr[0][1],\
                   ) for st in sorted(dirs)]
        print acc_bond[0].shape
        acc_count_traj = []
        for ab in acc_bond:
            print ab.shape
            acc_count_traj.append(ab.shape[0])
        min_val = min(acc_count_traj)
        print type(min_val),'min_val',min_val
        acc_count_trim_traj = []
        for ab in acc_bond:
            print ab.shape
            b = ab[:min_val]
            print b.shape
            acc_count_trim_traj.append(b)
        bnd2 = np.concatenate(acc_bond,axis=0)
        y_d = len(bnd2)/100
        ext2 = np.linspace(spos,spos+dist,bnd2.shape[0])
        x_d = len(ext2)/100
        plt.plot(ext2[::x_d],bnd2[::y_d],'k-',linewidth=1.0,label='hydrogen bonds')
    ####
    acc_domain=[]
    acc_bond  =[]
    sel = 'ihb'
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
        y_d = bnd.shape[2]/100
        ext = np.linspace(spos,spos+dist,bnd.shape[2])
        x_d = len(ext)/100
        print bnd.shape,'x',x_d,'y',y_d
        plt.plot(ext[::x_d],bnd[0,::,::y_d].mean(axis=0),'r-',linewidth=1.0, \
                                              label=r"i$\rightarrow$i+3")
        plt.plot(ext[::x_d],bnd[1,::,::y_d].mean(axis=0),'b-',linewidth=1.0, \
                                              label=r"i$\rightarrow$i+4")
        plt.plot(ext[::x_d],bnd[2,::,::y_d].mean(axis=0),'g-',linewidth=1.0, \
                                              label=r"i$\rightarrow$i+5")
    sel = 'hb'
    # LEGEND
    plt.legend(loc='upper right',prop={'size':12})
    plt.gca().set_ylim(ymin=-0.3)
    #plt.gca().set_ylim(ymax=+0.3) #doesn't seem to work, may18
    leg = plt.gca().get_legend()
    leg.draw_frame(False)

    plt.subplots_adjust(left=0.11,right=0.99,top=0.97,bottom=0.16)
    fig.set_size_inches(6.48,4)
    # DRAW
    plt.draw()

    '''
    def tex_pic(num,name):
        # SAVE:  ../../tex_workdir/fig_bond/plotname.png|.eps
        texdir = os.path.join(('/'.join(my_dir.split('/')[0:-3])), \
                       'tex_%s/fig_bond' % my_dir.split('/')[-4])
        if not os.path.exists(texdir): os.makedirs(texdir)
        plotname = 'danvtnamdexpvv1000asmd_%s' % sel
        plt.savefig('%s/%s.png' % (texdir,plotname))
        plt.savefig('%s/%s.eps' % (texdir,plotname))
    # matplotlib - end
    '''

    def tex_pic(num,name):
        texdir = os.path.join(('/'.join(my_dir.split('/')[0:-3])), \
                'tex_%s/fig_bond') % my_dir.split('/')[-4]
        if not os.path.exists(texdir): os.makedirs(texdir)
        plt.savefig('%s/%s.png' %(texdir,name))
        plt.savefig('%s/%s.eps' %(texdir,name))
    def cwd_pic(num,name):
        plt.savefig('%s.png' % name)
        plt.savefig('%s.eps' % name)
    def final_pic(num,name):
        plt.savefig('../../../../fig/%s.png' % name)
        plt.savefig('../../../../fig/%s.eps' % name)

    lst_name=['','','','','','']
    velcode=my_dir.split('/')[-1]
    dct_vel={'01':'1000','02':'100','03':'10','04':'1','05':'p1'}
    dct_case={100:'100',150:'150',200:'200',250:'250',400:'400',600:'600',800:'800'}
    print quota
    if quota[0] > 99:
        n = dct_case.get(quota[0],dct_case[min(dct_case.keys(), \
                                 key=lambda k:abs(k-quota[0]))])
    else:
        n = str(quota[0])
    # TITLE
    #plt.title('xxmoleculexx - xxngnxx - ASMD \n xxenvironxx xxvelxx $\AA$/ns')
    lst_name[0]=('bond_xxngnxx_xxmoleculexx_')
    lst_name[1]=my_dir.split('/')[-2].split('.')[1]+'_'
    lst_name[2]='v'+dct_vel[velcode]+'_'
    lst_name[3]='n'+n+'_'
    if int(num) > 1:
        lst_name[4]='asmd_'
    else:
        lst_name[4]='smd_'
    lst_name[5]=sel
    print ''.join(lst_name)
    name = ''.join(lst_name)

    tex_pic(num,name)
    #cwd_pic(num,name)
    #final_pic(num,name)

#___main_call_'hb','wp','ihb'_________________________________________________
main_bond('hb')
main_bond('ihb')
if my_dir.split('/')[-2].split('.')[1]=='exp':
    main_bond('wp')
