#!/usr/bin/env python
import sys,os,pickle,shutil,fnmatch,itertools
import os.path,datetime,time
from glob import glob
from sys import argv
import numpy as np
from random import *

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
kb=-0.001987  
temp=xxtempxx
beta=kb*temp  #beta=-0.5961 @ 300
quota=xxquotaxx*xxhowmanyxx

def calc_work(data,st,w_c,pmf_c):
    data[::,3] = np.cumsum(data[::,3]*path_v_aps[int(st)-1]*dt)
    wf = data[::,3][-1]
    return data,wf
def calc_pmf(data,st,w_c,pmf_c):
    data[::,::,3] = np.cumsum(data[::,::,3]*path_v_aps[int(st)-1]*dt,axis=1)
    deltaf = np.log(np.exp(data[::,::,3]*beta).mean(axis=0))*(1/beta)
    JA = deltaf[-1]
    return JA

class asmd_calcs:
    def __init__(self,wrk_pkl,w_c,pmf_c,d_cp):
        self.wrk = wrk_pkl
        self.w_c = w_c
        self.pmf_c = pmf_c
        self.d_cp = d_cp
    def create_pkl(self,st):
        acc = []
        self.wrk[st]={}
        self.d_cp[st]={}
        stdir = os.path.join(my_dir,st)
        folds=[os.path.join(stdir,f) for f in os.listdir(stdir) \
                if os.path.isdir(os.path.join(stdir,f))]
        seeds=[p.split('/')[-1].split('.')[2] for f in folds \
                 for p in glob(os.path.join(f,'*tef.dat*'))]
        for path in glob(os.path.join(my_dir,'%s/*/*tef.dat*' % st)):
            folder = path.split('/')[-2]
            seed = path.split('/')[-1].split('.')[2]
            sample_i = np.loadtxt(path)
            # errcheck
            if len(sample_i)==xxlenarrayxx:
                self.wrk[st][seed]={}
                acc.append(sample_i)
                data_1=np.array(sample_i)
                tew,wf=calc_work(data_1,st,self.w_c,self.pmf_c)#sample_i/data => tew
                self.wrk[st][seed]=folder,tew,wf
                #os.remove(path)
            else:
                pass
                #os.remove(path)
        data=np.array(acc)
        JA=calc_pmf(data,st,self.w_c,self.pmf_c)               # get JA
        wf_sd=dict([(self.wrk[st][s][2],s) for s in seeds])
        sel_seed=wf_sd.get(JA, wf_sd[min(wf_sd.keys(), key=lambda k: abs(k-JA))])
        work_ss=self.wrk[st][sel_seed][2]
        self.d_cp[st][sel_seed]=self.wrk[st][sel_seed][0]
        print wf_sd
        print 'JA:',JA
        self.pmf_c[int(st)]=JA
        print 'selected seed',sel_seed
        print 'work for selected seed:',work_ss
        self.w_c[int(st)]=work_ss
        return seeds

def cp_file(f_dir,f,d_dir,d):
    shutil.copy(os.path.join(f_dir,f),os.path.join(d_dir,d))

def main_call(st,w_c,pmf_c,d_cp):
    print 'prelims'
    print w_c
    print pmf_c
    print d_cp
    wrk_pkl={}
    c1 = asmd_calcs(wrk_pkl,w_c,pmf_c,d_cp)
    seed_l = c1.create_pkl(st)
    nextnum=str(int(st)+1).zfill(2)
    seed_folder = d_cp[st]
    if st == num:
        for s,f in seed_folder.items():
            if not os.path.exists(os.path.join(my_dir,nextnum)):
                os.makedirs(os.path.join(my_dir,nextnum))
            cp_file(os.path.join(my_dir,st,f),'daOut.coor.%s' % s, \
                    os.path.join(my_dir,nextnum),'00.coor')
            cp_file(os.path.join(my_dir,st,f),'daOut.vel.%s' % s, \
                    os.path.join(my_dir,nextnum),'00.vel')
            if os.path.isfile(os.path.join(my_dir,nextnum,'00.vel'))==True:
                print 'COOR VEL copied into %s' % nextnum
                for path in glob(os.path.join(my_dir,'%s/*/*daOut*' % st)):
                    any_seed = path.split('/')[-1].split('.')[-1]
                    #print any_seed
                    if any_seed == s:
                        pass
                    else:
                        #os.remove(path)
                        pass
    pickle.dump(w_c,open('%s-wc.pkl' % st,'w'))
    pickle.dump(pmf_c,open('%s-pmfc.pkl' % st,'w'))
    pickle.dump(wrk_pkl,open('%s-sfwf.pkl' % st,'w'))
    print 'post-evaluations'
    print w_c
    print pmf_c
    print d_cp

def load_work(st,correction):
    dct_corrects = {}
    if st == '00':
        dct_corrects[0]=0
    else:
        dct_corrects = pickle.load(open('%s-%s.pkl' % (st,correction), 'rb'))
    return dct_corrects

# main call
w_c    = {}
w_c[0] = 0
w_c = load_work(str(int(num)-1).zfill(2),'wc')
pmf_c  = {}
pmf_c[0]=0
pmf_c  = load_work(str(int(num)-1).zfill(2),'pmfc')
d_cp   = {} # dict for stage,directory,and seed copied

main_call(num,w_c,pmf_c,d_cp)
