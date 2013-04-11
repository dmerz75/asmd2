#!/usr/bin/env python
import sys,os,pickle,shutil,fnmatch,itertools
import os.path,datetime,time
from glob import glob
from sys import argv
import numpy as np
from random import *
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

def main():
    def load_pmf(solvent,v,method,f_older,select):
        config = pickle.load(open('%s/config%s_%s.pkl'%(f_older,v,method),'rb'))
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

        def plot_work(data,st):
            rnd = np.random.RandomState(0x2913)
            indices = np.arange(data.shape[0])
            rnd.shuffle(indices)
            plot_indices = indices[1:data.shape[0]:1]   # plot this many
            phase = int(st)-1
            if phase == 0:
                d = np.linspace(spos,spos+domain[phase],data.shape[1])
            else:
                d = np.linspace(spos+domain[phase-1],spos+domain[phase],data.shape[1])
            for index in plot_indices:
                w_i = data[index,::,3]
                plt.plot(d,w_i,'#000000',linewidth=0.3)
        def plot_pmf(data,st,c_lin):
            if st=='01':
                print 'STAGE 01',data.shape[0]
                phase = int(st)-1
                deltaf= np.log(np.exp(data[::,::,3]*beta).mean(axis=0))*(1/beta)
                if phase == 0:
                    d = np.linspace(spos,spos+domain[phase],deltaf.shape[0])
                else:
                    d = np.linspace(spos+domain[phase-1],spos+domain[phase],deltaf.shape[0])
                lb = str(data.shape[0])+' '+method
                plt.plot(d,deltaf,'%s' % c_lin,linewidth=4.0,label=lb)
                #plt.plot(d,deltaf,'k--',linewidth=1.4)
            else:
                print data.shape[0]
                phase = int(st)-1
                deltaf= np.log(np.exp(data[::,::,3]*beta).mean(axis=0))*(1/beta)
                if phase == 0:
                    d = np.linspace(spos,spos+domain[phase],deltaf.shape[0])
                else:
                    d = np.linspace(spos+domain[phase-1],spos+domain[phase],deltaf.shape[0])
                lb = st+' '+str(data.shape[0])
                plt.plot(d,deltaf,'%s' % c_lin,linewidth=4.0)
                #plt.plot(d,deltaf,'k--',linewidth=1.4)
    #_____________________________________________________________________
        def pkl2(tp,st,c_lin):
            print tp,st,c_lin
            wpkl=fsvs_dct[tp][st]
            acc = []
            wrk_pkl={}
            wrk_pkl=wpkl
            print len(wrk_pkl[st])
            seeds = wrk_pkl[st].keys()
            for s in seeds:
                sample_i = wrk_pkl[st][s][1]
                acc.append(sample_i)
            data = np.array(acc)
            print data.shape
            #plot_work(data,st)
            plot_pmf(data,st,c_lin)

        def pkl_func(tp,clin):
            [pkl2(tp,s,clin) for s in stg_list]
    #_____________________________________________________________________
        fsvs_list=[]
        fsvs_dct={}
        stg_list=[]

        def build_list(list,item):
            while item not in list:
                list.append(item)
            return list

        for path in glob(os.path.join(my_dir,'%s/*/*/*-sfwf.pkl' % f_older)):
            fold = path.split('/')[-4]
            solv = path.split('/')[-3]
            vel  = path.split('/')[-2]
            stg  = path.split('/')[-1].split('-')[0]
            fsvs_dct[(fold,solv,vel)]={}
            while stg not in stg_list:
                stg_list.append(stg)
        for path in glob(os.path.join(my_dir,'%s/*/*/*-sfwf.pkl' % f_older)):
            fold = path.split('/')[-4]
            solv = path.split('/')[-3]
            vel  = path.split('/')[-2]
            stg  = path.split('/')[-1].split('-')[0]
            sfwf_pkl = pickle.load(open(path,'rb'))
            fsvs_list =build_list(fsvs_list ,(fold,solv,vel))
            fsvs_dct[(fold,solv,vel)][stg]=sfwf_pkl

        colors=['k','b','r','c','g','y','m']
        line_sty=['-','--','-.',':',' ']
        #[pkl_func(tup,colors[select]) for tup in fsvs_list if (tup[1]==solvent) and (tup[2]==v.zfill(2))]
        lin = colors[select]+line_sty[select]
        [pkl_func(tup,lin) for tup in fsvs_list if (tup[1]==solvent) and (tup[2]==v.zfill(2))]

    # matplotlib begin______________________________________
    fig=plt.figure()
    plt.clf()
    v_dct ={'1':'1000','2':'100','3':'10','4':'1','5':'0.1'}
    plt.title('xxmoleculexx - xxngnxx - MULTI \n %s, %s $\AA$/ns' % \
            (solv_p.split('.')[1].upper(),v_dct[vel_p]))
    # FONTSIZE    xx-small,x-small,small,medium,large,x-large,xx-large
    fpropxxl=matplotlib.font_manager.FontProperties(size='xx-large')
    fpropxl=matplotlib.font_manager.FontProperties(size='x-large')
    # AXES labels
    plt.xlabel('end-to-end distance ($\AA$)',fontproperties=fpropxxl)
    plt.ylabel('PMF (kcal/mol)',fontproperties=fpropxxl)
    # TICKS
    #list = [0,10,20,30]
    #plt.yticks(list,fontproperties=fpropxl)
    plt.yticks((0,10,20,30),fontproperties=fpropxl)
    plt.xticks((15,20,25,30,),fontproperties=fpropxl)
    plt.xlim([13,33])                 # manually define
    # LEGEND
    #plt.legend(loc='lower right')
    #leg = plt.gca().get_legend()
    #leg.draw_frame(False)
    # pmf_plot____CALLS_____________________________________
    vel_p = '1'
    solv_p= '01.vac'
    #vel_p = sys.argv[1]
    #solv_p= sys.argv[2]
    # if smd, only 01 stages, use 0 in 5th spot
    load_pmf(solv_p,vel_p,'01','cwd',0)
    load_pmf(solv_p,vel_p,'10','cwd',1)
    load_pmf(solv_p,vel_p,'10','cwd',2)
    load_pmf(solv_p,vel_p,'10','cwd',3)
    # DRAW__________________________________________________
    plt.draw()

    def tex_pic(v,s):
        texdir = os.path.join(my_dir,'tex_%s/fig_pmf' % my_dir.split('/')[-1])
        if not os.path.exists(texdir): os.makedirs(texdir)
        plt.savefig('%s/danamd%s%s.png' % (texdir,s,v))
        plt.savefig('%s/danamd%s%s.eps' % (texdir,s,v))
    def continue_pic():
        plt.savefig('danamd.png') # % (texdir,n))
        plt.savefig('danamd.eps') # % (texdir,n))
    #tex_pic(vel_p,solv_p.split('.')[1])
    continue_pic()

if __name__=="__main__":
    main()
