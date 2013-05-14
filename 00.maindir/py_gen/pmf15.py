#!/usr/bin/env python
import sys,os,pickle,shutil,fnmatch,itertools
import os.path,datetime,time
from glob import glob
from sys import argv
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

def main():
    def load_pmf(solvent,v,stg_cnt,f_older,line_c,line_sty,lb,p_work='no'):
        config = pickle.load(open('%s/config%s_%s.pkl'%(f_older,v,stg_cnt),'rb'))
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
        quota=56*15
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
                #plt.plot(d,w_i,'#000000',linewidth=0.3)
        def plot_pmf(data,st,c_lin):
            if st=='01':
                print data.shape[0]
                phase = int(st)-1
                deltaf= np.log(np.exp(data[::,::,3]*beta).mean(axis=0))*(1/beta)
                if phase == 0:
                    d = np.linspace(spos,spos+domain[phase],deltaf.shape[0])
                else:
                    d = np.linspace(spos+domain[phase-1],spos+domain[phase],deltaf.shape[0])
                lb = str(data.shape[0])+' '+stg_cnt
                ##plt.plot(d,deltaf,'%s' % c_lin,linewidth=4.0,label=lb)
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
                ##plt.plot(d,deltaf,'%s' % c_lin,linewidth=4.0)
                #plt.plot(d,deltaf,'k--',linewidth=1.4)
            acc_df = []
            acc_df.append(d)
            acc_df.append(deltaf)
            stck_df = np.transpose(np.vstack(acc_df))
            return stck_df
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
            if p_work == 'yes':
                plot_work(data,st)      # 4change
            d_deltaf = plot_pmf(data,st,c_lin)
            total_pmf.append(d_deltaf)
        def pkl_func(tp,clin):
            [pkl2(tp,s,clin) for s in sort_stg]
    #_____________________________________________________________________
        fsvs_list=[]
        fsvs_dct={}
        stg_list=[]
        sort_stg=[]
        total_pmf=[]
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
        sort_stg = sorted(stg_list)
        for path in glob(os.path.join(my_dir,'%s/*/*/*-sfwf.pkl' % f_older)):
            fold = path.split('/')[-4]
            solv = path.split('/')[-3]
            vel  = path.split('/')[-2]
            stg  = path.split('/')[-1].split('-')[0]
            sfwf_pkl = pickle.load(open(path,'rb'))
            fsvs_list =build_list(fsvs_list ,(fold,solv,vel))
            fsvs_dct[(fold,solv,vel)][stg]=sfwf_pkl
        # RESOURCES
        # colors=['k','b','r','c','m','y','g','w']
        # gray = '0.75'
        # line_sty=['-','--','-.',':','.']
        # line_sty=['.','^','d','p','+','-.',':','--']
        lst_colors=['k','y','c','b','r','y','g','w','0.75']
        # COLOR CODE   ^^^   vvv  LINE STYLE
        lst_line_sty=['-','.',':','-.','--']
        lin = lst_colors[line_c]+lst_line_sty[line_sty]
        print fsvs_list
        [pkl_func(tup,lin) for tup in fsvs_list if (tup[1]==solvent) and \
                                                   (tup[2]==v.zfill(2))]
        print len(total_pmf)
        pmf = np.vstack(total_pmf)
        print pmf.shape
        print len(pmf[0]),len(pmf[1])
        print len(pmf[::,0])
        print len(pmf[::,1])
        num_pts = 48
        x_len = len(pmf[::,1])/num_pts
        counter = -1
        dom = np.linspace(13,33,num_pts+counter)
        while len(dom)!=len(pmf[::x_len,1]):
            dom = np.linspace(13,33,num_pts+counter)
            counter += 1
        print len(dom)
        print len(pmf[::x_len,1])
        '''
        plt.plot(pmf[::x_len,0],pmf[::x_len,1],'%s' % lin,linewidth=2.0, \
                                                      label=lb)
        '''
        plt.plot(dom,pmf[::x_len,1],'%s' % lin,linewidth=2.0, \
                                                      label=lb)
    #_MATPLOTLIB_begin_____________________________________________________
    fig,(ax1,ax2,ax3)=plt.subplots(3, sharex=True,sharey=True)
    plt.clf()
    v_dct ={'1':'1000','2':'100','3':'10','4':'1','5':'0.1'}
    #plt.title('DA - NAMD - MULTI')
    # FONTSIZE    xx-small,x-small,small,medium,large,x-large,xx-large
    fpropxxl=matplotlib.font_manager.FontProperties(size='xx-large')
    fpropxl=matplotlib.font_manager.FontProperties(size='x-large')
    #_pmf_plot____CALLS____________________________________________________
    v_2_s = '800' #'nda_smd_fgate_vac'
    v_2_1 = '100'
    v_2_2 = '200'
    v_2_4 = '400'
    v_2_8 = '800'
    #__
    v_3_s = 'nda_smd_v3smd_cg'
    v_3_1 = 'nda10_steelevac'
    v_3_2 = 'nda20_steelevac'
    v_3_4 = 'nda_vacfg1248'
    v_3_8 = 'nda_v3asmd2'
    #__
    i_2_s = 'ndanvt_smd_allvie'
    i_2_1 = 'i2_100'
    i_2_2 = 'i2_200'
    i_2_4 = 'i2_400'
    i_2_8 = 'i2_800'
    #__
    i_3_s = 'nda_r3_smd_r3i_smd1000'
    i_3_1 = 'nda_r310_fgate3_128'
    i_3_2 = 'nda_r320_fgate3_128'
    i_3_4 = 'nda_r3_r3i_all4g'
    i_3_8 = 'nda_r380_fgate3_128'
    #__
    e_2_s = 'ndanvt_smd_allvie'
    e_2_1 = 'ndanvt15_exp1248'
    e_2_2 = 'ndanvt29_exp1248'
    e_2_4 = 'ndanvt_e2asmdcg2'
    e_2_8 = 'ndanvt_e2asmd'
    #__
    e_3_s = 'ndanvt_smd_allvie'
    e_3_1 = 'ndanvt25_e3_1248'
    e_3_2 = 'ndanvt_stexp24'
    e_3_4 = 'ndanvt100_e3_1248'
    e_3_8 = 'ndanvt200_e3_1248'
    subplot(311)
    vel_p = '3'
    solv_p= '01.vac'
    #plt.xlabel('end-to-end distance ($\AA$)',fontproperties=fpropxxl)
    #plt.ylabel('PMF (kcal/mol)',fontproperties=fpropxxl)
    plt.yticks((0,5,10,15,20,25),fontproperties=fpropxl)
    plt.xticks((15,20,25,30),fontproperties=fpropxl)
    plt.ylim([-3,26])                 # manually define
    plt.xlim([12.0,34.0])             # manually define

    if vel_p =='2':
        # 2_100
        load_pmf(solv_p,vel_p,'01',v_2_s,0,0,'SMD 800','no')
        load_pmf(solv_p,vel_p,'10',v_2_1,1,1,'ASMD 100','no')
        load_pmf(solv_p,vel_p,'10',v_2_2,2,2,'ASMD 200','no')
        load_pmf(solv_p,vel_p,'10',v_2_4,3,3,'ASMD 400','no')
        load_pmf(solv_p,vel_p,'10',v_2_8,4,4,'ASMD 800','no')
    elif vel_p =='3':
        # 3_10
        load_pmf(solv_p,vel_p,'01',v_3_s,0,0,'SMD 800','no')
        load_pmf(solv_p,vel_p,'10',v_3_1,1,1,'ASMD 100','no')
        load_pmf(solv_p,vel_p,'10',v_3_2,2,2,'ASMD 200','no')
        load_pmf(solv_p,vel_p,'10',v_3_4,3,3,'ASMD 400','no')
        load_pmf(solv_p,vel_p,'10',v_3_8,4,4,'ASMD 800','no')
    #
    vel_p = '3'
    solv_p= '02.imp'
    load_pmf(solv_p,vel_p,'10',i_3_8,4,0,'ASMD 800 IMP','no')
    #
    plt.legend(loc='upper left',prop={'size':10})
    leg = plt.gca().get_legend()
    leg.draw_frame(False)
    axis = plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
    plt.annotate('(a)', xy=(-10,10),
                xycoords='axes points',
                horizontalalignment='right', verticalalignment='bottom',
                fontsize=12)
    #_____________________________________
    subplot(312)
    vel_p = '3'
    solv_p= '02.imp'
    #plt.xlabel('end-to-end distance ($\AA$)',fontproperties=fpropxxl)
    #plt.ylabel('PMF (kcal/mol)',fontproperties=fpropxxl)
    plt.yticks((0,5,10,15,20,25),fontproperties=fpropxl)
    plt.xticks((15,20,25,30),fontproperties=fpropxl)
    plt.ylim([-3,26.0])             # manually define
    plt.xlim([12.0,34.0])           # manually define

    if vel_p =='2':
        # 2_100
        load_pmf(solv_p,vel_p,'01',i_2_s,0,0,'SMD 800','no')
        load_pmf(solv_p,vel_p,'10',i_2_1,1,1,'ASMD 100','no')
        load_pmf(solv_p,vel_p,'10',i_2_2,2,2,'ASMD 200','no')
        load_pmf(solv_p,vel_p,'10',i_2_4,3,3,'ASMD 400','no')
        load_pmf(solv_p,vel_p,'10',i_2_8,4,4,'ASMD 800','no')
    elif vel_p =='3':
        # 3_10
        load_pmf(solv_p,vel_p,'01',i_3_s,0,0,'SMD 800','no')
        load_pmf(solv_p,vel_p,'10',i_3_1,1,1,'ASMD 100','no')
        load_pmf(solv_p,vel_p,'10',i_3_2,2,2,'ASMD 200','no')
        load_pmf(solv_p,vel_p,'10',i_3_4,3,3,'ASMD 400','no')
        load_pmf(solv_p,vel_p,'10',i_3_8,4,4,'ASMD 800','no')

    plt.legend(loc='upper left',prop={'size':10})
    leg = plt.gca().get_legend()
    leg.draw_frame(False)
    axis = plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
    plt.annotate('(b)', xy=(-10, 10),
                xycoords='axes points',
                horizontalalignment='right', verticalalignment='bottom',
                fontsize=12)
    #_______________________________________________________
    subplot(313)
    vel_p = '3'
    solv_p= '03.exp'
    #plt.xlabel('end-to-end distance ($\AA$)',fontproperties=fpropxxl)
    #plt.ylabel('PMF (kcal/mol)',fontproperties=fpropxxl)
    plt.yticks((0,5,10,15,20,25),fontproperties=fpropxl)
    plt.xticks((15,20,25,30),fontproperties=fpropxl)
    plt.ylim([-3,26.0])             # manually define
    plt.xlim([12.0,34.0])           # manually define

    if vel_p =='2':
        # 2_100
        load_pmf(solv_p,vel_p,'01',e_2_s,0,0,'SMD 800','no')
        load_pmf(solv_p,vel_p,'10',e_2_1,1,1,'ASMD 100','no')
        load_pmf(solv_p,vel_p,'10',e_2_2,2,2,'ASMD 200','no')
        load_pmf(solv_p,vel_p,'10',e_2_4,3,3,'ASMD 400','no')
        load_pmf(solv_p,vel_p,'10',e_2_8,4,4,'ASMD 800','no')
    elif vel_p =='3':
        # 3_10
        load_pmf(solv_p,vel_p,'01',e_3_s,0,0,'SMD 800','no')
        load_pmf(solv_p,vel_p,'10',e_3_1,1,1,'ASMD 100','no')
        load_pmf(solv_p,vel_p,'10',e_3_2,2,2,'ASMD 200','no')
        load_pmf(solv_p,vel_p,'10',e_3_4,3,3,'ASMD 400','no')
        load_pmf(solv_p,vel_p,'10',e_3_8,4,4,'ASMD 800','no')
    #
    vel_p = '3'
    solv_p= '02.imp'
    load_pmf(solv_p,vel_p,'10',i_3_8,4,0,'ASMD 800 IMP','no')
    #
    plt.legend(loc='upper left',prop={'size':10})
    leg = plt.gca().get_legend()
    leg.draw_frame(False)
    axis = plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
    plt.annotate('(c)', xy=(-10, 10),
                xycoords='axes points',
                horizontalalignment='right', verticalalignment='bottom',
                fontsize=12)
    #_END_SUBPLOT___________________________________________
    #_______________________________________________________
    # EXTRA
    fig.subplots_adjust(hspace=0)
    fig.set_size_inches(5.18,9.6)
    #plt.xlabel('end-to-end distance ($\AA$)',fontproperties=fpropxxl)
    #plt.ylabel('PMF (kcal/mol)',fontproperties=fpropxxl)
    fig.text(0.55,0.04,'end-to-end distance ($\AA$)',ha='center',va='center',fontproperties=fpropxxl)
    fig.text(0.06,0.5,'PMF (kcal/mol)',ha='center',va='center', \
                     rotation='vertical',fontproperties=fpropxxl)
    plt.subplots_adjust(left=0.175,right=0.99,top=0.99,bottom=0.085)
    # DRAW__________________________________________________
    plt.draw()
    #plt.show()
    def tex_pic():
        texdir = os.path.join(my_dir,'tex_%s/fig_pmf' % my_dir.split('/')[-1])
        if not os.path.exists(texdir): os.makedirs(texdir)
        plt.savefig('%s/pmf15.png' % (texdir))
        plt.savefig('%s/pmf15.eps' % (texdir))
    def continue_pic():
        plt.savefig('pmf15_%s.png' % v_dct[vel_p])
        plt.savefig('pmf15_%s.eps' % v_dct[vel_p])
    tex_pic()
    #continue_pic()

if __name__=="__main__":
    main()
