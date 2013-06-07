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

count = 0
for path in glob(os.path.join(my_dir,'*-sfwf.pkl*')):
    count +=1
num = str(count).zfill(2)

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

highest_work = []

spos=xxsposxx
kb  =-0.001987
temp=xxtempxx
beta=1/(kb*temp) # 1/kb*T
quota=xxquotaxx*xxhowmanyxx

def plot_work(data,st,w_c):
    rnd = np.random.RandomState(0x2913)
    indices = np.arange(data.shape[0])
    rnd.shuffle(indices)
    plot_indices = indices[1:data.shape[0]:1]   # plot this many
    phase = int(st)-1
    if phase == 0:
        d = np.linspace(spos,spos+domain[phase],data.shape[1])
    else:
        d = np.linspace(spos+domain[phase-1],spos+domain[phase],data.shape[1])
    data[::,::,3] = data[::,::,3] + w_c[-2]
    if st=='01':
        plt.plot(d,data[0,::,3],'k',linewidth=0.4,label='Work(i)')
    for index in plot_indices:
        w_i = data[index,::,3]
        plt.plot(d,w_i,'k',linewidth=0.3)
        highest_work.append(w_i[-1])
        
def plot_pmf(data,st,w_c):
    phase = int(st)-1
    print data.shape[0]
    deltaf= np.log(np.exp(data[::,::,3]*beta).mean(axis=0))*(1/beta)
    if st=='01':
        d = np.linspace(spos,spos+domain[phase],deltaf.shape[0])
    else:
        d = np.linspace(spos+domain[phase-1],spos+domain[phase],deltaf.shape[0])
    return (d,deltaf)

quota = []
collect_deltaf = []
def main_call(st,w_c,d_cp):
    phase = int(st)-1
    acc = []
    wrk_pkl={}
    wrk_pkl= pickle.load(open('%s-sfwf.pkl' % st,'rb'))
    pcorr = pickle.load(open('%s-pmfc.pkl' % st,'rb'))
    p_arr = np.array(pcorr.values()) 
    w_c   = np.cumsum(p_arr)
    print w_c
    print len(wrk_pkl[st])
    seeds = wrk_pkl[st].keys()
    for s in seeds:
        sample_i = wrk_pkl[st][s][1]
        acc.append(sample_i)
    data = np.array(acc)
    print data.shape
    quota.append(data.shape[0])
    plot_work(data,st,w_c)
    tup_pmf = plot_pmf(data,st,w_c)
    collect_deltaf.append(tup_pmf)

w_c={}
w_c[0]=0
d_cp={}

# matplotlib begin
fig,(ax1)=plt.subplots(1)
plt.clf()
subplot(111)

# main call
dirs = []
for i in range(1,int(num)+1):
    dirs.append(str(i).zfill(2))
[main_call(st,w_c,d_cp) for st in sorted(dirs)]

pmf_2d = np.hstack(collect_deltaf)
#plt.plot(pmf_2d[0,::50],pmf_2d[1,::50],'bs')
plt.plot(pmf_2d[0,::],pmf_2d[1,::],'r-',linewidth=4.0,label='PMF')
plt.plot(pmf_2d[0,::],pmf_2d[1,::],'k--',linewidth=1.4)
pickle.dump(pmf_2d,open('pmf_2d.pkl','w'))

# matplotlib end
####  plt.title('xxmoleculexx - xxngnxx - ASMD \n xxenvironxx xxvelxx $\AA$/ns')
# FONTSIZE    xx-small,x-small,small,medium,large,x-large,xx-large
fpropxxl=matplotlib.font_manager.FontProperties(size='xx-large')
fpropxl=matplotlib.font_manager.FontProperties(size='x-large')
fpropl=matplotlib.font_manager.FontProperties(size='large')
fpropm=matplotlib.font_manager.FontProperties(size='medium')
# AXES labels
plt.xlabel('end-to-end distance ($\AA$)',fontproperties=fpropxl)
plt.ylabel('Work & PMF (kcal/mol)',fontproperties=fpropxl)
plt.xlim([spos,spos+dist])
# TICKS
#list = [0,10,20,30]
#plt.yticks(list,fontproperties=fpropxl)

print max(highest_work)
print 'highest_work'

if max(highest_work)<=15:
    plt.yticks((0,5,10,15),fontproperties=fpropl)
    plt.ylim([-4,16])
elif max(highest_work)<=23:
    plt.yticks((0,5,10,15,20),fontproperties=fpropl)
    plt.ylim([-5,23])
elif max(highest_work)<=28:
    plt.yticks((0,10,20,30),fontproperties=fpropl)
    plt.ylim([-5,32])
elif max(highest_work)<=33:
    plt.yticks((0,10,20,30),fontproperties=fpropl)
    plt.ylim([-4,35])
elif max(highest_work)<=39:
    plt.yticks((0,10,20,30,40),fontproperties=fpropl)
    plt.ylim([-4,42])
elif max(highest_work)<=47:
    plt.yticks((0,15,30,45),fontproperties=fpropl)
    plt.ylim([-4,48])
elif max(highest_work)<=54:
    plt.yticks((0,10,20,30,40,50),fontproperties=fpropl)
    plt.ylim([-4,56])
elif max(highest_work)<=64:
    plt.yticks((0,15,30,45,60),fontproperties=fpropl)
    plt.ylim([-4,66])
elif max(highest_work)>64:
    pass

plt.xticks((15,20,25,30),fontproperties=fpropl)

# LEGEND
plt.legend(loc='lower right',prop={'size':14})
leg = plt.gca().get_legend()
leg.draw_frame(False)

axis = plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
axis = plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.subplots_adjust(left=0.11,right=0.99,top=0.97,bottom=0.13)
fig.set_size_inches(7.12,4.4)

# DRAW
plt.draw()

def tex_pic(num,name):
    texdir = os.path.join(('/'.join(my_dir.split('/')[0:-3])), \
             'tex_%s/fig_pmf') % my_dir.split('/')[-4]
    if not os.path.exists(texdir): os.makedirs(texdir)
    plt.savefig('%s/%s.png' %(texdir,name))
    plt.savefig('%s/%s.eps' %(texdir,name))
def cwd_pic(num,name):
    plt.savefig('%s.png' % name)
    plt.savefig('%s.eps' % name)
def final_pic(num,name):
    plt.savefig('../../../../fig/%s.png' % name)
    plt.savefig('../../../../fig/%s.eps' % name)

lst_name=['','','','','']
velcode=my_dir.split('/')[-1]
dct_vel={'01':'1000','02':'100','03':'10','04':'1','05':'p1'}
dct_case={100:'100',150:'150',200:'200',250:'250',400:'400',600:'600',800:'800'}
if quota[0] > 99:
    n = dct_case.get(quota[0],dct_case[min(dct_case.keys(), \
                                 key=lambda k:abs(k-quota[0]))])
else:
    n = str(quota[0])

lst_name[0]=('pmf_xxngnxx_xxmoleculexx_')
lst_name[1]=my_dir.split('/')[-2].split('.')[1]+'_'
lst_name[2]='v'+dct_vel[velcode]+'_'
lst_name[3]='n'+n+'_'
if int(num) > 1:
    lst_name[4]='asmd'
else:
    lst_name[4]='smd'
print ''.join(lst_name)
name = ''.join(lst_name)

#tex_pic(num,name)
cwd_pic(num,name)
#final_pic(num,name)

