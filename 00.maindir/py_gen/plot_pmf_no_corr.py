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
if '00' == num:
    print num,'no data acquired'
    sys.exit()

# load AsmdMethod_solv_vel_stage.pkl
# ex.: AsmdMethod_vac_02_10.pkl
solvent = my_dir.split('/')[-2].split('.')[1]
vel_dir = my_dir.split('/')[-1]
total_stages = 'xxtot_stagesxx'
asmd_pkl_name = 'AsmdMethod_%s_%s_%s.pkl' % (solvent,vel_dir,total_stages)
dir_loc_AsmdMethod_pkl = '/'.join(my_dir.split('/')[0:-2])
asmd_pkl = os.path.join(dir_loc_AsmdMethod_pkl,asmd_pkl_name)
sys.path.append(dir_loc_AsmdMethod_pkl)
from asmd.asmdwork import *
c_asmd = pickle.load(open(asmd_pkl,'r'))

print dir(c_asmd)

vel  = c_asmd.v
dist = c_asmd.dist
ts   = c_asmd.ts
path_seg   = c_asmd.path_seg
path_svel  = c_asmd.path_svel
path_vel   = c_asmd.path_vel
path_steps = c_asmd.path_steps
dt         = c_asmd.dt
path_v_aps = c_asmd.pv_aps
domain     = np.cumsum(((path_steps*ts)/1000)*path_v_aps)

print vel
print dist
print ts
print path_seg
print path_svel
print path_vel
print path_steps
print dt
print path_v_aps
print domain

highest_work = []

spos=xxsposxx
kb  =-0.001987
temp=xxtempxx
beta=1/(kb*temp) # 1/kb*T
quota=xxquotaxx*xxhowmanyxx
quota = []

#_____________________________________________________________________________
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

pmf_2d = np.transpose(np.hstack(collect_deltaf))
#plt.plot(pmf_2d[0,::50],pmf_2d[1,::50],'bs')
plt.plot(pmf_2d[::,0],pmf_2d[::,1],'r-',linewidth=4.0,label='PMF')
plt.plot(pmf_2d[::,0],pmf_2d[::,1],'k--',linewidth=1.4)
#pickle.dump(pmf_2d,open('pmf_2d.pkl','w'))

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

def save_pic_data(i,subdir,fname):
    content_dir = os.path.join('/'.join(my_dir.split('/')[0:i]),subdir)
    if not os.path.exists(content_dir): os.makedirs(content_dir)
    abs_file_name = os.path.join(content_dir,fname)
    plt.savefig('%s.png' % abs_file_name)
    plt.savefig('%s.eps' % abs_file_name)
    os.chdir(content_dir)
    pickle.dump(pmf_2d,open('%s.pkl' % fname,'w'))
    np.savetxt('%s.dat' % fname,pmf_2d,fmt=['%3.4f','%3.11f'],delimiter=' ')


lst_name=['','','','','']
velcode=my_dir.split('/')[-1]
dct_vel={'01':'1000','02':'100','03':'10','04':'1','05':'p1'}
dct_case={100:'100',150:'150',200:'200',250:'250',400:'400',600:'600',800:'800'}
if quota[0] > 91:
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
# levels back, -4:beyond,-3:default,-2:count,-1:env,'':cwd
# save_pic_data(levels_back,subdir,name)
# example: save_pic_data(-4,'fig',name)
# example: save_pic_data(-3,'',name)
save_pic_data(-3,'fig',name)
