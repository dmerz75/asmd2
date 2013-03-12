#!/usr/bin/env python
import sys,os,itertools,shutil,re,pickle,time
my_dir=os.path.abspath(os.path.dirname(__file__))
from asmd.asmdwork import *
import numpy as np

#_____MOLECULE___configurations________________________________________________
ngn    =['namd']                           # 'namd','amb,'gro'
mlist  =['da','rda','ee','le','el','oo','ti']  # da,rda  ti:42,270,50-0.02
molec  =[mlist[0]]                         # can use [0],[1] ... [n]
zcrd   = 13.0                              # z constraint:  13,33,4, start pos.
envdist={'01.vac':zcrd,'02.imp':zcrd,'03.exp':zcrd} # i.e. '01.vac':zc7...
dist   = 20.0                              # declare a float dist: 6.0, 28.0
ts     = 2.0                               # 0.5, 1.0, 2.0
n      =[2.]                            # [1.,2.] | [4.,5.]
#xcopy  ={'1':2,'2':2,'3':3,'4':5,'5':5}    # duplication for smd() only
environ=['01.vac','02.imp','03.exp']       # ['01.vac'] | ['01.vac','03.exp']
langevD='5'                                # langevin Damping: 0.2, 1, 5
direct = 1                     # untrusted # direction
#sf     = 1                     # untrusted # scale factor: -1, 1, or 5 if el
#_____GATE_______configurations________________________________________________
gate ='fgatecpu2'   # namd                 # steele2,fgatecpu2,ggatecpu2/gpu2
                    # amb                  # multisndr2,fgatecpu2
cn   ='1'                                  # ppn request
ppn_env={'01.vac':cn,'02.imp':cn,'03.exp':cn}
comp ='cpu'                                # gpu or cpu !TESLA: always 1
wallt='mwt'                    # asmd      # swt=1.5h,mwt:4h,lwt:72h,dwt:15d
wt_env={'01.vac':wallt,'02.imp':wallt,'03.exp':wallt}
queue='standby'                            # tg_'short'72 'workq'720
q_env={'01.vac':queue,'02.imp':queue,'03.exp':queue}
                                           # 'standby-8','standby','debug'
#_____ASMD_____________________________________________________________________
path_seg  =np.array([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1])# <--Sum to 1
path_svel =np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])# <--STAGES
#path_seg =np.linspace(0.02,0.02,50)
#path_svel=np.linspace(1,1,50)
#_________pickle_______________________________________________________________
def super_pickle(nset):
    config = mdict()
    def construct(n):
        vel = 1/(100*10**n)
        path_vel = (np.linspace(vel,vel,len(path_seg)))*ts*path_svel
        return path_vel
    path_vel  = construct(float(nset))
    path_steps=np.rint(path_seg*dist/path_vel)
    config[float(nset)]=dist,ts,path_seg,path_svel,path_vel, \
            path_steps,setup[nset]
    return config
#_____mdict____________________________________________________________________
class mdict(dict):
    def __setitem__(self,key,value):
        self.setdefault(key,[]).append(value)
def  print_dict(dt):
    for key,value in dt.items():
        print key,value
        print ''
#_____CODE_____________________________________________________________________
def asmd(dircount):
    def a_work_dir():
        w = a_make_JobDirSmd(ngn[0],molec[0],zcrd,workdir,jobdir,pack_dir)
        subdir = w.a_makeJobDir()
        w.reg_exp(subdir)
    def call_a_Struc(ng,mol,env,workdir,jobdir,pack_dir):
        s = a_Struc_Dirs(ng,mol,env,workdir,jobdir,pack_dir)
        s.a_makeStrucDir()
    def call_a_Smd(ng,mol,env,v,zc,workdir,jobdir,pack_dir):
        config = super_pickle(int(v))
        stages = len(config[int(v)][0][5])
        print 'ASMD is ready with',ng,'for',mol,'in',env,'at velocity',v, \
              'assembled inside',jobdir+'.'
        f=a_Smd_Method(ng,mol,env,v,ts,zc,langevD,workdir,jobdir,pack_dir,\
              gate,ppn_env[env],comp,wt_env[env],q_env[env],dircount,stages, \
              direct,dist,config)
        f.a_makeEnvDir()
        f.a_makeVelDir()
        f.a_makeContainDir()
        f.a_savepickle()
        f.a_makeSubDir()
        f.a_steering_control()
    # asmd():
    pack_dir=ngn[0]+'_'+molec[0]+'_'+jobid
    workdir=os.path.abspath(os.path.dirname(__file__))
    jobdir =ngn[0]+'_'+molec[0]+'_'+str(int(dircount)* \
               setup[2]['howmany'])+'_'+jobid
    a_work_dir()
    [call_a_Struc(ng,mol,env,workdir,jobdir,pack_dir) for ng in ngn for mol \
         in molec for env in environ]
    [call_a_Smd(ng,mol,env,v,envdist[env],workdir,jobdir,pack_dir) for ng \
        in ngn for mol in molec for env in environ for v in n]
    os.chdir(my_dir)
    return pack_dir

def tar_ball(t_dir):
    import tarfile
    tar = tarfile.open('%s.tar.gz' % t_dir,'w:gz')
    tar.add(t_dir)
    tar.close()

#___main_calls___
setup  ={1:{'howmany':100,'freq':50},
         2:{'howmany':2,'freq':50},    # 45*18t = 810, 29*28t = 812,
         3:{'howmany':2,'freq':50},    # 20*40t = 800, 25*32t = 800
         4:{'howmany':2,'freq':50},
         5:{'howmany':1,'freq':50}}
# --> FINISH HERE         # total trajectories = dircount*setup[n]['howmany']
dircounts=['5','10','15']    # ['1','2'], jobid=12
jobid=('').join(dircounts)+'_st'+str(len(path_seg))+'_a00' # <<<~~~~~!!!
#[asmd(str(d)) for d in dircounts]

# for tarball
for d in dircounts:
    pd=asmd(str(d))
tar_ball(pd)
