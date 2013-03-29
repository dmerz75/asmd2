#!/usr/bin/env python
import sys,os,itertools,shutil,re,pickle,time
my_dir=os.path.abspath(os.path.dirname(__file__))
from asmd.asmdwork import *
import numpy as np
import ConfigParser

dircounts=['15','30']    # 1. modify if necessary, ['1','2','5','10']
jobid=('').join(dircounts)+'_st15'       # 2. choose name

setup= {1:{'howmany':7,'freq':50},     # 3. change if necessary
        2:{'howmany':7, 'freq':50},
        3:{'howmany':7, 'freq':50},
        4:{'howmany':2,  'freq':50},
        5:{'howmany':1,  'freq':50}, }

# run ex1.  >>> ./gen.py namd da         # 4. issue command
# run ex2.  >>> ./gen.py namd da_smd
# run ex3.  >>> ./gen.py namd ee
#_____________________________________________________________________________
ngn  =[sys.argv[1]]
molc =[sys.argv[2]]

config = ConfigParser.ConfigParser()
config.read('%s.gconf' % ngn[0])

molec= molc[0]  # now it's da or da_smd
mol     = config.get(molec,'mol')  # now it's da. w/ molec = da_smd or da
zcrd    = float(config.get(molec,'zcrd'))
envdists= config.get(molec,'envdist')
envdist = {}
for entry in envdists.split(','):
    env = entry.split(':')[0]
    if entry.split(':')[1] =='zcrd':
        envdist[env]=zcrd
    else:
        envdist[env]=float(entry.split(':')[1])
dist    = float(config.get(molec,'dist'))
ts      = float(config.get(molec,'ts'))
n_conf  = config.get(molec,'n')
n       = []
for ni in range(len(n_conf.split(','))):
    n.append(float(n_conf.split(',')[ni]))
env_cnf = config.get(mol,'environ')
environ = []
for ei in range(len(env_cnf.split(','))):
    environ.append(str(env_cnf.split(',')[ei]))
langevD = config.get(molec,'langevD')
direct  = config.get(molec,'direct')
gate    = config.get(molec,'gate')
cn      = config.get(molec,'cn')
ppn_envs= config.get(molec,'ppn_env')
ppn_env = {}
for entry in ppn_envs.split(','):
    env = entry.split(':')[0]
    if entry.split(':')[1]=='cn':
        ppn_env[env]=cn
    else:
        ppn_env[env]=entry.split(':')[1]
comp    = config.get(molec,'comp')
wallt   = config.get(molec,'wallt')
wt_envs = config.get(molec,'wt_env')
wt_env  = {}
for entry in wt_envs.split(','):
    env = entry.split(':')[0]
    if entry.split(':')[1]=='wallt':
        wt_env[env]=wallt
    else:
        wt_env[env]=entry.split(':')[1]
queue   = config.get(molec,'queue')
q_envs  = config.get(molec,'q_env')
q_env   = {}
for entry in q_envs.split(','):
    env = entry.split(':')[0]
    if entry.split(':')[1]=='queue':
        q_env[env]=queue
    else:
        q_env[env]=entry.split(':')[1]
# constructing path_seg, path_svel
p_seg   = config.get(molec,'path_seg')
p_svel  = config.get(molec,'path_svel')
ln_spc  = config.get(molec,'lnspc')
path_seg = []
path_svel= []
if ln_spc == 'False':
    for i in range(len(p_seg.split(','))):
        path_seg.append(float(p_seg.split(',')[i]))
        path_svel.append(float(p_svel.split(',')[i]))
    path_seg = np.array(path_seg)
    path_svel= np.array(path_svel)
else:
    pseg = float(p_seg)
    psvel= float(p_svel)
    parts= int(1/pseg)
    path_seg = np.linspace(pseg,pseg,parts)
    path_svel= np.linspace(psvel,psvel,parts)

#_________pickle_______________________________________________________________
def super_pickle(nset):
    config = mdict()
    def construct(n):
        vel = 1/(100*10**n)
        path_vel=np.linspace(vel,vel,len(path_seg))*ts*path_svel
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
        w = a_make_JobDirSmd(ngn[0],mol,zcrd,workdir,jobdir,pack_dir)
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
    pack_dir=ngn[0][0]+molec+'_'+jobid
        # molec[0]
    workdir=os.path.abspath(os.path.dirname(__file__))
    jobdir =ngn[0][0]+molec+str(dircount)+'_'+jobid
        # molec[0]
    a_work_dir()
    [call_a_Struc(ng,mol,env,workdir,jobdir,pack_dir) for ng in ngn \
         for env in environ]
         # for mol in molec
    [call_a_Smd(ng,mol,env,v,envdist[env],workdir,jobdir,pack_dir) for ng \
         in ngn for env in environ for v in n]
         # for mol in molec
    os.chdir(my_dir)
    return pack_dir

def tar_ball(t_dir):
    import tarfile
    tar = tarfile.open('%s.tar.gz' % t_dir,'w:gz')
    tar.add(t_dir)
    tar.close()
#______________________________________________________________________________
for d in dircounts:
    pd=asmd(str(d))
tar_ball(pd)                            # tarball
#[asmd(str(d)) for d in dircounts]
