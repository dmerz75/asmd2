#!/usr/bin/env python
import sys,os,itertools,shutil,re,pickle,time
my_dir=os.path.abspath(os.path.dirname(__file__))
from asmd.asmdwork import *
import numpy as np
import ConfigParser

# EXAMPLE

ngn  =[sys.argv[1]] # >>> ./gen.py namd da
molc =[sys.argv[2]] # >>> ./gen.py ngn  molecule

#_____________________________________________________________________________
config = ConfigParser.ConfigParser()
config.read('%s.gconf' % ngn[0])   # get configurations from ngn.gconf

molec= molc[0]                     # now it's da or da_smd
mol     = config.get(molec,'mol')  # now it's da. w/ molec = da_smd or da
print 'molec',molec,'mol',mol
zcrd    = float(config.get(molec,'zcrd')) # z-coordinate
envdists= config.get(molec,'envdist') # distance per environment
envdist = {}
for entry in envdists.split(','):
    env = entry.split(':')[0]
    if entry.split(':')[1] =='zcrd':
        envdist[env]=zcrd
    else:
        envdist[env]=float(entry.split(':')[1])
dist    = float(config.get(molec,'dist')) # total distance pulled
ts      = float(config.get(molec,'ts'))   # timestep: ~ 1,2 fs
n_conf  = config.get(molec,'n')
n       = []                              # 1,2,3,4,5
for ni in range(len(n_conf.split(','))):  # 1000,100,10,1,0.1
    n.append(float(n_conf.split(',')[ni]))
#env_cnf = config.get(mol,'environ')
env_cnf = config.get(molec,'environ')
environ = []                              # 01.vac,02.imp,03.exp
for ei in range(len(env_cnf.split(','))):
    environ.append(str(env_cnf.split(',')[ei]))
langevD = config.get(molec,'langevD')     # langevin Damping: 0.2,1,5
temp    = config.get(molec,'temp')        # temperature: 300
direct  = config.get(molec,'direct')      # direction: !! may not work yet !!
gate    = config.get(molec,'gate')        # gate: cluster name
cn      = config.get(molec,'cn')          # cn: processors per node
ppn_envs= config.get(molec,'ppn_env')     # procs per node, by environment
ppn_env = {}
for entry in ppn_envs.split(','):
    env = entry.split(':')[0]
    if entry.split(':')[1]=='cn':
        ppn_env[env]=cn
    else:
        ppn_env[env]=entry.split(':')[1]
comp    = config.get(molec,'comp')        # computation: cpu or gpu
wallt   = config.get(molec,'wallt')       # walltime: see asmdwork.py
wt_envs = config.get(molec,'wt_env')      # walltime by environment
wt_env  = {}
for entry in wt_envs.split(','):
    env = entry.split(':')[0]
    if entry.split(':')[1]=='wallt':
        wt_env[env]=wallt
    else:
        wt_env[env]=entry.split(':')[1]
queue   = config.get(molec,'queue')       # queue: in case job.sh requires queue
q_envs  = config.get(molec,'q_env')
q_env   = {}
for entry in q_envs.split(','):
    env = entry.split(':')[0]
    if entry.split(':')[1]=='queue':
        q_env[env]=queue
    else:
        q_env[env]=entry.split(':')[1]
# constructing path_seg, path_svel
p_seg   = config.get(molec,'path_seg')    # path_segmentation: decimal, sum to 1
p_svel  = config.get(molec,'path_svel')   # path scaled velocity: scaling factor
ln_spc  = config.get(molec,'lnspc')       # linspace: 'evenly discretized'
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

# print path_seg
# print path_svel

vels    = config.get(molec,'vels')          # 1,2,3,4,5
vels_l  = [int(v) for v in vels.split(',')] 
hows    = config.get(molec,'howmany')       # 5,10,10,2,1
hows_l  = [(h) for h in hows.split(',')]    # tpd: traj / directory
freqs   = config.get(molec,'freq')
freq_l  = [int(f)for f in freqs.split(',')] # tcl freq, controls dt interval

setup={}
for i in range(len(vels_l)):
    setup[vels_l[i]]={}
    setup[vels_l[i]]['howmany']=hows_l[i]
    setup[vels_l[i]]['freq']=freq_l[i]

jobid   = config.get(molec,'jobid')        # jobdir specific name
dirs    = config.get(molec,'dircounts')    # #num of directories, in which
dircounts=[d for d in dirs.split(',')]     #   tpd will be acquired

'''
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
'''

#_____CODE_____________________________________________________________________
def asmd(dps,tps,v):
    ''' asmd method arguments: 1 dircount, tps:traj/stage, v(velocity, ~100)
    '''

    def work_dir():
        ''' makes the jobdir inside pack_dir
        '''
        w = est_JobDir(ngn[0],mol,zcrd,workdir,jobdir,\
                pack_dir)
        subdir = w.makeJobDir()
        w.reg_exp(subdir)

    def make_struc(ng,mol,env,workdir,jobdir,pack_dir):
        ''' copies the appropriate starting structure according to the
            environment selected
        '''
        s = est_StrucDir(ng,mol,env,workdir,jobdir,pack_dir)
        s.makeStrucDir()

    def make_asmd(ng,mol,env,v,zc,workdir,jobdir,\
            pack_dir,tpd):
        ''' establishes the class AsmdMethod from asmdwork.py
              subcalls:
              makeEnvDir - make directory for environment - 01.vac
              makeVelDir - make directory for velocity - 02 or 03
              makeSubDir - stage dirs created - 01,02,03,04
              steering_control - tcl file or appropriate steering file
                                 placed
        '''
        #config = super_pickle(int(v))
        # umpire
        vel = 1/(100*10**v)
        path_vel=np.linspace(vel,vel,len(path_seg))*ts*path_svel
        path_steps=np.rint(path_seg*dist/path_vel)
        stages = len(path_steps)
        print 'ASMD is ready with',ng,'for',mol,'in',env,'at velocity',v, \
              'assembled inside',jobdir+'.'
        f=AsmdMethod(ng,mol,env,v,ts,zc,langevD,workdir,jobdir,pack_dir,gate, \
            ppn_env[env],comp,wt_env[env],q_env[env],stages,direct, \
            dist,path_seg,path_svel,path_vel,path_steps,freq_l[int(v)-1],dps, \
            tpd,tps,temp)
            #dist,path_seg,path_svel,path_vel,path_steps,,tpd,temp)
        jdir = os.path.join(workdir,pack_dir,jobdir)
        asmd_mod = os.path.join(workdir,'asmd')
        if not os.path.exists(os.path.join(jdir,'asmd')):
            shutil.copytree(asmd_mod,os.path.join(jdir,'asmd'))
        pickle.dump(f,open('%s/AsmdMethod_%s_%s_%s.pkl' % (jdir, \
                env.split('.')[1],str(int(v)).zfill(2),stages),'w'))
        f.makeEnvDir()
        f.makeVelDir()
        f.makeContainDir()
        #f.savePickle()  remove june 7
        f.makeSubDir()
        f.steering_control()
        
    # begin def asmd():
    pack_dir=ngn[0][0]+molec+'_'+jobid
    workdir=os.path.abspath(os.path.dirname(__file__))
    # formerly, jobdir =ngn[0][0]+molec+str(dircount)+'_'+jobid
    # declare in asmdwork.py
    jobdir = tps
    tpd = int(tps)/int(d)
    work_dir()
    [make_struc(ng,mol,env,workdir,jobdir,pack_dir) \
            for ng in ngn for env in environ]
    [make_asmd(ng,mol,env,v,envdist[env],workdir,jobdir, \
            pack_dir,tpd) for ng in ngn for env in environ]
    os.chdir(my_dir)
    return pack_dir

def tar_ball(t_dir):
    ''' pack_dir is tarred
    '''
    import tarfile
    tar = tarfile.open('%s.tar.gz' % t_dir,'w:gz')
    tar.add(t_dir)
    tar.close()
#______________________________________________________________________________
pd=[asmd(str(d),str(int(hows_l[int(v)-1])*int(d)),v) for d in dircounts for v in n]
tar_ball(pd[0]) # tarball
#[asmd(str(d)) for d in dircounts]

