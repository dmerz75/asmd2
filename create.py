#!/usr/bin/env python
import sys,os,itertools,shutil,re,pickle,time
my_dir=os.path.abspath(os.path.dirname(__file__))
from asmd.asmdwork import *
import numpy as np
import ConfigParser

# EXAMPLE RUN COMMAND
# >>> ./create.py namd create

def parse_gconf(ngn,selection):
    ''' This function uses the Python ConfigParser module to
        read an ngn.gconf configuration file. From those configurations,
        multiple 'asmdmethods' will be constructed and placed in a working
        directory as a pickled class. 
        Each asmdmethod pickle can be expanded in its working directory to
        populate, and therefore permit acquisition, and finally analysis of
        any number of trajectories in the form of Steered or Adaptive Steered
        Molecular Dynamics.
        ngn: namd, amber, gromacs
        selection: [create],[da]
    '''
    config = ConfigParser.ConfigParser()
    # get configurations from ngn.gconf
    config.read('%s.gconf' % ngn)
    '''
    Example from namd.gconf:
    
    [create]
    molecule    = Decaalanine                   (variable _ string)
    molecule_id = da                            (variable _ string)
    job_id      = jobid                         (variable _ string)
    solvents    = vac,imp,exp                   (list of strings)
    start_coord = 13.0                          (variable _ float)
    end_to_end  = 33.0                          (variable _ float)
    extension   = 20.0                          (variable _ float)
    dct_solv_extension = vac:zcrd,imp:zcrd,exp:zcrd     (dct)
    velocities  = 100.0,10.0                    (list of floats)
    total_trajs = 24                            (list of integers)
    traj_per_dir= 3,6,12,24                     (list of integers)
    timestep    = 2.0                           (variable _ float)
    cluster     = fgatecpu2                     (variable _ string)
    ppn         = 1                             (variable _ integer)
    dct_solv_ppn= vac:1,imp:cn,exp:cn           (dct)
    walltime    = 072:00:00                     (variable _ string)
    dct_solv_walltime = vac:walltime,imp:walltime,exp:400:00:00   (dct)
    queue       = standby                       (variable _ string)
    dct_solv_queue = vac:queue,imp:queue,exp:queue     (dct)
    path_seg    = 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1    (list of floats)
    path_svel   = 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0    (list of floats)
    langevD     = 5.0                           (variable _ float)
    temperature = 300                           (variable _ integer)
    '''
    # begin config.get
    molecule    = config.get(selection,'molecule')
    molecule_id = config.get(selection,'molecule_id')
    job_id      = config.get(selection,'job_id')
    solvents    = config.get(selection,'solvents')
    start_coord = float(config.get(selection,'start_coord'))
    end_to_end  = float(config.get(selection,'end_to_end'))
    extension   = float(config.get(selection,'extension'))
    dct_solv_extension = config.get(selection,'dct_solv_extension')
    velocities  = config.get(selection,'velocities')
    total_trajs = config.get(selection,'total_trajs')
    traj_per_dir= config.get(selection,'traj_per_dir')
    timestep    = float(config.get(selection,'timestep'))
    cluster     = config.get(selection,'cluster')
    ppn         = int(config.get(selection,'ppn'))
    dct_solv_ppn= config.get(selection,'dct_solv_ppn')
    walltime    = config.get(selection,'walltime')
    dct_solv_walltime = config.get(selection,'dct_solv_walltime')
    queue       = config.get(selection,'queue')
    dct_solv_queue = config.get(selection,'dct_solv_queue')
    path_seg    = config.get(selection,'path_seg')
    path_svel   = config.get(selection,'path_svel')
    langevD     = float(config.get(selection,'langevD'))
    temperature = int(config.get(selection,'temperature'))
    # debugging by type
    def print_type(*args):
        for i,obj in enumerate(*args):
            print i,obj,type(obj)
    # CHECK: strings, floats, ints
    print_type([molecule,molecule_id,job_id,queue,\
                start_coord,end_to_end,extension,
                timestep,langevD,ppn,temperature])
    
    # convert to float,int
    def conversion(list,type='float'):
        if type == 'float':
            return_list = [float(i) for i in list]
        elif type == 'int':
            return_list = [int(i) for i in list]
        return return_list
    # list of ints
    traj_per_dir = conversion(traj_per_dir.split(','),'int')
    total_trajs = conversion(total_trajs.split(','),'int')
    # list of floats
    velocities = conversion(velocities.split(','))
    path_seg   = conversion(path_seg.split(','))
    path_svel  = conversion(path_svel.split(','))
    # list of strings
    solvents = solvents.split(',')
    # CHECK: lists
    print_type([traj_per_dir,total_trajs,velocities,path_seg,path_svel,\
                solvents])

    # build 3 dct's; solvent by extension,walltime,queue
    def build_dct(dct_parse,default_assign,place_holder,end=':00:00'):
        return_dct = {}
        for entry in dct_parse.split(','):
            key = entry.split(':')[0]
            value = entry.split(':')[1]
            if value==place_holder:
                return_dct[key]=default_assign
            else:
                if place_holder == 'walltime':
                    return_dct[key]=value + end
                elif place_holder == 'ext':
                    return_dct[key]=float(value)
                else:
                    return_dct[key]=value
        return return_dct

    dct_solv_extension = build_dct(dct_solv_extension,extension,'ext')
    dct_solv_walltime  = build_dct(dct_solv_walltime,walltime,'walltime')
    dct_solv_queue     = build_dct(dct_solv_queue,queue,'queue')
    # CHECK: dct's
    print_type([dct_solv_extension,dct_solv_walltime,dct_solv_queue])

    ''' All parameters acquired, now construct a # of asmdmethods
        asmd = project_directory/
                      working directory + asmdmethod.pkl/
                      working directory + asmdmethod.pkl/
                      working directory + asmdmethod.pkl/
                      ...
        Total number of working directories / created will be:
        number of velocities * number of solvents * number of total_trajs *
        number traj_per_dir,   2 * 3 * 1 * 4 = 24;
    '''
    def build_svtt_dct(*args):
        keys = ['solvent','velocity','total_traj','traj_per_dir']
        dct = dict(zip(keys,*args))
        print dct
        return dct
    dct_list = [build_svtt_dct([s,v,t,tpd]) for s in solvents for v in \
                velocities for t in total_trajs for tpd in \
                traj_per_dir]
    print dct_list

    # Now create the class
    class asmdmethod
    
# parse_gconf(namd,create)
parse_gconf(sys.argv[1],sys.argv[2])

'''

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
