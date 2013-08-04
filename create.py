#!/usr/bin/env python
import sys,os,itertools,shutil,re,time,subprocess
import cPickle as pickle
my_dir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(my_dir)
from asmd.asmdmethod import *
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
    dct_solv_extension = vac:24.0,imp:ext,exp:ext       (dct)
    velocities  = 100.0,10.0                    (list of floats)
    total_trajs = 24                            (list of integers)
    traj_per_dir= 3,6,12,24                     (list of integers)
    timestep    = 2.0                           (variable _ float)
    cluster     = fgatecpu2                     (variable _ string)
    ppn         = 1                             (variable _ integer)
    dct_solv_ppn= vac:1,imp:ppn,exp:ppn         (dct)
    walltime    = 072:00:00                     (variable _ string)
    dct_solv_walltime = vac:walltime,imp:walltime,exp:400   (dct)
    queue       = standby                       (variable _ string)
    dct_solv_queue = vac:queue,imp:queue,exp:tg_work    (dct)
    path_seg    = 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1  (list of floats)
    path_svel   = 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0  (list of floats)
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
    hbpkl_length= float(config.get(selection,'hbpkl_length'))
    freq_force_pr= int(config.get(selection,'freq_force_pr'))
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
    total_trajs  = conversion(total_trajs.split(','),'int')
    ### freq_force_pr= conversion(freq_force_pr.split(','),'int')
    # list of floats
    velocities = conversion(velocities.split(','))
    path_seg   = conversion(path_seg.split(','))
    path_svel  = conversion(path_svel.split(','))
    ### hbpkl_length = conversion(hbpkl_length.split(','))
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
                elif place_holder == 'ppn':
                    return_dct[key]=int(value)
                else:
                    return_dct[key]=value
        return return_dct

    dct_solv_extension = build_dct(dct_solv_extension,extension,'ext')
    dct_solv_walltime  = build_dct(dct_solv_walltime,walltime,'walltime')
    dct_solv_queue     = build_dct(dct_solv_queue,queue,'queue')
    dct_solv_ppn       = build_dct(dct_solv_ppn,ppn,'ppn')
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
        ''' Creates a dictionary with the following keys:
        '''
        keys = ['solvent','velocity','total_traj','traj_per_dir']
        dct = dict(zip(keys,*args))
        return dct

    def build_params_dct(*args):
        ''' Creates a dictionary with the following keys:
            ngn,molecule,molecule_id,job_id, \
            start_coord,end_to_end,extension,dct_solv_extension, \
            timestep,cluster,ppn,dct_solv_ppn,walltime, \
            dct_solv_walltime,queue,dct_solv_queue,path_seg,path_svel, \
            langevD,temperature
        '''
        keys = ['ngn','my_dir','molecule','molecule_id','job_id', \
                'start_coord','end_to_end','extension', \
                'dct_solv_extension','timestep','cluster','ppn', \
                'dct_solv_ppn','walltime','dct_solv_walltime', \
                'queue','dct_solv_queue','path_seg','path_svel',\
                'langevD','temperature','hbpkl_length','freq_force_pr']
        dct = dict(zip(keys,*args))
        return dct

    def combine_all_asmd_params(ssvt_params,all_params):
        ''' gconf has all the parameters for a varied range of datasets.
            the dict returned from combine_all_asmd_params
            represents all the essential but with zero excess parameters stored
            as a python dictionary required to build an AsmdMethod class.
        '''
        dct_combined = dict(all_params.items() + ssvt_params.items())
        return dct_combined

    def create_AsmdMethod(dct):
        ''' AsmdMethod called with dct which contains all essential elements
            for a single asmd implementation.
        '''
        f = AsmdMethod(dct)
        print dir(f)
        attrs = vars(f)
        print ' \n'.join("%s: %s" % item for item in attrs.items())
        f.make_proj_and_work_dir()
        name_pkl = '%s/AsmdMethod_%s_%s_%s.pkl' % (f.fpath_work_dir, \
                    f.solvent,str(f.velocity),str(len(f.path_seg)))
        pickle.dump(f,open(name_pkl,'wb'))
        f.place_templates_and_structures()
        return f.fpath_proj_dir

    # call build_svtt_dct()
    svtt_list = [build_svtt_dct([s,v,t,tpd]) for s in solvents for v in \
                 velocities for t in total_trajs for tpd in \
                 traj_per_dir]
    # call_build_params_dct()
    dct_all_params = build_params_dct([ngn,my_dir,molecule,molecule_id,job_id, \
              start_coord,end_to_end,extension,dct_solv_extension, \
              timestep,cluster,ppn,dct_solv_ppn,walltime, \
              dct_solv_walltime,queue,dct_solv_queue,path_seg,path_svel, \
              langevD,temperature,hbpkl_length,freq_force_pr])
    # call combine_all_asmd_params()
    asmd_list = [combine_all_asmd_params(svtt,dct_all_params) for svtt \
                 in svtt_list]
    # call create_AsmdMethod
    proj_dirs = [create_AsmdMethod(asmd_dct) for asmd_dct in asmd_list]
    return proj_dirs
    
# parse_gconf(namd,create)
proj_dirs = parse_gconf(sys.argv[1],sys.argv[2])

def expand_pickle(script):
    print script
    pipe=subprocess.Popen(['python',script],stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdout,stderr = pipe.communicate()
    print stdout
    print 'stderr >> ',stderr

expand_pickle(os.path.join(proj_dirs[0],'build_%s.py' % sys.argv[1]))

              
















