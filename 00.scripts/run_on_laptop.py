#!/usr/bin/env python
import sys,os,glob,subprocess,time
import os.path
from glob import glob
import fnmatch
import itertools

def main():
    my_dir = os.path.abspath(os.path.dirname(__file__))
    acc=[]
    def find_job(f,vel,solv,stages):
        def job_stage(st):
            print st
            traj_dirs = []
            for path in glob(os.path.join(my_dir,'%s/*.%s/%s/%s/*/go.py'%(f,solv,vel,st))):
                traj_dirs.append(path.split('/')[-2])
                print path
            for traj_dir in sorted(traj_dirs):
                for path in glob(os.path.join(my_dir,'%s/*.%s/%s/%s/%s/go.py' \
                                                     %(f,solv,vel,st,traj_dir))):
                    fd =(path.split('/')[-6])
                    num=(path.split('/')[-4])
                    sol=(path.split('/')[-5]).split('.')[1]
                    stg=(path.split('/')[-3])
                    jtype=num+sol+stg
                    acc.append(jtype)
                    if path.split('/')[1]=='export':
                        root='/'+'/'.join(path.split('/')[2:-1])
                        path='/'+'/'.join(path.split('/')[2:])
                    else:
                        root='/'.join(path.split('/')[:-1])
                    print root
                    print path,'submitted'
                    os.chdir(root)
                    pipe=subprocess.call(['/usr/bin/python',path],stdin=subprocess.PIPE, \
                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                    #qsub_job(stg,path)
            for path in glob(os.path.join(my_dir,'%s/*.%s/%s/%s-continue.py'%(f,solv,vel,st))):
                fd =(path.split('/')[-4])
                num=(path.split('/')[-2])
                sol=(path.split('/')[-3]).split('.')[1]
                stg=(path.split('/')[-1]).split('-')[0]
                if path.split('/')[1]=='export':
                    root='/'+'/'.join(path.split('/')[2:-1])
                    path='/'+'/'.join(path.split('/')[2:])
                else:
                    root='/'.join(path.split('/')[:-1])
                print root
                print path
                os.chdir(root)
                pipe=subprocess.call(['/usr/bin/python',path,st],stdin=subprocess.PIPE, \
                        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                #qsub_jobc(stg,path)
            for path in glob(os.path.join(my_dir,'%s/*.%s/%s/00-hb_pkl.py'%(f,solv,vel))):
                fd =(path.split('/')[-4])
                num=(path.split('/')[-2])
                sol=(path.split('/')[-3]).split('.')[1]
                stg=(path.split('/')[-1]).split('-')[0]
                if path.split('/')[1]=='export':
                    root='/'+'/'.join(path.split('/')[2:-1])
                    path='/'+'/'.join(path.split('/')[2:])
                else:
                    root='/'.join(path.split('/')[:-1])
                print root
                print path
                os.chdir(root)
                pipe=subprocess.call(['/usr/bin/python',path,st],stdin=subprocess.PIPE, \
                        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                #qsub_jobh(stg,path)
        [job_stage(st) for st in stages]

    #__________________________________________________________________________
    # submitted 02:
    # submitted 03:
    #__________________________________________________________________________
    dirs = [str(d) for d in sorted([int(f) for f in os.listdir(my_dir) \
                         if os.path.isdir(f)])]
    velocities = ['02']
    #velocities = ['02','03','04','05']
    solvents   = ['vac']
    #solvents   = ['vac','imp','exp']
    stages     = [str(x).zfill(2) for x in range(1,11)]
    # stages    = ['01','02','03','04','05','06','07','08','09','10']
    # alternatively, limit stages to ['01','02','03']
    # MAIN SUBMISSION CALL
    # alternatively, qsub_job('01','vac')
    [find_job(dirs[0],v,s,stages) for s in solvents for v in velocities]
    #[find_job(f,v,s,stages) for f in dirs for v in velocities for s in solvents]
    #__________________________________________________________________________

if __name__ == "__main__":
    main()
