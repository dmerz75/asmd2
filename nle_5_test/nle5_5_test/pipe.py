#!/usr/bin/env python
import sys,os,glob,subprocess,time
import os.path
from glob import glob
import fnmatch
import itertools

class mdict(dict):
    def __setitem__(self,key,value):
        self.setdefault(key,[]).append(value)

my_dir = os.path.abspath(os.path.dirname(__file__))
acc=[]
qsub_path='qsub'
jdict = mdict()
cdict = mdict()

def qsub_job(stage,job_path):
    if stage==stages[0]:
        dep_args = []
        pipe=subprocess.Popen([qsub_path] + dep_args +
                          [job_path],stdin=subprocess.PIPE, \
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        stout, stderr = pipe.communicate()
        print stout
        print stout.split('.')[0]
        print 'stderr >> ',stderr
        jdict[stage]=stout.split('.')[0]
    elif stage!=stages[0]:
        entry=str(int(stage)-1).zfill(2)
        print cdict[entry]
        job_deps = ':'.join(cdict[entry])
        dep_args = ['-W depend=afterany:%s' % job_deps]
        print 'dep_args',dep_args
        pipe=subprocess.Popen([qsub_path] + dep_args +
                          [job_path],stdin=subprocess.PIPE, \
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        stout, stderr = pipe.communicate()
        print stout
        print stout.split('.')[0]
        print 'stderr >> ',stderr
        jdict[stage]=stout.split('.')[0]
def qsub_jobc(stage,job_path):
    print jdict[stage]
    job_deps = ':'.join(jdict[stage])
    dep_args = ['-W depend=afterany:%s' % job_deps]
    print 'dep_args',dep_args
    pipe=subprocess.Popen([qsub_path] + dep_args +
                      [job_path],stdin=subprocess.PIPE, \
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stout, stderr = pipe.communicate()
    print stout
    print stout.split('.')[0]
    print 'stderr >> ',stderr
    cdict[stage]=stout.split('.')[0]

def find_job(vel,solv,st):
    for path in glob(os.path.join(my_dir,'*.%s/%s/%s/*/job.sh'%(solv,vel,st))):
        num=(path.split('/')[-4])
        sol=(path.split('/')[-5]).split('.')[1]
        stg=(path.split('/')[-3])
        jtype=num+sol+stg
        acc.append(jtype)
        root='/'.join(path.split('/')[:-1])
        #root='/'+'/'.join(path.split('/')[2:-1])
        #path='/'+'/'.join(path.split('/')[2:])
        print root
        print path
        os.chdir(root)
        qsub_job(stg,path)
    for path in glob(os.path.join(my_dir,'*.%s/%s/%s-job.sh'%(solv,vel,st))):
        num=(path.split('/')[-2])
        sol=(path.split('/')[-3]).split('.')[1]
        stg=(path.split('/')[-1]).split('-')[0]
        root='/'.join(path.split('/')[:-1])
        #root='/'+'/'.join(path.split('/')[2:-1])
        #path='/'+'/'.join(path.split('/')[2:])
        print root
        print path
        os.chdir(root)
        qsub_jobc(stg,path)

# submitted 02:
# submitted 03:
velocities = ['02']
solvents   = ['vac']   #solvents = ['imp'] #solvents = ['exp'] #one at a time
stages     = ['01','02','03','04','05','06','07','08','09','10']
# MAIN SUBMISSION CALL
# alternatively, qsub_job('01','vac')
[find_job(v,s,st) for v in velocities for s in solvents for st in stages]

# informational purposes only
print len(acc)
result={}
def count():
    for cond in acc:
        if cond not in result:
            result[cond]= 0
        result[cond] += 1
count()
os.chdir(my_dir)
for key in result:
    print key + '  ' + str(result[key])
    t = key + '  ' + str(result[key])
print jdict
