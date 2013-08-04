import shutil,os,itertools,re,random,fnmatch,time,pickle,sys
import numpy as np

# my_dir=os.path.abspath(os.path.dirname(__file__))

# get pylib
from pylib.cp import *
from pylib.regex import *

confign={'1':{'gpu':'nodes=1:ppn=1:gpus=1:TESLA','cpu':'nodes=1:ppn=1'},
         '2':{'gpu':'nodes=1:ppn=2:gpus=1:TESLA','cpu':'nodes=1:ppn=2'},
         '3':{'gpu':'nodes=1:ppn=3:gpus=1:TESLA','cpu':'nodes=1:ppn=3'},
         '4':{'gpu':'nodes=1:ppn=4:gpus=1:TESLA','cpu':'nodes=1:ppn=4'},
         '5':{'gpu':'nodes=1:ppn=5:gpus=1:TESLA','cpu':'nodes=1:ppn=5'},
         '6':{'gpu':'nodes=1:ppn=6:gpus=1:TESLA','cpu':'nodes=1:ppn=6'},
         '7':{'gpu':'nodes=1:ppn=7:gpus=1:TESLA','cpu':'nodes=1:ppn=7'},
         '8':{'gpu':'nodes=1:ppn=8:gpus=1:TESLA','cpu':'nodes=1:ppn=8'},
       '16':{'gpu':'nodes=1:ppn=16:gpus=1:TESLA','cpu':'nodes=1:ppn=16'}}

class AsmdMethod:
    def __init__(self,dictionary):
        for k,v in dictionary.items():
            setattr(self,k,v)
        self.dirs_per_stage = self.total_traj / self.traj_per_dir
        self.proj_dir = 'asmd_%s_%s_%s' % (self.ngn,self.molecule_id, \
                                           self.job_id)
        self.work_dir = '%s_%s_%s_tot%s_dirs%d_x%d_st%s' % (self.molecule_id, \
                        self.solvent,str(self.velocity),str(self.total_traj), \
                        self.dirs_per_stage,self.traj_per_dir, \
                        str(len(self.path_seg)))
        self.dt_ps = self.freq_force_pr*self.timestep/1000
        self.dist_seg = np.array(self.path_seg)*self.extension
        self.pv_ans = np.array(self.path_svel)*self.velocity
        self.pv_aps = self.pv_ans / 1000
        self.t_seg_ps= self.dist_seg / self.pv_aps
        self.steps_seg= self.t_seg_ps / (self.timestep/1000)
        self.v_A_timestep = self.dist_seg / self.steps_seg
        self.fpath_proj_dir = os.path.join(self.my_dir,self.proj_dir)
        self.fpath_work_dir = os.path.join(self.fpath_proj_dir,self.work_dir)
    def make_proj_and_work_dir(self):
        ''' Create a project directory:
            Create a working directory:
        '''
        if not os.path.exists(self.fpath_proj_dir):
            os.makedirs(self.fpath_proj_dir)
        if not os.path.exists(self.fpath_work_dir):
            os.makedirs(self.fpath_work_dir)
    def dump_asmdmethod_pkl():
        os.chdir(fpath_work_dir)
        # pickle.dump(
        # os.chdir(self.jdir)    # pdir = nda_jobid
        # pickle.dump(self.cfg,open('config%s_%s.pkl' %(self.vel,self.st),'w'))
        # pickle.dump(f,open('%s/AsmdMethod_%s_%s_%s.pkl' % (jdir, \
        #         env.split('.')[1],str(int(v)).zfill(2),stages),'w'))
        
'''
        self.pdir = os.path.join(self.workdir,self.packdir)
        self.jobdir = jobdir      
        # switching this to total trajs   # 100, 200
        # ^^^ from gen.py: jobdir=ngn[0][0]+molec+str(dircount)+'_'+jobid
        self.jdir = os.path.join(self.pdir,self.jobdir)
    def makeJobDir(self):
        if not os.path.exists(self.jdir):os.makedirs(self.jdir)
        texdir=os.path.join(self.jdir,'tex_%s' % self.jdir.split('/')[-1])
        if not os.path.exists(texdir): os.makedirs(texdir)
        cp_file(self.workdir,'gen.py',self.pdir,'.gen_%s.py' % \
                  self.jdir.split('/')[-1])
        cp_file(self.pydir,'pipe.py',self.pdir,'pipe.py')  # formerly jdir
        cp_file(self.mdir,'tex/tm.tex',texdir,'tm.tex')
        cp_file(self.mdir,'tex/pdflatex.sh',texdir,'pdflatex.sh')
        #cp_file(self.pydir,'del.py',self.jdir,'del.py')
        cp_file(self.pydir,'del.py',self.pdir,'del.py')
        cp_file(self.workdir,'%s.gconf' % self.ngn,self.pdir,'.%s.gconf' % \
                self.ngn)
        cp_file(os.path.join(self.workdir,'00.scripts'),'run_on_laptop.py',\
                self.pdir,'')
        if self.ngn == 'namd':    # NAMD
            if not os.path.exists(os.path.join(self.jdir,'toppar')):
                cp_tree(self.ndir,'toppar',self.jdir,'toppar')
        if not os.path.exists(os.path.join(self.jdir,'00.scripts')):
            cp_tree(self.workdir,'00.scripts',self.jdir,'00.scripts')
        return texdir
'''
