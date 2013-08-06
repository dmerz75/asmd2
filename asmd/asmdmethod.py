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

def reg_express(script,stage,i):
    phase = i-1
    reg_ex(proj_dir,script,'xxcstagexx',stage)

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
        self.fpath_script_dir = os.path.join(self.my_dir,'00.scripts')
        self.fpath_main_dir = os.path.join(self.my_dir,'00.maindir')
        self.fpath_py_gen   = os.path.join(self.fpath_main_dir,'py_gen')
        self.temp           = os.path.join(self.fpath_proj_dir,'templates')
        self.st = len(self.path_seg)
    def make_proj_and_work_dir(self):
        ''' Create a project directory:
            Create a working directory:
        '''
        if not os.path.exists(self.fpath_proj_dir):
            os.makedirs(self.fpath_proj_dir)
        if not os.path.exists(self.fpath_work_dir):
            os.makedirs(self.fpath_work_dir)
        if not os.path.exists(self.temp):
            os.makedirs(self.temp)
    def place_templates_and_structures(self):
        ''' Copy asmd running scripts, pbs resource managers,
            and template directories
        '''
        cp_file(self.fpath_py_gen,'pipe.py',self.fpath_proj_dir,'pipe.py')
        cp_file(self.fpath_py_gen,'del.py',self.fpath_proj_dir,'del.py')
        cp_file(self.my_dir,'%s.gconf' % self.ngn, \
                self.fpath_proj_dir,'.%s.gconf' % self.ngn)
        cp_file(self.fpath_script_dir,'run_on_laptop.py',\
                self.fpath_proj_dir,'')
        if self.ngn == 'namd':
            if not os.path.exists(os.path.join(self.fpath_proj_dir,'toppar')):
                cp_tree(os.path.join(self.fpath_main_dir,self.ngn),'toppar', \
                        self.fpath_proj_dir,'toppar')
                cp_file(self.fpath_script_dir,'build_%s.py' % self.ngn,\
                        self.fpath_proj_dir,'build_%s.py' % self.ngn)
        if not os.path.exists(os.path.join(self.fpath_proj_dir,'00.scripts')):
            cp_tree(self.fpath_script_dir,'',self.fpath_proj_dir,'00.scripts')
        if not os.path.exists(os.path.join(self.fpath_proj_dir,'00.struc')):
            cp_tree(os.path.join(self.fpath_main_dir,self.ngn),'struc/%s' \
                    % self.molecule_id,self.fpath_proj_dir,'00.struc')
        if not os.path.exists(os.path.join(self.fpath_proj_dir,'asmd')):
            cp_tree(self.my_dir,'asmd',self.fpath_proj_dir,'asmd')
        if not os.path.exists(os.path.join(self.fpath_proj_dir,'pylib')):
            cp_tree(self.my_dir,'pylib',self.fpath_proj_dir,'pylib')
        # build template dir
        cp_file('%s/%s/%s' % (self.fpath_main_dir,self.ngn,'continue'), \
                'pmf_no_corr.py',self.temp,'continue.py')
        cp_file('%s/%s/%s' % (self.fpath_main_dir,self.ngn,'go'), \
                'go-%s.py' % self.cluster,self.temp,'go.py')
        cp_file('%s/%s/%s' % (self.fpath_main_dir,self.ngn,'hb'), \
                'hb.py',self.temp,'hb.py')
        cp_file('%s/%s/%s' % (self.fpath_main_dir,self.ngn,'hb_pkl'), \
                'hb_pkl.py',self.temp,'hb_pkl.py')
        cp_file('%s/%s/%s' % (self.fpath_main_dir,self.ngn,'job'), \
                'job-%s.sh' % self.cluster,self.temp,'job.sh')
        cp_file('%s/%s/%s' % (self.fpath_main_dir,self.ngn,'jobc'), \
                'job-%s.sh' % self.cluster,self.temp,'jobc.sh')
        cp_file('%s/%s/%s' % (self.fpath_main_dir,self.ngn,'jobhb'), \
                'job-%s.sh' % self.cluster,self.temp,'jobhb.sh')
        cp_file(self.fpath_py_gen,'plothb.py',self.temp,'plothb.py')
        cp_file(self.fpath_py_gen,'plot_pmf_no_corr.py',self.temp,\
                'plotpkl.py')
        if self.ngn == 'namd':
            if not os.path.isdir(os.path.join(self.temp,self.solvent)):
                cp_tree('%s/%s/%s/%s/%s' % (self.fpath_main_dir,self.ngn, \
                    'mol.conf.tcl',self.molecule_id,self.solvent),
                    '',self.temp,self.solvent)
    # called from build_namd.py
    # run build_namd.py to continue from here
    def populate_work_dir(self):
        ''' da_exp... /   03   /   001    /   ...
            work_level:   ^  (inside work)
            stage_level:           ^  (inside a stage, 04)
            traj_level:      (inside traj folder) ^             
        '''
        # def reg_express(script,stage,i):
        #     phase = i-1
        #     reg_ex(proj_dir,script,'xxcstagexx',stage)
        cwd = os.getcwd()
        proj_dir = ('/').join(cwd.split('/')[0:-1])
        dir_temp = os.path.join(proj_dir,'templates')
        dir_scripts = os.path.join(proj_dir,'00.scripts')
        # print cwd
        # print proj_dir
        # print dir_temp
        # print dir_scripts
        def traj_level(tdir,stage):
            ''' go.py,job.sh,smdforce.tcl,smd.namd
            '''
            # tdir = os.path.join(stage_dir,traj_dirname)
            # print dir_temp
            # print tdir
            dir_temp_solv = os.path.join(dir_temp,self.solvent)
            cp_file(dir_temp,'go.py',tdir,'go.py')
            reg_ex(tdir,'go.py','xxhowmanyxx','BLANK')
            reg_ex(tdir,'go.py','xxquotaxx','BLANK')
            reg_ex(tdir,'go.py','xxnodecountxx','BLANK')
            cp_file(dir_temp,'job.sh',tdir,'job.sh')
            reg_ex(tdir,'job.sh','xxjobnamexx','BLANK')
            reg_ex(tdir,'job.sh','xxwalltimexx','BLANK')
            reg_ex(tdir,'job.sh','xxnodesxx','BLANK')
            cp_file(dir_temp_solv,'smd_force.tcl',tdir,'smdforce.tcl')
            reg_ex(tdir,'smdforce.tcl','xxfreqxx','BLANK')
            reg_ex(tdir,'smdforce.tcl','xxzcoordxx','BLANK')
            reg_ex(tdir,'smdforce.tcl','xxcur_zxx','BLANK')
            reg_ex(tdir,'smdforce.tcl','xxvelocityxx','BLANK')
            reg_ex(tdir,'smdforce.tcl','xxtsxx','BLANK')
            if stage == '01':
                cp_file(dir_temp_solv,'smd_initial.namd',tdir,'smd.namd')
            else:
                cp_file(dir_temp_solv,'smd_continue.namd',tdir,'smd.namd')
            reg_ex(tdir,'smd.namd','xxtsxx','BLANK')
            reg_ex(tdir,'smd.namd','xxlDxx','BLANK')
            reg_ex(tdir,'smd.namd','xxtempxx','BLANK')
            reg_ex(tdir,'smd.namd','xxdcdxx','BLANK')
            reg_ex(tdir,'smd.namd','xxstepsxx','BLANK')
        def stage_level(stage_dir,i,seeds):
            ''' 000, 001, 002 directories
                01.txt
                hb.py
            '''
            stage = str(i).zfill(2)
            # print stage_dir
            # print cwd
            # print proj_dir
            # print dir_temp
            # print dir_scripts
            # print seeds.shape
            os.chdir(stage_dir)
            with file('%s.txt' % stage ,'w') as outfile:
                if seeds.shape[1]==1:
                    np.savetxt(outfile,seeds,fmt='%5.0d')
                else:
                    np.savetxt(outfile,np.transpose(seeds),fmt='%5.0d')
            cp_file(dir_temp,'hb.py',stage_dir,'hb.py')
            reg_ex(stage_dir,'hb.py','xxenvironxx','BLANK')
            for t in range(0,self.dirs_per_stage):
                tdir = os.path.join(stage_dir,str(t).zfill(3))
                if not os.path.exists(tdir): os.makedirs(tdir)
                traj_level(tdir,stage)
        def work_level():
            ''' 01,02,03,04,05 directories
                01,02,03 ... - continue.py, jobhb.sh, job.sh
                00-hb_pkl.py,plothb.py,plotpkl.py,seeds.txt
            '''
            def gen_all_seeds():
                sds = np.random.randint(10000,high=99999, \
                         size=(self.st,self.dirs_per_stage,self.traj_per_dir))
                fname = 'seeds.txt'
                with file('seeds.txt','w') as outfile:
                    outfile.write('# %s stages' % sds.shape[0])
                    outfile.write(" %s directories per stage" % sds.shape[1])
                    outfile.write(" %s trajectories per directory\n" % \
                                  sds.shape[2])
                    outfile.write('# {0}\n'.format(sds.shape))
                    for stg_slice in sds:
                        ss = np.transpose(stg_slice)
                        outfile.write('# stage\n')
                        np.savetxt(outfile,ss,fmt='%5.0d')
                return sds
            # 00-hb_pkl.py,plothb.py,plotpkl.py,seeds.txt
            cp_file(dir_temp,'hb_pkl.py',cwd,'00-hb_pkl.py')
            reg_ex(cwd,'00-hb_pkl.py','xxlenbpklxx','BLANK')
            cp_file(dir_temp,'plothb.py',cwd,'plothb.py')
            reg_ex(cwd,'plothb.py','xxtot_stagesxx','BLANK')
            reg_ex(cwd,'plothb.py','xxsposxx','BLANK')
            reg_ex(cwd,'plothb.py','xxtempxx','BLANK')
            reg_ex(cwd,'plothb.py','xxquotaxx','BLANK')
            reg_ex(cwd,'plothb.py','xxhowmanyxx','BLANK')
            reg_ex(cwd,'plothb.py','xxmoleculexx','BLANK')
            reg_ex(cwd,'plothb.py','xxngnxx','BLANK')
            reg_ex(cwd,'plothb.py','xxenvironxx','BLANK')
            reg_ex(cwd,'plothb.py','xxvelxx','BLANK')
            cp_file(dir_temp,'plotpkl.py',cwd,'plotpkl.py')
            reg_ex(cwd,'plotpkl.py','xxtot_stagesxx','BLANK')
            reg_ex(cwd,'plotpkl.py','xxsposxx','BLANK')
            reg_ex(cwd,'plotpkl.py','xxtempxx','BLANK')
            reg_ex(cwd,'plotpkl.py','xxquotaxx','BLANK')
            reg_ex(cwd,'plotpkl.py','xxhowmanyxx','BLANK')
            reg_ex(cwd,'plotpkl.py','xxmoleculexx','BLANK')
            reg_ex(cwd,'plotpkl.py','xxngnxx','BLANK')
            reg_ex(cwd,'plotpkl.py','xxenvironxx','BLANK')
            reg_ex(cwd,'plotpkl.py','xxvelxx','BLANK')
            seeds = gen_all_seeds()
            # now make directories 01,02,03 ...
            # and copy 01,02,03 ... - continue.py, jobhb.sh, job.sh
            for i in range(1,self.st+1):
                continue_script = '%s-continue.py' % str(i).zfill(2)
                cp_file(dir_temp,'continue.py',cwd,\
                        continue_script)
                reg_ex(cwd,continue_script,'xxsposxx','BLANK')
                reg_ex(cwd,continue_script,'xxtempxx','BLANK')
                reg_ex(cwd,continue_script,'xxquotaxx','BLANK')
                reg_ex(cwd,continue_script,'xxhowmanyxx','BLANK')
                reg_ex(cwd,continue_script,'xxlenarrayxx','BLANK')
                jobhb_script = '%s-jobhb.sh' % str(i).zfill(2)
                cp_file(dir_temp,'jobhb.sh',cwd,\
                        jobhb_script)
                reg_ex(cwd,jobhb_script,'xxjobnamexx','BLANK')
                reg_ex(cwd,jobhb_script,'xxnumxx','BLANK')
                job_script = '%s-job.sh' % str(i).zfill(2)
                cp_file(dir_temp,'job.sh',cwd,\
                        job_script)
                reg_ex(cwd,job_script,'xxjobnamexx','BLANK')
                reg_ex(cwd,job_script,'xxwalltimexx','BLANK')
                reg_ex(cwd,job_script,'xxnodesxx','BLANK')
                stage_dir = os.path.join(cwd,str(i).zfill(2))
                if not os.path.exists(stage_dir): os.makedirs(stage_dir)
                stage_level(stage_dir,i,seeds[i-1])
        # call work_level, cascade into stage_level, traj_level
        work_level()

'''
            reg_ex(script,'xxcstagexx',stage)
            reg_ex(script,'xxtot_stagesxx',str(self.st))
            reg_ex(script,'xxquotaxx',str(self.dps))
            #reg_ex(script,'xxquotaxx',str(self.hm))
            #reg_ex(script,'xxhowmanyxx',str(self.dct['howmany']))
            reg_ex(script,'xxhowmanyxx',str(self.tpd))
            reg_ex(script,'xxenvironxx',self.e)
            reg_ex(script,'xxvelxx',str(self.pv_ans[0]))
            lenarray=self.path_steps[phase]/self.freq+1
            reg_ex(script,'xxlenarrayxx',str(int(lenarray)))
            #len_hb_pkl=500 # length of the hydrogen bond pkl
            #print str(int(self.ps[phase])/len_hb_pkl)  #xxdcdxx
            reg_ex(script,'xxlenbpklxx',str(self.hb_l))
            reg_ex(script,'xxdtxx',str(self.dt))
            reg_ex(script,'xxmoleculexx',self.mol.upper())
            reg_ex(script,'xxngnxx',self.ngn.upper())
            reg_ex(script,'xxstartconstraintxx',str(self.spos))
            reg_ex(script,'xxsposxx',str(self.spos))
            reg_ex(script,'xxtempxx',str(self.temp))
            plotname=self.mol+self.ngn+self.e \
                    +'v'+confige[self.vel]+'asmd'
            reg_ex(script,'xxplotnamexx',plotname)
            reg_ex(script,'xxnumxx',stage)
            reg_ex(script,'xxjobnamexx',self.mol+self.ngn[0]+'cjob'+stage)
        for i in range(1,int(self.st)+1):
            # write job.sh, jobh.sh, continue.py + stage
            os.makedirs(os.path.join(self.vdir,str(i).zfill(2)))
            #cp_file(os.path.join(self.ndir,'continue'),'continue.py', \
                            #self.vdir,str(i).zfill(2)+'-continue.py')

            #cp_file(os.path.join(self.ndir,'continue'),'perpetuate.py', \
            #                self.vdir,str(i).zfill(2)+'-continue.py')
            cp_file(os.path.join(self.ndir,'continue'),'pmf_no_corr.py', \
                            self.vdir,str(i).zfill(2)+'-continue.py')
            cp_file(os.path.join(self.ndir,'jobc'),'job-'+self.gate+'.sh',\
                            self.vdir,str(i).zfill(2)+'-job.sh')
            cp_file(os.path.join(self.ndir,'jobhb'),'job-'+self.gate+'.sh',\
                            self.vdir,str(i).zfill(2)+'-jobh.sh')
            stage=str(i).zfill(2)
            reg_exp_contd(os.path.join(self.vdir,str(i).zfill(2)+ \
                  '-continue.py'),stage,i)
            #reg_exp_contd(os.path.join(self.vdir,str(i).zfill(2)+ \
            #      '-perpetuate.py'),stage,i)

            reg_exp_contd(os.path.join(self.vdir,str(i).zfill(2)+ \
                  '-job.sh'),stage,i)
            reg_exp_contd(os.path.join(self.vdir,str(i).zfill(2)+ \
                  '-jobh.sh'),stage,i)
        cp_file(os.path.join(self.ndir,'hb_pkl'),'hb_pkl.py', \
                self.vdir,'00-hb_pkl.py')
        reg_exp_contd(os.path.join(self.vdir,'00-hb_pkl.py'),stage,i)
        #        '-hb_pkl.py'),stage,i)
        # ^ changing hb_pkl / continue
        #cp_file(self.pydir,'plotpkl.py',self.vdir,'plotpkl.py')
        #cp_file(self.pydir,'plotpkl2.py',self.vdir,'plotpkl.py')
        cp_file(self.pydir,'plot_pmf_no_corr.py',self.vdir,'plotpkl.py')
        cp_file(self.pydir,'plothb.py',self.vdir,'plothb.py')
        cp_file(self.pydir,'mpmf.py',self.pdir,'mpmf.py')
        cp_file(self.pydir,'pmf15.py',self.pdir,'pmf15.py')
        cp_file(self.pydir,'weighthb.py',self.pdir,'weighthb.py')
        # DESCRIPTION
        # regular expressions called
        # d_vis = i.e. .../01.vac/02/plotpkl.py
        # from above:   def reg_exp_contd(script,stage,i):
        d_vis=os.path.join(self.vdir,'plotpkl.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.vdir,'plothb.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.pdir,'mpmf.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.pdir,'pmf15.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.pdir,'weighthb.py')
        reg_exp_contd(d_vis,stage,1)

            def call_expavg(script):
                tefdir='0'+self.vel+'.*/*-tef.dat*'
                reg_ex(script,'xxtefdirxx',tefdir)
                reg_ex(script,'xxnumxx',self.v0)
                # gone - june 6: reg_ex(script,'xxpfxx',str(dictpf[self.vel]))
                plotname=self.mol+self.ngn+self.e+\
                         str(confige[self.vel])+str(self.spos)
                reg_ex(script,'xxplotnamexx',plotname)
                reg_ex(script,'xxsposxx',str(self.spos))
                reg_ex(script,'xxstartconstraintxx',str(self.spos))
                reg_ex(script,'xxmoleculexx',self.mol.upper())
                reg_ex(script,'xxenvironxx',self.e)
                reg_ex(script,'xxdtxx',str(self.dt))
            def call_smd(script):
                phase = int(script.split('/')[-3])-1
                reg_ex(script,'xxstepsxx',str(int(self.path_steps[phase])))
                # len_hb_pkl, self.hb_l
                reg_ex(script,'xxdcdxx',str(int(self.path_steps[phase]) \
                           /self.hb_l))
                if self.ngn=='amb':
                    reg_ex(script,'xxtsxx',str(self.ts/1000))
                elif self.ngn=='namd':
                    reg_ex(script,'xxtsxx',str(self.ts))
                reg_ex(script,'xxlDxx',self.lD)
                reg_ex(script,'xxtempxx',self.temp)
                reg_ex(script,'xxfreqxx',str(self.freq))
            for root, dirnames, filenames in os.walk(subdir):
                for f in filenames:
                    if len(f.split('-'))==1:
                        idn=f
                    elif len(f.split('-'))==2:
                        idn=f.split('-')[1]
                    script=os.path.join(root,f)
                    if idn=='job.sh':
                        bashjobname=self.mol+self.ngn[0]+\
                            str(int(self.path_vel[0]/2*(10**6)))+self.e[0]+ds
                        reg_ex(script,'xxjobnamexx',bashjobname)
                        reg_ex(script,'xxqueuexx',configq[self.q])
                        reg_ex(script,'xxnodesxx',confign[self.cn][self.comp])
                        reg_ex(script,'xxwalltimexx',configw[self.wt])
                    elif idn=='go.py':
                        #hw=str(self.dct['howmany'])
                        reg_ex(script,'xxhowmanyxx',str(self.tpd))
                        reg_ex(script,'xxnodecountxx',self.cn)
                        reg_ex(script,'xxstartconstraintxx',str(self.spos))
                        reg_ex(script,'xxquotaxx',str(self.dps))
                        reg_ex(script,'xxstrucequilxx',self.env)
                    elif idn=='smd.namd':
                        call_smd(script)
                    elif idn=='smd.in':
                        call_smd(script)
                    elif idn=='tm.tex':
                        call_expavg(script)
                    elif idn=='smdforce.tcl':       # NAMD STEERING FILE
                        phase = int(script.split('/')[-3])-1
                        reg_ex(script,'xxvelocityxx',str(self.path_vel[phase]))
                        reg_ex(script,'xxzcoordxx',str(self.spos))
                        #zdist = (self.path_vel*self.path_seg).cumsum()
                        # june 7
                        zdist = (self.path_vel*self.path_steps).cumsum()
                        if phase == 0:
                            reg_ex(script,'xxcur_zxx',str(0))
                        else:
                            reg_ex(script,'xxcur_zxx',str(zdist[phase-1]))
                        reg_ex(script,'xxtsxx',str(self.ts))
                        reg_ex(script,'xxfreqxx',str(self.freq))
                        reg_ex(script,'yyyyy',str(int(script.split('/')[-3])-1))
                    elif idn=='dist.RST':           # AMB STEERING FILE
                        phase = int(script.split('/')[-3])-1
                        #zdist = (self.path_vel*self.path_seg).cumsum()
                        zdist = (self.path_vel*self.path_steps).cumsum()
                        print zdist
                        if phase == 0:
                            reg_ex(script,'xxsposxx',str(self.spos))
                            reg_ex(script,'xxeposxx',str(self.spos+zdist[phase]))
                        else:
                            reg_ex(script,'xxsposxx',str(self.spos+zdist[phase-1]))
                            reg_ex(script,'xxeposxx',str(self.spos+zdist[phase]))
        def how_many(tmpdir):
            #rlist=[]
            #rlist=[str(x).zfill(3) for x in range(int(self.hm))]
            # change v (below) on june 7
            #for i in range(int(self.hm)):
            for i in range(int(self.dps)):


                def reg_seed(subdir,seed):
                    if self.ngn=='namd':
                        f1=os.path.join(subdir,'smd.namd')
                    elif self.ngn=='amb':
                        f1=os.path.join(subdir,'smd.in')
                r=random.randint(10000,99999)
                while r in rlist:
                    r=random.randint(10000,99999)
                rlist.append(r)


                tdi='/'.join(tmpdir.split('/')[0:-1])
                td=os.path.join(tdi,str(i).zfill(3))
                shutil.copytree(tmpdir,td)
                #reg_seed(td,r)
            shutil.rmtree(tmpdir)
        dls=[]
        [dls.append(d) for d in os.listdir(self.vdir) if \
             os.path.isdir(os.path.join(self.vdir,d))]
        gs = gen_all_seeds()
        for ds in dls:
            ddir=os.path.join(self.vdir,ds,'tmp')
            os.makedirs(ddir)
            cp_file(os.path.join(self.ndir,'job'),'job-'+self.gate+'.sh',\
                   ddir,'job.sh')
            cp_file(os.path.join(self.ndir,'go'),'go-'+self.gate+ \
                   '.py',ddir,'go.py')
            if ds=='01':
                if self.ngn=='namd':    # NAMD
                    cp_file(os.path.join(self.ndir,'mol.conf.tcl',self.mol, \
                            self.env),'smd_initial.namd',ddir,'smd.namd')
                    cp_file(os.path.join(self.ndir,'mol.conf.tcl',self.mol, \
                            self.env),'smd_force.tcl',ddir,'smdforce.tcl')
                elif self.ngn=='amb':   # AMB
                    cp_file(os.path.join(self.ndir,'mol.conf',self.mol,self.env),\
                            'smd.in',ddir,'smd.in')
                    cp_file(os.path.join(self.ndir,'mol.conf',self.mol,self.env),\
                            'dist.RST',ddir,'dist.RST')
                elif self.ngn=='gro':   # GRO
                    cp_file(os.path.join(self.ndir,'mol.conf',self.mol,self.env),\
                            'smd.in',ddir,'smd.in')
                    cp_file(os.path.join(self.ndir,'mol.conf',self.mol,self.env),\
                            'dist.RST',ddir,'dist.RST')
            elif ds!='01':
                if self.ngn=='namd':    # NAMD
                    cp_file(os.path.join(self.ndir,'mol.conf.tcl',self.mol, \
                            self.env),'smd_continue.namd',ddir,'smd.namd')
                    cp_file(os.path.join(self.ndir,'mol.conf.tcl',self.mol, \
                            self.env),'smd_force.tcl',ddir,'smdforce.tcl')
                elif self.ngn=='amb':   # AMB - significantly developed, not ready
                    cp_file(os.path.join(self.ndir,'mol.conf',self.mol,self.env),\
                            'smd_r.in',ddir,'smd.in')
                    cp_file(os.path.join(self.ndir,'mol.conf',self.mol,self.env),\
                            'dist.RST',ddir,'dist.RST')
                elif self.ngn=='gro':   # GRO - minimally developed
                    cp_file(os.path.join(self.ndir,'mol.conf',self.mol,self.env),\
                            'smd_r.in',ddir,'smd.in')
                    cp_file(os.path.join(self.ndir,'mol.conf',self.mol,self.env),\
                            'dist.RST',ddir,'dist.RST')
            reg_exp(ddir,ds,gs[int(ds)-1])
            how_many(ddir)
    def steering_control(self):
        def reg_exp(dir_loc,idn):
            script=os.path.join(dir_loc,idn)
            if idn=='smdforce.tcl':
                reg_ex(script,'xxzcoordxx',str(self.spos))
                reg_ex(script,'xxtsxx',str(self.ts))
                reg_ex(script,'xxfreqxx',str(setup[self.vel]['freq']))
                reg_ex(script,'yyyyy',str(int(dir_loc.split('/')[-1])-1))
            elif idn=='dist.RST':
                scon=str(self.zc+(int(dir_loc.split('/')[-1])-1)*2)
                econ=str(self.zc+int(dir_loc.split('/')[-1])*2)
                reg_ex(script,'xxstartconstraintxx',scon)
                reg_ex(script,'xxendconstraintxx',econ)
            elif idn=='hb.py':
                reg_ex(script,'xxenvironxx',self.env)
            #____END OF FUNCTION_____________
        dls=[]
        [dls.append(d) for d in os.listdir(self.vdir) if \
             os.path.isdir(os.path.join(self.vdir,d))]
        for ds in dls:
            ddir=os.path.join(self.vdir,ds)
            cp_file(os.path.join(self.ndir,'hb'),'hb.py',ddir,'hb.py')
            reg_exp(ddir,'hb.py')
'''
