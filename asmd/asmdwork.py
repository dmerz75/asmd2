import shutil,os,itertools,re,random,fnmatch,time,pickle,sys
import numpy as np
#_______DICTIONARY_____________________________________________________________
'''
    confign: regular expressions for substituting the xxnodesxx from job.sh
    configw: regular expressions for substituting the xxwalltimexx from job.sh
    configq: regular expressions for substituting the xxqueuexx from job.sh
    confige: almost unnecessary
'''
confign={'1':{'gpu':'nodes=1:ppn=1:gpus=1:TESLA','cpu':'nodes=1:ppn=1'},
         '2':{'gpu':'nodes=1:ppn=2:gpus=1:TESLA','cpu':'nodes=1:ppn=2'},
         '3':{'gpu':'nodes=1:ppn=3:gpus=1:TESLA','cpu':'nodes=1:ppn=3'},
         '4':{'gpu':'nodes=1:ppn=4:gpus=1:TESLA','cpu':'nodes=1:ppn=4'},
         '5':{'gpu':'nodes=1:ppn=5:gpus=1:TESLA','cpu':'nodes=1:ppn=5'},
         '6':{'gpu':'nodes=1:ppn=6:gpus=1:TESLA','cpu':'nodes=1:ppn=6'},
         '7':{'gpu':'nodes=1:ppn=7:gpus=1:TESLA','cpu':'nodes=1:ppn=7'},
         '8':{'gpu':'nodes=1:ppn=8:gpus=1:TESLA','cpu':'nodes=1:ppn=8'},
       '16':{'gpu':'nodes=1:ppn=16:gpus=1:TESLA','cpu':'nodes=1:ppn=16'}}
configw={'sst':'walltime=15:00','swt':'walltime=00:90:00','mwt': \
         'walltime=03:59:00',
        'lwt':'walltime=072:00:00','dwt':'walltime=15:00:00:00'}
configq={'short':'tg_short','workq':'tg_workq','standby':'standby',
         'standby-8':'standby-8','debug':'tg_debug'}
confige={'1':'v1000','2':'v100','3':'v10','4':'v1','5':'vp1'}
# dictpf ={'1':1,'2':1,'3':50,'4':100,'5':500}

#__________global_____use______________________________________________________
'''
class mdict(dict):
    def __setitem__(self,key,value):
        self.setdefault(key,[]).append(value)

def print_dict(dct):
    for key,value in dct.items():
        print key,value
        print ''
'''

def cp_file(f_dir,f,d_dir,d):
    ''' used extensively for copying a single template file into a temporary dir
    '''
    shutil.copy(os.path.join(f_dir,f),os.path.join(d_dir,d))

def cp_tree(f_dir,f,d_dir,d):
    ''' used extensively for copying a template dir into a created directory
    '''
    shutil.copytree(os.path.join(f_dir,f),os.path.join(d_dir,d))

def reg_ex(script,subout,subin):
    ''' performs all regular expressions when a text file is opened
    '''
    o=open(script,'r+')
    text=o.read()
    text=re.sub(subout,subin,text)
    o.close()
    o=open(script,'w+')
    o.write(text)
    o.close()
    
#__class_a_make_JobDirSmd______________________________________________________
class est_JobDir:
    def __init__(self,ngn,mol,zc,workdir,jobdir,pack_dir):
        self.ngn  = ngn
        self.mol  = mol
        self.zc   = zc
        self.workdir = workdir
        # ^^^ from gen.py: workdir=os.path.abspath(os.path.dirname(__file__))
        self.mdir    = os.path.join(workdir,'00.maindir')
        self.ndir    = os.path.join(workdir,'00.maindir',ngn)
        self.pydir= os.path.join(workdir,'00.maindir','py_gen')
        self.packdir = pack_dir
        # ^^^ from gen.py: pack_dir=ngn[0][0]+molec+'_'+jobid
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
    def reg_exp(self,subdir):
        for root, dirnames, filenames in os.walk(subdir):
            for f in filenames:
                script=os.path.join(root,f)
                if f=='tm.tex':
                    plotname=self.mol+self.ngn+'env'+confige['2']+'asmd'
                    reg_ex(script,'xxplotnamexx',plotname)
#__class_a_Struc_Dirs__________________________________________________________
class est_StrucDir:
    def __init__(self,ngn,mol,env,workdir,jobdir,pack_dir):
        self.ngn  = ngn
        self.mol  = mol
        self.env  = env
        self.workdir = workdir
        self.mdir = os.path.join(workdir,'00.maindir')
        self.ndir = os.path.join(workdir,'00.maindir',ngn)
        self.pydir= os.path.join(workdir,'00.maindir','py_gen')
        self.packdir = pack_dir
        self.pdir = os.path.join(self.workdir,self.packdir)
        self.jobdir = jobdir
        self.jdir = os.path.join(self.pdir,self.jobdir)
    def makeStrucDir(self):
        if not os.path.exists(os.path.join(self.jdir,'00.struc')):
            os.makedirs(os.path.join(self.jdir,'00.struc'))
        if not os.path.exists(os.path.join(self.jdir,'00.struc',self.env)):
            cp_tree(os.path.join(self.ndir,'struc',self.mol,self.env),'',\
                os.path.join(self.jdir,'00.struc'),self.env)
#__class_AsmdMethod____________________________________________________________
class AsmdMethod:
    def __init__(self,ngn,mol,env,v,ts,zc,lD,workdir,jobdir,pack_dir,\
          gate,cn,comp,wallt,queue,stages,direct,\
          dist,path_seg,path_svel,path_vel,path_steps,freq,dps,tpd,tps,temp):
        self.ngn  = ngn
        self.mol  = mol
        self.env  = env                     # 01.vac
        self.e    = self.env.split('.')[1]  #    vac
        self.v0   = str(int(v)).zfill(2)    # 02
        self.vel  = str(int(v))             #  2
        self.v    = v                       #  2.
        self.ts   = ts
        self.spos = zc                      # zc ~ zcrd
        self.lD   = lD
        self.temp = temp
        self.workdir = workdir
        self.mdir = os.path.join(workdir,'00.maindir')   # 00.maindir
        self.ndir = os.path.join(workdir,'00.maindir',ngn) # 00.maindir/namd
        self.pydir= os.path.join(workdir,'00.maindir','py_gen')
        self.packdir = pack_dir                          # ^ 00.maindir/py_gen
        self.pdir = os.path.join(self.workdir,self.packdir)
        self.jobdir = jobdir
        self.jdir = os.path.join(self.pdir,self.jobdir)
        self.gate = gate
        self.cn   = cn
        self.comp = comp
        self.wt   = wallt
        self.q    = queue
        self.edir = os.path.join(self.jdir,self.env)  # i.e. 02.imp
        self.vdir = os.path.join(self.edir,self.v0)   # i.e. 02,03,...
        self.tpd  = tpd  # T   D   T
        self.dps  = dps  # / = / * / ; traj/stage = dir/stage * traj/dir
        self.tps  = tps  # S   S   D
        #self.hm   = howmany
        #self.hmtpd= tpd
        self.st   = stages
        # replacement
        self.dist = dist
        self.path_seg = path_seg
        self.path_svel= path_svel
        self.path_vel = path_vel
        self.path_steps=path_steps
        self.freq = freq
        self.dt   = freq*ts/1000
        self.pv_aps=self.path_vel/ts*1000
        self.pv_ans=self.path_vel/ts*(10**6)
        '''
        self.cfg  = config   # << print_dict(self.cfg)
        #print_dict(self.cfg)
        # DESCRIPTION: writing the config.pkl
        #___CONFIG_SECTION_BEGIN
        #self.d  = self.cfg[v][0][0]    # total distance
        #self.ts = self.cfg[v][0][1]    # timestep
        #path.seg=self.cfg[v][0][2]    # path_seg
        #path.sv =self.cfg[v][0][3]    # path_svel
        self.pv = self.cfg[v][0][4]    # path_vel
        self.ps = self.cfg[v][0][5]    # path_steps
        self.dct= self.cfg[v][0][6]    # dict: 'howmany','freq'
        self.dt   = self.dct['freq']*ts/1000
        self.pv_aps=self.pv/ts*1000
        self.pv_ans=self.pv/ts*(10**6)
        '''
        # back to normal
        zdist_c = format((self.path_vel*self.path_steps).cumsum()[-1],'.2f')
        dist_c  = format(self.dist,'.2f')
        if zdist_c!=dist_c:
            print '?@#$! total distance doesn\'t match path_seg dist'
        #len_hb_pkl=500 # length of the hydrogen bond pkl
        self.hb_l = 100 # length of the hydrogen bond pkl
        #___CONFIG_SECTION_END
    def makeEnvDir(self):            # i.e.   02.imp or 03.exp
        if not os.path.exists(self.edir):os.makedirs(self.edir)
    def makeVelDir(self):            # i.e.   01, 02, 03, 04, 05
        os.makedirs(self.vdir)
    def makeContainDir(self):
        def reg_exp_contd(script,stage,i):
            phase = i-1
            reg_ex(script,'xxcstagexx',stage)
            reg_ex(script,'xxtot_stagesxx',str(self.st))
            # before tps = tpd * dps
            #print self.dps,self.tpd
            # june 7
            reg_ex(script,'xxquotaxx',str(self.dps))
            #reg_ex(script,'xxquotaxx',str(self.hm))
            #reg_ex(script,'xxhowmanyxx',str(self.dct['howmany']))
            reg_ex(script,'xxhowmanyxx',str(self.tpd))
            reg_ex(script,'xxenvironxx',self.e)
            reg_ex(script,'xxvelxx',str(self.pv_ans[0]))
            lenarray=self.path_steps[phase]/self.freq+1
            #print 'lenarray',lenarray
            #print 'self.ps[phase]',self.ps[phase]
            #print 'self.dct[freq]',self.dct['freq']+1
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
            ''' switching from work-corrected to pmf-corrected
            '''
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
            ''' because perpetuate.py is now 01-continue.py
            '''
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
        #d_vis=os.path.join(self.vdir,'plotpkl2.py')
        #reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.vdir,'plothb.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.pdir,'mpmf.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.pdir,'pmf15.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.pdir,'weighthb.py')
        reg_exp_contd(d_vis,stage,1)
        '''  june 7
    def savePickle(self):      # dumps config.pkl in
        os.chdir(self.jdir)    # pdir = nda_jobid
        pickle.dump(self.cfg,open('config%s_%s.pkl' %(self.vel,self.st),'w'))
        os.chdir(self.vdir)
        pickle.dump(self.cfg,open('config.pkl','w'))
        '''
#_____________________________________________________________________________
    def makeSubDir(self):
        def gen_all_seeds():
            os.chdir(self.vdir)
            sds = np.random.randint(10000,high=99999, \
                              size=(self.st,int(self.dps),self.tpd))
            fname = 'seeds.txt'
            with file('seeds.txt','w') as outfile:
                outfile.write('# %s stages' % sds.shape[0])
                outfile.write(" %s directories per stage" % sds.shape[1])
                outfile.write(" %s trajectories per directory\n" % sds.shape[2])
                outfile.write('# {0}\n'.format(sds.shape))
                for stg_slice in sds:
                    ss = np.transpose(stg_slice)
                    outfile.write('# stage\n')
                    np.savetxt(outfile,ss,fmt='%5.0d')
            return sds
        def reg_exp(subdir,ds,seed_list):
            os.chdir(self.vdir)
            with file('%s/%s.txt' %(ds,ds) ,'w') as outfile:
                if seed_list.shape[1]==1:
                    np.savetxt(outfile,seed_list,fmt='%5.0d')
                else:
                    np.savetxt(outfile,np.transpose(seed_list),fmt='%5.0d')
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
                '''
                def reg_seed(subdir,seed):
                    if self.ngn=='namd':
                        f1=os.path.join(subdir,'smd.namd')
                    elif self.ngn=='amb':
                        f1=os.path.join(subdir,'smd.in')
                r=random.randint(10000,99999)
                while r in rlist:
                    r=random.randint(10000,99999)
                rlist.append(r)
                '''

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
