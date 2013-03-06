import shutil,os,itertools,re,random,fnmatch,time,pickle,sys
#_______DICTIONARY_____________________________________________________________
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
strdir ={'01.vac':'00.struc-equil.v','02.imp':'00.struc-equil.i',
         '03.exp':'00.struc-equil.e'}
dictpf ={'1':1,'2':1,'3':50,'4':100,'5':500}
setup0 ={'1':{'vel':0.002,'steps':1000,'dcd':100,'howmany':50,'freq':50},
       '2':{'vel':0.0002,'steps':1000,'dcd':1000,'howmany':25,'freq':50},
       '3':{'vel':0.00002,'steps':1000,'dcd':10000,'howmany':25,'freq':50},
       '4':{'vel':0.000002,'steps':1000,'dcd':100000,'howmany':5,'freq':50},
       '5':{'vel':0.0000002,'steps':1000,'dcd':1000000,'howmany':5,'freq':50}}
#__________global_____use______________________________________________________
class mdict(dict):
    def __setitem__(self,key,value):
        self.setdefault(key,[]).append(value)
def print_dict(dct):
    for key,value in dct.items():
        print key,value
        print ''
def cp_file(f_dir,f,d_dir,d):
    shutil.copy(os.path.join(f_dir,f),os.path.join(d_dir,d))
def cp_tree(f_dir,f,d_dir,d):
    shutil.copytree(os.path.join(f_dir,f),os.path.join(d_dir,d))
def reg_ex(script,subout,subin):
    o=open(script,'r+')
    text=o.read()
    text=re.sub(subout,subin,text)
    o.close()
    o=open(script,'w+')
    o.write(text)
    o.close()
#______________________________________________________________________________
class a_make_JobDirSmd:
    def __init__(self,ngn,mol,zc,workdir,jobdir):
        self.ngn  = ngn
        self.mol  = mol
        self.zc   = zc
        self.workdir = workdir
        self.jobdir = jobdir
        self.jdir = os.path.join(self.workdir,self.jobdir)
    def a_makeJobDir(self):
        if not os.path.exists(self.jdir):os.makedirs(self.jdir)
        texdir=os.path.join(self.jobdir,'tex_%s' % self.jdir.split('/')[-1])
        if not os.path.exists(texdir): os.makedirs(texdir)
        cp_file(self.workdir,'gen.py',self.jobdir,'.gen_%s.py' % \
                  self.jdir.split('/')[-1])
        cp_file(self.workdir,'04.scripts/pipe.py',self.jobdir,'pipe.py')
        cp_file(self.workdir,'04.scripts/tm.tex',texdir,'tm.tex')
        cp_file(self.workdir,'04.scripts/pdflatex.sh',texdir,'pdflatex.sh')
        cp_file(self.workdir,'04.scripts/del.py',self.jobdir,'del.py')
        if self.ngn == 'namd':
            cp_tree(self.workdir,'04.toppar',self.jobdir,'toppar')
        cp_tree(self.workdir,'04.scripts',self.jobdir,'04.scripts')
        return texdir
    def reg_exp(self,subdir):
        for root, dirnames, filenames in os.walk(subdir):
            for f in filenames:
                script=os.path.join(root,f)
                if f=='tm.tex':
                    plotname=self.mol+self.ngn+'env'+confige['2']+'asmd'
                    reg_ex(script,'xxplotnamexx',plotname)
#______________________________________________________________________________
class a_Struc_Dirs:
    def __init__(self,ngn,mol,env,workdir,jobdir):
        self.ngn  = ngn
        self.mol  = mol
        self.env  = env
        self.workdir = workdir
        self.jobdir = jobdir
        self.jdir = os.path.join(self.workdir,self.jobdir)
    def a_makeStrucDir(self):
        cp_tree(os.path.join(self.workdir,'01.struc-equil.'+self.ngn, \
                self.mol,self.env),'',self.jdir,strdir[self.env])

class a_Smd_Method:
    def __init__(self,ngn,mol,env,v,ts,zc,lD,sf,workdir,jobdir, \
          gate,cn,comp,wallt,queue,howmany,stages,direct,dist,config):
        self.ngn  = ngn
        self.mol  = mol
        self.env  = env                     # 01.vac
        self.e    = self.env.split('.')[1]  #    vac
        self.v0   = str(int(v)).zfill(2)    # 02
        self.vel  = str(int(v))             #  2
        self.v    = v                       #  2.
        self.ts   = ts
        self.spos = zc
        self.lD   = lD
        self.sf   = sf
        self.workdir = workdir
        self.jobdir = jobdir
        self.gate = gate
        self.cn   = cn
        self.comp = comp
        self.wt   = wallt
        self.q    = queue
        self.jdir = os.path.join(self.workdir,self.jobdir)
        self.edir = os.path.join(self.jdir,self.env)
        self.vdir = os.path.join(self.edir,self.v0)
        self.hm   = howmany
        self.st   = stages
        self.cfg  = config   # << print_dict(self.cfg)
        print_dict(self.cfg)
        #___CONFIG_SECTION_BEGIN
        '''
        print self.cfg[v][0][0],'dist'
        print self.cfg[v][0][1],'ts'
        print self.cfg[v][0][2],'path_seg'
        print self.cfg[v][0][3],'path_svel'
        print self.cfg[v][0][4],'path_vel'
        print self.cfg[v][0][5],'path_steps'
        print self.cfg[v][0][6],'setup'
        '''
        self.d  = self.cfg[v][0][0]    # total distance
        self.ts = self.cfg[v][0][1]    # timestep
        self.pv = self.cfg[v][0][4]    # path_vel
        self.ps = self.cfg[v][0][5]    # path_steps
        self.dct= self.cfg[v][0][6]    # dict: 'howmany','freq'
        self.dt   = self.dct['freq']*ts/1000
	self.pv_aps=self.pv/ts*1000
	self.pv_ans=self.pv/ts*(10**6)
        zdist_c = format((self.pv*self.ps).cumsum()[-1],'.2f')
        dist_c  = format(self.d,'.2f')
        print 'zdist',zdist_c,'dist',dist_c
        if zdist_c!=dist_c:
            print '?@#$!_correct total distance!'
            sys.exit()
        #___CONFIG_SECTION_END
    def a_makeEnvDir(self):
        if not os.path.exists(self.edir):os.makedirs(self.edir)
    def a_makeVelDir(self):
        os.makedirs(self.vdir)
    def a_makeContainDir(self):
        def reg_exp_contd(script,stage,i):
            phase = i-1
            reg_ex(script,'xxcstagexx',stage)
            reg_ex(script,'xxquotaxx',str(self.hm))
            reg_ex(script,'xxhowmanyxx',str(self.dct['howmany']))
            reg_ex(script,'xxenvironxx',self.e)
            reg_ex(script,'xxvelxx',str(self.pv_ans[0]))
            lenarray=self.ps[phase]/self.dct['freq']+1
            reg_ex(script,'xxlenarrayxx',str(lenarray))
            reg_ex(script,'xxdtxx',str(self.dt))
            reg_ex(script,'xxmoleculexx',self.mol.upper())
            reg_ex(script,'xxngnxx',self.ngn.upper())
            reg_ex(script,'xxstartconstraintxx',str(self.spos))
            reg_ex(script,'xxsposxx',str(self.spos))
            plotname=self.mol+self.ngn+self.e \
                    +'v'+confige[self.vel]+'asmd'
            reg_ex(script,'xxplotnamexx',plotname)
        for i in range(1,int(self.st)+1):
            os.makedirs(os.path.join(self.vdir,str(i).zfill(2)))
            cp_file(os.path.join(self.workdir,'00.maindir.%s.asmd' % self.ngn \
                    ),'continue.py',self.vdir,str(i).zfill(2)+'-continue.py')
            stage=str(i).zfill(2)
            reg_exp_contd(os.path.join(self.vdir,str(i).zfill(2)+ \
                        '-continue.py'),stage,i)
        cp_file(os.path.join(self.workdir,'04.scripts'),'env_allhb.py', \
                self.vdir,'env_allhb.py')
        cp_file(os.path.join(self.workdir,'04.scripts'),'env_allwp.py', \
                self.vdir,'env_allwp.py')
        cp_file(os.path.join(self.workdir,'04.scripts'),'env_ihbond.py', \
                self.vdir,'env_ihbond.py')
        cp_file(os.path.join(self.workdir,'04.scripts'),'discrete.py', \
                self.vdir,'discrete.py')
        cp_file(os.path.join(self.workdir,'04.scripts'),'plotpkl.py', \
                self.vdir,'plotpkl.py')
        d_vis=os.path.join(self.vdir,'env_allhb.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.vdir,'env_allwp.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.vdir,'env_ihbond.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.vdir,'discrete.py')
        reg_exp_contd(d_vis,stage,1)
        d_vis=os.path.join(self.vdir,'plotpkl.py')
        reg_exp_contd(d_vis,stage,1)
    def a_savepickle(self):
        os.chdir(self.vdir)
        pickle.dump(self.cfg,open('config.pkl','w'))
#_____________________________________________________________________________
    def a_makeSubDir(self):
        def reg_exp(subdir,ds):
            def call_expavg(script):
                tefdir='0'+self.vel+'.*/*-tef.dat*'
                reg_ex(script,'xxtefdirxx',tefdir)
                reg_ex(script,'xxnumxx',self.v0)
                reg_ex(script,'xxpfxx',str(dictpf[self.vel]))
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
                reg_ex(script,'xxstepsxx',str(int(self.ps[phase])))
                reg_ex(script,'xxdcdxx',str(int(self.ps[phase])/100))
                if self.ngn=='amb':
                    reg_ex(script,'xxtsxx',str(self.ts/1000))
                elif self.ngn=='namd':
                    reg_ex(script,'xxtsxx',str(self.ts))
                reg_ex(script,'xxlDxx',self.lD)
                reg_ex(script,'xxfreqxx',str(self.dct['freq']))
            for root, dirnames, filenames in os.walk(subdir):
                for f in filenames:
                    if len(f.split('-'))==1:
                        idn=f
                    elif len(f.split('-'))==2:
                        idn=f.split('-')[1]
                    script=os.path.join(root,f)
                    if idn=='job.sh':
                        bashjobname=self.mol+self.ngn+\
                            str(self.pv[0]/2*(10**6))+self.e+ds
                        reg_ex(script,'xxjobnamexx',bashjobname)
                        reg_ex(script,'xxqueuexx',configq[self.q])
                        reg_ex(script,'xxnodesxx',confign[self.cn][self.comp])
                        reg_ex(script,'xxwalltimexx',configw[self.wt])
                    elif idn=='go.py':
                        hw=str(self.dct['howmany'])
                        reg_ex(script,'xxhowmanyxx',hw)
                        reg_ex(script,'xxnodecountxx',self.cn)
                        reg_ex(script,'xxstartconstraintxx',str(self.spos))
                        reg_ex(script,'xxquotaxx',self.hm)
                        reg_ex(script,'xxstrucequilxx',strdir[self.env])
                    elif idn=='smd.namd':
                        call_smd(script)
                    elif idn=='smd.in':
                        call_smd(script)
                    elif idn =='expavg.py':
                        call_expavg(script)
                    elif idn=='tm.tex':
                        call_expavg(script)
                    elif idn=='dualplot.py':
                        call_expavg(script)
                    elif idn=='npy.py':
                        call_expavg(script)
                    elif idn=='allhb.py':
                        call_expavg(script)
                    elif idn=='allwp.py':
                        call_expavg(script)
                    elif idn=='ihbond.py':
                        call_expavg(script)
                    elif idn=='smdforce.tcl':
                        phase = int(script.split('/')[-3])-1
                        reg_ex(script,'xxvelocityxx',str(self.pv[phase]))
                        reg_ex(script,'xxzcoordxx',str(self.spos))
                        zdist = (self.pv*self.ps).cumsum()
                        if phase == 0:
                            reg_ex(script,'xxcur_zxx',str(0))
                        else:
                            reg_ex(script,'xxcur_zxx',str(zdist[phase-1]))
                        reg_ex(script,'xxtsxx',str(self.ts))
                        reg_ex(script,'xxfreqxx',str(self.dct['freq']))
                        reg_ex(script,'yyyyy',str(int(script.split('/')[-3])-1))
                    elif idn=='dist.RST':
                        phase = int(script.split('/')[-3])-1
                        zdist = (self.pv*self.ps).cumsum()
                        #print zdist
                        #time.sleep(1)
                        if phase == 0:
                            reg_ex(script,'xxsposxx',str(self.spos))
                            reg_ex(script,'xxeposxx',str(self.spos+zdist[phase]))
                        else:
                            reg_ex(script,'xxsposxx',str(self.spos+zdist[phase-1]))
                            reg_ex(script,'xxeposxx',str(self.spos+zdist[phase]))
                        # need mod
                        #scon=str(self.zc+(int(dir_loc.split('/')[-1])-1)*2)
                        #econ=str(self.zc+int(dir_loc.split('/')[-1])*2)
                        #reg_ex(script,'xxstartconstraintxx','scon')
                        #reg_ex(script,'xxendconstraintxx','econ')
        def how_many(tmpdir):
            rlist=[]
            for i in range(int(self.hm)):
                def reg_seed(subdir,seed):
                    if self.ngn=='namd':
                        f1=os.path.join(subdir,'smd.namd')
                        # ### reg_ex(f1,'xxxxx',str(seed))
                    elif self.ngn=='amb':
                        f1=os.path.join(subdir,'smd.in')
                        # ### reg_ex(f1,'xxxxx',str(seed))
                r=random.randint(10000,99999)
                while r in rlist:
                    r=random.randint(10000,99999)
                rlist.append(r)
                tdi='/'.join(tmpdir.split('/')[0:-1])
                td=os.path.join(tdi,str(r))
                shutil.copytree(tmpdir,td)
                reg_seed(td,r)
            shutil.rmtree(tmpdir)
        mdir='00.maindir.'+self.ngn+'.asmd'
        pre=os.path.join(self.workdir,mdir)
        sdir=os.path.join(pre,self.mol,self.env)
        dls=[]
        [dls.append(d) for d in os.listdir(self.vdir) if \
             os.path.isdir(os.path.join(self.vdir,d))]
        for ds in dls:
            ddir=os.path.join(self.vdir,ds,'tmp')
            os.makedirs(ddir)
            # switching tcl and dist back into subdir
            if ds=='01':
                if self.ngn=='namd':
                    cp_file(sdir,'smd_initial.namd',ddir,'smd.namd')
                    cp_file(sdir,'smd_force.tcl',ddir,'smdforce.tcl')
                    cp_file(os.path.join(pre,'VDIR'),'job-'+self.gate+ \
                            '.sh',ddir,'job.sh')
                    cp_file(os.path.join(pre,'VDIR'),'go-'+self.gate+ \
                            '.py',ddir,'go.py')
                if self.ngn=='amb':
                    cp_file(sdir,'smd.in',ddir,'smd.in')
                    # back into subdir
                    cp_file(sdir,'dist.RST',ddir,'dist.RST')
                    cp_file(os.path.join(pre,'VDIR'),'job-'+self.gate+'.sh',\
                            ddir,'job.sh')
                    cp_file(os.path.join(pre,'VDIR'),'go-'+self.gate+'.py',\
                            ddir,'go.py')
            elif ds!='01':
                if self.ngn=='namd':
                    cp_file(sdir,'smd_continue.namd',ddir,'smd.namd')
                    cp_file(sdir,'smd_force.tcl',ddir,'smdforce.tcl')
                    cp_file(os.path.join(pre,'VDIR'),'job-'+self.gate+ \
                            '.sh',ddir,'job.sh')
                    cp_file(os.path.join(pre,'VDIR'),'go-'+self.gate+ \
                            '.py',ddir,'go.py')
                if self.ngn=='amb':
                    cp_file(sdir,'smd_r.in',ddir,'smd.in')
                    # back into subdir
                    cp_file(sdir,'dist.RST',ddir,'dist.RST')
                    cp_file(os.path.join(pre,'VDIR'),'job-'+self.gate+'.sh',\
                            ddir,'job.sh')
                    cp_file(os.path.join(pre,'VDIR'),'go-'+self.gate+ \
                            '_r.py',ddir,'go.py')
            reg_exp(ddir,ds)
            how_many(ddir)
    def a_steering_control(self):
        def reg_exp(dir_loc,idn):
            script=os.path.join(dir_loc,idn)
            if idn=='smdforce.tcl':
                # removed self.sf
                #reg_ex(script,'xxvelocityxx',str(self.v_ats))
                #reg_ex(script,'xxvelocityxx',self.pv[
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
                reg_ex(script,'xxenvironxx',strdir[self.env])
            #____END OF FUNCTION_____________
        mdir='00.maindir.'+self.ngn+'.asmd'
        pre=os.path.join(self.workdir,mdir)
        sdir=os.path.join(pre,self.mol,self.env)
        svdir=os.path.join(pre,'VDIR')
        dls=[]
        [dls.append(d) for d in os.listdir(self.vdir) if \
             os.path.isdir(os.path.join(self.vdir,d))]
        if self.ngn=='namd':
            for ds in dls:
                ddir=os.path.join(self.vdir,ds)
                #cp_file(sdir,'smd_force.tcl',ddir,'smdforce.tcl')
                #reg_exp(ddir,'smdforce.tcl')
                cp_file(svdir,'hb.py',ddir,'hb.py')
                reg_exp(ddir,'hb.py')
                # dist and tcl back into subdir
        if self.ngn=='amb':
            for ds in dls:
                ddir=os.path.join(self.vdir,ds)
                #cp_file(sdir,'dist.RST',ddir,'dist.RST')
                #reg_exp(ddir,'dist.RST')
                cp_file(svdir,'hb.py',ddir,'hb.py')
                reg_exp(ddir,'hb.py')
