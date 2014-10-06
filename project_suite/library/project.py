
import os

from generic import obj
from project_base import Pobj,modes

from project_manager import ProjectManager
from machines import Job,Machine,Supercomputer,get_machine

from structure import Structure,generate_structure,generate_cell
from physical_system import PhysicalSystem,generate_physical_system
from pseudopotential import Pseudopotential,Pseudopotentials
from bundle import bundle
from opium       import Opium      , OpiumInput      , OpiumAnalyzer
from sqd         import Sqd        , SqdInput        , SqdAnalyzer        , generate_sqd_input        , generate_sqd, hunds_rule_filling
from pwscf       import Pwscf      , PwscfInput      , PwscfAnalyzer      , generate_pwscf_input      , generate_pwscf
from pw2casino   import Pw2casino  , Pw2casinoInput  , Pw2casinoAnalyzer  , generate_pw2casino_input  , generate_pw2casino
from pw2qmcpack  import Pw2qmcpack , Pw2qmcpackInput , Pw2qmcpackAnalyzer , generate_pw2qmcpack_input , generate_pw2qmcpack
from wfconvert   import Wfconvert  , WfconvertInput  , WfconvertAnalyzer  , generate_wfconvert_input  , generate_wfconvert
from gamess      import Gamess     , GamessInput     , GamessAnalyzer     , generate_gamess_input     , generate_gamess, FormattedGroup
from convert4qmc import Convert4qmc, Convert4qmcInput, Convert4qmcAnalyzer, generate_convert4qmc_input, generate_convert4qmc
from qmcpack     import Qmcpack    , QmcpackInput    , QmcpackAnalyzer    , generate_qmcpack_input    , generate_qmcpack
from vasp        import Vasp       , VaspInput       , VaspAnalyzer       , generate_vasp_input       , generate_vasp
from qmcpack import loop,linear,cslinear,vmc,dmc
from qmcpack import generate_jastrows,generate_jastrow,generate_jastrow1,generate_jastrow2,generate_jastrow3,generate_opt,generate_opts
from auxiliary import *


#set the machine if known, otherwise user will provide
hostmachine = Machine.get_hostname()
if Machine.exists(hostmachine):
    Job.machine = hostmachine
    ProjectManager.machine = Machine.get(hostmachine)
#end if



class Settings(Pobj):
    singleton = None

    project_vars = set([
            'machine','account','generate_only','verbose','debug','mode',
            'local_directory','remote_directory','runs','results','sleep',
            'file_locations','load_images','trace','machine_mode','stages',
            'pseudo_dir','skip_submit','interactive_cores','monitor',
            'status_only','machine_info',
            ])

    gamess_vars = set('ericfmt mcppath'.split())

    all_vars = project_vars | gamess_vars

    @staticmethod
    def kw_set(vars,source=None):
        kw = obj()
        if source!=None:
            for n in vars:
                if n in source:
                    kw[n]=source[n]
                    del source[n]
                #end if
            #end for
        #end if
        return kw
    #end def null_kw_set


    def __init__(self):
        if Settings.singleton is None:
            Settings.singleton = self
        else:
            self.error('attempted to create a second Settings object\n  please just use the original')
        #end if
    #end def __init__


    def error(self,message,header='settings',exit=True,trace=True):
        Pobj.error(self,message,header,exit,trace)
    #end def error


    def __call__(self,**kwargs):

        # guard against invalid settings
        not_allowed = set(kwargs.keys()) - self.all_vars
        if len(not_allowed)>0:
            self.error('unrecognized variables provided.\n  You provided: '+str(list(not_allowed))+'\n  Allowed variables are: '+str(list(vars)))
        #end if

        # extract settings based on keyword groups
        kw = Settings.kw_set(self.project_vars,kwargs)       # project settings
        gamess_kw = Settings.kw_set(self.gamess_vars,kwargs) # gamess settings
        if len(kwargs)>0:
            self.error('some settings keywords have not been accounted for\nleftover keywords: {0}\nthis is a developer error'.format(sorted(kwargs.keys())))
        #end if

        # transfer project variables to project base class
        for name,value in kw.iteritems():
            Pobj.__dict__[name] = value
        #end for

        # process project manager settings
        if 'debug' in kw and kw.debug:
            Pobj.verbose = True
        #end if
        if 'mode' in kw:
            Pobj.set_mode(kw.mode)
        #end if

        # process machine settings
        if 'machine_info' in kw:
            machine_info = Pobj.machine_info
            del Pobj.machine_info
            if isinstance(machine_info,dict) or isinstance(machine_info,obj):
                for machine_name,minfo in machine_info.iteritems():
                    mname = machine_name.lower()
                    if Machine.exists(mname):
                        machine = Machine.get(mname)
                        machine.incorporate_user_info(minfo)
                    else:
                        self.error('machine {0} is unknown\n  cannot set machine_info'.format(machine_name))
                    #end if
                #end for
            else:
                self.error('machine_info must be a dict or obj\n  you provided type '+machine_info.__class__.__name__)
            #end if
        #end if
        if 'machine' in kw:
            machine_name = kw.machine
            if not Machine.exists(machine_name):
                Pobj.class_error('machine {0} is unknown'.format(machine_name))
            #end if
            Job.machine = machine_name
            ProjectManager.machine = Machine.get(machine_name)
            #del Pobj.machine
            if 'account' in kw:
                account = Pobj.account
                #del Pobj.account
                if not isinstance(account,str):
                    self.error('account for {0} must be a string\n  you provided: {1}'.format(machine_name,account))
                #end if
                ProjectManager.machine.account = account
            #end if
            if 'machine_mode' in kw:
                machine_mode = kw.machine_mode
                if machine_mode in Machine.modes:
                    machine_mode = Machine.modes[machine_mode]
                #end if
                if machine_mode==Machine.modes.interactive:
                    if ProjectManager.machine==None:
                        ProjectManager.class_error('no machine specified for interactive mode')
                    #end if
                    if not isinstance(ProjectManager.machine,Supercomputer):
                        self.error('interactive mode is not supported for machine type '+ProjectManager.machine.__class__.__name__)
                    #end if
                    if not 'interactive_cores' in kw:
                        self.error('interactive mode requested, but interactive_cores not set')
                    #end if
                    ProjectManager.machine = ProjectManager.machine.interactive_representation(Pobj.interactive_cores)
                    del Pobj.interactive_cores
                #end if
                del Pobj.machine_mode
            #end if
        #end if

        # process simulation settings
        if 'local_directory' in kw:
            Pobj.file_locations.append(kw.local_directory)
        #end if
        if 'skip_submit' in kw:
            Pobj.skip_submission = Pobj.skip_submit
            del Pobj.skip_submit
        #end if
        if 'file_locations' in kw:
            fl = kw.file_locations
            if isinstance(fl,str):
                Pobj.file_locations.extend([fl])
            else:
                Pobj.file_locations.extend(list(fl))
            #end if
        #end if
        if not 'pseudo_dir' in kw:
            Pobj.pseudopotentials = Pseudopotentials()
        else:
            pseudo_dir = kw.pseudo_dir
            Pobj.file_locations.append(pseudo_dir)
            if not os.path.exists(pseudo_dir):
                self.error('pseudo_dir "{0}" does not exist'.format(pseudo_dir),trace=False)
            #end if
            files = os.listdir(pseudo_dir)
            ppfiles = []
            for f in files:
                pf = os.path.join(pseudo_dir,f)
                if os.path.isfile(pf):
                    ppfiles.append(pf)
                #end if
            #end for
            Pobj.pseudopotentials = Pseudopotentials(ppfiles)        
        #end if

        # more simulation/project manager settings processing
        mode = Pobj.mode
        modes = Pobj.modes
        if mode==modes.stages:
            stages = Pobj.stages
        elif mode==modes.all:
            stages = list(Pobj.primary_modes)
        else:
            stages = [kw.mode]
        #end if
        allowed_stages = set(Pobj.primary_modes)
        if isinstance(stages,str):
            stages = [stages]
        #end if
        if len(stages)==0:
            stages = list(Pobj.primary_modes)
            #self.error('variable stages must be a list of primary modes.\n  Options are '+str(list(allowed_stages)))
        elif 'all' in stages:
            stages = list(Pobj.primary_modes)
        else:
            forbidden = set(Pobj.stages)-allowed_stages
            if len(forbidden)>0:
                self.error('some stages provided are not primary stages.\n  You provided '+str(list(forbidden))+'\n  Options are '+str(list(allowed_stages)))
            #end if
        #end if
        Pobj.mode = modes.stages
        Pobj.stages     = stages
        Pobj.stages_set = set(Pobj.stages)

        # process gamess settings
        Gamess.settings(**gamess_kw)

        return
    #end def __call__
#end class Settings


settings = Settings()


def run_project(*args,**kwargs):
    pm = ProjectManager()
    pm.add_simulations(*args,**kwargs)
    pm.run_project()
#end def run_project
