
import os

from generic import obj
from project_base import Pobj,modes

from project_manager import ProjectManager
from machines import Job,Machine,Supercomputer

from structure import Structure,generate_structure,generate_cell
from physical_system import PhysicalSystem,generate_physical_system
from pseudopotential import Pseudopotential,Pseudopotentials
from bundle import bundle
from sqd import Sqd, SqdInput, SqdAnalyzer, generate_sqd_input, generate_sqd, hunds_rule_filling
from pwscf      import Pwscf     , PwscfInput     , PwscfAnalyzer     , generate_pwscf_input, generate_pwscf
from pw2casino import Pw2casino, Pw2casinoInput, Pw2casinoAnalyzer, generate_pw2casino_input, generate_pw2casino
from pw2qmcpack import Pw2qmcpack, Pw2qmcpackInput, Pw2qmcpackAnalyzer, generate_pw2qmcpack_input, generate_pw2qmcpack
from wfconvert  import Wfconvert , WfconvertInput , WfconvertAnalyzer , generate_wfconvert_input, generate_wfconvert
from qmcpack    import Qmcpack   , QmcpackInput   , QmcpackAnalyzer   , generate_qmcpack_input, generate_qmcpack
from qmcpack import loop,linear,cslinear,vmc,dmc
from qmcpack import generate_jastrows,generate_jastrow,generate_jastrow1,generate_jastrow2,generate_jastrow3,generate_opt,generate_opts


#set the machine if known, otherwise user will provide
hostmachine = Machine.get_hostname()
if Machine.exists(hostmachine):
    Job.machine = hostmachine
    ProjectManager.machine = Machine.get(hostmachine)
#end if



class Settings(Pobj):
    singleton = None

    def __init__(self):
        if Settings.singleton is None:
            Settings.singleton = self
        else:
            self.error('attempted to create a second Settings object\n  please just use the original')
        #end if
    #end def __init__


    def __call__(self,**kwargs):
        vars = set(['machine','account','generate_only','verbose','debug','mode',
                    'local_directory','remote_directory','runs','results','sleep',
                    'file_locations','load_images','trace','machine_mode','stages',
                    'pseudo_dir','skip_submit','interactive_cores','monitor',
                    'status_only','machine_info'])

        keys = set(kwargs.keys())
        not_allowed = keys - vars
        if len(not_allowed)>0:
            Pobj.class_error('unrecognized variables provided.\n  You provided: '+str(list(not_allowed))+'\n  Allowed variables are: '+str(list(vars)),'settings')
        #end if

        for name,value in kwargs.iteritems():
            Pobj.__dict__[name] = value
        #end for

        if 'debug' in kwargs and kwargs['debug']:
            Pobj.verbose = True
        #end if
        if 'mode' in kwargs:
            Pobj.set_mode(kwargs['mode'])
        #end if
        if 'machine_info' in kwargs:
            machine_info = Pobj.machine_info
            del Pobj.machine_info
            if isinstance(machine_info,dict) or isinstance(machine_info,obj):
                for machine_name,minfo in machine_info.iteritems():
                    mname = machine_name.lower()
                    if Machine.exists(mname):
                        machine = Machine.get(mname)
                        machine.incorporate_user_info(minfo)
                    else:
                        Pobj.class_error('machine {0} is unknown\n  cannot set machine_info'.format(machine_name),'settings')
                    #end if
                #end for
            else:
                Pobj.class_error('machine_info must be a dict or obj\n  you provided type '+machine_info.__class__.__name__,'settings')
            #end if
        #end if
        if 'machine' in kwargs:
            machine_name = kwargs['machine']
            if not Machine.exists(machine_name):
                Pobj.class_error('machine {0} is unknown'.format(machine_name),'settings')
            #end if
            Job.machine = machine_name
            ProjectManager.machine = Machine.get(machine_name)
            #del Pobj.machine
            if 'account' in kwargs:
                account = Pobj.account
                #del Pobj.account
                if not isinstance(account,str):
                    Pobj.class_error('account for {0} must be a string\n  you provided: {1}'.format(machine_name,account),'settings')
                #end if
                ProjectManager.machine.account = account
            #end if
            if 'machine_mode' in kwargs:
                machine_mode = kwargs['machine_mode']
                if machine_mode in Machine.modes:
                    machine_mode = Machine.modes[machine_mode]
                #end if
                if machine_mode==Machine.modes.interactive:
                    if ProjectManager.machine==None:
                        ProjectManager.class_error('no machine specified for interactive mode','settings')
                    #end if
                    if not isinstance(ProjectManager.machine,Supercomputer):
                        Pobj.class_error('interactive mode is not supported for machine type '+ProjectManager.machine.__class__.__name__,'settings')
                    #end if
                    if not 'interactive_cores' in kwargs:
                        Pobj.class_error('interactive mode requested, but interactive_cores not set','settings')
                    #end if
                    ProjectManager.machine = ProjectManager.machine.interactive_representation(Pobj.interactive_cores)
                    del Pobj.interactive_cores
                #end if
                del Pobj.machine_mode
            #end if
        #end if
        if 'local_directory' in kwargs:
            Pobj.file_locations.append(kwargs['local_directory'])
        #end if
        if 'skip_submit' in kwargs:
            Pobj.skip_submission = Pobj.skip_submit
            del Pobj.skip_submit
        #end if
        if 'file_locations' in kwargs:
            fl = kwargs['file_locations']
            if isinstance(fl,str):
                Pobj.file_locations.extend([fl])
            else:
                Pobj.file_locations.extend(list(fl))
            #end if
        #end if
        if not 'pseudo_dir' in kwargs:
            Pobj.pseudopotentials = Pseudopotentials()
        else:
            pseudo_dir = kwargs['pseudo_dir']
            Pobj.file_locations.append(pseudo_dir)
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

        mode = Pobj.mode
        modes = Pobj.modes
        if mode==modes.stages:
            stages = Pobj.stages
        elif mode==modes.all:
            stages = list(Pobj.primary_modes)
        else:
            stages = [kwargs['mode']]
        #end if
        allowed_stages = set(Pobj.primary_modes)
        if isinstance(stages,str):
            stages = [stages]
        #end if
        if len(stages)==0:
            stages = list(Pobj.primary_modes)
            #Pobj.class_error('variable stages must be a list of primary modes.\n  Options are '+str(list(allowed_stages)),'settings')
        elif 'all' in stages:
            stages = list(Pobj.primary_modes)
        else:
            forbidden = set(Pobj.stages)-allowed_stages
            if len(forbidden)>0:
                Pobj.class_error('some stages provided are not primary stages.\n  You provided '+str(list(forbidden))+'\n  Options are '+str(list(allowed_stages)),'settings')
            #end if
        #end if
        Pobj.mode = modes.stages
        Pobj.stages     = stages
        Pobj.stages_set = set(Pobj.stages)

        return
    #end def __call__
#end class Settings


settings = Settings()

