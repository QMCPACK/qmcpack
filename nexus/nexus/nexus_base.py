##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  nexus_base.py                                                     #
#    Provides base class functionality and access to 'global' data   #
#    for core Nexus classes.                                         #
#                                                                    #
#  Content summary:                                                  #
#    NexusCore                                                       #
#      Base class for core Nexus classes.                            #
#      Data intended to be 'global' among core classes is assigned   #
#      by the Settings class.                                        #
#                                                                    #
#    nexus_core                                                      #
#      Namespace to be accessed by core Nexus classes.               #
#      These are classes that inherit from NexusCore directly.       #
#                                                                    #
#    nexus_noncore                                                   #
#      Namespace to be accessed in read-only fashion by non-core     #
#      classes.                                                      #
#                                                                    #
#====================================================================#


import os
import gc as garbage_collector
from os import PathLike
from copy import deepcopy
import pickle
from pickle import UnpicklingError
from pathlib import Path
from .utilities import path_string
from .nexus_version import nexus_version
from .memory import resident
from .developer import DevBase, obj, log


# Nexus namespaces
#  nexus_core:   to be used by NexusCore classes only
#  nexus_noncore: allows read only access to some nexus_core data to non-core classes
nexus_core    = obj()
nexus_noncore = obj()
nexus_core_noncore = obj()

status_modes = obj(
    none     = 0,
    standard = 1,
    active   = 2,
    failed   = 3,
    ready    = 4,
    )

modes = obj(
    none       = 0,
    setup      = 1,
    send_files = 2,
    submit     = 3,
    get_output = 4,
    analyze    = 5,
    stages     = 6,
    all        = 7
    )

garbage_collector.enable()


nexus_noncore_defaults = obj(
    basis_dir         = None,
    basissets         = None,
    )

# core namespace elements that can be accessed by noncore classes
nexus_core_noncore_defaults = obj(
    pseudo_dir = None, # used by: Settings, VaspInput
    )

nexus_core_defaults = obj(
    status_only       = False,             # used by: ProjectManager
    generate_only     = False,             # used by: Simulation,Machine
    sleep             = 3,                 # used by: ProjectManager
    timeout           = 5*60,              # used by: Simulation
    runs              = 'runs',            # used by: Simulation,Machine
    results           = '',                # used by: Simulation
    local_directory   = './',              # used by: Simulation,Machine
    remote_directory  = './',              # used by: Simulation
    file_locations    = ['./'],            # used by: Simulation
    monitor           = True,              # used by: ProjectManager,Simulation,Machine
    skip_submit       = False,             # used by: Simulation
    load_images       = True,              # used by: ProjectManager
    modes             = modes,             # used by: ProjectManager,Simulation
    mode              = modes.stages,      # used by: Simulation
    stages_set        = set(),             # used by: ProjectManager,Simulation
    stages            = [],                # used by: Simulation
    primary_modes     = ['setup','send_files','submit','get_output','analyze'], # used by: Settings
    dependent_modes   = set(['submit']),   # used by: ProjectManager,Simulation
    verbose           = True,              # used by: NexusCore
    debug             = False,             # used by: NexusCore
    trace             = False,             # used by: NexusCore
    indent            = '  ',              # used by: NexusCore
    status_modes      = status_modes,      # used by: ProjectManager
    status            = status_modes.none, # used by: ProjectManager
    emulate           = False,             # unused
    progress_tty      = False,             # used by: ProjectManager
    graph_sims        = False,             # used by: ProjectManager
    command_line      = True,              # used by: Settings
    dynamic           = False,             # used by: DynamicWorkflowManager
                                           #          Simulation
    **nexus_core_noncore_defaults
    )

def restore_nexus_core_defaults():
    nexus_core.clear()
    nexus_noncore.clear()
    nexus_core_noncore.clear()

    nexus_core.update(**deepcopy(nexus_core_defaults))
    nexus_noncore.update(**deepcopy(nexus_noncore_defaults))
    for k in nexus_core_noncore_defaults.keys():
        nexus_core_noncore[k] = nexus_core[k]
#end def restore_nexus_core_defaults

restore_nexus_core_defaults()


nexus_core_no_process = {'status_only', 'generate_only', 'sleep', 'timeout'}

nexus_modules = [mod.stem for mod in Path(__file__).parent.iterdir() if mod.suffix == ".py"]

class NexusUnpickler(pickle.Unpickler):
    """This class is designed for backwards compatibility with pickles generated
    before Nexus was packaged (PR #5700, December 20, 2025). 
    It shouldn't touch anything but old Nexus pickles.
    """
    def find_class(self, module, name):
        if module in nexus_modules and "nexus." not in module:
            module = "nexus." + module
        if module == "nexus.generic":
            if name == "obj":
                module = "nexus.developer_tools"
            elif name == "DevBase":
                module = "nexus.developer"

        return super().find_class(module, name)


class NexusCore(DevBase):

    # garbage collector
    gc = garbage_collector

    # mutable/dynamic nexus core data
    wrote_something   = False # for pretty printing
    working_directory = None
    wrote_splash      = False

    @staticmethod
    def write_splash():
        if not NexusCore.wrote_splash:
            splash_text = '''
_____________________________________________________

                     Nexus {}.{}.{}

        (c) Copyright 2012-  Nexus developers

                     Please cite:
  J. T. Krogel Comput. Phys. Commun. 198 154 (2016)
     https://doi.org/10.1016/j.cpc.2015.08.012
_____________________________________________________
          
            '''.format(*nexus_version)
            log(splash_text)
            NexusCore.wrote_splash = True
        #end if
    #end def write_splash

    @staticmethod
    def write_end_splash():
        return # don't do this yet
        splash_text = '''
_____________________________________________________
_____________________________________________________
            '''
        print(splash_text)
    #end def write_end_splash

    def mem_usage(self):
        return int(resident()/1e6)
    #end def mem_usage

    def log(self,*texts,**kwargs):
        """Write output to log file.

        Parameters
        ----------
        *texts
            Strings that will be joined by newlines
        n : int, kwargs
            Spaces to indent
        progress : bool, kwargs
            If ``True`` and output is to a terminal, overwrite and update the
            last line, rather than scrolling.
        """
        if nexus_core.verbose:
            if len(kwargs)>0:
                n = kwargs['n']
            else:
                n=0
            #end if
            is_progress = kwargs.get('progress',False)
            text=''
            for t in texts:
                text+=str(t)+' '
            #end for
            pad = n*nexus_core.indent
            output_text = pad+text.replace('\n','\n'+pad)
            if nexus_core.progress_tty and is_progress and self._logfile.isatty():
                # spaces to ensure previous line is overwritten.  Need better solution.
                self._logfile.write(output_text+'        \r')
                self._logfile.flush()
            else:
                self._logfile.write(output_text+'\n')
        #end if
        NexusCore.wrote_something = True
    #end def log

    def dlog(self,*texts,**kwargs):
        if nexus_core.debug:
            #self.log('mem_usage',self.mem_usage(),n=5)
            self.log(*texts,**kwargs)
        #end if
    #end def dlog

    def tlog(self,*texts,**kwargs):
        if nexus_core.trace:
            self.log(*texts,**kwargs)
            w,s,j,f,g,a=int(self.setup),int(self.submitted),int(self.job.finished),int(self.finished),int(self.got_output),int(self.analyzed)
            self.log('w,s,j,f,g,a',w,s,j,f,g,a,n=kwargs['n']+1)
            #self.log('dependencies',self.dependencies.keys(),n=kwargs['n']+1)
            #self.log('dependents  ',self.dependents.keys(),n=kwargs['n']+1)
        #end if
    #end def tlog

    def enter(self, directory: PathLike, *, changedir: bool = True, msg: str = ''):
        """Have Nexus enter a directory and change its current working directory.
        
        Parameters
        ----------
        directory : PathLike
            Directory to enter. Can be a ``str`` or ``pathlib.Path``
            object.
        changedir : bool, default=True
            Default of ``True`` will change the CWD, setting to ``False``
            will not change the CWD.
        msg : str, optional
            Optional message to pass to the output log.
        """
        NexusCore.working_directory = os.getcwd()
        directory = path_string(directory)

        self.log('    Entering ' + directory, msg)
        if changedir:
            os.chdir(directory)
        #end if
        pad = '      '
        return pad
    #end def enter

    def leave(self):
        os.chdir(NexusCore.working_directory)
    #end def leave

    def load(self, fpath: PathLike | None = None):
        if fpath is None:
            fpath = f'./{type(self).__name__}.p'

        with open(fpath, 'rb') as fobj:
            try:
                tmp = pickle.load(fobj)
            except (ImportError, ModuleNotFoundError):
                fobj.seek(0)
                try:
                    # Old pickles from before Nexus was packaged (PR #5700, December 20 2025)
                    # won't have the correct module path. The custom unpickler will handle this by 
                    # prepending "nexus." to the module path
                    tmp = NexusUnpickler(fobj).load()
                except UnpicklingError:
                    # NumPy pickles can use latin1 encoding
                    # They will likely still fail from an underflow since they are not pickle-compliant
                    tmp = NexusUnpickler(fobj).load(encoding='latin1')

        d = self.__dict__
        d.clear()
        for k, v in tmp.__dict__.items():
            d[k] = v
    #end def load
#end class NexusCore


# support dynamic workflows
dynamic_storage = obj(
    simulations         = obj(), # all sims, in dyn proc or not
    simulation_ids      = set(),
    dynamic_processes   = obj(),
    dynamic_process_ids = set(),
    )
