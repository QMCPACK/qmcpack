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
import traceback
import gc as garbage_collector
from memory import resident
from generic import obj
from developer import DevBase


# Nexus namespaces
#  nexus_core:   to be used by NexusCore classes only
#  nexus_noncore: allows read only access to some nexus_core data to non-core classes
nexus_core    = obj()
nexus_noncore = obj()

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


nexus_noncore.set(
    basis_dir         = None,
    basissets         = None,
    )

# core namespace elements that can be accessed by noncore classes
nexus_core_noncore = obj(
    pseudo_dir        = None,              # used by: Settings, VaspInput
    pseudopotentials  = None,              # used by: Simulation, GamessInput
    )

nexus_core.set(
    status_only       = False,             # used by: ProjectManager
    generate_only     = False,             # used by: Simulation,Machine
    sleep             = 3,                 # used by: ProjectManager
    runs              = 'runs',            # used by: Simulation,Machine
    results           = 'results',         # used by: Simulation
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
    **nexus_core_noncore
    )

nexus_core_no_process = set('''
  status_only  generate_only  sleep
  '''.split())


class NexusCore(DevBase):

    # garbage collector
    gc = garbage_collector

    # mutable/dynamic nexus core data
    wrote_something   = False # for pretty printing
    working_directory = None
    wrote_splash      = False

    @staticmethod
    def write_splash():
        return # don't do this yet
        if not NexusCore.wrote_splash:
            splash_text = '''
_____________________________________________________
   _         _______                       _______   
  | \    /| |  ____ \ |\     /| |\     /| |  ____ \  
  |  \  | | | |    \/ | \   / | | |   | | | |    \/  
  |   \ | | | |__      \ \_/ /  | |   | | | |_____   
  | |\ \| | |  __|      | _ |   | |   | | |_____  |  
  | | \   | | |        / / \ \  | |   | |       | |  
  | |  \  | | |____/\ | /   \ | | |___| | /\____| |  
  |/    \_| |_______/ |/     \| |_______| \_______|  
_____________________________________________________
                                               
            Main author: Jaron T. Krogel            
            ____________________________
          
            '''
            print splash_text
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
        print splash_text
    #end def write_end_splash

    def mem_usage(self):
        return int(resident()/1e6)
    #end def mem_usage

    def log(self,*texts,**kwargs):
        """Write output to log file.
           Keyword arguments
            n - spaces to indent
            progress - if True and output is to a terminal, overwrite and
                       update the last line, rather than scrolling.
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

    def enter(self,directory,changedir=True,msg=''):
        NexusCore.working_directory = os.getcwd()
        self.log('    Entering '+directory,msg)
        if changedir:
            os.chdir(directory)
        #end if
        pad = '      '
        return pad
    #end def enter

    def leave(self):
        os.chdir(NexusCore.working_directory)
    #end def leave
#end class NexusCore
