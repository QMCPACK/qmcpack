
import os
import traceback
import gc as garbage_collector
from memory import resident
from generic import obj
from developer import DevBase

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


class Pobj(DevBase):
    gc = garbage_collector
    gc.enable()
    mode = modes.stages
    generate_only = False
    status_only   = False
    monitor       = True
    skip_submission = False
    emulate       = False #override application and use application_emulator
    verbose       = True
    debug         = False
    trace         = False
    load_images   = True

    sleep = 3
    local_directory  = './'
    remote_directory = local_directory
    file_locations   = [local_directory]
    runs     = 'runs'
    results  = 'results'
    #pseudo_dir = os.path.join(local_directory,'pseudopotentials')
    pseudo_dir = None
    pseudopotentials = None

    modes = modes

    primary_modes = ['setup','send_files','submit','get_output','analyze']
    dependent_modes = set(['submit'])
    stages = []
    stages_set = set(stages)


    wrote_something = False # for pretty printing

    @staticmethod
    def set_mode(mode):
        if mode in Pobj.modes:
            Pobj.mode = Pobj.modes[mode]
        else:
            print 'settings Error: invalid mode specified: '+mode+'\n  valid modes are '+str(Pobj.modes.keys())
        #end if
    #end def set_mode

    def mem_usage(self):
        return int(resident()/1e6)
    #end def mem_usage

    indent = '  '
    def log(self,*texts,**kwargs):
        if self.verbose:
            if len(kwargs)>0:
                n = kwargs['n']
            else:
                n=0
            #end if
            text=''
            for t in texts:
                text+=str(t)+' '
            #end for
            pad = n*self.indent
            self.logfile.write(pad+text.replace('\n','\n'+pad)+'\n')
        #end if
        Pobj.wrote_something = True
    #end def log

    #@classmethod
    #def class_error(cls,msg,source=None,n=0,trace=True):
    #    if source==None:
    #        source = cls.__name__
    #    #end if
    #    pad = n*cls.indent
    #    text=pad+source+' Error: '+msg
    #    text = '\n'+text.replace('\n','\n'+pad)+'\n\n'
    #    cls.logfile.write(text)
    #    if trace:
    #        traceback.print_stack()
    #    #end if
    #    exit()
    ##end def class_error
    #
    #@classmethod
    #def class_warn(cls,msg,source=None,n=0):
    #    if source==None:
    #        source = cls.__name__
    #    #end if
    #    pad = n*cls.indent
    #    text=pad+source+' Warning: '+msg
    #    cls.logfile.write(text.replace('\n','\n'+pad)+'\n')
    ##end def class_warn

    def dlog(self,*texts,**kwargs):
        if self.debug:
            #self.log('mem_usage',self.mem_usage(),n=5)
            self.log(*texts,**kwargs)
        #end if
    #end def dlog

    def tlog(self,*texts,**kwargs):
        if self.trace:
            self.log(*texts,**kwargs)
            w,s,j,f,g,a=int(self.setup),int(self.submitted),int(self.job.finished),int(self.finished),int(self.got_output),int(self.analyzed)
            self.log('w,s,j,f,g,a',w,s,j,f,g,a,n=kwargs['n']+1)
            #self.log('dependencies',self.dependencies.keys(),n=kwargs['n']+1)
            #self.log('dependents  ',self.dependents.keys(),n=kwargs['n']+1)
        #end if
    #end def tlog

    working_directory = None
    def enter(self,directory,changedir=True,msg=''):
        self.working_directory = os.getcwd()
        self.log('    Entering '+directory,msg)
        if changedir:
            os.chdir(directory)
        #end if
        pad = '      '
        return pad
    #end def enter

    def leave(self):
        os.chdir(self.working_directory)
    #end def leave
#end class Pobj


