
import os
import time
#from multiprocessing import cpu_count
from socket import gethostname
from subprocess import Popen,PIPE
from numpy import array,mod,floor,ceil,round,log
from generic import obj
from project_base import Pobj
from debug import *


import re,subprocess
def  cpu_count():
    """ Number of virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError,NotImplementedError):
        None

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError,ValueError):
        None
#end def cpu_count



class Options(Pobj):
    def __init__(self,**kwargs):
        self.add(**kwargs)
    #end def __init__


    def read(self,options):
        nopts = 0
        intext = False
        nstart = -2
        nend   = -2
        n = 0
        for c in options:
            if c=='-' and nstart!=n-1 and not intext:
                prevdash=True
                if nopts>0:
                    opt  = options[nstart:n].strip()
                    name = opt.replace('-','').replace('=',' ').split()[0]
                    self[name]=opt
                #end if
                nopts+=1
                nstart = n
            elif c=='"' or c=="'":
                intext=not intext
            #end if
            n+=1
        #end for
        if nopts>0:
            opt  = options[nstart:n].strip()
            name = opt.replace('-','').replace('=',' ').split()[0]
            self[name]=opt
        #end if
    #end def read


    def add(self,**kwargs):
        self.transfer_from(kwargs)
    #end def add


    def write(self):
        s = ''
        for o in self:
            s += ' '+str(o)
        #end for
        return s
    #end def write
#end class Options



class Job(Pobj):

    machine = None #default machine if none is specified in settings

    states = obj(
        none      = 0,
        waiting   = 1,
        submitted = 2,
        running   = 3,
        finished  = 4
        )
    state_names = states.inverse()

    job_count = 0


    intel_compilers = set(['intel','icc','icpc','ifort'])
    gnu_compilers   = set(['gnu','gcc','g++','gfortran'])
    pgi_compilers   = set(['pgi'])
    cray_compilers  = set(['cray'])


    @staticmethod
    def generate_jobid():
        Job.job_count += 1
        return Job.job_count
    #end def generate_jobid


    @classmethod
    def zero_time(cls):
        time = obj(days=0,hours=0,minutes=0,seconds=0)
        return time
    #end def zero_time


    def __init__(self,
                 name         = 'jobname',
                 directory    = None,
                 sub_dir      = None,
                 app_name     = None,
                 app_flags    = None,
                 app_command  = None,
                 app_props    = None,
                 app          = None,
                 env          = None,
                 user_env     = True,
                 presub       = '',
                 postsub      = '',
                 outfile      = None,
                 errfile      = None,
                 mode         = None,
                 machine      = None,
                 account      = None,
                 queue        = None,
                 bundled_jobs = None,
                 relative     = False,
                 cores        = None,
                 nodes        = None,
                 threads      = 1,
                 ppn          = None,
                 compiler     = None,
                 override     = None,
                 options      = None,
                 app_options  = None,
                 run_options  = None,
                 sub_options  = None,
                 serial       = False,
                 days         = 0,
                 hours        = 0,
                 minutes      = 0,
                 seconds      = 0,
                 subfile      = None,
                 grains       = None,
                 procs        = None,
                 processes    = None,
                 processes_per_proc = None,
                 processes_per_node = None
                 ):

        self.directory   = directory
        self.subdir      = sub_dir
        self.app_name    = app_name
        self.app_command = app_command
        self.app_props   = app_props
        self.outfile     = outfile
        self.errfile     = errfile
        self.user_env    = user_env
        self.presub      = presub
        self.postsub     = postsub
        self.name        = name
        self.queue       = queue
        self.bundled_jobs= bundled_jobs
        self.relative    = relative
        self.cores       = cores
        self.nodes       = nodes
        self.threads     = threads
        self.ppn         = ppn
        self.compiler    = compiler
        self.override    = override
        self.app_options = Options()
        self.run_options = Options()
        self.sub_options = Options()
        self.serial      = serial
        self.days        = days
        self.hours       = hours
        self.minutes     = minutes
        self.seconds     = seconds
        self.subfile     = subfile
        self.grains      = grains
        self.procs       = procs
        self.processes   = processes
        self.processes_per_proc = processes_per_proc
        self.processes_per_node = processes_per_node        
        self.account     = account
        self.internal_id = None
        self.system_id   = None
        self.tot_cores   = None
        self.identifier  = None
        self.submitted   = False
        self.status      = self.states.none
        self.crashed     = False
        self.overtime    = False
        self.successful  = False
        self.finished    = False

        if app != None:
            self.app_name = app
        #end if
        if app_options != None:
            self.app_options.read(app_options)
        #end if
        if run_options != None:
            self.run_options.read(run_options)
        #end if
        if sub_options != None:
            self.sub_options.read(sub_options)
        #end if
        if app_flags != None:
            self.app_options.read(app_flags)
        #end if
        if options != None:
            self.run_options.read(options)
        #end if
        if env is None:
            self.env = None
        else:
            self.set_environment(**env)
        #end if

        if app_props==None:
            self.app_props = []
        #end if
        if compiler!=None:
            if compiler in self.intel_compilers:
                self.compiler = 'intel'
            elif compiler in self.gnu_compilers:
                self.compiler = 'gnu'
            elif compiler in self.pgi_compilers:
                self.compiler = 'pgi'
            elif compiler in self.cray_compilers:
                self.compiler = 'cray'
            #end if
        #end if

        if machine!=None:
            self.machine = machine
        #end if
        #check that the machine exists and have it complete the job info
        machine = self.get_machine()
        machine.process_job(self)

        self.batch_mode = machine.in_batch_mode()

        if bundled_jobs!=None: 
            if not machine.batch_capable:
                self.error('running batched/bundled jobs on {0} is either not possible or not yet implemented, sorry.'.format(machine.name))
            #end if
        #end if

        self.normalize_time()
    #end def __init__


    def get_machine(self):
        return Machine.get(self.machine)
    #end def get_machine


    def initialize(self,sim):
        self.set_id()
        self.identifier = sim.identifier
        machine = self.get_machine()
        if machine.prefixed_output:
            sim.outfile = sim.identifier + machine.outfile_extension
            sim.errfile = sim.identifier + machine.errfile_extension
        #end if
        if self.directory is None:
            self.directory = sim.locdir
            self.abs_dir    = os.path.abspath(sim.locdir)
        elif self.abs_dir is None:
            self.abs_dir = os.path.abspath(self.directory)
        #end if
        if self.subdir is None:
            if machine.local_directory!=None:
                self.subdir = os.path.join(machine.local_directory,sim.runs,sim.path)
                self.abs_subdir = self.subdir
            else:
                self.subdir = self.directory
                self.abs_subdir = self.abs_dir
            #end if
        #end if
        if self.app_name is None:
            app_name = sim.app_name
        else:
            app_name = self.app_name
        #end if
        if app_name!=None and not '/' in app_name:
            ads = machine.app_directories
            ad  = machine.app_directory
            new_app_name = None
            if ads!=None and self.app_name in ads:
                new_app_name = os.path.join(ads[self.app_name],app_name)
            elif ad!=None:
                new_app_name = os.path.join(ad,app_name)
            #end if
            if new_app_name!=None:
                self.app_name = new_app_name
            #end if
        #end if
        sim.set_app_name(app_name)
        self.set(
            name        = sim.identifier,
            simid       = sim.simid,
            app_command = sim.app_command(),
            outfile     = sim.outfile,
            errfile     = sim.errfile
            )
        if self.app_props==None:
            self.app_props   = list(sim.app_props)
        #end if
    #end def initialize


    def set_id(self):
        self.internal_id = Job.generate_jobid()
    #end def set_id


    def set_processes(self):
        if self.processes is None:
            self.error('processes should have been set before now\n  contact the developers and have them fix this','Developer')
            self.processes = int(ceil(float(self.cores)/self.threads))
        #end if
    #end def set_processes

            
    def set_environment(self,limited_env=False,clear_env=False,**env):
        machine = self.get_machine()
        if isinstance(machine,Supercomputer):
            limited_env = True
        #end if
        if self.env is None:
            self.env = os.environ.copy()
            if limited_env:
                self.env.clear()
            #end if
        #end if
        if clear_env:
            self.env.clear()
        #end if
        for n,v in env.iteritems():
            self.env[n]=str(v)
        #end for
    #end def set_environment


    def get_time(self):
        time = obj(
            days = self.days,
            hours = self.hours,
            minutes = self.minutes,
            seconds = self.seconds
            )
        return time
    #end def get_time

    
    def max_time(self,time):
        t  = time.seconds + 60*(time.minutes+60*(time.hours+24*time.days))
        ts = self.seconds + 60*(self.minutes+60*(self.hours+24*self.days))
        if ts>t:
            time.days    = self.days
            time.hours   = self.hours
            time.minutes = self.minutes
            time.seconds = self.seconds                
        #end if
        return time
    #end def max_time


    def serial_only(self):
        return 'serial' in self.app_props and len(self.app_props)==1
    #end if


    def determine_end_status(self,status):
        if not self.generate_only:
            self.successful = False # not really implemented yet
            #errfile = os.path.join(self.directory,self.errfile)
            #if os.path.exists(errfile):
            #    fobj = open(errfile,'r')
            #    errors = fobj.read()
            #    fobj.close()
            #    if errors.strip()=='':
            #        self.successful=True
            #    #end if
            ##end if
        #end if
    #end def determine_end_status


    def write(self):
        machine = self.get_machine()
        return machine.write_job(self)
    #end def write


    def submit(self):
        machine = self.get_machine()
        machine.add_job(self)
        self.submitted = True
    #end def submit


    def reenter_queue(self):
        machine = self.get_machine()
        machine.requeue_job(self)
    #end def reenter_queue


    def run_command(self,launcher,redirect=False):
        c = ''
        if self.bundled_jobs is None:
            if isinstance(self.override,str):
                c = self.override
            else:
                if self.app_command is None:
                    self.error('app_command has not been provided')
                #end if
                if launcher=='runjob':
                    separator = ' : '
                else:
                    separator = ' '
                #end if
                if self.serial and self.processes==1:
                    c = ''
                else:
                    c = launcher + self.run_options.write() + separator
                #end if
                c+=self.app_command+self.app_options.write()
                if redirect:
                    c+=' >'+self.outfile+' 2>'+self.errfile+'&'
                #end if
            #end if
        elif self.relative:
            cdir = self.abs_subdir
            c+='\n'
            for job in self.bundled_jobs:
                c+='\ncd '+os.path.relpath(job.abs_subdir,cdir)+'\n'
                c+=job.run_command(launcher,redirect=True)+'\n'
                cdir = job.abs_subdir
            #end for
            c+='\nwait\n'
        else:
            c+='\n'
            for job in self.bundled_jobs:
                c+='\ncd '+job.abs_subdir+'\n'
                c+=job.run_command(launcher,redirect=True)+'\n'
            #end for
            c+='\nwait\n'
        #end if
        return c
    #end def run_command


    def pbs_walltime(self):
        walltime=\
            str(int(self.hours   )).zfill(2)+':'\
            +str(int(self.minutes)).zfill(2)+':'\
            +str(int(self.seconds)).zfill(2)
        if self.days!=0:
            walltime = str(self.days)+':'+walltime
        #end if
        return walltime
    #end def pbs_walltime


    def normalize_time(self):
        t = self.total_seconds()
        d = int(t/(24*3600))
        t -= d*24*3600
        h = int(t/3600)
        t -= h*3600
        m = int(t/60)
        t -= m*60
        s = int(t)
        self.days    = d
        self.hours   = h
        self.minutes = m
        self.seconds = s
    #end def normalize_time


    def total_seconds(self):
        return self.seconds+60*(self.minutes+60*(self.hours+24*self.days))
    #end def total_seconds

    
    def total_minutes(self):
        return int(self.total_seconds()/60)
    #end def total_minutes


    def total_hours(self):
        return int(self.total_seconds()/3600)
    #end def total_hours


    def total_days(self):
        return int(self.total_seconds()/(24*3600))
    #end def total_days
#end class Job




class Machine(Pobj):

    machines = obj()

    modes = obj(
        none        = 0,
        interactive = 1,
        batch       = 2
        )
    mode = modes.none

    batch_capable      = False
    requires_account   = False
    executable_subfile = False

    prefixed_output    = False
    outfile_extension  = None
    errfile_extension  = None

    @staticmethod
    def get_hostname():
        hostname = gethostname()
        if '.' in hostname:
            machine_name = hostname.split('.')[0]
        else:
            machine_name = hostname
        #end if
        return machine_name.lower()
    #end def get_hostname


    @staticmethod
    def exists(machine_name):
        return machine_name in Machine.machines
    #end def exists


    @staticmethod
    def is_unique(machine):
        return id(machine)==id(Machine.machines[machine.name])
    #end def is_valid


    @staticmethod
    def add(machine):
        if not isinstance(machine,Machine):
            Machine.class_error('attempted to add non-machine instance')
        #end if
        if not 'name' in machine:
            Machine.class_error('attempted to add a machine without a name')
        #end if
        name = machine.name
        if not name in Machine.machines:
            Machine.machines[name] = machine
        else:
            Machine.class_error('attempted to create machine {0}, but it already exists'.format(name))
        #end if
    #end def add


    @staticmethod
    def get(machine_name):
        if isinstance(machine_name,str):
            machine_name = machine_name.lower()
        else:
            Machine.class_error('machine name must be a string, you provided a '+machine_name.__class__.__name__)
        #end if
        if Machine.exists(machine_name):
            machine = Machine.machines[machine_name]
        else:
            machs = Machine.machines.keys()
            machs.sort()
            Machine.class_error('attempted to get machine '+machine_name+', but it is unknown\n  known options are '+str(machs))
        #end if
        return machine
    #end def get


    def validate(self):
        if Machine.exists(self.name):
            if not Machine.is_unique(self):
                self.error('duplicate instance of machine '+self.name+' encountered\n  this is either a developer error, or you have created a duplicate machine')
            #end if
        else:
            self.error('machine {0} id {1} was created without calling Machine.__init__() and is therefore invalid'.format(self.name,id(self)))
        #end if
    #end def validate
        

    def in_batch_mode(self):
        return self.mode==self.modes.batch
    #end def in_batch_mode


    def query_queue(self):
        self.not_implemented()
    #end def query_queue

    def submit_jobs(self):
        self.not_implemented()
    #end def submit_jobs

    def process_job(self,job):
        self.not_implemented()
    #end def process_job

    def write_job(self,job):
        self.not_implemented()
    #end def write_job

    def submit_job(self,job):
        self.not_implemented()
    #end def submit_job


    def __init__(self,name,queue_size=0):
        self.name = name
        self.queue_size = queue_size
        self.processes = obj()
        self.jobs = obj()
        self.waiting = set()
        self.running = set()
        self.finished= set()

        #user defined variables
        self.account         = None
        self.local_directory = None
        self.app_directory   = None
        self.app_directories = None

        if not isinstance(name,str):
            self.error('machine name must be a string\n  you provided '+str(name))
        #end if

        Machine.add(self)
    #end def __init__


    def add_job(self,job):
        if isinstance(job,Job):
            self.process_job(job)
            self.write_job(job)
            jid = job.internal_id
            self.jobs[jid] = job
            job.status = job.states.waiting
            self.waiting.add(jid)
            #self.write_job_states('add_job')
        else:
            self.error('add_job received non-Job instance '+job.__class__.__name__)
        #end if
    #end def add_job


    def requeue_job(self):
        None
    #end def requeue_job


    allowed_user_info = set(['account','local_directory','app_directory','app_directories'])
    def incorporate_user_info(self,infoin):
        info = obj(**infoin)
        vars = set(info.keys())
        invalid = vars-self.allowed_user_info
        if len(invalid)>0:
            self.error('invalid inputs encountered in incorporate_user_info\n  allowed inputs: {0}\n  invalid inputs: {1}'.format(list(self.allowed_user_info),list(invalid)))
        #end if
        if 'app_directories' in info:
            ad = info.app_directories
            if not isinstance(ad,dict) and not isinstance(ad,obj):
                self.error('app_directories must be of type dict or obj\n you provided '+ad.__class__.__name__)
            #end if
        #end if
        self.transfer_from(info)
    #end def incorporate_user_info
#end class Machine




class Workstation(Machine):

    mode = Machine.modes.interactive

    batch_capable = False

    def __init__(self,name='workstation',cores=None,app_launcher='mpirun',process_granularity=1):
        Machine.__init__(self,name)
        self.app_launcher = app_launcher
        if cores==None:
            self.cores = cpu_count()
        else:
            self.cores = cores
        #end if
        self.queue_size = cores
        self.process_granularity = process_granularity
    #end def __init__

    
    def process_job(self,job):
        if job.serial_only():
            job.cores=1
        elif job.cores==None:
            if job.processes!=None:
                job.cores = job.processes*job.threads
            else:
                job.cores = self.cores
            #end if
        #end if
        job.processes = max(1,int(floor(float(job.cores)/job.threads)))
        grains = int(ceil(float(job.cores)/self.process_granularity))
        if abs(grains-1-float(job.cores)/self.process_granularity)<1e-6:
            grains-=1
        #end if
        job.grains = grains
        job.cores = grains*self.process_granularity

        job.run_options.add(np='-np '+str(job.processes))
    #end def process_job


    def write_job_states(self,title=''):
        self.log(title,n=2)
        n=3
        self.log('{0} {1} {2} job states'.format(self.__class__.__name__,self.name,id(self)),n=n )
        self.log('processes',n=n+1)
        for process in self.processes:
            job = process.job
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('jobs',n=n+1)
        jobids = list(self.jobs.keys())
        jobids.sort()
        for jobid in jobids:
            job = self.jobs[jobid]
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('waiting',n=n+1)
        jobids = list(self.waiting)
        jobids.sort()
        for jobid in jobids:
            job = self.jobs[jobid]
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('running',n=n+1)
        jobids = list(self.running)
        jobids.sort()
        for jobid in jobids:
            job = self.jobs[jobid]
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('finished',n=n+1)
        jobids = list(self.finished)
        jobids.sort()
        for jobid in jobids:
            job = self.jobs[jobid]
            self.log('{0:>4} {1:>10} {2:>4} {3}'.format(job.internal_id,job.name,job.simid,job.directory),n=n+2)
        #end for
        self.log('end job states',n=1)
    #end def write_job_states

    def query_queue(self):
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~query queue',n=1)
        #self.write_job_states('query queue')
        self.validate()
        done = []
        self.dlog('qq',self.processes.keys(),n=2)
        for pid,process in self.processes.iteritems():
            if self.generate_only or not self.monitor:
                qpid,status = pid,0
            else:
                qpid,status = os.waitpid(pid,os.WNOHANG) 
            #end if
            if pid==qpid:
                job = process.job
                job.status = job.states.finished
                job.finished = True
                job.determine_end_status(status)
                iid = job.internal_id
                self.running.remove(iid)
                self.finished.add(iid)
                done.append(pid)
                if not self.generate_only:
                    job.out.close()
                    job.err.close()
                #end if
            #end if
        #end for
        for pid in done:
            del self.processes[pid]
        #end for
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end query queue',n=1)
    #end def query_queue



    def submit_jobs(self):
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~submit jobs',n=1)
        self.dlog('sj',self.waiting,n=2)
        cores_used = 0
        for process in self.processes:
            cores_used += process.job.cores
        #end for
        cores_available = self.cores-cores_used

        core_req = []
        job_req  = []
        for iid in self.waiting:
            job = self.jobs[iid]
            job_req.append(job)
            core_req.append(job.cores)
        #end for
        core_req = array(core_req)
        job_req  = array(job_req)
        order = core_req.argsort()
        job_req = job_req[order]

        self.dlog('cores available:',cores_available,n=2)
        for job in job_req:
            if job.cores>self.cores and not self.generate_only:                
                self.error('job '+str(job.internal_id)+' is too large to run on this machine\n  cores requested: '+str(job.cores)+'\n  machine cores: '+str(self.cores))
            #end if
            if job.cores<=cores_available:
                iid = job.internal_id
                self.waiting.remove(iid)
                self.running.add(iid)
                self.submit_job(job)
                cores_available-=job.cores
            elif job.cores>self.cores:
                self.error('job requested more cores than are present on '+self.name+'\n  cores requested: {0}\n  cores present: {1}'.format(job.cores,self.cores))
            else:
                break
            #end if
        #end for                
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end submit jobs',n=1)
    #end def submit_jobs


    def write_job(self,job):
        return None
    #end def write_job


    def submit_job(self,job):
        pad = self.enter(job.directory,msg=job.simid)
        command = 'export OMP_NUM_THREADS='+str(job.threads)+'\n'
        if len(job.presub)>0:
            command += job.presub+'\n'
        #end if
        command += job.run_command(self.app_launcher)
        if len(job.postsub)>0:
            command += job.postsub+'\n'
        #end if
        command = ('\n'+command).replace('\n','\n  '+pad)
        job.status = job.states.running
        process = obj()
        process.job = job
        if self.generate_only:
            self.log(pad+'Would have executed:  '+command)
            job.system_id = job.internal_id
        else:
            if self.monitor:
                self.log(pad+'Executing:  '+command)
                job.out = open(job.outfile,'w')
                job.err = open(job.errfile,'w')
                p = Popen(command,env=job.env,stdout=job.out,stderr=job.err,shell=True)
                process.popen = p
                job.system_id = p.pid
            else:
                command+=' >'+job.outfile+' 2>'+job.errfile+'&'
                self.log(pad+'Executing:  '+command)
                os.system(command)
                job.system_id = job.internal_id
            #end if
        #end if
        self.processes[job.system_id] = process
        self.leave()
    #end def submit_job
#end class Workstation




class InteractiveCluster(Workstation):
    def __init__(self,*args,**kwargs):
        
        if len(args)==0 or not isinstance(args[0],Supercomputer):
            self.init_from_args(*args,**kwargs)
        else:
            super = args[0]
            cores = args[1]
            self.init_from_supercomputer(super,cores)
        #end if
        Machine.__init__(self,self.name,self.queue_size)
    #end def __init__


    def init_from_args(self,name='icluster',nodes=None,procs_per_node=None,
                 cores_per_proc=None,process_granularity=None,ram_per_node=None,
                 app_launcher=None):
        self.name           = name
        self.nodes          = nodes
        self.procs_per_node = procs_per_node
        self.cores_per_proc = cores_per_proc
        self.process_granularity = process_granularity
        self.ram_per_node   = ram_per_node
        self.app_launcher   = app_launcher

        self.cores_per_node = self.cores_per_proc*self.procs_per_node
        if process_granularity is None:
            self.process_granularity = self.cores_per_node
        #end if

        self.procs = self.procs_per_node*self.nodes
        self.cores = self.cores_per_proc*self.procs
        self.ram   = self.ram_per_node*self.nodes

        self.queue_size = self.cores
    #end def init_from_args


    def init_from_supercomputer(self,super,cores):
        nodes = cores/super.cores_per_node
        if cores-nodes*super.cores_per_node!=0:
            self.error('interactive cores corresponds to a fractional number of nodes\n  cores '+str(cores)+'\n  cores per node '+str(super.cores_per_node))
        #end if
        self.init_from_args(super.name+'_interactive',nodes,super.procs_per_node,
                            super.cores_per_proc,super.cores_per_node,
                            super.ram_per_node,super.app_launcher)
    #end def init_from_supercomputer

    
    def process_job(self,job):
        job.cores = min(job.cores,self.cores)

        Workstation.process_job(self,job)

        job.nodes = job.grains
        job.procs = job.nodes*self.procs_per_node

        if mod(job.processes,job.nodes)!=0:
            job.processes_per_node = None
        else:
            job.processes_per_node = job.processes/job.nodes
        #end if

        if mod(job.processes,job.procs)!=0:
            job.processes_per_proc = None
        else:
            job.processes_per_proc = job.processes/job.procs
        #end if
    #end def process_job
#end class InteractiveCluster




class Supercomputer(Machine):
    mode = Machine.modes.batch
    name = 'supercomputer'

    batch_capable = False #only set to true for specific machines

    aprun_options = set(['n','d'])

    def __init__(self,nodes=None,procs_per_node=None,
                 cores_per_proc=None,ram_per_node=None,queue_size=0,
                 app_launcher=None,sub_launcher=None,queue_querier=None,
                 job_remover=None,name=None):
        if name is None:
            if self.name!=None:
                name = self.name
            else:
                name = self.__class__.__name__.lower()
            #end if
        #end if
        Machine.__init__(self,name)
        self.nodes          = nodes          #  # of nodes
        self.procs_per_node = procs_per_node #  # of processors/sockets on a node
        self.cores_per_proc = cores_per_proc #  # of cores on a node
        self.ram_per_node   = ram_per_node
        self.queue_size     = queue_size
        self.app_launcher   = app_launcher
        self.sub_launcher   = sub_launcher
        self.queue_querier  = queue_querier
        self.job_remover    = job_remover
     
        required = ['nodes','procs_per_node','cores_per_proc','ram_per_node','queue_size','app_launcher','sub_launcher','queue_querier','job_remover']
        for var in required:
            if self[var]==None:
                self.error('input variable '+var+' is required to initialize Supercomputer object.')
            #end if
        #end for

        self.cores_per_node = self.cores_per_proc*self.procs_per_node

        self.procs = self.procs_per_node*self.nodes
        self.cores = self.cores_per_proc*self.procs
        self.ram   = self.ram_per_node*self.nodes

        self.system_queue = obj()
        if self.queue_querier=='qstat':
            self.job_states=dict(E = 'exiting',
                                 H = 'held',
                                 Q = 'queued',
                                 R = 'running',
                                 S = 'suspended',
                                 T = 'transferring',
                                 W = 'waiting',
                                 C = 'complete'
                                 )
        elif self.queue_querier=='qstata':
            #already gives status as queued, running, etc.
            None
        elif  self.queue_querier=='squeue':
            self.job_states=dict(CG = 'exiting',
                                 TO = 'timeout',
                                 CA = 'failed',
                                 F = 'failed',
                                 NF = 'node_fail',
                                 PD = 'waiting',
                                 R = 'running',
                                 S = 'suspended',
                                 CD = 'complete'
                                 )
        elif  self.queue_querier=='sacct':
            self.job_states=dict(CANCELLED = 'failed',  #long form
                                 COMPLETED = 'complete',
                                 COMPLETING = 'exiting',
                                 CONFIGURING = 'waiting',
                                 FAILED = 'failed',
                                 PREMEEMPTED = 'failed',
                                 PENDING = 'waiting',
                                 NODE_FAIL = 'failed',
                                 RESIZING = 'resizing',
                                 RUNNING = 'running',
                                 SUSPENDED = 'suspended',
                                 TIMEOUT = 'failed',
                                 CA = 'failed',        #short form
                                 CD = 'complete',
                                 CG = 'exiting',
                                 CF = 'waiting',
                                 F = 'failed',
                                 PR = 'failed',
                                 PD = 'waiting',
                                 NF = 'failed',
                                 RS = 'resizing',
                                 R = 'running',
                                 S = 'suspended',
                                 TO = 'failed'
                                 )
        else:
            self.error('ability to query queue with '+self.queue_querier+' has not yet been implemented')
        #end if
            
    #end def __init__


    def interactive_representation(self,cores):
        return InteractiveCluster(self,cores)
    #end def interactive_representation


    def requeue_job(self,job):
        if isinstance(job,Job):
            jid = job.internal_id
            pid = job.system_id
            if pid is None:
                self.error('job {0} does not have a process id issued by the scheduler'.format(jid))
            #end if
            self.process_job(job)
            self.jobs[jid] = job
            job.status = job.states.running
            self.running.add(jid)
            process = obj(job=job)
            self.processes[pid] = process
        else:
            self.error('requeue_job received non-Job instance '+job.__class__.__name__)
        #end if
    #end def requeue_job


    def process_job(self,job):
        job.subfile = job.name+'.'+self.sub_launcher+'.in'
        no_cores = job.cores==None
        no_nodes = job.nodes==None
        if no_cores and no_nodes:
            self.error('job did not specify cores or nodes\n  At least one must be provided')
        elif no_cores:
            job.cores = self.cores_per_node*job.nodes
        elif no_nodes:
            job.nodes = int(ceil(float(job.cores)/self.cores_per_node))
            if abs(job.nodes-1-float(job.cores)/self.cores_per_node)<1e-6:
                job.nodes-=1
            #end if
        else:
            job.cores = min(job.cores,job.nodes*self.cores_per_node)
        #end if
        job.processes = max(1,int(float(job.cores)/job.threads))
        job.tot_cores = job.nodes*self.cores_per_node
        job.procs = job.nodes*self.procs_per_node

        if mod(job.processes,job.nodes)!=0:
            job.processes_per_node = None
        else:
            job.processes_per_node = job.processes/job.nodes
        #end if

        if mod(job.processes,job.procs)!=0:
            job.processes_per_proc = None
        else:
            job.processes_per_proc = job.processes/job.procs
        #end if

        if job.ppn is None:
            job.ppn = self.cores_per_node
        #end if

        if job.account is None and self.requires_account:
            if self.account is None:
                self.error('account not specified for job on '+self.name)
            #end if
            job.account = self.account
        #end if

        self.process_job_extra(job)

        job.set_environment(OMP_NUM_THREADS=job.threads)

        launcher = self.app_launcher
        if launcher=='mpirun':
            job.run_options.add(np='-np '+str(job.processes))
        elif launcher=='aprun':
            if 'n' in self.aprun_options:
                job.run_options.add(n='-n '+str(job.processes))
            #end if
            if 'd' in self.aprun_options and job.threads>1:
                job.run_options.add('-d '+str(job.threads))
            #end if
            if 'N' in self.aprun_options and job.processes_per_node!=None:
                job.run_options.add(N='-N '+str(job.processes_per_node))
            #end if
            if 'S' in self.aprun_options and job.processes_per_proc!=None:
                job.run_options.add(S='-S '+str(job.processes_per_proc))
            #end if
        elif launcher=='runjob':
            #bypass setup_environment
            envs='--envs'
            for name,value in job.env.iteritems():
                envs+=' {0}={1}'.format(name,value)
            #end for
            job.env = None
            job.run_options.add(
                np      = '--np '+str(job.processes),
                p       = '-p '+str(job.processes_per_node),
                locargs = '$LOCARGS',
                verbose = '--verbose=INFO',
                envs    = envs
                )
        elif launcher=='srun':  # Amos contribution from Ryan McAvoy
            None # anything needed here?
        elif launcher=='ibrun': # Lonestar contribution from Paul Young
	    job.run_options.add(
	    	np	= '-n '+str(job.processes),
	    	p	= '-o '+str(0),
	    	)
        else:
            self.error(launcher+' is not yet implemented as an application launcher')
        #end if
    #end def process_job


    def process_job_extra(self,job):
        None
    #end def process_job_extra


    def query_queue(self):

#        print
#        print 'in query queue'

        self.system_queue.clear()
        if self.queue_querier=='qstat':
            out,err = Popen('qstat -a',shell=True,stdout=PIPE,stderr=PIPE,close_fds=True).communicate()            
            lines = out.splitlines()
            for line in lines:

#                print line

                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[0]
                    else:
                        spid = tokens[0]
                    #endif
                    if spid.isdigit() and len(tokens)==11:
                        pid = int(spid)
                        jid,uname,queue,jname,sessid,nodes,tasks,mem,rtime,status,etime = tokens
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                        #end if
                    #end if
                #end if
            #end for
        elif self.queue_querier=='qstata':
            out,err = Popen('qstat',shell=True,stdout=PIPE,stderr=PIPE,close_fds=True).communicate()            
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[0]
                    else:
                        spid = tokens[0]
                    #endif
                    if spid.isdigit() and len(tokens)==6:
                        pid = int(spid)
                        jid,uname,wtime,nodes,status,loc = tokens
                        self.system_queue[pid] = status
                    #end if
                #end if
            #end for
        elif self.queue_querier=='squeue': # contributed by Ryan McAvoy
            out,err = Popen('squeue',shell=True,stdout=PIPE,stderr=PIPE,close_fds=True).communicate()
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[0]
                    else:
                        spid = tokens[0]
                    #endif
                    if spid.isdigit() and len(tokens)==8:
                        pid = int(spid)
                        jid,loc,name,uname,status,wtime,nodes,reason = tokens
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                        #end if
                    #end if
                #end if
            #end for
        elif self.queue_querier=='sacct': # contributed by Ryan McAvoy
            out,err = Popen('sacct',shell=True,stdout=PIPE,stderr=PIPE,close_fds=True).communicate()
            lines = out.splitlines()
            for line in lines:
                tokens=line.split()
                if len(tokens)>0:
                    if '.' in tokens[0]:
                        spid = tokens[0].split('.')[0]
                    else:
                        spid = tokens[0]
                    #endif
                    #print tokens
                    #print '----------'
                    if spid.isdigit() and len(tokens)==6:  #if account is empty, only 6 tokens.

                        pid = int(spid)
                        jid,name,loc,cores,status,exit_code = tokens
                        status = status.split('+')[0]  ## get rid of '+' in the end
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                    elif spid.isdigit() and len(tokens)==7:

                        pid = int(spid)
                        jid,name,loc,uname,cores,status,exit_code = tokens
                        status = status.split('+')[0]  ## get rid of '+' in the end
                        if status in self.job_states:
                            self.system_queue[pid] = self.job_states[status]
                        else:
                            self.error('job state '+status+' is unrecognized')
                        #end if
                    #end if
                #end if
            #end for
        else:
            self.error('ability to query queue with '+self.queue_querier+' has not yet been implemented')
        #end if

#        print
#        print 'system_queue'
#        print self.system_queue

        done = []
        for pid,process in self.processes.iteritems():
            if not pid in self.system_queue or self.system_queue[pid]=='complete' or self.generate_only:
                job = process.job
                job.status = job.states.finished
                job.finished = True
                iid = job.internal_id
                self.running.remove(iid)
                self.finished.add(iid)
                done.append(pid)
            #end if
        #end for
        for pid in done:
            del self.processes[pid]
        #end for

#        print
#        print 'processes'
#        print self.processes
    #end def query_queue


    def submit_jobs(self):
        nprocesses_running = len(self.processes)
        queue_slots_available = self.queue_size-nprocesses_running
        remove = []
        for iid in self.waiting:
            if queue_slots_available>0:
                remove.append(iid)
                self.running.add(iid)
                job = self.jobs[iid]
                self.submit_job(job)
                queue_slots_available -= 1
            else:
                break
            #end if
        #end for
        for iid in remove:
            self.waiting.remove(iid)
        #end for
    #end def submit_jobs

        
    def submit_job(self,job):
        pad = self.enter(job.directory,msg=job.internal_id)
        if job.subfile==None:
            self.error('submission file not specified for job')
        elif not os.path.exists(job.subfile):
            self.error('job submission file was not written prior to submission\n  submission file: '+os.path.join(job.directory,job.subfile))
        #end if
        command = self.sub_launcher+job.sub_options.write()+' '+job.subfile
        if self.generate_only:
            self.log(pad+'Would have executed:  '+command)
            job.status = job.states.running
            process = obj()
            process.job = job
            self.processes[job.internal_id] = process
        else:
            self.log(pad+'Executing:  '+command)
            job.status = job.states.running
            process = obj()
            process.job = job
            out,err = Popen(command,shell=True,stdout=PIPE,stderr=PIPE,close_fds=True).communicate()
            output=out+'\n'+err
            pid = self.read_process_id(output)
            if pid is None:
                self.error('process id could not be determined from submission output\n  output:\n'+output)
            else:
                self.log(pad+'  pid: {0}'.format(pid))
            #end if
            #pid = 'fakepid_'+str(job.internal_id)
            job.system_id = pid
            self.processes[pid] = process
        #end if
        self.leave()
    #end def submit_job
            

    def remove_job(self,job):
        if self.job_remover=='qdel':
            command = 'qdel '+str(job.system_id)
        elif self.job_remover=='scancel':
            command = 'scancel '+str(job.system_id)
        else:
            self.error('ability to remove job using '+self.job_remover+' has not yet been implemented')
        #endif
        os.system(command)
    #end def remove_job


    def setup_environment(self,job):
        env = ''
        if job.env!=None:
            for name,val in job.env.iteritems():
                env +='export {0}={1}\n'.format(name,val)
            #end for
        #end if
        return env
    #end def setup_environment


    def write_job(self,job):
        job.subfile = job.name+'.'+self.sub_launcher+'.in'
        env = self.setup_environment(job)
        command = job.run_command(self.app_launcher)
        
        c = self.write_job_header(job)+'\n'
        if len(job.presub)>0:
            c+=job.presub+'\n'
        #end if
        c+=env
        c+=command+'\n'
        if len(job.postsub)>0:
            c+=job.postsub+'\n'
        #end if

        filepath = os.path.join(job.directory,job.subfile)
        fobj = open(filepath,'w')
        fobj.write(c)
        fobj.close()
        if self.executable_subfile:
            os.system('chmod +x '+filepath)
        #end if
    #end def write_job


    def write_job_header(self,job):
        self.not_implemented()
    #end def write_job_header


    def read_process_id(self,output):
        pid = None
        lines = output.splitlines()
        for line in lines:
            ls = line.strip()
            if ls.isdigit():
                pid = int(ls)
                break
            elif '.' in line:
                spid = line.split('.')[0]
                if spid.isdigit():
                    pid = int(spid)
                    break
                #end if
            elif ' ' in line: # specialized for Amos?
                spid = line.split(' ')[-1]
                if spid.isdigit():
                    pid = int(spid)
                    break
                #end if    
            #end if
        #end for
        return pid
    #end def read_process_id

#end class Supercomputer




class Kraken(Supercomputer):

    name = 'kraken'

    requires_account = True

    def write_job_header(self,job):
        c='#!/bin/bash\n'
        c+='#PBS -A '+str(job.account)+'\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l size='+str(job.tot_cores)+'\n'
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
cd $PBS_O_WORKDIR
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_MAX_SHORT_MSG_SIZE=1024
export MPICH_PTL_UNEX_EVENTS=800000
export MPICH_UNEX_BUFFER_SIZE=16M
export MPI_MSGS_PER_PROC=32768
'''
        return c
    #end def write_job_header


    def read_process_id(self,output):
        self.not_implemented()
    #end def read_process_id
#end class Kraken




class Jaguar(Supercomputer):
    name = 'jaguar'

    requires_account = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue='batch'
        #end if
        c='#!/bin/bash\n'
        c+='#PBS -A '+str(job.account)+'\n'
        c+='#PBS -q '+job.queue+'\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l size='+str(job.tot_cores)+'\n'
        c+='#PBS -l gres=widow2%widow3\n'
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
export MPICH_PTL_SEND_CREDITS=-1
export MPICH_MAX_SHORT_MSG_SIZE=1024
export MPICH_PTL_UNEX_EVENTS=800000
export MPICH_UNEX_BUFFER_SIZE=16M
export MPI_MSGS_PER_PROC=32768
'''
        return c
    #end def write_job_header


    def read_process_id(self,output):
        self.not_implemented()
    #end def read_process_id
#end class Jaguar




class Taub(Supercomputer):

    name = 'taub'
    
    def write_job_header(self,job):
        if job.queue is None:
            job.queue='cse'
        #end if
        c=''
        c+='#PBS -q '+job.queue+'\n'
        c+='#PBS -N '+job.name+'\n'
        c+='#PBS -l nodes={0}:ppn={1}\n'.format(job.nodes,job.ppn)
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        c+='#PBS -o '+job.outfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
cd ${PBS_O_WORKDIR}

export MPICH_PTL_SEND_CREDITS=-1
export MPICH_MAX_SHORT_MSG_SIZE=1024
export MPICH_PTL_UNEX_EVENTS=800000
export MPICH_UNEX_BUFFER_SIZE=16M
export MPI_MSGS_PER_PROC=32768

module load mvapich2/1.6-intel
'''
        return c
    #end def write_job_header


    def read_process_id(self,output):
        lines = output.splitlines()
        found_pid = False
        for line in lines:
            tokens = line.split()
            for token in tokens:
                t = token.strip().split('.',2)
                if len(t)==3 and t[0].isdigit() and t[1].startswith('taub') and t[2]=='campuscluster.illinois.edu':
                    pid = int(t[0])
                    found_pid = True
                    return pid
                #end if
            #end for
        #end for
        if not found_pid:
            return None
        #end if
    #end def read_process_id
#end class Taub




class OIC5(Supercomputer):

    name = 'oic5'
    batch_capable = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'mstqmc13q'
        #end if

        ppn = job.processes_per_node
        #ppn = 32/job.threads
        #if ppn*job.threads!=32:
        #    self.error('ppn is not being set properly for OIC5\n  perhaps the number of threads requested does not evenly divide the 32 cores\n  you requested {0} threads'.format(job.threads))
        ##end if

        c='#!/bin/bash\n'
        c+='#PBS -q '+job.queue+'\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l nodes={0}:ppn={1}\n'.format(job.nodes,ppn)
        c+='#PBS -W x=\"NACCESSPOLICY:SINGLEJOB\"\n'
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header


    def read_process_id(self,output):
        pid = None
        lines = output.splitlines()
        for line in lines:
            if 'oic.ornl.gov' in line:
                spid = line.split('.')[0]
                if spid.isdigit():
                    pid = int(spid)
                #end if
            #end if
        #end for
        return pid
    #end def read_process_id
#end class OIC5




class NerscMachine(Supercomputer):
    batch_capable = True

    def process_job_extra(self,job):
        if job.threads>1:
            if job.options is None:
                if mod(8,job.threads)==0:
                    job.options = '-ss'
                #end if
            #end if
            if job.compiler!=None and job.compiler=='intel':
                if job.options is None:
                    job.options = '-cc numa_node'
                else:
                    job.options = '-cc numa_node '+job.options
                #end if
            #end if
        #end if
    #end def process_job_extra

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'regular'
        #end if

        c='#!/bin/bash\n'
        c+='#PBS -q '+job.queue+'\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l mppwidth={0}\n'.format(job.tot_cores)
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class NerscMachine


class Edison(NerscMachine):
    name = 'edison'
#end class Edison

class Hopper(NerscMachine):
    name = 'hopper'
#end class Hopper



class BlueWatersXK(Supercomputer):
    
    name = 'bluewaters_xk'
    requires_account = False
    batch_capable    = True

    def write_job_header(self,job):
        c='#!/bin/bash\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l nodes={0}:ppn={1}:xk\n'.format(job.nodes,job.ppn)
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class BlueWatersXK




class BlueWatersXE(Supercomputer):
    
    name = 'bluewaters_xe'
    requires_account = False
    batch_capable    = True

    def write_job_header(self,job):
        c='#!/bin/bash\n'
        c+='#PBS -N '+str(job.name)+'\n'
        c+='#PBS -l walltime='+job.pbs_walltime()+'\n'
        c+='#PBS -l nodes={0}:ppn={1}:xe\n'.format(job.nodes,job.ppn)
        c+='#PBS -o '+job.outfile+'\n'
        c+='#PBS -e '+job.errfile+'\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class BlueWatersXE




class Titan(Supercomputer):
    
    name = 'titan'
    requires_account = True
    batch_capable    = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'batch'
        #end if
        c= '#!/bin/bash\n'
        c+='#PBS -A {0}\n'.format(job.account)
        c+='#PBS -q {0}\n'.format(job.queue)
        c+='#PBS -N {0}\n'.format(job.name)
        c+='#PBS -o {0}\n'.format(job.outfile)
        c+='#PBS -e {0}\n'.format(job.errfile)
        c+='#PBS -l walltime={0}\n'.format(job.pbs_walltime())
        c+='#PBS -l nodes={0}\n'.format(job.nodes)
        #c+='#PBS -l gres=widow3\n'
        c+='#PBS -l gres=atlas1\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class Titan



class EOS(Supercomputer):
    
    name = 'eos'
    requires_account = True
    batch_capable    = True

    def process_job_extra(self,job):
        if job.threads>1:
            if job.threads<=8: 
                job.run_options.add(ss='-ss')
            #end if
            job.run_options.add(cc='-cc numa_node')
        #end if
    #end def process_job_extra


    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'batch'
        #end if
        c= '#!/bin/bash\n'
        c+='#PBS -A {0}\n'.format(job.account)
        c+='#PBS -q {0}\n'.format(job.queue)
        c+='#PBS -N {0}\n'.format(job.name)
        c+='#PBS -o {0}\n'.format(job.outfile)
        c+='#PBS -e {0}\n'.format(job.errfile)
        c+='#PBS -l walltime={0}\n'.format(job.pbs_walltime())
        c+='#PBS -l nodes={0}\n'.format(job.nodes)
        c+='#PBS -l gres=atlas1\n'
        if job.user_env:
            c+='#PBS -V\n'
        #end if
        c+='''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        return c
    #end def write_job_header
#end class EOS



class ALCF_Machine(Supercomputer):    
    requires_account   = True
    batch_capable      = True
    executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    base_partition = None

    def process_job_extra(self,job):
        job.sub_options.add(
            env  = '--env BG_SHAREDMEMSIZE=32',
            mode = '--mode script'
            )
        if job.nodes<self.base_partition:
            self.warn('!!! ATTENTION !!!\n  number of nodes on {0} should not be less than {1}\n  you requested: {2}'.format(self.name,self.base_partition,job.nodes))
        else:
            partition = log(float(job.nodes)/self.base_partition)/log(2.)
            if abs(partition-int(partition))>1e-6:
                self.warn('!!! ATTENTION !!!\n  number of nodes on {0} must be {1} times a power of two\n  you requested: {2}\n  nearby valid node count: {3}'.format(self.name,self.base_partition,job.nodes,self.base_partition*2**int(round(partition))))
            #end if
        #end if
    #end def process_job_extra

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'default'
        #end if
        c= '#!/bin/bash\n'
        c+='#COBALT -q {0}\n'.format(job.queue)
        c+='#COBALT -A {0}\n'.format(job.account)
        c+='#COBALT -n {0}\n'.format(job.nodes)
        c+='#COBALT -t {0}\n'.format(job.total_minutes())
        c+='#COBALT -O {0}\n'.format(job.identifier)
        c+='\nLOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"\n'
        c+='echo "Cobalt location args: $LOCARGS" >&2\n\n'
        return c
    #end def write_job_header
#end class ALCF_Machine


class Vesta(ALCF_Machine):
    name = 'vesta'
    base_partition = 32
#end class Vesta

class Cetus(ALCF_Machine):
    name = 'cetus'
    base_partition = 128
#end class Cetus

class Mira(ALCF_Machine):
    name = 'mira'
    base_partition = 512
#end class Mira



class Lonestar(Supercomputer):  # Lonestar contribution from Paul Young

    name = 'lonestar' # will be converted to lowercase anyway
    requires_account = False
    batch_capable    = True

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'batch'
        #end if
        c= '#!/bin/bash\n'
        #c+='#$ -A {0}\n'.format(job.account)
        c+='#$ -q {0}\n'.format(job.queue)
        c+='#$ -N {0}\n'.format(job.name)
        c+='#$ -o {0}\n'.format(job.outfile)
        c+='#$ -e {0}\n'.format(job.errfile)
        c+='#$ -l h_rt={0}\n'.format(job.pbs_walltime())
        c+='#$ -pe 12way {0}\n'.format(job.nodes*12)
	c+='#$ -cwd\n'
        if job.user_env:
            c+='#$ -V\n'
        #end if
        return c
    #end def write_job_header


    def read_process_id(self,output):
        pid = None
        lines = output.splitlines()
        
        for line in lines:
            if 'Your job' in line:
                spid = line.split(' ')[2]
                if spid.isdigit():
                    pid = int(spid)
                #end if
            #end if
        #end for
        return pid
    #end def read_process_id
#end class Lonestar


class ICMP_Machine(Supercomputer): # ICMP and Amos contributions from Ryan McAvoy
    batch_capable      = True
    executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'defq'
        #end if
        c= '#!/bin/bash -x\n'
        c+='#SBATCH --export=ALL\n'
        c+='#SBATCH -J {0}\n'.format(job.identifier)
        c+='#SBATCH -p {0}\n'.format(job.queue)
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        c+='#SBATCH --nodes {0}\n'.format(job.nodes)
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        return c
    #end def write_job_header
#end class ICMP_Machine


class Komodo(ICMP_Machine):
    name = 'komodo'
#end class Komodo

class Matisse(ICMP_Machine):
    name = 'matisse'
#end class Matisse


class Amos(Supercomputer):
    name = 'amos'

    #requires_account   = True
    batch_capable      = True
    executable_subfile = True

    prefixed_output    = True
    outfile_extension  = '.output'
    errfile_extension  = '.error'

    def write_job_header(self,job):
        if job.queue is None:
            job.queue = 'debug'
        #end if
        if job.queue == 'debug':
            base_partition = 1
            max_partition = 32
            max_time =1
        elif job.queue == 'small':
            base_partition = 1
            max_partition = 64
            max_time =24
        elif job.queue == 'medium':
            base_partition = 128
            max_partition = 512
            max_time =12
        elif job.queue == 'large':
            base_partition = 1024
            max_partition = 2048
            max_time =6
        elif job.queue == 'verylarge':
            base_partition = 3072
            max_partition = 4096
            max_time =6
        #end if
        job.total_hours = job.days*24 + job.hours + job.minutes/60.0 + job.seconds/3600.0
        if job.total_hours > max_time:
            self.warn('!!! ATTENTION !!!\n  the maximum runtime on {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_time,job.total_hours))
            job.hours   = max_time
            job.minutes =0
            job.seconds =0
        #end if
        if job.nodes<base_partition:
            self.warn('!!! ATTENTION !!!\n  number of nodes in {0} should not be less than {1}\n  you requested: {2}'.format(job.queue,self.base_partition,job.nodes))
        elif job.nodes>max_partition:
            self.warn('!!! ATTENTION !!!\n  number of nodes in {0} should not be more than {1}\n  you requested: {2}'.format(job.queue,max_partition,job.nodes))
        else:
            if job.queue != 'verylarge':
                partition = log(float(job.nodes)/base_partition)/log(2.)
                if abs(partition-int(partition))>1e-6:
                    self.warn('!!! ATTENTION !!!\n  number of nodes on {0} must be {1} times a power of two\n  you requested: {2}\n  nearby valid node count: {3}'.format(self.name,self.base_partition,job.nodes,self.base_partition*2**int(round(partition))))
            elif job.nodes != 3072 and job.nodes != 4096:
                self.warn('!!! ATTENTION !!!\n  number of nodes on {0} must be 3072 or 4096 you requested {1}'.format(self.name,job.nodes))
            #end if
        #end if

        c= '#!/bin/bash -x\n'
        c+='#SBATCH --export=ALL\n'
        #c+=#SBATCH -D /gpfs/sb/data/<project>/<user>/
        c+='#SBATCH -J {0}\n'.format(job.identifier)
        c+='#SBATCH -p {0}\n'.format(job.queue)
        c+='#SBATCH -o {0}\n'.format(job.outfile)
        c+='#SBATCH -e {0}\n'.format(job.errfile)
        c+='#SBATCH --nodes {0}\n'.format(job.nodes)
        c+='#SBATCH --ntasks-per-node={0}\n'.format(job.processes_per_node)
        c+='#SBATCH --cpus-per-task={0}\n'.format(job.threads)
        c+='#SBATCH -t {0}:{1}:{2}\n'.format(str(job.hours+24*job.days).zfill(2),str(job.minutes).zfill(2),str(job.seconds).zfill(2))
        # c+='#SBATCH --mail-type=ALL'
        # c+='#SBATCH --mail-user=<{0}>'

        return c
    #end def write_job_header
#end class Amos



#Known machines
#  workstations
for cores in range(1,128+1):
    Workstation('ws'+str(cores),cores,'mpirun'),
    Workstation('node'+str(cores),cores,'mpirun'),
#end for
#  supercomputers and clusters
#            nodes sockets cores ram qslots  qlaunch  qsubmit     qstatus   qdelete
Jaguar(      18688,   2,     8,   32,  100,  'aprun',   'qsub',   'qstat',    'qdel')
Kraken(       9408,   2,     6,   16,  100,  'aprun',   'qsub',   'qstat',    'qdel')
Taub(          400,   2,     6,   24,   50, 'mpirun',   'qsub',   'qstat',    'qdel')
OIC5(           28,   2,    16,  128, 1000, 'mpirun',   'qsub',   'qstat',    'qdel')
Hopper(       6384,   2,    12,   64, 1000,  'aprun',   'qsub',   'qstat',    'qdel')
Edison(        664,   2,     8,   64,  100,  'aprun',   'qsub',   'qstat',    'qdel')
BlueWatersXK( 3072,   1,    16,   32,  100,  'aprun',   'qsub',   'qstat',    'qdel')
BlueWatersXE(22640,   2,    16,   64,  100,  'aprun',   'qsub',   'qstat',    'qdel')
Titan(       18688,   1,    16,   32,  100,  'aprun',   'qsub',   'qstat',    'qdel')
EOS(           744,   2,     8,   64, 1000,  'aprun',   'qsub',   'qstat',    'qdel')
Vesta(        2048,   1,    16,   16,   10, 'runjob',   'qsub',  'qstata',    'qdel')
Cetus(        1024,   1,    16,   16,   10, 'runjob',   'qsub',  'qstata',    'qdel')
Mira(        49152,   1,    16,   16,   10, 'runjob',   'qsub',  'qstata',    'qdel')
Lonestar(    22656,   2,     6,   12,  128,  'ibrun',   'qsub',   'qstat',    'qdel')
Matisse(        20,   2,     8,   64,    2, 'mpirun', 'sbatch',   'sacct', 'scancel')
Komodo(         24,   2,     6,   48,    2, 'mpirun', 'sbatch',   'sacct', 'scancel')
Amos(         5120,   1,    16,   16,  128,   'srun', 'sbatch',   'sacct', 'scancel')

#machine accessor functions
get_machine_name = Machine.get_hostname
get_machine      = Machine.get


