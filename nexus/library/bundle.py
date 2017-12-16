##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  bundle.py                                                         #
#    Enables bundled simulations.  Groups of bundled simulations     #
#    share a single batch submission file on supercomputers and      #
#    will execute concurrently.                                      #
#                                                                    #
#  Content summary:                                                  #
#    SimulationBundle                                                #
#      Simulation class to represent a bundle of simulations.        #
#      Contains a list of simulations and inherits dependencies      #
#        from all component simulations.                             #
#      Merges individual Job information (e.g. nodes requested) into #
#        a bundled job.                                              #
#                                                                    #
#    bundle                                                          #
#      Function to bundle together many simulations.                 #
#      Returns a SimulationBundle object.                            #
#      Syntax:                                                       #
#        sim_bundle = bundle(sim1,sim2,sim3,...)                     #
#                                                                    #                                        
#====================================================================#


from machines import Workstation,Job
from simulation import Simulation,NullSimulationInput,NullSimulationAnalyzer

class SimulationBundleInput(NullSimulationInput):
    None
#end class SimulationBundleInput

class SimulationBundleAnalyzer(NullSimulationAnalyzer):
    None
#end class SimulationBundleAnalyzer



class SimulationBundle(Simulation):

    input_type = SimulationBundleInput
    analyzer_type = SimulationBundleAnalyzer
    generic_identifier = 'bundle'
    image_directory    = 'bundle'

    is_bundle = True

    def __init__(self,*sims,**kwargs):
        if len(sims)==1 and isinstance(sims[0],list):
            sims = sims[0]
        #end if
        sims = list(sims) # make a copy
        if len(sims)==0:
            self.error('attempted to bundle 0 simulations\n  at least one simulation must be provided to bundle')
        #end if
        for sim in sims:
            if not isinstance(sim,Simulation):
                self.error('attempted to bundle non-simulation object: '+sim.__class__.__name__)
            #end if
            sim.bundled = True
            sim.bundler = self
        #end for
        relative_paths = False
        if 'relative' in kwargs:
            relative_paths = kwargs['relative']
            del kwargs['relative']
        #end if

        self.sims = sims
        self.bundle_jobs(relative=relative_paths)
        self.system = None

        if not 'path' in kwargs:
            kwargs['path'] = self.sims[0].path
        #end if
        if not 'job' in kwargs:
            kwargs['job'] = self.job
        #end if
        Simulation.__init__(self,**kwargs)
        self.infile = None
        if isinstance(self.job.get_machine(),Workstation):
            self.outfile = None
            self.errfile = None
        #end if
        self.bundle_dependencies()

        # constituent sims do not actually submit to the job queue
        for sim in self.sims:
            sim.skip_submit = True
        #end for

        # restrict bundle from getting ahead of any constituent sims
        self.allow_create_directories = False
        self.allow_get_dependencies   = False
        self.allow_write_inputs       = False
        self.allow_send_files         = False
        self.allow_submit             = False
        self.allow_get_output         = False
        self.allow_analyze            = False
    #end def __init__

            
    #def init_job(self):
    #    None # this is to override the default behavior of Simulation
    ##end def init_job


    def bundle_dependencies(self):
        deps = []
        sim_ids = set()
        depsim_ids = set()
        for sim in self.sims:
            sim_ids.add(sim.simid)
            depsim_ids |= sim.downstream_simids()
            for d in sim.dependencies:
                deps.append((d.sim,'other'))
            #end for
        #end for
        # guard against dependencies within the bundle
        internal_deps = sim_ids & depsim_ids
        if len(internal_deps)>0:
            sims = dict()
            for sim in self.sims:
                sims[sim.simid]=sim
            #end for
            msg = 'attempted to bundle simulations that depend on each other\nsimulations can only be bundled if they can be executed simultaneously\nbundle identifier, simid, and directory:\n  {0:<8} {1:>4} {2}\n'.format(self.identifier,self.simid,self.locdir)
            msg+='sims in the bundle that can remain:\n'
            for simid in sorted(sim_ids-internal_deps):
                sim = sims[simid]
                msg +='  {0:<8} {1:>4} {2}\n'.format(sim.identifier,sim.simid,sim.locdir)
            #end for
            msg+='sims in the bundle that need to be removed:\n'
            for simid in sorted(internal_deps):
                sim = sims[simid]
                msg +='  {0:<8} {1:>4} {2}\n'.format(sim.identifier,sim.simid,sim.locdir)
            #end for
            msg+='please remove the necessary sims from the bundle and try again\nthe excluded sims can likely be bundled separately'
            self.error(msg,'bundle')
        #end if
        self.depends(*deps)
    #end def bundle_dependencies


    def bundle_jobs(self,relative=False):
        jobs = []
        job0 = self.sims[0].job
        time    = Job.zero_time()
        nodes   = 0
        cores   = 0
        thread_set = set()
        queue_set  = set()
        presub_set = set()
        machine_set = set()
        for sim in self.sims:
            job = sim.job
            nodes += job.nodes
            cores += job.cores
            time    = job.max_time(time)
            machine = job.get_machine()
            machine_set.add(machine.name)
            thread_set.add(job.threads)
            queue_set.add(job.queue)
            presub_set.add(job.presub)
            jobs.append(job)
        #end for
        if len(thread_set)>1:
            self.error('bundling jobs with different numbers of threads is not yet supported\nthread inputs provided: {0}'.format(sorted(thread_set),trace=False))
        #end if
        if len(queue_set)>1:
            self.error('bundling jobs with different queues is not allowed\nqueue inputs provided: {0}'.format(sorted(queue_set)),trace=False)
        #end if
        if len(presub_set)>1:
            ps = ''
            for psub in sorted(presub_set):
                ps+=psub+'\n\n'
            #end for
            self.error('bundling jobs with different pre-submission commands is not allowed\npresub inputs provided: \n{0}'.format(ps),trace=False)
        #end if
        if len(machine_set)>1:
            self.error('attempted to bundle jobs across these machines: {0}\n  jobs may only be bundled on the same machine'.format(sorted(machine_set)),trace=False)
        #end if
        threads = list(thread_set)[0]
        queue   = list(queue_set)[0]
        presub  = list(presub_set)[0]
        machine = list(machine_set)[0]
        self.job = Job(
            bundled_jobs = jobs,
            relative     = relative,
            queue        = queue,
            nodes        = nodes,
            cores        = cores,
            threads      = threads,
            machine      = machine,
            presub       = presub,
            **time
            )
    #end def bundle_jobs


    def completed(self):
        bsims_comp = True
        for sim in self.sims:
            bsims_comp &= sim.completed()
        #end for
        return Simulation.completed(self) & bsims_comp
    #end def completed


    def check_allowed(self,indicator):
        allowed = True
        for sim in self.sims:
            allowed &= sim[indicator]
        #end for
        return allowed
    #end def check_allowed


    def progress(self,dependency_id=None):
        if dependency_id!=None and dependency_id in self.wait_ids:
            self.wait_ids.remove(dependency_id)
        #end if
        if len(self.wait_ids)==0 and not self.block and not self.failed:
            # allow all bundled sims to progress
            for sim in self.sims:
                if len(sim.wait_ids)>0:
                    sim.wait_ids = set()
                #end if
                sim.progress()
            #end for

            # restrict bundle from getting ahead of any constituent sims
            if not self.allow_create_directories:
                self.allow_create_directories = self.check_allowed('created_directories')
            #end if
            if not self.allow_get_dependencies:
                self.allow_get_dependencies = self.check_allowed('got_dependencies')
            #end if
            if not self.allow_write_inputs:
                self.allow_write_inputs = self.check_allowed('setup')
            #end if
            if not self.allow_send_files:
                self.allow_send_files = self.check_allowed('sent_files')
            #end if
            if not self.allow_submit:
                self.allow_submit = self.check_allowed('submitted')
            #end if
            if not self.allow_get_output:
                self.allow_get_output = self.check_allowed('got_output')
            #end if
            if not self.allow_analyze:
                self.allow_analyze = self.check_allowed('analyzed')
            #end if
            
            # progress the bundle itself
            Simulation.progress(self)
        #end if
    #end def progress


    def create_directories(self,*args,**kwargs):
        if self.allow_create_directories:
            Simulation.create_directories(self,*args,**kwargs)
        #end if
    #end def create_directories

    def get_dependencies(self,*args,**kwargs):
        if self.allow_get_dependencies:
            Simulation.get_dependencies(self,*args,**kwargs)
        #end if
    #end def get_dependencies

    def write_inputs(self,*args,**kwargs):
        if self.allow_write_inputs:
            Simulation.write_inputs(self,*args,**kwargs)
        #end if
    #end def write_inputs

    def send_files(self,*args,**kwargs):
        if self.allow_send_files:
            Simulation.send_files(self,*args,**kwargs)
        #end if
    #end def send_files

    def submit(self,*args,**kwargs):
        if self.allow_submit:
            Simulation.submit(self,*args,**kwargs)
            if self.job.finished:
                for sim in self.sims:
                    sim.job.finished = True
                #end for
            #end if
        #end if
    #end def submit

    def check_sim_status(self):
        finished = True
        for sim in self.sims:
            finished = finished and sim.finished
        #end for
        self.finished = finished
    #end def check_sim_status

    def get_output(self,*args,**kwargs):
        if self.allow_get_output:
            Simulation.get_output(self,*args,**kwargs)
        #end if
    #end def get_output

    def analyze(self,*args,**kwargs):
        if self.allow_analyze:
            Simulation.analyze(self,*args,**kwargs)
        #end if
    #end def analyze


    def get_output_files(self):
        return list()
    #end def get_output_files


    def app_command(self):
        return None
    #end def app_command
#end class SimulationBundle



def bundle(*sims,**kwargs):
    return SimulationBundle(*sims,**kwargs)
#end def bundle
