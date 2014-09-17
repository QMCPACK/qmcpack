
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

    preserve = Simulation.preserve & set(['sims'])


    def __init__(self,*sims,**kwargs):
        if len(sims)==1 and isinstance(sims[0],list):
            sims = sims[0]
        #end if
        if len(sims)==0:
            self.error('attempted to bundle 0 simulations\n  at least one simulation must be provided to bundle')
        #end if
        for sim in sims:
            if not isinstance(sim,Simulation):
                self.error('attempted to bundle non-simulation object: '+sim.__class__.__name__)
            #end if
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
        #sims in bundle should not submit jobs
        for sim in sims:
            sim.skip_submit = True
        #end for
    #end def __init__

            
    #def init_job(self):
    #    None # this is to override the default behavior of Simulation
    ##end def init_job


    def bundle_dependencies(self):
        deps = []
        for sim in self.sims:
            for d in sim.dependencies:
                deps.append((d.sim,'other'))
            #end for
        #end for
        self.depends(*deps)
    #end def bundle_dependencies


    def bundle_jobs(self,relative=False):
        jobs = []
        job0 = self.sims[0].job
        time    = Job.zero_time()
        nodes   = 0
        cores   = 0
        threads = job0.threads
        queue   = job0.queue
        same_threads = True
        same_queue   = True
        machine_names = set()
        for sim in self.sims:
            job = sim.job
            nodes += job.nodes
            cores += job.cores
            same_threads = same_threads and threads==job.threads
            same_queue   = same_queue   and queue==job.queue
            time    = job.max_time(time)
            machine = job.get_machine()
            machine_names.add(machine.name)
            jobs.append(job)
        #end for
        if not same_threads:
            self.error('bundling jobs with different numbers of threads is not yet supported',trace=False)
        #end if
        if not same_queue:
            self.error('bundling jobs with different queues is not allowed',trace=False)
        #end if
        machine_names = list(machine_names)
        if len(machine_names)!=1:
            self.error('attempted to bundle jobs across these machines: '+str(machine_names)+'\n  jobs may only be bundled on the same machine',trace=False)
        #end if
        self.job = Job(
            bundled_jobs = jobs,
            relative     = relative,
            queue        = queue,
            nodes        = nodes,
            cores        = cores,
            threads      = threads,
            machine      = machine_names[0],
            **time
            )
    #end def bundle_jobs

        
    def pre_write_inputs(self,save_image):
        for sim in self.sims:
            if not sim.setup:
                sim.write_inputs(save_image)
            #end if
        #end for
    #end def pre_write_inputs


    def pre_send_files(self,enter):
        for sim in self.sims:
            if not sim.sent_files:
                sim.send_files(enter)
            #end if
        #end for
    #end def pre_send_files


    def post_submit(self):
        for sim in self.sims:
            sim.submitted = True
        #end for
    #end def post_submit


    def pre_check_status(self):
        if self.job.finished:
            for sim in self.sims:
                sim.job.finished = True
            #end for
        #end if
        for sim in self.sims:
            sim.check_status()
        #end for
    #end def pre_check_status


    def check_sim_status(self):
        finished = True
        for sim in self.sims:
            finished = finished and sim.finished
        #end for
        self.finished = finished
    #end def check_sim_status


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
