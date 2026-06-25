##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  project_manager.py                                                #
#    Supports the active management of simulation cascades           #
#    (workflows) at the heart of Nexus functionality.                #
#                                                                    #
#  Content summary:                                                  #
#    ProjectManager                                                  #
#      Class actively manages simulation workflows.                  #
#      Works closely with Simulation, Machine, and Job objects.      #
#                                                                    #
#====================================================================#


import time
from . import memory
from .developer import obj, error
from .nexus_base import NexusCore, nexus_core, dynamic_storage
from .simulation import Simulation
from .machines import Machine,Job


def trivial(sim,*args,**kwargs):
    None
#end def trivial


class ProjectManager(NexusCore):

    machine = None

    @staticmethod
    def restore_default_settings():
        ProjectManager.machine = None
    #end def restore_default_settings

    def __init__(self):
        modes = nexus_core.modes
        self.persistent_modes = set([modes.submit,modes.all])
        self.simulations = obj()
        self.cascades = obj()
        self.progressing_cascades = obj()
    #end def __init__


    def add_simulations(self,*simulations):
        if len(simulations)==0:
            self.add_simulations(Simulation.all_sims)
        #end if
        if len(simulations)>0 and not isinstance(simulations[0],Simulation):
            simulations = simulations[0]
        #end if
        for sim in simulations:
            if not sim.fake():
                if len(sim.dependencies)==0:
                    self.add_cascade(sim)
                #end if
                self.simulations[sim.simid]=sim
            #end if
        #end for
    #end def add_simulations


    def add_cascade(self,cascade):
        cid = cascade.simid
        self.cascades[cid]=cascade
        self.progressing_cascades[cid]=cascade
    #end def add_cascade


    def run_project(self,status=False,status_only=False):
        self.log('\nProject starting',n=0)
        self.init_cascades()
        status_only = status_only or nexus_core.status_only
        status = status or status_only or nexus_core.status!=nexus_core.status_modes.none
        if status:
            self.write_simulation_status()
            if status_only:
                NexusCore.write_end_splash()
                return
            #end if
        #end if
        self.log('\nstarting runs:\n'+30*'~',n=1)
        if nexus_core.dependent_modes <= nexus_core.stages_set:
            if nexus_core.monitor:
                start_time = time.time()
                ipoll = 0
                while len(self.progressing_cascades)>0:
                    elapsed_time = time.time() - start_time
                    self.log('elapsed time %.1f s'%elapsed_time,
                             ' memory %3.2f MB'%(memory.resident(children=True)/1e6),
                             n=1,progress=True)
                    NexusCore.wrote_something = False
                    ipoll+=1
                    self.machine.query_queue()
                    self.progress_cascades()
                    self.machine.submit_jobs()
                    self.update_process_ids()
                    time.sleep(nexus_core.sleep)
                    if NexusCore.wrote_something:
                        self.log()
                    #end if
                #end while
            elif len(self.progressing_cascades)>0:
                self.machine.query_queue()
                self.progress_cascades()
                self.machine.submit_jobs()
                self.update_process_ids()
            #end if
        else:
            self.progress_cascades()
        #end if
        self.log('Project finished\n')
        NexusCore.write_end_splash()
    #end def run_project


    def init_cascades(self):
        self.screen_fake_sims()
        self.resolve_file_collisions()
        self.propagate_blockages()
        self.log('loading cascade images',n=1)
        if nexus_core.load_images:
            self.load_cascades()
        else:
            self.log('cascades',n=1)
        #end if
        for c in self.progressing_cascades:
            self.log('cascade',c.simid,'checking in',n=2)
        #end for
        self.check_dependencies()
    #end def init_cascades


    def screen_fake_sims(self):
        def collect_fake(sim,fake):
            if sim.fake():
                fake.append(sim)
            #end if
        #end def collect_fake
        fake = []
        self.traverse_cascades(collect_fake,fake)
        if len(fake)>0:
            msg = 'fake/temporary simulation objects detected in cascade\nthis is a developer error\nlist of fake sims and directories:\n'
            for sim in fake:
                msg +='  {0:>8}  {1}\n'.format(sim.simid,sim.locdir)
            #end for
            self.error(msg)
        #end if
    #end def screen_fake_sims


    def resolve_file_collisions(self):
        self.log('checking for file collisions',n=1)
        entry_order = obj()
        def set_entry_order(sim,entry_order):
            locdir = sim.locdir
            if locdir not in entry_order:
                entry_order[locdir] = [sim]
            else:
                entry_order[locdir].append(sim)
            #end if
        #end def set_entry_order
        self.traverse_cascades(set_entry_order,entry_order)        
        any_collisions = False
        collpath = ''
        for path,simlist in entry_order.items():
            if len(simlist)>1:
                #raise an error if any in/out/err files will collide
                filespace = dict()
                for sim in simlist:
                    if not sim.allow_overlapping_files:
                        files = sim.list('infile','outfile','errfile')
                        for f in files:
                            if f not in filespace:
                                filespace[f] = [sim]
                            else:
                                filespace[f].append(sim)
                            #end if
                        #end for
                    #end if
                #end for
                for f,sims in filespace.items():
                    if len(sims)>1 and f is not None:
                        any_collisions = True
                        msg = 'collision: file '+f+' is overwritten by '
                        for sim in sims:
                            msg +=str(sim.identifier)+' '+str(sim.simid)+','
                        #end for
                        self.log(msg[:-1],n=2)
                        collpath = path
                    #end if
                #end for
            #end if
        #end for
        if any_collisions:
            self.error('file collisions found in directory\n  '+path+'\n  set a unique identifier for each simulation')
        #end if
    #end def resolve_file_collisions


    def propagate_blockages(self):
        def collect_blocked(sim,blocked):
            if sim.block or sim.block_subcascade:
                blocked.append(sim)
            #end if
        #end def collect_blocks
        blocked=[]
        self.traverse_cascades(collect_blocked,blocked)
        for sim in blocked:
            sim.block_dependents(block_self=False)
        #end for
    #end def propagate_blockages

    
    def load_cascades(self):
        cascades = obj()
        progressing_cascades = obj()
        for cid,cascade in self.cascades.items():
            rc = cascade.reconstruct_cascade()
            cascades[rc.simid] = rc 
            progressing_cascades[rc.simid] = rc
        #end for
        self.cascades = cascades
        self.progressing_cascades = progressing_cascades
    #end def load_cascades


    def check_dependencies(self):
        self.log('checking cascade dependencies',n=1)
        result = obj()
        result.dependencies_satisfied = True
        self.traverse_cascades(Simulation.check_dependencies,result)
        if result.dependencies_satisfied:
            self.log('all simulation dependencies satisfied',n=2)
        else:
            self.error('some simulation dependecies are not satisfied')
        #end if
    #end def check_dependencies

                    
    def traverse_cascades(self,operation=trivial,*args,**kwargs):
        for cascade in self.cascades:
            cascade.reset_wait_ids()
        #end for
        for cascade in self.cascades:
            cascade.traverse_cascade(operation,*args,**kwargs)
        #end for
        return
    #end def traverse_cascades


    def write_simulation_status(self):
        status       = nexus_core.status
        status_modes = nexus_core.status_modes
        self.log('\ncascade status',n=1)
        self.log('setup, sent_files, submitted, finished, got_output, analyzed, failed',n=2)
        all_sids = set()
        for sim in self.simulations:
            add = False
            if status==status_modes.active:
                add = sim.active()
            elif status==status_modes.ready:
                add = sim.ready()
            elif status==status_modes.failed:
                add = sim.failed
            else:
                add = True
            #end if
            if add:
                all_sids.add(sim.simid)
            #end if
        #end for
        sids = set()
        for isim in sorted(all_sids):
            sim = self.simulations[isim]
            if not sim.bundled:
                if status==status_modes.active and not sim.active():
                    continue
                #end if
                self.status_line(sim)
                sids.add(sim.simid)
                if sim.is_bundle:
                    for s in sim.sims:
                        self.status_line(s)
                        sids.add(s.simid)
                    #end for
                #end if
            #end if
        #end for
        sids = all_sids-sids
        if len(sids)>0:
            self.log('==== sims missed, part of bundles? ====',n=2)
            for isim in sorted(sids):
                sim = self.simulations[isim]
                self.status_line(sim)
            #end for
        #end if
        self.log('setup, sent_files, submitted, finished, got_output, analyzed, failed',n=2)
    #end def write_simulation_status

        
    def status_line(self,sim,extra=''):
        indicators = ('setup','sent_files','submitted','finished','got_output','analyzed')
        stats = sim.tuple(*indicators)
        status = ''
        for stat in stats:
            status+=str(int(stat))
        #end for
        if sim.process_id is None:
            pid = '------'
        else:
            pid = sim.process_id
        #end if
        sline = '{0}  {1}  {2:<8}  {3:<6}  {4}'.format(status,str(int(sim.failed)),pid,sim.identifier,sim.locdir)
        self.log(sline,extra,n=2)
    #end def status_line


    def progress_cascades(self):
        NexusCore.gc.collect()
        finished = []
        progressing_cascades = self.progressing_cascades
        for cid,cascade in progressing_cascades.items():
            cascade.reset_wait_ids()
        #end for
        for cid,cascade in progressing_cascades.items():
            if not cascade.bundled or cascade.bundler.finished:
                cascade.progress()
            #end if
            cascade.check_subcascade()
            if cascade.subcascade_finished:
                finished.append(cid)
            #end if
        #end for
        for cid in finished:
            del progressing_cascades[cid]
        #end for
    #end def progress_cascades


    def update_process_ids(self):
        for sim in self.simulations:
            sim.update_process_id()
        #end for
    #end def update_process_ids


    # test needed
    def write_sim_dependencies(self,idkey=None):
        for simid in sorted(self.simulations.keys()):
            sim = self.simulations[simid]
            if idkey is None or sim.identifier==idkey:
                self.log('\n{0} {1} {2}'.format(sim.identifier,simid,sim.locdir))
                for did in sorted(sim.dependencies.keys()):
                    dep = sim.dependencies[did]
                    dsim  = dep.sim
                    names = dep.result_names
                    self.log('  {0} {1} {2} {3}'.format(dsim.identifier,dsim.simid,names,dsim.locdir))
                #end for
            #end if
        #end for
    #end def write_sim_dependencies


    # test needed
    def write_cascade_dependents(self):
        self.log('cascade dependents',n=1)
        for cascade in self.cascades:
            cascade.reset_wait_ids()
        #end for
        for cascade in self.cascades:
            self.log(cascade.__class__.__name__+' '+str(cascade.simid),n=2)
            cascade.write_dependents(n=2)
        #end for
        return
    #end def write_cascade_dependents
#end class ProjectManager




class DynamicWorkflowManager(NexusCore):
    '''Replacement for ProjectManager for dynamic workflows'''

    # all dynamic processes encountered
    all_dp         = set() # all dp found via add_new_dyn_procs
    all_sims       = set() # all sims associated w/ dp

    def __init__(self):
        # dynamic processes by status
        self.untouched_dp = obj() # newly minted, perform initial setup
        self.blocked_dp   = obj() # blocked by the user, do nothing
        self.active_dp    = obj() # actively progressing through run stages
        self.finished_dp  = obj() # system or batch process no longer running
        self.succeeded_dp = obj() # may have failed if detection is faulty
        self.failed_dp    = obj() # known to have failed
        self.machine = Machine.get(Job.machine)
        self.add_new_dyn_procs()
    #end def __init__


    def add_new_dyn_procs(self):
        new_dps = dynamic_storage.dynamic_process_ids-self.all_dp
        if len(new_dps)>0:
            for dpid in new_dps:
                dp = dynamic_storage.dynamic_processes[dpid]
                self.all_dp.add(dpid)
                self.untouched_dp[dpid] = dp
                sim = dp.sim
                self.all_sims.add(sim.simid)
                if len(sim.dependencies)>0 or len(sim.dependents)>0:
                    self.error('encountered simulation with explicit dependencies, but these are not allowed in dynamic workflows.\nSimulation id: {}\nSimulation directory: {}'.format(sim.simid,sim.locdir))
        # screen for simulations not associated with dyn process
        all_sims = dynamic_storage.simulation_ids
        rogue_sims = all_sims-self.all_sims
        if len(rogue_sims)>0:
            msg = 'encountered simulations not associated with any dynamic process.\nSimulation ids: {}\nSimulation directories:'
            paths = [all_sims[simid].locdir for simid in rogue_sims]
            for p in sorted(paths):
                msg += '\n'+p
            msg += '\nThis is likely a developer error.\nPlease contact the developers.'
            self.error(msg)
    #end def add_new_dyn_procs


    def poll(self,sleep=None):
        if sleep is None:
            sleep = nexus_core.sleep

        # find and add newly created dynamic process objects
        self.add_new_dyn_procs()

        # update the machine queue for currently running jobs
        self.machine.query_queue()

        # setup new sims
        rem = []
        for dp in self.untouched_dp.values():
            sim = dp.sim
            assert len(sim.dependencies)==0
            assert len(sim.dependents)==0
            if sim.block:
                self.blocked_dp[dp.dpid] = dp
                del self.untouched_dp[dp.dpid]
                continue
            if not sim.created_directories:
                sim.create_directories()
            if not sim.setup:
                sim.write_inputs()
            if not sim.sent_files:
                sim.send_files()
            self.active_dp[dp.dpid] = dp
            rem.append(dp.dpid)
        for dpid in rem:
            del self.untouched_dp[dpid]

        # manage active sims
        rem = []
        for dp in self.active_dp.values():
            #if not dp.requirements_met():
            #    self.error('Requirements not met for simulation execution.\nSimulation type:      {}\nSimulation id:        {}\nSimulation directory: {}\nDynamic process id:   {}\nUnmet requirements:   {}'.format(sim.__class__.__name__,sim.simid,sim.locdir,dp.dpid,dp.unmet_reqs))
            sim = dp.sim
            assert len(sim.dependencies)==0
            assert len(sim.dependents)==0
            if len(dp.unmet_reqs)>0:
                continue
            if not sim.finished:
                sim.submit()
            if sim.finished:
                if not sim.got_output:
                    sim.get_output()
                if not sim.analyzed:
                    sim.analyze()
                self.finished_dp[dp.dpid] = dp
                if sim.failed:
                    self.failed_dp[dp.dpid] = dp
                else:
                    self.succeeded_dp[dp.dpid] = dp
                rem.append(dp.dpid)
        for dpid in rem:
            del self.active_dp[dpid]

        # submit jobs on the machine and store machine pid in sim
        self.machine.submit_jobs()
        for dp in self.active_dp.values():
            dp.sim.update_process_id()

        # wait until next poll
        time.sleep(sleep)
    #end def poll

#end class DynamicWorkflowManager




def workflow_manager(**kw):
    if not hasattr(workflow_manager,'first'):
        workflow_manager.first = True
    else:
        workflow_manager.first = False
    if not nexus_core.dynamic:
        error('workflow_manager is only compatible with dynamic workflows.\nIf you intend to use dynamic workflows, please set dynamic=True in settings.')
    elif not workflow_manager.first:
        error('function "workflow_manager" should only be called once')
    wm = DynamicWorkflowManager(**kw)
    return wm
#end def workflow_manager
