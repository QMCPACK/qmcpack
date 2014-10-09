
import os
import time
import memory
from generic import obj
from project_base import Pobj
from simulation import Simulation


def trivial(sim,*args,**kwargs):
    None
#end def trivial


class ProjectManager(Pobj):
    def __init__(self):
        #variables determined by self
        modes = self.modes
        self.persistent_modes = set([modes.submit,modes.all])
        self.simulations = obj()
        self.cascades = obj()
        self.progressing_cascades = obj()
        self.operations = []
    #end def __init__

    def add_simulations(self,*simulations):
        if len(simulations)>0 and not isinstance(simulations[0],Simulation):
            simulations = simulations[0]
        #end if
        for sim in simulations:
            if len(sim.dependencies)==0:
                self.add_cascade(sim)
            #end if
            self.simulations[sim.simid]=sim
        #end for
    #end def add_simulations

    def add_cascade(self,cascade):
        cid = cascade.simid
        self.cascades[cid]=cascade
        self.progressing_cascades[cid]=cascade
    #end def add_cascade


    def init_cascades(self):
        self.resolve_file_collisions()
        self.propagate_blockages()
        self.log('loading cascade images',n=1)
        if self.load_images:
            self.load_cascades()
        else:
            self.log('cascades',n=1)
        #end if
        for c in self.progressing_cascades:
            self.log('cascade',c.simid,'checking in',n=2)
        #end for
        self.perform_operations()
        #self.write_cascade_status()
        self.check_dependencies()
    #end def init_cascades


    def run_project(self,status=False,status_only=False):
        self.log('\nProject starting',n=0)
        self.init_cascades()
        status_only = status_only or self.status_only
        status = status or status_only
        if status:
            self.write_simulation_status()
            if status_only:
                return
            #end if
        #end if
        self.log('\nstarting runs:\n'+30*'~',n=1)
        if self.dependent_modes <= self.stages_set:
            if self.monitor:
                ipoll = 0
                while len(self.progressing_cascades)>0:
                    self.dlog('\n\n\n'+70*'=',n=1)
                    self.log('poll',ipoll,' memory %3.2f MB'%(memory.resident()/1e6),n=1)
                    Pobj.wrote_something = False
                    self.dlog('cascades',self.progressing_cascades.keys(),n=2)
                    ipoll+=1
                    self.machine.query_queue()
                    self.progress_cascades()
                    self.machine.submit_jobs()
                    self.update_process_ids()
                    self.dlog('sleeping',self.sleep,n=2)
                    time.sleep(self.sleep)
                    self.dlog('awake',n=2)
                    if Pobj.wrote_something:
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
    #end def run_project

    
    def load_cascades(self):
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~load cascades',n=1)
        cascades = obj()
        progressing_cascades = obj()
        for cid,cascade in self.cascades.iteritems():
            self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~reconstruct cascade',n=1)
            rc = cascade.reconstruct_cascade()
            self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end reconstruct cascade',n=1)
            cascades[rc.simid] = rc 
            progressing_cascades[rc.simid] = rc
        #end for
        self.cascades = cascades
        self.progressing_cascades = progressing_cascades
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end load cascades',n=1)
    #end def load_cascades


    def perform_operations(self):
        for op in self.operations:
            operation = op.operation
            sims      = op.sims
            for sim in sims:
                operation(sim)
            #end for
        #end for
        self.operations = []
    #end def perform_operations


    def reset_indicators(self,sims):
        op = obj(
            operation = Simulation.reset_indicators,
            sims = sims
            )
        self.operations.append(op)
    #end def reset_indicators

                    
    def traverse_cascades(self,operation=trivial,*args,**kwargs):
        for cascade in self.cascades:
            cascade.reset_wait_ids()
        #end for
        for cascade in self.cascades:
            cascade.traverse_cascade(operation,*args,**kwargs)
        #end for
        return
    #end def traverse_cascades


    def save_cascades(self):
        def save(sim):
            sim.save_image()
        #end def save
        self.traverse_cascades(save)
    #end def save_cascades


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

            
    def propagate_values(self,**values):
        def set_values(sim,**values):
            for name,value in values.iteritems():
                sim[name] = value
            #end for
        #end def set_values
        self.traverse_cascades(set_values,**values)
    #end def propagate_values

        
    def write_simulation_status(self):
        self.log('cascade status',n=1)
        self.log('setup, sent_files, submitted, finished, got_output, analyzed',n=2)
        indicators = ('setup','sent_files','submitted','finished','got_output','analyzed')
        for isim in self.simulations.keys():
            sim = self.simulations[isim]
            stats = sim.tuple(*indicators)
            status = ''
            for stat in stats:
                status+=str(int(stat))
            #end for
            self.log('{0}  {1}  {2}'.format(status,sim.identifier,sim.locdir),n=2)
        #end for
        self.log('setup, sent_files, submitted, finished, got_output, analyzed',n=2)
    #end def write_simulation_status


    def write_cascade_status(self):
        self.log('cascade status',n=1)

        self.log('setup, sent_files, submitted, finished, got_output, analyzed',n=2)
        def write_status(sim):
            indicators = ('setup','sent_files','submitted','finished','got_output','analyzed')
            stats = sim.tuple(*indicators)
            status = ''
            for stat in stats:
                status+=str(int(stat))
            #end for
            self.log('{0}  {1}  {2}'.format(status,sim.identifier,sim.locdir),n=2)
            #self.log(str(sim.simid)+' '+str(sim.identifier),n=2)
            #self.log('setup      = '+str(sim.setup     ),n=4)
            #self.log('sent_files = '+str(sim.sent_files),n=4)
            #self.log('submitted  = '+str(sim.submitted ),n=4)
            #self.log('finished   = '+str(sim.finished  ),n=4)
            #self.log('got_output = '+str(sim.got_output),n=4)
            #self.log('analyzed   = '+str(sim.analyzed  ),n=4)
        #end def write_status
        self.traverse_cascades(write_status)
        self.log('setup, sent_files, submitted, finished, got_output, analyzed',n=2)
    #end def write_cascade_status


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


    def resolve_file_collisions(self):
        self.log('checking for file collisions',n=1)
        entry_order = obj()
        def set_entry_order(sim,entry_order):
            locdir = sim.locdir
            if not locdir in entry_order:
                entry_order[locdir] = [sim]
            else:
                entry_order[locdir].append(sim)
            #end if
        #end def set_entry_order
        self.traverse_cascades(set_entry_order,entry_order)        
        any_collisions = False
        collpath = ''
        for path,simlist in entry_order.iteritems():
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
                for f,sims in filespace.iteritems():
                    if len(sims)>1 and f!=None:
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

    def progress_cascades(self):
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~progress cascades',n=1)
        self.gc.collect()
        finished = []
        progressing_cascades = self.progressing_cascades
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~reset wait ids',n=1)
        for cid,cascade in progressing_cascades.iteritems():
            cascade.reset_wait_ids()
        #end for
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end reset wait ids',n=1)
        for cid,cascade in progressing_cascades.iteritems():
            self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~progress',n=1)
            cascade.progress()
            self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end progress',n=1)
            self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~check subcascade',n=1)
            cascade.check_subcascade()
            self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end check subcascade',n=1)
            if cascade.subcascade_finished:
                finished.append(cid)
            #end if
        #end for
        for cid in finished:
            self.dlog('removing cascade:',cid,n=1)
            del progressing_cascades[cid]
        #end for
        self.dlog('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end progress cascades',n=1)
    #end def progress_cascades


    def update_process_ids(self):
        for sim in self.simulations:
            sim.update_process_id()
        #end for
    #end def update_process_ids
#end class ProjectManager


