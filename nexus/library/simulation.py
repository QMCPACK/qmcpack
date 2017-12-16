##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  simulation.py                                                     #
#    Provides base classes for simulation objects, including input   #
#    and analysis.  The Simulation base class enables a large amount #
#    of the functionality of Nexus, including workflow construction  #
#    and monitoring, in tandem with the ProjectManager class.        #
#                                                                    #
#  Content summary:                                                  #
#    SimulationInput                                                 #
#      Abstract base class for simulation input.                     #
#                                                                    #
#    SimulationAnalyzer                                              #
#      Abstract base class for simulation data analysis.             #
#                                                                    #
#    Simulation                                                      #
#      Major Nexus class representing a simulation prior to, during, #
#      and after execution.  Checks dependencies between simulations #
#      connected in workflows, manages input file writing,           #
#      participates in job submission, checks for successful         #
#      simulation completion, and analyzes output data.  Saves state #
#      image of simulation progress regularly.  Contains             #
#      SimulationInput, SimulationAnalyzer, and Job objects (also    #
#      optionally contains a PhysicalSystem object).  Derived        #
#      classes tailor specific functions such as passing dependency  #
#      data and checking simulation state to a target simulation     #
#      code.  Derived classes include Qmcpack, Pwscf, Vasp, Gamess,  #
#      Opium, Sqd, Convert4qmc, Wfconvert, Pw2qmcpack, Pw2casino,    #
#      SimulationBundle, and TemplateSimulation.                     #
#                                                                    #
#    NullSimulationInput                                             #
#      Simulation input class intended for codes that do not use an  #
#      input file.                                                   #
#                                                                    #
#    NullSimulationAnalyzer                                          #
#      Simulation input class intended for codes that do not produce #
#      or need to analyze output data.                               #
#                                                                    #
#    SimulationInputTemplate                                         #
#      Supports template input files.  A template input file is a    #
#      standard text input file provided by the user that optionally #
#      has specially marked keywords.  Using find and replace        #
#      operations, Nexus can produce variations on the template      #
#      input file (e.g. to scan over a parameter).  In this way      #
#      Nexus can drive codes that do not have specialized classes    #
#      derived from Simulation, SimulationInput, or                  #
#      SimulationAnalyzer.                                           #
#                                                                    #
#    SimulationInputMultiTemplate                                    #
#      Supports templated input files for codes that take many       #
#      different files as input (VASP is an example of this).  The   #
#      multi-template is essentially a collection of individual      #
#      template input files.                                         #
#                                                                    #
#    input_template                                                  #
#      User-facing function to create SimulationInputTemplate's.     #
#                                                                    #
#    multi_input_template                                            #
#      User-facing function to create SimulationInputMultiTemplate's.#
#                                                                    #
#====================================================================#


import os
import sys
import shutil
import string
from subprocess import Popen,PIPE
from developer import unavailable,ci
from generic import obj
from periodic_table import is_element
from physical_system import PhysicalSystem
from machines import Job
from pseudopotential import ppset
from nexus_base import NexusCore,nexus_core

 
class SimulationInput(NexusCore):
    def is_valid(self):
        self.not_implemented()
    #end def is_valid

    def read_file_text(self,filepath):
        if not os.path.exists(filepath):
            self.error('file does not exist:  '+filepath)
        #end if
        fobj = open(filepath,'r')
        text = fobj.read()
        fobj.close()
        return text
    #end def read_file_text

    def write_file_text(self,filepath,text):
        fobj = open(filepath,'w')
        fobj.write(text)
        fobj.flush()
        fobj.close()
    #end def write_file_text

    def read(self,filepath):
        tokens = filepath.split(None,1)
        if len(tokens)>1:
            text = filepath
            self.read_text(text)
        else:
            text = self.read_file_text(filepath)
            self.read_text(text,filepath)
        #end if
    #end def read

    def write(self,filepath=None):
        text = self.write_text(filepath)
        if filepath!=None:
            self.write_file_text(filepath,text)
        #end if
        return text
    #end def write

    def read_text(self,text,filepath=None):
        self.not_implemented()
    #end def read_text

    def write_text(self,filepath=None):
        self.not_implemented()
    #end def write_text

    def incorporate_system(self,system):
        #take information from a physical system object and fill in input file
        self.not_implemented()
    #end def incorporate_system

    def return_system(self):
        #create a physical system object from input file information
        self.not_implemented()
    #end def return_system

    def incorporate_descriptor(self,descriptor):
        #take information from a generic simulation descriptor object (ie DFT or QMC)
        #and fill in input file
        self.not_implemented()
    #end def incorporate_descriptor

    def return_descriptor(self):
        #create a general simulation descriptor object from input file information
        self.not_implemented()
    #end def return_descriptor
#end class SimulationInput




class SimulationAnalyzer(NexusCore):
    def __init__(self,sim):
        self.not_implemented()
    #end def __init__

    def analyze(self):
        self.not_implemented()
    #end def analyze
#end class SimulationAnalyzer




class SimulationEmulator(NexusCore):
    def run(self):
        self.not_implemented()
    #end def run
#end class SimulationEmulator



class SimulationImage(NexusCore):
    save_only_fields = set([
            # user block (temporary) of (sim+) subcascade
            'block',
            'block_subcascade',
            # local/remote/results directories
            'locdir',
            'remdir',
            'resdir',
            # image directories
            'imlocdir',
            'imremdir',
            'imresdir',
            ])

    load_fields = set([
            # important sim variables
            'identifier',
            'path',
            'process_id',
            # properties of the executable
            'app_name',
            'app_props',
            # names of in/out/err files
            'infile',
            'outfile',
            'errfile',
            # directory and image file names for sim/input/analyzer
            'image_dir',
            'sim_image',
            'input_image',
            'analyzer_image',
            # files copied in/out before/after run
            'files',
            'outputs',
            # simulation status flags
            'setup',
            'sent_files',
            'submitted',
            'finished',
            'failed',
            'got_output',
            'analyzed',
            # cascade status flag
            'subcascade_finished',
            ])

    save_fields = load_fields | save_only_fields

    def __init__(self):
        None
    #end def __init__

    def save_image(self,sim,imagefile):
        self.clear()
        self.transfer_from(sim,SimulationImage.save_fields)
        self.save(imagefile)
        self.clear()
    #end def save_image

    def load_image(self,sim,imagefile):
        self.clear()
        self.load(imagefile)
        self.transfer_to(sim,SimulationImage.load_fields)
        self.clear()
    #end def load_image

#end class SimulationImage



class Simulation(NexusCore):
    input_type    = SimulationInput
    analyzer_type = SimulationAnalyzer
    generic_identifier = 'sim'
    infile_extension   = '.in'
    outfile_extension  = '.out'
    errfile_extension  = '.err'
    application   = 'simapp'
    application_properties = set(['serial'])
    application_results    = set()
    allow_overlapping_files = False
    allowed_inputs = set(['identifier','path','infile','outfile','errfile','imagefile',
                          'input','job','files','dependencies','analysis_request',
                          'block','block_subcascade','app_name','app_props','system',
                          'skip_submit','force_write','simlabel','fake_sim',
                          'restartable','force_restart'])
    sim_imagefile      = 'sim.p'
    input_imagefile    = 'input.p'
    analyzer_imagefile = 'analyzer.p'
    image_directory    = 'sim'
    supports_restarts  = False

    is_bundle = False

    sim_count = 0
    creating_fake_sims = False

    sim_directories = dict()

    @classmethod
    def code_name(cls):
        return cls.generic_identifier
    #end def code_name

    @classmethod
    def separate_inputs(cls,kwargs,overlapping_kw=-1,copy_pseudos=True):
        if overlapping_kw==-1:
            overlapping_kw = set(['system'])
        elif overlapping_kw is None:
            overlapping_kw = set()
        #end if
        kw       = set(kwargs.keys())
        sim_kw   = kw & Simulation.allowed_inputs
        inp_kw   = (kw - sim_kw) | (kw & overlapping_kw)    
        sim_args = obj()
        inp_args = obj()
        sim_args.transfer_from(kwargs,sim_kw)
        inp_args.transfer_from(kwargs,inp_kw)
        if 'system' in inp_args:
            system = inp_args.system
            if not isinstance(system,PhysicalSystem):
                extra=''
                if not isinstance(extra,obj):
                    extra = '\nwith value: {0}'.format(system)
                #end if
                cls.class_error('invalid input for variable "system"\nsystem object must be of type PhysicalSystem\nyou provided type: {0}'.format(system.__class__.__name__)+extra)
            #end if
        #end if
        if 'pseudos' in inp_args and inp_args.pseudos!=None:
            pseudos = inp_args.pseudos
            # support ppset labels
            if isinstance(pseudos,str):
                code = cls.code_name()
                if not ppset.supports_code(code):
                    cls.class_error('ppset labeled pseudopotential groups are not supported for code "{0}"'.format(code))
                #end if
                if 'system' not in inp_args:
                    cls.class_error('system must be provided when using a ppset label')
                #end if
                system = inp_args.system
                pseudos = ppset.get(pseudos,code,system)
                if 'pseudos' in sim_args:
                    sim_args.pseudos = pseudos
                #end if
                inp_args.pseudos = pseudos
            #end if
            if copy_pseudos:
                if 'files' in sim_args:
                    sim_args.files = list(sim_args.files)
                else:
                    sim_args.files = list()
                #end if
                sim_args.files.extend(list(pseudos))
            #end if
            if 'system' in inp_args:
                system = inp_args.system
                species_labels,species = system.structure.species(symbol=True)
                pseudopotentials = nexus_core.pseudopotentials
                for ppfile in pseudos:
                    if not ppfile in pseudopotentials:
                        cls.class_error('pseudopotential file {0} cannot be found'.format(ppfile))
                    #end if
                    pp = pseudopotentials[ppfile]
                    if not pp.element_label in species_labels and not pp.element in species:
                        cls.class_error('the element {0} for pseudopotential file {1} is not in the physical system provided'.format(pp.element,ppfile))
                    #end if
                #end for
            #end if
        #end if
        # this is already done in Simulation.__init__()
        #if 'system' in inp_args and isinstance(inp_args.system,PhysicalSystem):
        #    inp_args.system = inp_args.system.copy()
        ##end if
        return sim_args,inp_args
    #end def separate_inputs


    def __init__(self,**kwargs):
        #user specified variables
        self.path          = None   #directory where sim will be run
        self.job           = None   #Job object for machine
        self.dependencies  = obj()  #Simulation results on which sim serially depends
        self.restartable   = False  #if True, job can be automatically restarted as deemed appropriate
        self.force_restart = False  #force a restart of the run

        #variables determined by self
        self.identifier     = self.generic_identifier
        self.simid          = Simulation.sim_count
        self.simlabel       = None
        Simulation.sim_count+=1
        self.files          = set()
        self.app_name       = self.application
        self.app_props      = list(self.application_properties)
        self.sim_image      = self.sim_imagefile
        self.input_image    = self.input_imagefile
        self.analyzer_image = self.analyzer_imagefile
        self.image_dir      = self.image_directory
        self.input          = self.input_type() 
        self.system         = None
        self.dependents     = obj()
        self.created_directories = False
        self.got_dependencies = False
        self.setup          = False
        self.sent_files     = False
        self.submitted      = False
        self.finished       = False
        self.failed         = False
        self.got_output     = False
        self.analyzed       = False
        self.subcascade_finished = False
        self.dependency_ids = set()
        self.wait_ids       = set()
        self.block          = False
        self.block_subcascade = False
        self.skip_submit    = nexus_core.skip_submit
        self.force_write    = False
        self.loaded         = False
        self.ordered_dependencies = []
        self.process_id     = None
        self.infile         = None
        self.outfile        = None
        self.errfile        = None
        self.bundled        = False
        self.bundler        = None
        self.fake_sim       = Simulation.creating_fake_sims

        #variables determined by derived classes
        self.outputs = None  #object representing output data 
                             # accessed by dependents when calling get_dependencies

        self.set(**kwargs)
        self.set_directories()
        self.set_files()
        self.propagate_identifier()
        if len(kwargs)>0:
            self.init_job()
        #end if
        self.post_init()

    #end def __init__


    def fake(self):
        return self.fake_sim
    #end def fake


    def init_job(self):
        if self.job==None:
            self.error('job not provided.  Input field job must be set to a Job object.')
        elif not isinstance(self.job,Job):
            self.error('Input field job must be set to a Job object\nyou provided an object of type: {0}\nwith value: {1}'.format(self.job.__class__.__name__,self.job))
        #end if
        self.job = self.job.copy()
        self.job.initialize(self)
    #end def init_job


    def set_app_name(self,app_name):
        self.app_name = app_name
    #end def set_app_name


    def set(self,**kw):
        cls = self.__class__
        if 'dependencies' in kw:
            self.depends(*kw['dependencies'])
            del kw['dependencies']
        #end if
        kwset = set(kw.keys())
        invalid = kwset - self.allowed_inputs
        if len(invalid)>0:
            self.error('received invalid inputs\ninvalid inputs: {0}\nallowed inputs are: {1}'.format(sorted(invalid),sorted(self.allowed_inputs)))
        #end if
        allowed =  kwset & self.allowed_inputs
        for name in allowed:
            self[name] = kw[name]
        #end for
        if 'path' in allowed:
            p = self.path
            if not isinstance(p,str):
                self.error('path must be a string, you provided {0} (type {1})'.format(p,p.__class__.__name__))
            #end if
            if p.startswith('./'):
                p = p[2:]
            #end if
            ld = nexus_core.local_directory
            if p.startswith(ld):
                p = p.split(ld)[1].lstrip('/')
            #end if
            self.path = p
        #end if
        if 'files' in allowed:
            self.files = set(self.files)
        #end if
        if not isinstance(self.input,(self.input_type,GenericSimulationInput)):
            self.error('input must be of type {0}\nreceived {1}\nplease provide input appropriate to {2}'.format(self.input_type.__name__,self.input.__class__.__name__,self.__class__.__name__))
        #end if
        if isinstance(self.system,PhysicalSystem):
            self.system = self.system.copy()
            consistent,msg = self.system.check_consistent(exit=False,message=True)
            if not consistent:
                locdir = os.path.join(nexus_core.local_directory,nexus_core.runs,self.path)
                self.error('user provided physical system is not internally consistent\nsimulation identifier: {0}\nlocal directory: {1}\nmore details on the user error are given below\n\n{2}'.format(self.identifier,locdir,msg))
            #end if
        elif self.system!=None:
            self.error('system must be a PhysicalSystem object\nyou provided an object of type: {0}'.format(self.system.__class__.__name__))
        #end if
        if self.restartable or self.force_restart:
            if not cls.supports_restarts:
                self.warn('restarts are not supported by {0}, request ignored'.format(cls.__name__))
            #end if
        #end if
    #end def set


    def set_directories(self):
        self.locdir = os.path.join(nexus_core.local_directory,nexus_core.runs,self.path)
        self.remdir = os.path.join(nexus_core.remote_directory,nexus_core.runs,self.path)
        self.resdir = os.path.join(nexus_core.local_directory,nexus_core.results,nexus_core.runs,self.path)
        
        if not self.fake():
            #print '  creating sim {0} in {1}'.format(self.simid,self.locdir)

            if not self.locdir in self.sim_directories:
                self.sim_directories[self.locdir] = set([self.identifier])
            else:
                idset = self.sim_directories[self.locdir]
                if not self.identifier in idset:
                    idset.add(self.identifier)
                else:
                    self.error('multiple simulations in a single directory have the same identifier\nplease assign unique identifiers to each simulation\nsimulation directory: {0}\nrepeated identifier: {1}\nother identifiers: {2}\nbetween the directory shown and the identifiers listed, it should be clear which simulations are involved\nmost likely, you described two simulations with identifier {3}'.format(self.locdir,self.identifier,sorted(idset),self.identifier))
                #end if
            #end if
        #end if

        self.image_dir = self.image_dir+'_'+self.identifier
        self.imlocdir = os.path.join(self.locdir,self.image_dir) 
        self.imremdir = os.path.join(self.remdir,self.image_dir) 
        self.imresdir = os.path.join(self.resdir,self.image_dir) 
    #end def set_directories


    def set_files(self):
        if self.infile is None:
            self.infile  = self.identifier + self.infile_extension
        #end if
        if self.outfile is None:
            self.outfile = self.identifier + self.outfile_extension
        #end if
        if self.errfile is None:
            self.errfile = self.identifier + self.errfile_extension
        #end if
    #end def set_files


    def reset_indicators(self):
        #this is now needed to enable restart support
        #self.error('remove this error call if you really want to use reset_indicators')
        self.got_dependencies = False
        self.setup          = False
        self.sent_files     = False
        self.submitted      = False
        self.finished       = False
        self.failed         = False
        self.got_output     = False
        self.analyzed       = False
    #end def reset_indicators

    def completed(self):
        completed  = self.setup
        completed &= self.sent_files 
        completed &= self.submitted  
        completed &= self.finished   
        completed &= self.got_output 
        completed &= self.analyzed   
        completed &= not self.failed
        return completed
    #end def completed

    def active(self):
        deps_completed = True
        for dep in self.dependencies:
            deps_completed &= dep.sim.completed()
        #end for
        active = deps_completed and not self.completed()
        return active
    #end def active

    def ready(self):
        ready = self.active()
        ready &= not self.submitted
        ready &= not self.finished
        ready &= not self.got_output
        ready &= not self.analyzed
        ready &= not self.failed
        return ready
    #end def ready

    def check_result(self,result_name):
        self.not_implemented()
    #end def check_result

    def get_result(self,result_name):
        self.not_implemented()
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        self.not_implemented()
    #end def incorporate_result

    def app_command(self):
        self.not_implemented()
    #end def app_command

    def check_sim_status(self):
        self.not_implemented()
    #end def check_sim_status

    def get_output_files(self): # returns list of output files to save
        self.not_implemented()
    #end def get_output_files


    def propagate_identifier(self):
        None
    #end def propagate_identifier

    def post_init(self):
        None
    #end def post_init

    def pre_create_directories(self):
        None
    #end def pre_create_directories

    def write_prep(self):
        None
    #end def write_prep

    def pre_write_inputs(self,save_image):
        None
    #end def pre_write_inputs

    def pre_send_files(self,enter):
        None
    #end def pre_send_files

    def post_submit(self):
        None
    #end def post_submit

    def pre_check_status(self):
        None
    #end def pre_check_status

    def post_analyze(self,analyzer):
        None
    #end def post_analyze


    def condense_name(self,name):
        return name.strip().lower().replace('-','_').replace(' ','_')
    #end def condense_name


    def has_generic_input(self):
        return isinstance(self.input,GenericSimulationInput)
    #end def has_generic_input

    def outfile_text(self):
        return self._file_text('outfile')
    #end def outfile_text

    def errfile_text(self):
        return self._file_text('errfile')
    #end def errfile_text

    def _file_text(self,filename):
        filepath = os.path.join(self.locdir,self[filename])
        fobj = open(filepath,'r')
        text = fobj.read()
        fobj.close()
        return text
    #end def _file_text


    def _create_dir(self,dir):
        if not os.path.exists(dir):
            os.makedirs(dir)
        elif os.path.isfile(dir):
            self.error('cannot create directory {0}\na file exists at this location'.format(dir))
        #end if
    #end def _create_dir

    def create_directories(self):
        self.pre_create_directories()
        self._create_dir(self.locdir)
        self._create_dir(self.imlocdir)
        self.created_directories = True
    #end def create_directories
            

    def depends(self,*dependencies):
        if len(dependencies)==0:
            return
        #end if
        if isinstance(dependencies[0],Simulation):
            dependencies = [dependencies]
        #end if
        for d in dependencies:
            sim = d[0]
            if not isinstance(sim,Simulation):
                self.error('first element in a dependency tuple must be a Simulation object\n  you provided a '+sim.__class__.__name__)
            #end if
            dep = obj()
            dep.sim = sim
            rn = []
            unrecognized_names = False
            app_results = sim.application_results | set(['other'])
            for name in d[1:]:
                result_name = self.condense_name(name)
                if result_name in app_results:
                    rn.append(result_name)
                else:
                    unrecognized_names = True
                    self.error(name+' is not known to be a result of '+sim.generic_identifier,exit=False)
                #end if
            #end for
            if unrecognized_names:
                self.error('unrecognized dependencies specified for simulation '+self.identifier)
            #end if
            dep.result_names = rn
            dep.results = obj()
            if not sim.simid in self.dependencies:
                self.ordered_dependencies.append(dep)
                self.dependencies[sim.simid]=dep
                sim.dependents[self.simid]=self
                self.dependency_ids.add(sim.simid)
                self.wait_ids.add(sim.simid)
            else:
                self.dependencies[sim.simid].result_names.extend(dep.result_names)
            #end if
        #end for
    #end def depends


    def undo_depends(self,sim):
        i=0
        for dep in self.ordered_dependencies:
            if dep.sim.simid==sim.simid:
                break
            #end if
            i+=1
        #end for
        self.ordered_dependencies.pop(i)
        del self.dependencies[sim.simid]
        del sim.dependents[self.simid]
        self.dependency_ids.remove(sim.simid)
        if sim.simid in self.wait_ids:
            self.wait_ids.remove(sim.simid)
        #end if
    #end def undo_depends


    def acquire_dependents(self,sim):
        # acquire the dependents from the other simulation
        dsims = obj(sim.dependents)
        for dsim in dsims:
            dep = dsim.dependencies[sim.simid]
            dsim.depends(self,*dep.result_names)
        #end for
        # eliminate the other simulation
        #   this renders it void (fake) and removes all dependency relationships
        sim.eliminate()
    #end def acquire_dependents


    def eliminate(self):
        # reverse relationship of dependents (downstream)
        dsims = obj(self.dependents)
        for dsim in dsims:
            dsim.undo_depends(self)
        #end for
        # reverse relationship of dependencies (upstream)
        deps = obj(self.dependencies)
        for dep in deps:
            self.undo_depends(dep.sim)
        #end for
        # mark sim to be ignored in all future interactions
        self.fake_sim = True
    #end def eliminate


    def check_dependencies(self,result):
        dep_satisfied = result.dependencies_satisfied
        for dep in self.dependencies:
            sim = dep.sim
            for result_name in dep.result_names:
                if result_name!='other':
                    if sim.has_generic_input():
                        calculating_result = False
                        cls = self.__class__
                        self.warn('a simulation result cannot be inferred from generic formatted or template input\nplease use {0} instead of {1}\nsee error below for information identifying this simulation instance'.format(cls.input_type.__class__.__name__,sim.input.__class__.__name__))
                    else:
                        calculating_result = sim.check_result(result_name,self)
                    #end if
                    if not calculating_result:
                        self.error('simulation {0} id {1} is not calculating result {2}\nrequired by simulation {3} id {4}\n{5} {6} directory: {7}\n{8} {9} directory: {10}'.format(sim.identifier,sim.simid,result_name,self.identifier,self.simid,sim.identifier,sim.simid,sim.locdir,self.identifier,self.simid,self.locdir),exit=False)
                    #end if
                else:
                    calculating_result = True
                #end if
                dep_satisfied = dep_satisfied and calculating_result
            #end for
        #end for
        result.dependencies_satisfied = dep_satisfied
    #end def check_dependencies


    def get_dependencies(self):
        if nexus_core.generate_only or self.finished:
            for dep in self.dependencies:
                for result_name in dep.result_names:
                    dep.results[result_name] = result_name
                #end for
            #end for
        else:
            for dep in self.dependencies:
                sim = dep.sim
                for result_name in dep.result_names:
                    if result_name!='other':
                        if sim.has_generic_input():
                            self.error('a simulation result cannot be inferred from generic formatted or template input\nplease use {0} instead of {1}\nsim id: {2}\ndirectory: {3}\nresult: {4}'.format(cls.input_type.__class__.__name__,sim.input.__class__.__name__,sim.id,sim.locdir,result_name))
                        #end if
                        dep.results[result_name] = sim.get_result(result_name,sim)
                    else:
                        dep.results['other'] = obj()
                    #end if
                #end for
            #end for
            if not self.got_dependencies:
                for dep in self.ordered_dependencies:
                    sim = dep.sim
                    for result_name,result in dep.results.iteritems():
                        if result_name!='other':
                            if self.has_generic_input():
                                self.error('a simulation result cannot be incorporated into generic formatted or template input\nplease use {0} instead of {1}\nsim id: {2}\ndirectory: {3}\nresult: {4}'.format(cls.input_type.__class__.__name__,self.input.__class__.__name__,self.id,self.locdir,result_name))
                            #end if
                            self.incorporate_result(result_name,result,sim)
                        #end if
                    #end for
                #end for
            #end if
        #end if
        self.got_dependencies = True
    #end def get_dependencies
        

    def downstream_simids(self,simids=None):
        if simids is None:
            simids = set()
        #end if
        for sim in self.dependents:
            simids.add(sim.simid)
            sim.downstream_simids(simids)
        #end for
        return simids
    #end def downstream_simids


    def copy_file(self,sourcefile,dest):
        src = os.path.dirname(os.path.abspath(sourcefile))
        dst = os.path.abspath(dest)
        if src!=dst:
            shutil.copy2(sourcefile,dest)
        #end if
    #end def copy_file


    def save_image(self,all=False):
        imagefile = os.path.join(self.imlocdir,self.sim_image)
        if os.path.exists(imagefile):
            os.system('rm '+imagefile)
        #end if
        if not all:
            sim_image = SimulationImage()
            sim_image.save_image(self,imagefile)
        else:
            self.error('attempting to save full object!')
            self.save(imagefile)
        #end if
    #end def save_image


    def load_image(self,imagepath=None,all=False):
        if imagepath==None:
            imagepath=os.path.join(self.imlocdir,self.sim_image)
        #end if
        if not all:
            sim_image = SimulationImage()
            sim_image.load_image(self,imagepath)
        else:
            self.load(imagepath)
        #end if
        # update process id for backwards compatibility
        if 'process_id' not in self:
            self.process_id = self.job.system_id
        #end if
    #end def load_image


    def load_analyzer_image(self,imagepath=None):
        if imagepath==None:
            imagepath = os.path.join(self.imresdir,self.analyzer_image)
        #end if
        analyzer = self.analyzer_type(self)
        analyzer.load(imagepath)
        return analyzer
    #end def load_analyzer_image


    def save_attempt(self):
        local = [self.infile,self.outfile,self.errfile]
        filepaths = []
        for file in local:
            filepath = os.path.join(self.locdir,file)
            if os.path.exists(filepath):
                filepaths.append(filepath)
            #end if
        #end for
        if len(filepaths)>0:
            prefix = self.identifier+'_attempt'
            n=0
            for dir in os.listdir(self.locdir):
                if dir.startswith(prefix):
                    n=max(n,int(dir.replace(prefix,'')))
                #end if
            #end for
            n+=1
            attempt_dir = os.path.join(self.locdir,prefix+str(n))
            os.makedirs(attempt_dir)
            for filepath in filepaths:
                os.system('mv {0} {1}'.format(filepath,attempt_dir))
            #end for
            #print self.locdir
            #os.system('ls '+self.locdir)
            #print attempt_dir
            #os.system('ls '+attempt_dir)
            #exit()
        #end if
        #self.error('save_attempt')
    #end def save_attempt


    def idstr(self):
        return '  '+str(self.simid)+' '+str(self.identifier)
    #end def idstr


    def write_inputs(self,save_image=True):
        self.pre_write_inputs(save_image)
        self.enter(self.locdir,False,self.simid)
        self.log('writing input files'+self.idstr(),n=3)
        self.write_prep()
        if self.infile!=None:
            infile = os.path.join(self.locdir,self.infile)
            self.input.write(infile)
        #end if
        self.job.write(file=True)
        self.setup = True
        if save_image:
            self.save_image()
            self.input.save(os.path.join(self.imlocdir,self.input_image))
        #end if
        #try to also write structure information
        if self.system!=None:
            filebase = os.path.join(self.locdir,self.identifier+'.struct')
            try:
                self.system.structure.write(filebase+'.xyz')
            except:
                None
            #end try
            try:
                if self.system.structure.has_axes():
                    self.system.structure.write(filebase+'.xsf')
                #end if
            except:
                None
            #end try
        #end if
    #end def write_inputs


    def send_files(self,enter=True):
        self.pre_send_files(enter)
        if enter:
            self.enter(self.locdir,False,self.simid)
        #end if
        self.log('sending required files'+self.idstr(),n=3)
        if not os.path.exists(self.remdir):
            os.makedirs(self.remdir)
        #end if
        if not os.path.exists(self.imremdir):
            os.makedirs(self.imremdir)
        #end if
        if self.infile!=None:
            self.files.add(self.infile)
        #end if
        send_files = self.files
        file_locations = [self.locdir]+nexus_core.file_locations
        remote = self.remdir
        for file in send_files:
            found_file = False
            for location in file_locations:                
                local = os.path.join(location,file)
                found_file = os.path.exists(local)
                if found_file:
                    break
                #end if
            #end if
            if found_file:
                self.copy_file(local,remote)
            else:
                self.error('file {0} not found\n  locations checked: {1}'.format(file,file_locations))
            #end if
        #end for
        self.sent_files = True
        self.save_image()
        send_imfiles=[self.sim_image,self.input_image]
        remote = self.imremdir
        for imfile in send_imfiles:
            local = os.path.join(self.imlocdir,imfile)
            if os.path.exists(local):
                self.copy_file(local,remote)
            #end if
        #end for
    #end def send_files


    def submit(self):
        if not self.submitted:
            self.log('submitting job'+self.idstr(),n=3)
            if not self.skip_submit:
                if not self.job.local:
                    self.job.submit()
                else:
                    self.execute() # execute local job immediately
                #end if
            #end if
            self.submitted = True
            if (self.job.batch_mode or not nexus_core.monitor) and not nexus_core.generate_only:
                self.save_image()
            #end if
        elif not self.finished:
            self.check_status()
        #end if
        self.post_submit()
    #end def submit


    def update_process_id(self):
        if self.process_id is None and self.job.system_id!=None:
            self.process_id = self.job.system_id
            self.save_image()
        #end if
    #end def update_process_id


    def check_status(self):
        self.pre_check_status()
        if nexus_core.generate_only: 
            self.finished = self.job.finished
        elif self.job.finished:
            should_check = True
            if self.outfile!=None:
                outfile = os.path.join(self.locdir,self.outfile)
                should_check &= os.path.exists(outfile)
            #end if
            if self.errfile!=None:
                errfile = os.path.join(self.locdir,self.errfile)
                should_check &= os.path.exists(errfile)
            #end if
            if not self.finished and should_check:
                self.check_sim_status()
            #end if
            if self.failed:
                self.finished = True
                # commented out block dependents 15/09/30
                # try to rely on persistent failed flag instead
                #self.block_dependents() 
            #end if
        #end if
        if self.finished:
            self.save_image()
        #end if
    #end def check_status


    def get_output(self):
        if not os.path.exists(self.resdir):
            os.makedirs(self.resdir)
        #end if
        if not os.path.exists(self.imresdir):
            os.makedirs(self.imresdir)
        #end if
        images = [self.sim_image,self.input_image]
        for image in images:
            remote_image = os.path.join(self.imremdir,image)
            if os.path.exists(remote_image):
                self.copy_file(remote_image,self.imresdir)
            #end if
        #end for
        results_image = os.path.join(self.imresdir,self.sim_image)
        if os.path.exists(results_image):
            self.load_image(results_image)
        #end if
        if self.finished:
            self.enter(self.locdir,False,self.simid)
            self.log('copying results'+self.idstr(),n=3)
            if not nexus_core.generate_only:
                output_files = self.get_output_files()
                if self.infile!=None:
                    output_files.append(self.infile)
                #end if
                if self.outfile!=None:
                    output_files.append(self.outfile)
                #end if
                if self.errfile!=None:
                    output_files.append(self.errfile)
                #end if
                files_missing = []
                for file in output_files:
                    remfile = os.path.join(self.remdir,file)
                    if os.path.exists(remfile):
                        self.copy_file(remfile,self.resdir)
                    else:
                        files_missing.append(file)
                    #end if
                #end for
                if len(files_missing)>0:
                    self.log('warning: the following files were missing',n=4)
                    for file in files_missing:
                        self.log(file,n=5)
                    #end for
                #end if
            #end if
            self.got_output = True
            self.save_image()
        #end if
    #end def get_output

        
    def analyze(self):
        if self.finished:
            self.enter(self.locdir,False,self.simid)
            self.log('analyzing'+self.idstr(),n=3)
            if not nexus_core.generate_only:
                analyzer = self.analyzer_type(self)
                analyzer.analyze()
                self.post_analyze(analyzer)
                analyzer.save(os.path.join(self.imresdir,self.analyzer_image))
                del analyzer
            #end if
            self.analyzed = True
            self.save_image()
        #end if
    #end def analyze


    def reset_wait_ids(self):
        self.wait_ids = set(self.dependency_ids)
        for sim in self.dependents:
            sim.reset_wait_ids()
        #end for
    #end def reset_wait_ids


    def check_subcascade(self):
        finished = self.finished or self.block
        if not self.block and not self.block_subcascade and not self.failed:
            for sim in self.dependents:
                finished &= sim.check_subcascade()
            #end for
        #end if
        self.subcascade_finished = finished
        return finished
    #end def check_subcascade


    def block_dependents(self,block_self=True):
        if block_self:
            self.block = True
        #end if
        self.block_subcascade = True
        for sim in self.dependents:
            sim.block_dependents()
        #end for
    #end def block_dependents


    def progress(self,dependency_id=None):
        if dependency_id!=None:
            self.wait_ids.remove(dependency_id)
        #end if
        if len(self.wait_ids)==0 and not self.block and not self.failed:
            modes = nexus_core.modes
            mode  = nexus_core.mode
            progress = True
            if mode==modes.none:
                return
            elif mode==modes.setup:
                self.write_inputs()
            elif mode==modes.send_files:
                self.send_files()
            elif mode==modes.submit:
                self.submit()
                progress = self.finished
            elif mode==modes.get_output:
                self.get_output()
                progress = self.finished
            elif mode==modes.analyze:
                self.analyze()
                progress = self.finished
            elif mode==modes.stages:
                if not self.created_directories:
                    self.create_directories()
                #end if
                if not self.got_dependencies:
                    self.get_dependencies()
                #end if
                if not self.setup and 'setup' in nexus_core.stages:
                    self.write_inputs()
                #end if
                if not self.sent_files and 'send_files' in nexus_core.stages:
                    self.send_files()
                #end if
                if not self.finished and 'submit' in nexus_core.stages:
                    self.submit()
                #end if
                if nexus_core.dependent_modes <= nexus_core.stages_set:
                    progress_post = self.finished
                    progress = self.finished and self.analyzed
                else:
                    progress_post = progress
                #end if
                if progress_post:
                    if not self.got_output and 'get_output' in nexus_core.stages:
                        self.get_output()
                    #end if
                    if not self.analyzed and 'analyze' in nexus_core.stages:
                        self.analyze()
                    #end if
                #end if
            elif mode==modes.all:
                if not self.setup:
                    self.write_inputs()
                    self.send_files(False)
                #end if
                if not self.finished:
                    self.submit()
                #end if
                if self.finished:
                    if not self.got_output:
                        self.get_output()
                    #end if
                    if not self.analyzed:
                        self.analyze()
                    #end if
                #end if
                progress = self.finished
            #end if
            if progress and not self.block_subcascade and not self.failed:
                for sim in self.dependents:
                    if not sim.bundled:
                        sim.progress(self.simid)
                    #end if
                #end for
            #end if
        elif len(self.wait_ids)==0 and self.force_write:
            modes = nexus_core.modes
            mode  = nexus_core.mode
            if mode==modes.stages:
                if not self.got_dependencies:
                    self.get_dependencies()
                #end if
                if 'setup' in nexus_core.stages:
                    self.write_inputs()
                #end if
                if not self.sent_files and 'send_files' in nexus_core.stages:
                    self.send_files()
                #end if
            #end if
        #end if
    #end def progress


    def reconstruct_cascade(self):
        imagefile = os.path.join(self.imlocdir,self.sim_image)
        if os.path.exists(imagefile) and not self.loaded:
            self.load_image()
            # continue from interruption
            if self.submitted and not self.finished and self.process_id!=None:
                self.job.system_id = self.process_id # load process id of job
                self.job.reenter_queue()
            #end if
            self.loaded = True
        #end if
        for sim in self.dependents:
            sim.reconstruct_cascade()
        #end for
        return self
    #end def reconstruct_cascade


    def traverse_cascade(self,operation,*args,**kwargs):
        if 'dependency_id' in kwargs:
            self.wait_ids.remove(kwargs['dependency_id'])
            del kwargs['dependency_id']
        #end if
        if len(self.wait_ids)==0:
            operation(self,*args,**kwargs)
            for sim in self.dependents:
                kwargs['dependency_id'] = self.simid
                sim.traverse_cascade(operation,*args,**kwargs)
            #end for
        #end if
    #end def traverse_cascade


#    def write_dependents(self,n=0,location=False):
#        n+=1
#        for sim in self.dependents:
#            if not location:
#                self.log(sim.__class__.__name__,sim.identifier,sim.simid,list(sim.dependency_ids),n=n)
#            else:
#                self.log(sim.__class__.__name__,sim.identifier,sim.simid,sim.locdir,list(sim.dependency_ids),n=n)
#            #end if
#            sim.write_dependents(n=n,location=location)
#        #end for
#    #end def write_dependents

    def write_dependents(self,n=0,location=False,block_status=False):
        outs = [self.__class__.__name__,self.identifier,self.simid]
        if location:
            outs.append(self.locdir)
        #end if
        if block_status:
            if self.block:
                outs.append('blocked')
            else:
                outs.append('unblocked')
            #end if
        #end if
        outs.append(list(self.dependency_ids))
        self.log(*outs,n=n)
        n+=1
        for sim in self.dependents:
            sim.write_dependents(n=n,location=location,block_status=block_status)
        #end for
    #end def write_dependents


    def execute(self,run_command=None):
        pad = self.enter(self.locdir)
        if run_command is None:
            job = self.job
            command = 'export OMP_NUM_THREADS='+str(job.threads)+'\n'
            if len(job.presub)>0:
                command += job.presub+'\n'
            #end if
            machine = job.get_machine()
            command += job.run_command(machine.app_launcher)
            if len(job.postsub)>0:
                command += job.postsub+'\n'
            #end if
            command = ('\n'+command).replace('\n','\n  '+pad)
            run_command = command
        #end if
        if self.job is None:
            env = os.environ.copy()
        else:
            env = job.env
        #end if
        if nexus_core.generate_only:
            self.log(pad+'Would have executed:  '+command)
        else:
            self.log(pad+'Executing:  '+command)
            fout = open(self.outfile,'w')
            ferr = open(self.errfile,'w')
            out,err = Popen(command,env=env,stdout=fout,stderr=ferr,shell=True,close_fds=True).communicate()
        #end if
        self.leave()
        self.submitted = True
        if self.job!=None:
            job.status = job.states.finished
            self.job.finished = True
        #end if
    #end def execute

#end class Simulation










 
class NullSimulationInput(SimulationInput):
    def is_valid(self):
        return True
    #end def is_valid

    def read(self,filepath):
        None
    #end def read

    def write(self,filepath=None):
        None
    #end def write

    def read_text(self,text,filepath=None):
        None
    #end def read_text

    def write_text(self,filepath=None):
        None
    #end def write_text

    def incorporate_system(self,system):
        None
    #end def incorporate_system

    def return_system(self):
        self.not_implemented()
    #end def return_system

    def incorporate_descriptor(self,descriptor):
        None
    #end def incorporate_descriptor

    def return_descriptor(self):
        self.not_implemented()
    #end def return_descriptor
#end class NullSimulationInput




class NullSimulationAnalyzer(SimulationAnalyzer):
    def __init__(self,sim):
        None
    #end def __init__

    def analyze(self):
        None
    #end def analyze
#end class NullSimulationAnalyzer


class GenericSimulationInput: # marker class for generic user input
    None
#end class GenericSimulationInput


class GenericSimulation(Simulation):
    def __init__(self,**kwargs):
        self.input_type    = NullSimulationInput
        self.analyzer_type = NullSimulationAnalyzer
        if 'input_type' in kwargs:
            self.input_type = kwargs['input_type']
            del kwargs['input_type']
        #end if
        if 'analyzer_type' in kwargs:
            self.input_type = kwargs['analyzer_type']
            del kwargs['analyzer_type']
        #end if
        if 'input' in kwargs:
            self.input_type = kwargs['input'].__class__
        #end if
        if 'analyzer' in kwargs:
            self.input_type = kwargs['analyzer'].__class__
        #end if
        Simulation.__init__(self,**kwargs)
    #end def __init__

    def check_sim_status(self):
        self.finished = True
    #end def check_sim_status

    def get_output_files(self):
        return []
    #end def get_output_files
#end class GenericSimulation



from string import Template
class SimulationInputTemplateDev(SimulationInput):
    def __init__(self,filepath=None,text=None):
        self.reset()
        if filepath!=None:
            self.read(filepath)
        elif text!=None:
            self.read_text(text)
        #end if
    #end def __init__
            
    def reset(self):
        self.template = None
        self.keywords = set()
        self.values   = obj()
    #end def reset
        
    def clear(self):
        self.values.clear()
    #end def clear

    def assign(self,**values):
        if self.template is None:
            self.error('cannot assign values prior to reading template')
        #end if
        invalid = set(values.keys()) - self.keywords
        if len(invalid)>0:
            self.error('attempted to assign invalid keywords\ninvalid keywords: {0}\nvalid options are: {1}'.format(sorted(invalid),sorted(self.keywords)))
        #end if
        self.values.set(**values)
    #end def assign

    def read_text(self,text,filepath=None):
        text = self.preprocess(text,filepath) # for derived class intervention
        try:
            template   = Template(text)
            key_tuples = Template.pattern.findall(text)
        except Exception,e:
            self.error('exception encountered during read\nfile: {0}\nexception: {1}'.format(filepath,e))
        #end try
        for ktup in key_tuples:
            if len(ktup[1])>0:   # normal keyword, e.g. $key
                self.keywords.add(ktup[1])
            elif len(ktup[2])>0: # braced keyword, e.g. ${key}
                self.keywords.add(ktup[2])
            #end if
        #end for
        self.template = template
    #end def read_text

    def write_text(self,filepath=None):
        kw_rem = self.keywords-set(self.values.keys())
        if len(kw_rem)>0:
            self.error('not all keywords for this template have been assigned\nkeywords remaining: {0}'.format(sorted(kw_rem)))
        #end if
        try:
            text = self.template.substitute(**self.values)
        except Exception,e:
            self.error('exception encountered during write:\n'+str(e))
        #end try
        return text
    #end def write_text

    def preprocess(self,text,filepath):
        return text # derived classes can modify text prior to template creation
    #end def preprocess
#end class SimulationInputTemplateDev




class SimulationInputMultiTemplateDev(SimulationInput):
    def __init__(self,delimiter='|',conditionals=None,defaults=None,**file_templates):
        self.filenames = obj()
        if len(file_templates)>0:
            self.set_templates(delimiter,conditionals,defaults,**file_templates)
        #end if
    #end def __init_


    def set_filenames(self,**filenames):
        self.filenames.set(**filenames)
    #end def set_filenames
        

    def set_templates(self,delimiter='|',conditionals=None,defaults=None,**file_templates):
        for name,val in file_templates.iteritems():
            if isinstance(val,str):
                if ' ' in val:
                    self.error('filename cannot have any spaces\nbad filename provided with keyword '+name)
                #end if
                self.filenames[name] = val
            else:
                filename,template_path = val
                self[name] = SimulationInputTemplate(
                    filepath     = template_path,
                    delimiter    = delimiter,
                    conditionals = conditionals,
                    defaults     = defaults
                    )
                self.filenames[name] = filename
            #end if
        #end for
        if len(self)>1 and len(self.filenames)!=len(self)-1:
            self.error('keyword inputs must either be all filenames or all filename/filepath pairs')
        #end if
    #end def set_templates


    def read(self,filepath):
        if len(self.filenames)==0:
            self.error('cannot perform read, filenames are not set')
        #end if
        base,filename = os.path.split(filepath)
        filenames = self.filenames
        self.clear()
        templates = obj()
        for name,filename in filenames.iteritems():
            templates[name] = filename, os.path.join(base,filename)
        #end for
        self.set_templates(**templates)
    #end def read


    def write(self,filepath=None):
        if filepath is None:
            contents = obj()
            for name in self.filenames.keys():
                contents[name] = self[name].write()
            #end for
            return contents
        else:
            base,filename = os.path.split(filepath)
            for name,filename in self.filenames.iteritems():
                self[name].write(os.path.join(base,filename))
            #end for
        #end if
    #end def write
#end class SimulationInputMultiTemplateDev


# these are for user access, *Dev are for development
class SimulationInputTemplate(SimulationInputTemplateDev,GenericSimulationInput):
    None
#end class SimulationInputTemplate

class SimulationInputMultiTemplate(SimulationInputMultiTemplateDev,GenericSimulationInput):
    None
#end class SimulationInputMultiTemplate



# developer functions

def input_template_dev(*args,**kwargs):
    return SimulationInputTemplateDev(*args,**kwargs)
#end def input_template_dev


def multi_input_template_dev(*args,**kwargs):
    return SimulationInputMultiTemplateDev(*args,**kwargs)
#end def multi_input_template_dev






# user functions

def input_template(*args,**kwargs):
    return SimulationInputTemplate(*args,**kwargs)
#end def input_template


def multi_input_template(*args,**kwargs):
    return SimulationInputMultiTemplate(*args,**kwargs)
#end def multi_input_template


def generate_template_input(*args,**kwargs):
    return SimulationInputTemplate(*args,**kwargs)
#end def generate_template_input


def generate_multi_template_input(*args,**kwargs):
    return SimulationInputMultiTemplate(*args,**kwargs)
#end def generate_multi_template_input


def generate_simulation(**kwargs):
    sim_type='generic'
    if 'sim_type' in kwargs:
        sim_type = kwargs['sim_type']
        del kwargs['sim_type']
    #end if
    if sim_type=='generic':
        return GenericSimulation(**kwargs)
    else:
        Simulation.class_error('sim_type {0} is unrecognized'.format(sim_type),'generate_simulation')
    #end if
#end def generate_simulation




# ability to graph simulation workflows
try:
    from pydot import Dot,Node,Edge
except:
    Dot,Node,Edge = unavailable('pydot','Dot','Node','Edge')
#end try
try:
    import Image
except:
    Image = unavailable('Image')
#end try
import tempfile
exit_call = sys.exit
def graph_sims(sims,useid=False,exit=True,quants=True):
    graph = Dot(graph_type='digraph')
    graph.set_label('simulation workflows')
    graph.set_labelloc('t')
    nodes = obj()
    for sim in sims:
        if 'fake_sim' in sim and sim.fake_sim:
            continue
        #end if
        if sim.simlabel is not None and not useid:
            nlabel = sim.simlabel+' '+str(sim.simid)
        else:
            nlabel = sim.identifier+' '+str(sim.simid)
        #end if
        node = obj(
            id    = sim.simid,
            sim   = sim,
            node  = Node(nlabel,style='filled',shape='Mrecord'),
            edges = obj(),
            )
        nodes[node.id] = node
        graph.add_node(node.node)
    #end for
    for node in nodes:
        for simid,dep in node.sim.dependencies.iteritems():
            other = nodes[simid].node
            if quants:
                for quantity in dep.result_names:
                    edge = Edge(other,node.node,label=quantity,fontsize='10.0')
                    graph.add_edge(edge)
                #end for
            else:
                edge = Edge(other,node.node)
                graph.add_edge(edge)
            #end if
        #end for
    #end for

    fout = tempfile.NamedTemporaryFile(suffix='png')
    savefile = fout.name
    graph.write(savefile,format='png',prog='dot')

    image = Image.open(savefile)
    image.show()

    if exit:
        exit_call()
    #end if
#end def graph_sims



#def write_sims(sims,
