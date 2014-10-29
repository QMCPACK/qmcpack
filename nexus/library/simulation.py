import os
import shutil
import string
from generic import obj
from physical_system import PhysicalSystem
from machines import Job
from project_base import Pobj

 
class SimulationInput(Pobj):
    @classmethod
    def templates(cls,template_name):
        print 'Developer Error: templates function has not been implemented in '+cls.__name__
        exit()
        template = None
        return template
    #end def

    def is_valid(self):
        self.not_implemented()
    #end def is_valid

    def read(self,filepath):
        tokens = filepath.split(None,1)
        if len(tokens)>1:
            contents = filepath
            self.read_contents(contents)
        else:
            if os.path.exists(filepath):
                self.read_contents(open(filepath,'r').read())
            else:
                self.error('failed read for filepath below\n    '+filepath)
            #end if
        #end if
    #end def read

    def write(self,filepath=None):
        contents = self.write_contents()
        if filepath is None:
            return contents
        else:
            fobj = open(filepath,'w')
            fobj.write(contents)
            fobj.flush()
            fobj.close()
        #end if
    #end def write

    def read_contents(self,contents):
        self.not_implemented()
    #end def read_contents

    def write_contents(self):
        self.not_implemented()
    #end def write_contents

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




class SimulationAnalyzer(Pobj):
    def __init__(self,sim):
        self.not_implemented()
    #end def __init__

    def analyze(self):
        self.not_implemented()
    #end def analyze
#end class SimulationAnalyzer




class SimulationEmulator(Pobj):
    def run(self):
        self.not_implemented()
    #end def run
#end class SimulationEmulator







class Simulation(Pobj):
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
                          'skip_submit','force_write'])
    sim_imagefile      = 'sim.p'
    input_imagefile    = 'input.p'
    analyzer_imagefile = 'analyzer.p'
    image_directory    = 'sim'

    updatable = set(['block','block_subcascade','force_write'])
    preserve  = set(['simid','got_dependencies','dependencies',
                     'ordered_dependencies','dependents','dependency_ids',
                     'wait_ids','input','locdir','remdir','resdir',
                     'imlocdir','imremdir','imresdir',
                     'skip_submit','job','system'])
    sim_count = 0


    sim_directories = dict()


    @staticmethod
    def separate_inputs(kwargs,overlapping_kw=-1,copy_pseudos=True):
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
        if 'pseudos' in inp_args:
            pseudos = inp_args.pseudos
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
                if not isinstance(system,PhysicalSystem):
                    Simulation.class_error('system object must be of type PhysicalSystem')
                #end if
                pseudopotentials = Simulation.pseudopotentials
                for ppfile in pseudos:
                    if not ppfile in pseudopotentials:
                        Simulation.class_error('pseudopotential file {0} cannot be found'.format(ppfile))
                    #end if
                    pp = pseudopotentials[ppfile]
                    if not pp.element in system.structure.elem:
                        Simulation.class_error('the element {0} for pseudopotential file {1} is not in the physical system provided'.format(pp.element,ppfile))
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
        self.path         = None   #directory where sim will be run
        self.job          = None   #Job object for machine
        self.dependencies = obj()  #Simulation results on which sim serially depends

        #variables determined by self
        self.identifier     = self.generic_identifier
        self.simid          = Simulation.sim_count
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
        self.got_dependencies = False
        self.setup          = False
        self.sent_files     = False
        self.submitted      = False
        self.finished       = False
        self.failed         = False
        self.got_output     = False
        self.analyzed       = False
        self.succeeded      = False
        self.subcascade_finished = False
        self.dependency_ids = set()
        self.wait_ids       = set()
        self.block          = False
        self.block_subcascade = False
        self.skip_submit    = self.skip_submission
        self.force_write    = False
        self.loaded         = False
        self.ordered_dependencies = []
        self.process_id     = None

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
        if self.system!=None:
            self.system.check_folded_system()
        #end if
    #end def __init__


    def init_job(self):
        if self.job==None:
            self.error('job not provided.  Input field job must be set to a Job object.')
        else:
            self.job = self.job.copy()
        #end if
        self.job.initialize(self)
    #end def init_job


    def set_app_name(self,app_name):
        self.app_name = app_name
    #end def set_app_name


    def set(self,**kw):
        if 'dependencies' in kw:
            self.depends(*kw['dependencies'])
            del kw['dependencies']
        #end if
        allowed = set(kw.keys()) & self.allowed_inputs
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
            ld = self.local_directory
            if p.startswith(ld):
                p = p.split(ld)[1].lstrip('/')
            #end if
            self.path = p
        #end if
        if 'files' in allowed:
            self.files = set(self.files)
        #end if
        if not isinstance(self.input,self.input_type):
            self.error('input must be of type {0}\nreceived {1}\nplease provide input appropriate to {2}'.format(self.input_type.__name__,self.input.__class__.__name__,self.__class__.__name__))
        #end if
        if isinstance(self.system,PhysicalSystem):
            self.system = self.system.copy()
        #end if
    #end def set


    def set_directories(self):
        self.locdir = os.path.join(self.local_directory,self.runs,self.path)
        self.remdir = os.path.join(self.remote_directory,self.runs,self.path)
        self.resdir = os.path.join(self.local_directory,self.results,self.runs,self.path)
        
        if not self.locdir in self.sim_directories:
            self.sim_directories[self.locdir] = set([self.identifier])
        else:
            idset = self.sim_directories[self.locdir]
            if not self.identifier in idset:
                idset.add(self.identifier)
            else:
                self.error('multiple simulations in a single directory have the same identifier\n  please assign unique identifiers to each simulation\n  simulation directory: {0}\n  repeated identifier: {1}\n  other identifiers: {2}\n  between the directory shown and the identifiers listed, it should be clear which simulations are involved\n  most likely, you described two simulations with identifier {3}'.format(self.locdir,self.identifier,list(idset),self.identifier))
            #end if
        #end if

        self.image_dir = self.image_dir+'_'+self.identifier
        self.imlocdir = os.path.join(self.locdir,self.image_dir) 
        self.imremdir = os.path.join(self.remdir,self.image_dir) 
        self.imresdir = os.path.join(self.resdir,self.image_dir) 
    #end def set_directories


    def set_files(self):
        self.infile  = self.identifier + self.infile_extension
        self.outfile = self.identifier + self.outfile_extension
        self.errfile = self.identifier + self.errfile_extension
    #end def set_files


    def reset_indicators(self):
        self.error('remove this error call if you really want to use reset_indicators')
        self.got_dependencies = False
        self.setup          = False
        self.sent_files     = False
        self.submitted      = False
        self.finished       = False
        self.failed         = False
        self.got_output     = False
        self.analyzed       = False
    #end def reset_indicators

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

    def get_output_files(self):
        self.not_implemented()
    #end def get_output_files


    def propagate_identifier(self):
        None
    #end def propagate_identifier

    def post_init(self):
        None
    #end def post_init

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


    def condense_name(self,name):
        return name.strip().lower().replace('-','_').replace(' ','_')
    #end def condense_name


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


    def check_dependencies(self,result):
        dep_satisfied = result.dependencies_satisfied
        for dep in self.dependencies:
            sim = dep.sim
            for result_name in dep.result_names:
                if result_name!='other':
                    calculating_result = sim.check_result(result_name,self)
                    if not calculating_result:
                        self.error('simulation {0} id {1} is not calculating result {2}\nrequired by simulation {3} id {4}\n{5} {6} directory: {7}\n{8} {9} directory: {10}'.format(sim.identifier,sim.simid,result_name,self.identifier,self.simid,sim.identifier,sim.simid,sim.locdir,self.identifier,self.simid,self.locdir),exit=False)
                else:
                    calculating_result = True
                #end if
                dep_satisfied = dep_satisfied and calculating_result
            #end for
        #end for
        result.dependencies_satisfied = dep_satisfied
    #end def check_dependencies


    def get_dependencies(self):
        if self.generate_only or self.finished:
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
                            self.incorporate_result(result_name,result,sim)
                        #end if
                    #end for
                #end for
            #end if
        #end if
        self.got_dependencies = True
    #end def get_dependencies
        

    def copy_file(self,sourcefile,dest):
        src = os.path.dirname(os.path.abspath(sourcefile))
        dst = os.path.abspath(dest)
        if src!=dst:
            shutil.copy2(sourcefile,dest)
        #end if
    #end def copy_file


    def save_image(self,all=False):
        self.tlog('save image',self.simid,n=5)
        imagefile = os.path.join(self.imlocdir,self.sim_image)
        if os.path.exists(imagefile):
            os.system('rm '+imagefile)
        #end if
        if not all:
            pres = dict()
            for v in self.preserve:
                pres[v] = self[v]
                del self[v]
            #end for
        #end if
        self.save(imagefile)
        if not all:
            for k,v in pres.iteritems():
                self[k] = v
            #end for
            del pres
        #end if
        self.tlog('end save image',self.simid,n=5)
    #end def save_image


    def load_image(self,imagepath=None,all=False):
        self.tlog('load image',self.simid,n=5)
        if imagepath==None:
            imagepath=os.path.join(self.imlocdir,self.sim_image)
        #end if
        if not all:
            pres = dict()
            for v in self.preserve:
                pres[v] = self[v]
            #end for
            for v in self.updatable:
                pres[v] = self[v]
            #end for
        #end if
        self.load(imagepath)
        if not all:
            for k,v in pres.iteritems():
                self[k] = v
            #end for
            del pres
        #end if
        # update process id for backwards compatibility
        if not 'process_id' in self:
            self.process_id = self.job.system_id
        #end if
        self.tlog('end load image',self.simid,n=5)
        return
    #end def load_image


    def load_analyzer_image(self,imagepath=None):
        if imagepath==None:
            imagepath = os.path.join(self.imresdir,self.analyzer_image)
        #end if
        analyzer = self.analyzer_type(self)
        analyzer.load(imagepath)
        return analyzer
    #end def load_analyzer_image


    def idstr(self):
        return '  '+str(self.simid)+' '+str(self.identifier)
    #end def idstr


    def write_inputs(self,save_image=True):
        self.tlog('write inputs',self.simid,n=4)
        self.pre_write_inputs(save_image)
        self.enter(self.locdir,False,self.simid)
        self.log('writing input files'+self.idstr(),n=3)
        if not os.path.exists(self.locdir):
            os.makedirs(self.locdir)
        #end if
        if not os.path.exists(self.imlocdir):
            os.makedirs(self.imlocdir)
        #end if
        self.write_prep()
        if self.infile!=None:
            infile = os.path.join(self.locdir,self.infile)
            self.input.write(infile)
        #end if
        job_file = self.job.write()
        if job_file!=None:
            self.files.add(job_file)
        #end if
        self.setup = True
        if save_image:
            self.save_image()
            self.input.save(os.path.join(self.imlocdir,self.input_image))
        #end if
        self.tlog('end write inputs',self.simid,n=4)
    #end def write_inputs


    def send_files(self,enter=True):
        self.tlog('send files',self.simid,n=4)
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
        file_locations = [self.locdir]+self.file_locations
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
        self.tlog('end send files',self.simid,n=4)
    #end def send_files


    def submit(self):
        self.tlog('submit',self.simid,n=4)
        if not self.submitted:
            self.log('submitting job'+self.idstr(),n=3)
            if not self.skip_submit:
                self.job.submit()
            #end if
            self.submitted = True
            if (self.job.batch_mode or not self.monitor) and not self.generate_only:
                self.save_image()
            #end if
        elif not self.finished:
            self.check_status()
        #end if
        self.post_submit()
        self.tlog('end submit',self.simid,n=4)
    #end def submit


    def update_process_id(self):
        if self.process_id is None and self.job.system_id!=None:
            self.process_id = self.job.system_id
            self.save_image()
        #end if
    #end def update_process_id


    def check_status(self):
        self.pre_check_status()
        self.tlog('check status',self.simid,n=5)
        if self.generate_only: 
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
                self.block_dependents()
            #end if
        #end if
        if self.finished:
            self.save_image()
        #end if
        self.tlog('end check status',self.simid,n=5)
    #end def check_status


    def get_output(self):
        self.tlog('get output',self.simid,n=4)
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
            if not self.generate_only:
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
        self.tlog('end get output',self.simid,n=4)
    #end def get_output

        
    def analyze(self):
        self.tlog('analyze',self.simid,n=4)
        if self.finished:
            self.enter(self.locdir,False,self.simid)
            self.log('analyzing'+self.idstr(),n=3)
            analyzer = self.analyzer_type(self)
            if not self.generate_only:
                analyzer.analyze()
            #end if
            analyzer.save(os.path.join(self.imresdir,self.analyzer_image))
            del analyzer
            self.analyzed = True
            self.save_image()
        #end if
        self.tlog('end analyze',self.simid,n=4)
    #end def analyze


    def reset_wait_ids(self):
        self.wait_ids = set(self.dependency_ids)
        self.dlog(self.simid,'wids',self.wait_ids,n=2)
        for sim in self.dependents:
            sim.reset_wait_ids()
        #end for
    #end def reset_wait_ids


    def check_subcascade(self):
        finished = self.finished or self.block
        if not finished:
            self.dlog(self.simid,'is not finished',n=2)
        else:
            self.dlog(self.simid,'is finished',n=2)
        #end if
        if not self.block and not self.block_subcascade:
            for sim in self.dependents:
                dfin = sim.check_subcascade()
                finished = finished and dfin
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
        self.tlog('progress',self.simid,n=2)
        if dependency_id!=None:
            self.dlog('wid',self.wait_ids,n=3)
            self.wait_ids.remove(dependency_id)
        #end if
        if len(self.wait_ids)==0 and not self.block:
            self.dlog('running',self.simid,n=3)
            modes = self.modes
            mode  = self.mode
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
                if not self.got_dependencies:
                    self.get_dependencies()
                #end if
                if not self.setup and 'setup' in self.stages:
                    self.write_inputs()
                #end if
                if not self.sent_files and 'send_files' in self.stages:
                    self.send_files()
                #end if
                if not self.finished and 'submit' in self.stages:
                    self.submit()
                #end if
                if self.dependent_modes <= self.stages_set:
                    progress_post = self.finished
                    progress = self.finished and self.analyzed
                else:
                    progress_post = progress
                #end if
                if progress_post:
                    if not self.got_output and 'get_output' in self.stages:
                        self.get_output()
                    #end if
                    if not self.analyzed and 'analyze' in self.stages:
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
            self.dlog('progress sub',self.simid,n=4)
            w,s,j,f,g,a=int(self.setup),int(self.submitted),int(self.job.finished),int(self.finished),int(self.got_output),int(self.analyzed)
            self.dlog('w,s,j,f,g,a',w,s,j,f,g,a,n=4)
            #self.dlog('f,jf,p,bs,nd',self.finished,self.job.finished,progress,self.block_subcascade,len(self.dependents),n=4)
            if progress and not self.block_subcascade:
                for sim in self.dependents:
                    self.dlog(self.simid,'is progressing',sim.simid,n=2)
                    sim.progress(self.simid)
                #end for
            #end if
        elif len(self.wait_ids)==0 and self.force_write:
            self.dlog('running',self.simid,n=3)
            modes = self.modes
            mode  = self.mode
            if mode==modes.stages:
                if not self.got_dependencies:
                    self.get_dependencies()
                #end if
                if 'setup' in self.stages:
                    self.write_inputs()
                #end if
                if not self.sent_files and 'send_files' in self.stages:
                    self.send_files()
                #end if
            #end if
        #end if
        self.tlog('end progress',self.simid,n=2)
    #end def progress


    def reconstruct_cascade(self):
        imagefile = os.path.join(self.imlocdir,self.sim_image)
        self.dlog(self.simid,'reconstructing, dids',self.dependency_ids,n=2)
        self.dlog(os.path.exists(imagefile),self.loaded,n=3)
        if os.path.exists(imagefile) and not self.loaded:
            self.load_image()
            # continue from interruption
            if self.submitted and not self.finished and self.process_id!=None:
                self.job.system_id = self.process_id # load process id of job
                self.job.reenter_queue()
            #end if
            self.loaded = True
        #end if
        self.dlog(self.simid,'reconstructed, dids',self.dependency_ids,n=3)
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

    def write_dependents(self,n=0,location=False):
        if not location:
            self.log(self.__class__.__name__,self.identifier,self.simid,list(self.dependency_ids),n=n)
        else:
            self.log(self.__class__.__name__,self.identifier,self.simid,self.locdir,list(self.dependency_ids),n=n)
        #end if
        n+=1
        for sim in self.dependents:
            sim.write_dependents(n=n,location=location)
        #end for
    #end def write_dependents

#end class Simulation

Simulation.set_user_interface(
    member_variables=['path','infile','outfile','input','job','files','dependencies','app_name'],
    member_functions=['set','depends','progress']                       
    )
Simulation.set_dev_instruction(
    'writing derived class',
    class_variables=['input_type','analyzer_type','infile_name'],
    member_functions=['check_result','get_result','app_command','check_sim_status','get_output_files']
    )










 
class NullSimulationInput(SimulationInput):
    @classmethod
    def templates(cls,template_name):
        print 'Developer Error: templates function has not been implemented in '+cls.__name__
        exit()
    #end def

    def is_valid(self):
        return True
    #end def is_valid

    def read(self,filepath):
        None
    #end def read

    def write(self,filepath=None):
        None
    #end def write

    def read_contents(self,contents):
        None
    #end def read_contents

    def write_contents(self):
        None
    #end def write_contents

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




class SimulationInputTemplate(SimulationInput):
    name_chars = string.ascii_letters+string.digits+'_'
    name_characters = set(name_chars)

    comparators = '== != <= >= < >'.split()

    def __init__(self,filepath=None,delimiter='|',conditionals=None,defaults=None):
        self.delimiter = delimiter
        self.keywords  = obj()
        self.values    = obj()
        self.use_default_value = obj()
        self.template  = ''
        if defaults is None:
            self.default_values = obj()
        else:
            self.default_values = obj(**defaults)
        #end if
        if conditionals is None:
            self.conditional_variables = obj()
        else:
            self.conditional_variables = obj(**conditionals)
        #end if
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__


    def set_delimiter(self,delimiter):
        if not isinstance(delimiter,str):
            self.error('delimiter must be a character, received '+str(delimiter))
        elif len(delimiter)!=1:
            self.error('delimiter must be a single character, received'+str([delimiter]))
        #end if
        self.delimiter = delimiter
    #end def set_delimiter


    def read(self,filepath):
        tokens = filepath.split(None,1)
        if len(tokens)>1:
            contents = filepath
            self.read_contents(contents)
        else:
            if os.path.exists(filepath):
                self.read_contents(open(filepath,'r').read(),filepath)
            else:
                self.error('failed read for filepath below\n    '+filepath)
            #end if
        #end if
    #end def read


    def preprocess(self,contents,filepath=None):
        return contents
    #end def preprocess

    def read_contents(self,contents,filepath=None):
        contents = self.preprocess(contents,filepath)
        delimiter = self.delimiter
        dlocs = []
        i=0
        for c in contents:
            if c==delimiter:
                dlocs.append(i)
            #end if
            i+=1
        #end for
        klocs=[]
        keywords = []
        nconditionals = len(self.conditional_variables)
        conditional_expressions = obj()
        conditional_words = []
        for i in range(len(dlocs)-1):
            d1 = dlocs[i]
            d2 = dlocs[i+1]
            word = contents[d1+1:d2]
            expressions = None
            if nconditionals>0 and ':' in word:
                old_word = word
                word,expressions = word.split(':')
                conditional_words.append((old_word,word))
            #end if
            if set(word)<=self.name_characters:
                keywords.append(word)
                klocs.append(i)
                if expressions!=None:
                    self.screen_conditional_expressions(word,expressions)
                    conditional_expressions[word] = expressions
                #end if
            #end if
        #end for
        for i in range(len(klocs)-1):
            if klocs[i+1]-klocs[i]==1:
                self.error('keywords cannot be adjacent\n  offending text: {0}{1}{0}{2}{0}'.format(delimiter,keywords[i],keywords[i+1]))
            #end if
        #end for
        keywords = set(keywords)
        for keyword in keywords:
            self.use_default_value[keyword] = False
        #end for
        for keyword in keywords:
            kw = keyword.lower()
            self.keywords[kw] = keyword
        #end for
        defaults = self.default_values
        for keyword in list(defaults.keys()):
            kw = keyword.lower()
            if kw in self.keywords:
                value = defaults[keyword]
                del defaults[keyword]
                defaults[self.keywords[kw]] = value
            else:
                self.error('default value provided for non-existent template keyword\n  non-existent keyword: {0}\n  default value provided: {1}\n  valid keywords are:\n    {2}'.format(keyword,value,self.keywords.keys()))
            #end if
        #end for
        self.evaluate_conditional_expressions(conditional_expressions)
        for old_word,word in conditional_words:
            old_word = delimiter+old_word+delimiter
            word = delimiter+word+delimiter
            contents = contents.replace(old_word,word)
        #end for
        self.template = contents
    #end def read_contents


    def write_contents(self):
        kw_rem = self.required_keywords_remaining()
        if len(kw_rem)>0:
            kw_rem = list(kw_rem)
            kw_rem.sort()
            self.error('not all keywords for this template have been assigned\n  keywords remaining:\n  '+str(kw_rem))
        #end if
        contents = self.fillout()
        return contents
    #end def write_contents


    def keywords_remaining(self):
        return set(self.keywords.values())-set(self.values.keys())
    #end def keywords_remaining


    def required_keywords_remaining(self):
        return self.keywords_remaining()-set(self.default_values.keys())
    #end def required_keywords_remaining


    def assign(self,**values):
        for keyword,value in values.iteritems():
            kw = keyword.lower()
            if not kw in self.keywords:
                self.error('cannot assign {0} because it is not a valid keyword for this template\n  valid options are: {1}'.format(kw,self.keywords.keys()))
#end class SimulationInputTemplate
            #end if
            kw = self.keywords[kw]
            self.values[kw] = value
        #end for
    #end def assign


    def fillout(self):
        d = self.delimiter
        template = str(self.template)
        for keyword in self.values.keys():
            if self.use_default_value[keyword]:
                value = self.default_values[keyword]
            else:
                value = self.values[keyword]
            #end if
            template = template.replace(d+keyword+d,str(value))
        #end for
        for keyword,value in self.default_values.iteritems():
            template = template.replace(d+keyword+d,str(value))
        #end for
        return template
    #end def fillout


    def screen_conditional_expressions(self,word,expressions):
        es = str(expressions)
        for comp in self.comparators:
            es = es.replace(comp,'@')
        #end for
        es = es.replace('not ','').replace(' or ',':').replace(' and ',':')
        es = es.split(':')
        for e in es:
            tokens = e.split('@')
            if not len(tokens)==2:
                self.error('invalid sentinel expression encountered for template keyword {0}\n  invalid expression: {1}\n  this expression must be of the form "name [comparator] value" or a boolean combination of them "expr1 and expr2 or expr3 ..."'.format(word,expressions))
            #end if
            name,value = tokens
            if not name.strip() in self.conditional_variables:
                self.error('conditional variable {0} has not been provided\n  valid options are: {1}'.format(name,self.conditional_variables.keys()))
            #end if
            #if not comparator in self.comparators:
            #    self.error('invalid comparator encountered in conditional expression for template keyword {0}\n  invalid comparator encountered: {1}\n  valid options are: {2}'.format(word,comparator,list(self.comparators)))
            ##end if
        #end for
    #end def screen_conditional_expressions


    def evaluate_conditional_expressions(self,conditional_expressions):
        if len(conditional_expressions)>0:
            for _name,_value in self.conditional_variables.iteritems():
                try:
                    exec('{0}=_value'.format(_name,_value))
                except Exception,e:
                    self.error('assignment of conditional variable failed\n  conditional variable: {0}\n  error message: {1}'.format(_name,e))
                #end try
            #end for
            for _keyword,_expression in conditional_expressions.iteritems():
                if not _keyword in self.default_values:
                    self.error('a default value must be provided for template keyword {0} because it has a conditional expression:\n  {1}'.format(_keyword,_expression))
                #end if
                try:
                    exec('_passes = '+_expression)
                except Exception,e:
                    self.error('evaluation of conditional expression for template keyword failed\n  template_keyword: {0}\n  conditional expression: {1}\n  error message: {2}'.format(_keyword,_expression,e))
                #end try
                self.use_default_value[_keyword] = not _passes
            #end for
        #end for
    #end if
#end class SimulationInputTemplate




class SimulationInputMultiTemplate(SimulationInput):
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
#end class SimulationInputMultiTemplate




def input_template(*args,**kwargs):
    return SimulationInputTemplate(*args,**kwargs)
#end def input_template


def multi_input_template(*args,**kwargs):
    return SimulationInputMultiTemplate(*args,**kwargs)
#end def input_template

