
import os
from numpy import array,dot,pi
from numpy.linalg import inv,norm
from generic import obj
from periodic_table import periodic_table
from physical_system import PhysicalSystem
from simulation import Simulation
from qmcpack_input import QmcpackInput,generate_qmcpack_input
from qmcpack_input import BundledQmcpackInput,TracedQmcpackInput
from qmcpack_input import QmcpackInputTemplate
from qmcpack_input import loop,linear,cslinear,vmc,dmc,collection,determinantset,hamiltonian,init,pairpot,bspline_builder
from qmcpack_input import generate_jastrows,generate_jastrow,generate_jastrow1,generate_jastrow2,generate_jastrow3
from qmcpack_input import generate_opt,generate_opts
from qmcpack_analyzer import QmcpackAnalyzer
from converters import Pw2qmcpack,Wfconvert
from convert4qmc import Convert4qmc
from sqd import Sqd
from debug import ci,ls,gs
from developer import unavailable
try:
    import h5py
except ImportError:
    h5py = unavailable('h5py')
#end try



class Qmcpack(Simulation):
    input_type    = QmcpackInput
    analyzer_type = QmcpackAnalyzer
    generic_identifier = 'qmcpack'
    infile_extension   = '.in.xml'
    #application   = 'qmcapp'
    application   = 'qmcapp_complex' # always use complex version until kpoint handling is fixed
    application_properties = set(['serial','omp','mpi'])
    application_results    = set(['jastrow'])
    preserve = Simulation.preserve | set(['should_twist_average'])


    def post_init(self):
        #jtk mark
        #  may need to put this back
        #  removed because particleset is not required by qmcpack
        #   and thus the system cannot be determined without access to the h5file
        #if self.system is None:
        #    self.error('system must be specified to determine type of run')
        ##end if
        if self.system is None:
            if not isinstance(self.input,QmcpackInputTemplate):
                self.warn('system must be specified to determine whether to twist average\n  proceeding under the assumption of no twist averaging')
            #end if
            self.should_twist_average = False
        else:
            self.system.group_atoms()
            self.system.change_units('B')
            twh = self.input.get_host('twist')
            tnh = self.input.get_host('twistnum')
            htypes = bspline_builder,determinantset
            user_twist_given  = isinstance(twh,htypes) and twh.twist!=None
            user_twist_given |= isinstance(tnh,htypes) and tnh.twistnum!=None
            many_kpoints = len(self.system.structure.kpoints)>1
            self.should_twist_average = many_kpoints and not user_twist_given
        #end if
    #end def post_init


    def propagate_identifier(self):
        if not isinstance(self.input,QmcpackInputTemplate):
            self.input.simulation.project.id = self.identifier
        #end if
    #end def propagate_identifier


    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='jastrow':
            calctypes = self.input.get_output_info('calctypes')
            calculating_result = 'opt' in calctypes
        else:
            self.error('ability to check for result '+result_name+' has not been implemented')
        #end if        
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        analyzer = self.load_analyzer_image()
        if result_name=='jastrow':
            if not 'results' in analyzer or not 'optimization' in analyzer.results:
                self.error('analyzer did not compute results required to determine jastrow')
            #end if
            opt_file = str(analyzer.results.optimization.optimal_file)
            result.opt_file = os.path.join(self.locdir,opt_file)
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        del analyzer
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        input = self.input
        system = self.system
        if result_name=='orbitals':
            if isinstance(sim,Pw2qmcpack) or isinstance(sim,Wfconvert):

                h5file = result.h5file

                wavefunction = input.get('wavefunction')
                if isinstance(wavefunction,collection):
                    wavefunction = wavefunction.get_single('psi0')
                #end if
                wf = wavefunction
                if 'sposet_builder' in wf and wf.sposet_builder.type=='bspline':
                    orb_elem = wf.sposet_builder
                elif 'sposet_builders' in wf and 'bspline' in wf.sposet_builders:
                    orb_elem = wf.sposet_builders.bspline
                elif 'sposet_builders' in wf and 'einspline' in wf.sposet_builders:
                    orb_elem = wf.sposet_builders.einspline
                elif 'determinantset' in wf and wf.determinantset.type in ('bspline','einspline'):
                    orb_elem = wf.determinantset
                else:
                    self.error('could not incorporate pw2qmcpack/wfconvert orbitals\n  bspline sposet_builder and determinantset are both missing')
                #end if
                orb_elem.href = os.path.relpath(h5file,self.locdir)
                if system.structure.folded_structure!=None:
                    orb_elem.tilematrix = array(system.structure.tmatrix)
                #end if
                defs = obj(
                    twistnum   = 0,
                    meshfactor = 1.0
                    )
                for var,val in defs.iteritems():
                    if not var in orb_elem:
                        orb_elem[var] = val
                    #end if
                #end for

                system = self.system
                structure = system.structure
                nkpoints = len(structure.kpoints)
                if nkpoints==0:
                    self.error('system must have kpoints to assign twistnums')
                #end if
                    
                if not os.path.exists(h5file):
                    self.error('wavefunction file not found:  \n'+h5file)
                #end if

                twistnums = range(len(structure.kpoints))
                if self.should_twist_average:
                    self.twist_average(twistnums)
                elif orb_elem.twistnum is None:
                    orb_elem.twistnum = twistnums[0]
                #end if

            elif isinstance(sim,Sqd):

                h5file  = os.path.join(result.dir,result.h5file)
                h5file  = os.path.relpath(h5file,self.locdir)

                sqdxml_loc = os.path.join(result.dir,result.qmcfile)
                sqdxml = QmcpackInput(sqdxml_loc)

                #sqd sometimes puts the wrong ionic charge
                #  rather than setting Z to the number of electrons
                #  set it to the actual atomic number
                g = sqdxml.qmcsystem.particlesets.atom.group
                elem = g.name
                if not elem in periodic_table.elements:
                    self.error(elem+' is not an element in the periodic table')
                #end if
                g.charge = periodic_table.elements[elem].atomic_number

                input = self.input
                s = input.simulation
                qsys_old = s.qmcsystem
                del s.qmcsystem
                s.qmcsystem = sqdxml.qmcsystem
                if 'jastrows' in qsys_old.wavefunction:
                    s.qmcsystem.wavefunction.jastrows = qsys_old.wavefunction.jastrows
                    for jastrow in s.qmcsystem.wavefunction.jastrows:
                        if 'type' in jastrow:
                            jtype = jastrow.type.lower().replace('-','_')
                            if jtype=='one_body':
                                jastrow.source = 'atom'
                            #end if
                        #end if
                    #end for
                #end if
                s.qmcsystem.hamiltonian = hamiltonian(
                    name='h0',type='generic',target='e',
                    pairpots = [
                        pairpot(name='ElecElec',type='coulomb',source='e',target='e'),
                        pairpot(name='Coulomb' ,type='coulomb',source='atom',target='e'),
                        ]
                    )
                s.init = init(source='atom',target='e')

                abset = input.get('atomicbasisset')
                abset.href = h5file

            elif isinstance(sim,Convert4qmc):

                res = QmcpackInput(result.location)
                wfn = res.qmcsystem.wavefunction
                qs = input.simulation.qmcsystem
                oldwfn = qs.wavefunction
                newwfn = wfn.copy()
                if 'jastrows' in newwfn:
                    del newwfn.jastrows
                #end if
                if 'jastrows' in oldwfn:
                    newwfn.jastrows = oldwfn.jastrows
                #end if
                qs.wavefunction = newwfn

            else:
                self.error('incorporating orbitals from '+sim.__class__.__name__+' has not been implemented')
            #end if
        elif result_name=='jastrow':
            if isinstance(sim,Qmcpack):
                opt_file = result.opt_file
                opt = QmcpackInput(opt_file)
                wavefunction = input.get('wavefunction')
                optwf = opt.qmcsystem.wavefunction
                def process_jastrow(wf):                
                    if 'jastrow' in wf:
                        js = [wf.jastrow]
                    elif 'jastrows' in wf:
                        js = wf.jastrows.values()
                    else:
                        js = []
                    #end if
                    jd = dict()
                    for j in js:
                        jtype = j.type.lower().replace('-','_').replace(' ','_')
                        jd[jtype] = j
                    #end for
                    return jd
                #end def process_jastrow
                if wavefunction==None:
                    qs = input.get('qmcsystem')
                    qs.wavefunction = optwf.copy()
                else:
                    jold = process_jastrow(wavefunction)
                    jopt = process_jastrow(optwf)
                    jnew = list(jopt.values())
                    for jtype in jold.keys():
                        if not jtype in jopt:
                            jnew.append(jold[jtype])
                        #end if
                    #end for
                    if len(jnew)==1:
                        wavefunction.jastrow = jnew[0].copy()
                    else:
                        wavefunction.jastrows = collection(jnew)
                    #end if
                #end if
                del optwf
            elif isinstance(sim,Sqd):
                wavefunction = input.get('wavefunction')
                jastrows = []
                if 'jastrows' in wavefunction:
                    for jastrow in wavefunction.jastrows:
                        jname = jastrow.name
                        if jname!='J1' and jname!='J2':
                            jastrows.append(jastrow)
                        #end if
                    #end for
                    del wavefunction.jastrows
                #end if

                ionps = input.get_ion_particlesets()
                if ionps is None or len(ionps)==0:
                    self.error('ion particleset does not seem to exist')
                elif len(ionps)==1:
                    ionps_name = list(ionps.keys())[0]
                else:
                    self.error('multiple ion species not supported for atomic calculations')
                #end if

                jastrows.extend([
                        generate_jastrow('J1','bspline',8,result.rcut,iname=ionps_name,system=self.system),
                        generate_jastrow('J2','pade',result.B)
                        ])

                wavefunction.jastrows = collection(jastrows)

            else:
                self.error('incorporating jastrow from '+sim.__class__.__name__+' has not been implemented')
            #end if
        elif result_name=='structure':
            structure = self.system.structure
            relstruct = result.structure
            structure.set(
                pos   = relstruct.positions,
                atoms = relstruct.atoms
                )
            self.input.incorporate_system(self.system)
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if        
    #end def incorporate_result


    def check_sim_status(self):
        outfile = os.path.join(self.locdir,self.outfile)
        errfile = os.path.join(self.locdir,self.errfile)
        fobj = open(outfile,'r')
        output = fobj.read()
        fobj.close()
        fobj = open(errfile,'r')
        errors = fobj.read()
        fobj.close()

        ran_to_end  = 'Total Execution' in output
        files_exist = True
        outfiles = self.input.get_output_info('outfiles')
        for file in outfiles:
            file_loc = os.path.join(self.locdir,file)
            files_exist = files_exist and os.path.exists(file_loc)
        #end for
            
        if ran_to_end and not files_exist:
            self.warn('run finished successfully, but output files do not seem to exist')
            print outfiles
            print os.listdir(self.locdir)
        #end if

        aborted = 'Fatal Error' in errors

        self.succeeded = ran_to_end
        self.failed    = aborted
        self.finished  = files_exist and (self.job.finished or ran_to_end) and not aborted 



        #print
        #print self.__class__.__name__
        #print 'identifier ',self.identifier
        #print 'ran_to_end ',ran_to_end
        #print 'files_exist',files_exist
        #print 'aborted    ',aborted
        #print 'job done   ',self.job.finished
        #print 'finished   ',self.finished
        #print

    #end def check_sim_status


    def get_output_files(self):
        if self.should_twist_average and not isinstance(self.input,TracedQmcpackInput):
            self.twist_average(range(len(self.system.structure.kpoints)))
            br = self.bundle_request
            input = self.input.trace(br.quantity,br.values)
            input.generate_filenames(self.infile)
            self.input = input
        #end if

        output_files = self.input.get_output_info('outfiles')

        return output_files
    #end def get_output_files


    def app_command(self):
        return self.app_name+' '+self.infile      
    #end def app_command


    def twist_average(self,twistnums):
        br = obj()
        br.quantity = 'twistnum'
        br.values   = list(twistnums)
        self.bundle_request = br
        #self.app_name = 'qmcapp_complex'
        #print 'twist_average'
        #print '  setting bundle request:'
        #print self.bundle_request
    #end def twist_average


    def write_prep(self):
        if self.got_dependencies:
            if 'bundle_request' in self and not isinstance(self.input,TracedQmcpackInput):
                br = self.bundle_request
                input = self.input.trace(br.quantity,br.values)
                input.generate_filenames(self.infile)
                if self.infile in self.files:
                    self.files.remove(self.infile)
                #end if
                for file in input.filenames:
                    self.files.add(file)
                #end for
                self.infile = input.filenames[-1]
                self.input  = input
                self.job.app_command = self.app_command()
            #end if
        #end if
    #end def write_prep
#end class Qmcpack




class BundledQmcpack(Qmcpack):
    infile_extension = '.in'
    application_results = set([])

    preserve = set(Simulation.preserve)
    preserve.add('sims')

    def __init__(self,**kwargs):
        if not 'sims' in kwargs:
            self.error('sims must be provided')
        #end if
        sims = kwargs['sims']
        self.sims = sims
        del kwargs['sims']
        files = set()
        for sim in sims:
            files = files | sim.files
        #end for
        kwargs['files'] = files

        inputs = []
        filenames = []
        for sim in sims:
            inputs.append(sim.input)
            filenames.append(sim.infile)
        #end for
        kwargs['input'] = BundledQmcpackInput(inputs=inputs,filenames=filenames)

        Simulation.__init__(self,**kwargs)
        deps = []
        for sim in sims:
            for dep in sim.dependencies:
                deps.append((dep.sim,'other'))
            #end for
        #end for
        self.depends(*deps)
    #end def __init__

    def propagate_identifier(self):
        for sim in self.sims:
            sim.propagate_identifier()
        #end for
    #end def propagate_identifier

    def check_result(self,result_name,sim):
        return False
    #end def check_result

    def get_result(self,result_name,sim):
        self.error(result_name+' is not calculated by BundledQmcpack')
    #end def get_result

    def check_dependencies(self,result):
        for sim in self.sims:
            sim.check_dependencies(results)
        #end for
        Simulation.check_dependencies(self,result)
    #end def check_dependencies

    def get_dependencies(self):
        for sim in self.sims:
            sim.get_dependencies()
        #end for
        Simulation.get_dependencies(self)
    #end def get_dependencies



    def check_sim_status(self):
        outfile = os.path.join(self.locdir,self.outfile)
        errfile = os.path.join(self.locdir,self.errfile)
        fobj = open(outfile,'r')
        output = fobj.read()
        fobj.close()
        fobj = open(errfile,'r')
        errors = fobj.read()
        fobj.close()

        ran_to_end  = 'Total Execution' in output
        files_exist = True
        outfiles = self.input.get_output_info('outfiles')
        for file in outfiles:
            file_loc = os.path.join(self.locdir,file)
            files_exist = files_exist and os.path.exists(file_loc)
        #end for
            
        if ran_to_end and not files_exist:
            self.warn('run finished successfully, but output files do not seem to exist')
            print outfiles
            print os.listdir(self.locdir)
        #end if

        aborted = 'Fatal Error' in errors

        self.failed   = aborted
        self.finished = files_exist and self.job.finished and not aborted 

        #print
        #print self.__class__.__name__
        #print 'identifier ',self.identifier
        #print 'ran_to_end ',ran_to_end
        #print 'files_exist',files_exist
        #print 'aborted    ',aborted
        #print 'job done   ',self.job.finished
        #print 'finished   ',self.finished
        #print
        #
        #import code
        #code.interact(local=dict(locals(),**globals()))

    #end def check_sim_status

#end class BundledQmcpack




def generate_qmcpack(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)

    if not 'input' in sim_args:
        input_type = inp_args.input_type
        del inp_args.input_type
        sim_args.input = generate_qmcpack_input(input_type,**inp_args)
    #end if
    qmcpack = Qmcpack(**sim_args)

    return qmcpack
#end def generate_qmcpack
