##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack.py                                                        #
#    Nexus interface with the QMCPACK simulation code.               #
#                                                                    #
#                                                                    #
#  Content summary:                                                  #
#    Qmcpack                                                         #
#      Simulation class for QMCPACK.                                 #
#      Handles incorporation of structure, orbital, and Jastrow      #
#        data from other completed simulations.                      #
#                                                                    #
#    generate_qmcpack                                                #
#      User-facing function to create QMCPACK simulation objects.    #
#                                                                    #
#    generate_cusp_correction                                        #
#      User-facing function to run QMCPACK as an intermediate tool   #
#        to add cusps to Gaussian orbitals coming from GAMESS.       #
#                                                                    #
#====================================================================#


import os
from numpy import array,dot,pi
from numpy.linalg import inv,norm
from generic import obj
from periodic_table import periodic_table
from physical_system import PhysicalSystem
from simulation import Simulation
from qmcpack_input import QmcpackInput,generate_qmcpack_input
from qmcpack_input import TracedQmcpackInput
from qmcpack_input import loop,linear,cslinear,vmc,dmc,collection,determinantset,hamiltonian,init,pairpot,bspline_builder
from qmcpack_input import generate_jastrows,generate_jastrow,generate_jastrow1,generate_jastrow2,generate_jastrow3
from qmcpack_input import generate_opt,generate_opts
from qmcpack_analyzer import QmcpackAnalyzer
from qmcpack_converters import Pw2qmcpack,Wfconvert,Convert4qmc
from sqd import Sqd
from debug import ci,ls,gs
from developer import unavailable
from nexus_base import nexus_core
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
    application   = 'qmcpack'
    application_properties = set(['serial','omp','mpi'])
    application_results    = set(['jastrow','cuspcorr','wavefunction'])


    def post_init(self):
        generic_input = self.has_generic_input()

        if self.system is None:
            if not generic_input:
                self.warn('system must be specified to determine whether to twist average\nproceeding under the assumption of no twist averaging')
            #end if
            self.should_twist_average = False
        else:
            if generic_input:
                cls = self.__class__
                self.error('cannot twist average generic or templated input\nplease provide {0} instead of {1} for input'.format(cls.input_type.__class__.__name__,self.input.__class__.__name__))
            #end if
            self.system.group_atoms()
            self.system.change_units('B')
            twh = self.input.get_host('twist')
            tnh = self.input.get_host('twistnum')
            htypes = bspline_builder,determinantset
            user_twist_given  = isinstance(twh,htypes) and twh.twist!=None
            user_twist_given |= isinstance(tnh,htypes) and tnh.twistnum!=None
            many_kpoints = len(self.system.structure.kpoints)>1
            self.should_twist_average = many_kpoints and not user_twist_given
            if self.should_twist_average:
                # correct the job app command to account for the change in input file name
                # this is necessary for twist averaged runs in bundles
                app_comm = self.app_command()
                prefix,ext = self.infile.split('.',1)
                self.infile = prefix+'.in'
                app_comm_new = self.app_command()
                if self.job.app_command==app_comm:
                    self.job.app_command=app_comm_new
                #end if
            #end if
        #end if
    #end def post_init


    def propagate_identifier(self):
        if not self.has_generic_input():
            self.input.simulation.project.id = self.identifier
        #end if
    #end def propagate_identifier


    def pre_write_inputs(self,save_image):
        # fix to make twist averaged input file under generate_only
        if self.system is None:
            self.should_twist_average = False
        elif nexus_core.generate_only:
            twistnums = range(len(self.system.structure.kpoints))
            if self.should_twist_average:
                self.twist_average(twistnums)
            #end if
        #end if
    #end def pre_write_inputs


    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='jastrow' or result_name=='wavefunction':
            calctypes = self.input.get_output_info('calctypes')
            calculating_result = 'opt' in calctypes
        elif result_name=='cuspcorr':
            calculating_result = self.input.cusp_correction()
        else:
            self.error('ability to check for result '+result_name+' has not been implemented')
        #end if        
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        analyzer = self.load_analyzer_image()
        if result_name=='jastrow' or result_name=='wavefunction':
            if not 'results' in analyzer or not 'optimization' in analyzer.results:
                self.error('analyzer did not compute results required to determine jastrow')
            #end if
            opt_file = analyzer.results.optimization.optimal_file
            opt_file = str(opt_file)
            result.opt_file = os.path.join(self.locdir,opt_file)
        elif result_name=='cuspcorr':
            result.spo_up_cusps = os.path.join(self.locdir,self.identifier+'.spo-up.cuspInfo.xml')
            result.spo_dn_cusps = os.path.join(self.locdir,self.identifier+'.spo-dn.cuspInfo.xml')
            result.updet_cusps = os.path.join(self.locdir,'updet.cuspInfo.xml')
            result.dndet_cusps = os.path.join(self.locdir,'downdet.cuspInfo.xml')
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
                if 'href' in orb_elem and isinstance(orb_elem.href,str) and os.path.exists(orb_elem.href):
                    # user specified h5 file for orbitals, bypass orbital dependency
                    orb_elem.href = os.path.relpath(orb_elem.href,self.locdir)
                else:
                    orb_elem.href = os.path.relpath(h5file,self.locdir)
                    if system.structure.folded_structure!=None:
                        orb_elem.tilematrix = array(system.structure.tmatrix)
                    #end if
                #end if
                defs = obj(
                    #twistnum   = 0,
                    meshfactor = 1.0
                    )
                for var,val in defs.iteritems():
                    if not var in orb_elem:
                        orb_elem[var] = val
                    #end if
                #end for
                has_twist    = 'twist' in orb_elem
                has_twistnum = 'twistnum' in orb_elem
                if  not has_twist and not has_twistnum:
                    orb_elem.twistnum = 0
                #end if

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
                elif not has_twist and orb_elem.twistnum is None:
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
                qs  = input.simulation.qmcsystem
                oldwfn = qs.wavefunction
                newwfn = res.qmcsystem.wavefunction
                if 'jastrows' in newwfn:
                    del newwfn.jastrows
                #end if
                if 'jastrows' in oldwfn:
                    newwfn.jastrows = oldwfn.jastrows
                #end if
                if input.cusp_correction():
                    newwfn.determinantset.cuspcorrection = True
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
        elif result_name=='particles':
            if isinstance(sim,Convert4qmc):
                ptcl_file = result.location
                qi = QmcpackInput(ptcl_file)
                self.input.simulation.qmcsystem.particlesets = qi.qmcsystem.particlesets
            else:
                self.error('incorporating particles from '+sim.__class__.__name__+' has not been implemented')
            # end if
        elif result_name=='structure':
            relstruct = result.structure.copy()
            relstruct.change_units('B')
            self.system.structure = relstruct
            self.system.remove_folded()
            self.input.incorporate_system(self.system)

        elif result_name=='cuspcorr':

            ds = self.input.get('determinantset')
            ds.cuspcorrection = True
            try: # multideterminant
              ds.sposets['spo-up'].cuspinfo = os.path.relpath(result.spo_up_cusps,self.locdir)
              ds.sposets['spo-dn'].cuspinfo = os.path.relpath(result.spo_dn_cusps,self.locdir)
            except: # single determinant
              sd = ds.slaterdeterminant
              sd.determinants['updet'].cuspinfo = os.path.relpath(result.updet_cusps,self.locdir)
              sd.determinants['downdet'].cuspinfo = os.path.relpath(result.dndet_cusps,self.locdir)
            # end try

        elif result_name=='wavefunction':
            if not isinstance(sim,Qmcpack):
                self.error('incorporating wavefunction from '+sim.__class__.__name__+' has not been implemented')
            #end if
            print '        getting optimal wavefunction from: '+result.opt_file
            opt = QmcpackInput(result.opt_file)
            qs = input.get('qmcsystem')
            qs.wavefunction = opt.qmcsystem.wavefunction.copy()
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if        
    #end def incorporate_result


    def check_sim_status(self):
        output = self.outfile_text()
        errors = self.errfile_text()

        ran_to_end  = 'Total Execution' in output
        aborted     = 'Fatal Error' in errors
        files_exist = True
        cusp_run    = False

        if not self.has_generic_input():
            if not isinstance(self.input,TracedQmcpackInput):
                cusp_run = self.input.cusp_correction()
            #end if
            if cusp_run:
                sd = self.input.get('slaterdeterminant')
                if sd!=None:
                    cuspfiles = []
                    for d in sd.determinants:
                        cuspfiles.append(d.id+'.cuspInfo.xml')
                    #end for
                else: # assume multideterminant sposet names
                    cuspfiles = ['spo-up.cuspInfo.xml','spo-dn.cuspInfo.xml']
                #end if
                outfiles   = cuspfiles
            else:
                outfiles = self.input.get_output_info('outfiles')
            #end if

            for file in outfiles:
                file_loc = os.path.join(self.locdir,file)
                files_exist = files_exist and os.path.exists(file_loc)
            #end for

            if ran_to_end and not files_exist:
                self.warn('run finished successfully, but output files do not seem to exist')
                print outfiles
                print os.listdir(self.locdir)
            #end if
        #end if


        self.succeeded = ran_to_end
        self.failed    = aborted
        self.finished  = files_exist and (self.job.finished or ran_to_end) and not aborted 

        if cusp_run and files_exist:
            for cuspfile in cuspfiles:
                cf_orig = os.path.join(self.locdir,cuspfile)
                cf_new  = os.path.join(self.locdir,self.identifier+'.'+cuspfile)
                os.system('cp {0} {1}'.format(cf_orig,cf_new))
            #end for
        #end if
    #end def check_sim_status


    def get_output_files(self):
        if self.has_generic_input():
            output_files = []
        else:
            if self.should_twist_average and not isinstance(self.input,TracedQmcpackInput):
                self.twist_average(range(len(self.system.structure.kpoints)))
                br = self.bundle_request
                input = self.input.trace(br.quantity,br.values)
                input.generate_filenames(self.infile)
                self.input = input
            #end if
            output_files = self.input.get_output_info('outfiles')
        #end if
        return output_files
    #end def get_output_files

    
    def post_analyze(self,analyzer):
        calctypes = self.input.get_output_info('calctypes')
        opt_run = calctypes!=None and 'opt' in calctypes
        if opt_run:
            opt_file = analyzer.results.optimization.optimal_file
            if opt_file is None:
                self.failed = True
            #end if
        #end if
    #end def post_analyze


    def app_command(self):
        return self.app_name+' '+self.infile      
    #end def app_command


    def twist_average(self,twistnums):
        br = obj()
        br.quantity = 'twistnum'
        br.values   = list(twistnums)
        self.bundle_request = br
    #end def twist_average


    def write_prep(self):
        if self.got_dependencies:
            traced_input  = isinstance(self.input,TracedQmcpackInput)
            generic_input = self.has_generic_input()
            if 'bundle_request' in self and not traced_input and not generic_input:
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



def generate_qmcpack(**kwargs):
    sim_args,inp_args = Qmcpack.separate_inputs(kwargs)

    if not 'input' in sim_args:
        input_type = inp_args.input_type
        del inp_args.input_type
        sim_args.input = generate_qmcpack_input(input_type,**inp_args)
    #end if
    qmcpack = Qmcpack(**sim_args)

    return qmcpack
#end def generate_qmcpack


def generate_cusp_correction(**kwargs):
    kwargs['input_type']   = 'basic'
    kwargs['bconds']       = 'nnn'
    kwargs['jastrows']     = []
    kwargs['corrections']  = []
    kwargs['calculations'] = []

    sim_args,inp_args = Simulation.separate_inputs(kwargs)

    input_type = inp_args.input_type
    del inp_args.input_type
    input = generate_qmcpack_input(input_type,**inp_args)

    wf = input.get('wavefunction')
    if not 'determinantset' in wf:
        Qmcpack.class_error('wavefunction does not have determinantset, cannot create cusp correction','generate_cusp_correction')
    #end if
    wf.determinantset.cuspcorrection = True

    sim_args.input = input
    qmcpack = Qmcpack(**sim_args)

    return qmcpack
#end def generate_cusp_correction
