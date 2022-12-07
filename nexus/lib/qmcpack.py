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
import numpy as np
from numpy import array,dot,pi
from numpy.linalg import inv,norm
from generic import obj
from periodic_table import periodic_table
from physical_system import PhysicalSystem
from simulation import Simulation,NullSimulationAnalyzer
from qmcpack_input import QmcpackInput,generate_qmcpack_input
from qmcpack_input import TracedQmcpackInput
from qmcpack_input import loop,linear,cslinear,vmc,dmc,collection,determinantset,hamiltonian,init,pairpot,bspline_builder
from qmcpack_input import generate_jastrows,generate_jastrow,generate_jastrow1,generate_jastrow2,generate_jastrow3
from qmcpack_input import generate_opt,generate_opts
from qmcpack_input import check_excitation_type
from qmcpack_analyzer import QmcpackAnalyzer
from qmcpack_converters import Pw2qmcpack,Convert4qmc,PyscfToAfqmc
from debug import ci,ls,gs
from developer import unavailable
from nexus_base import nexus_core
from copy import deepcopy
try:
    import h5py
except:
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


    def has_afqmc_input(self):
        afqmc_input = False
        if not self.has_generic_input():
            afqmc_input = self.input.is_afqmc_input()
        #end if
        return afqmc_input
    #end def has_afqmc_input


    def post_init(self):
        generic_input = self.has_generic_input()

        if self.has_afqmc_input():
            self.analyzer_type = NullSimulationAnalyzer
            self.should_twist_average = False
        elif self.system is None:
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
            twistnums = list(range(len(self.system.structure.kpoints)))
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
        #end if        
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        result = obj()
        if result_name=='jastrow' or result_name=='wavefunction':
            analyzer = self.load_analyzer_image()
            if not 'results' in analyzer or not 'optimization' in analyzer.results:
                if self.should_twist_average:
                    self.error('Wavefunction optimization was performed for each twist separately.\nCurrently, the transfer of per-twist wavefunction parameters from\none QMCPACK simulation to another is not supported.  Please either\nredo the optimization with a single twist (see "twist" or "twistnum"\noptions), or request that this feature be implemented.')
                else:
                    self.error('analyzer did not compute results required to determine jastrow')
                #end if
            #end if
            opt_file = analyzer.results.optimization.optimal_file
            opt_file = str(opt_file)
            result.opt_file = os.path.join(self.locdir,opt_file)
            del analyzer
        elif result_name=='cuspcorr':
            result.spo_up_cusps = os.path.join(self.locdir,self.identifier+'.spo-up.cuspInfo.xml')
            result.spo_dn_cusps = os.path.join(self.locdir,self.identifier+'.spo-dn.cuspInfo.xml')
            result.updet_cusps = os.path.join(self.locdir,'updet.cuspInfo.xml')
            result.dndet_cusps = os.path.join(self.locdir,'downdet.cuspInfo.xml')
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result


    def incorporate_result(self,result_name,result,sim):
        input = self.input
        system = self.system
        if result_name=='orbitals':
            if isinstance(sim,Pw2qmcpack):

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
                    self.error('could not incorporate pw2qmcpack orbitals\nbspline sposet_builder and determinantset are both missing')
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
                for var,val in defs.items():
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
                    self.error('wavefunction file not found:\n'+h5file)
                #end if

                twistnums = list(range(len(structure.kpoints)))
                if self.should_twist_average:
                    self.twist_average(twistnums)
                elif not has_twist and orb_elem.twistnum is None:
                    orb_elem.twistnum = twistnums[0]
                #end if

            elif isinstance(sim,Convert4qmc):

                res = QmcpackInput(result.location)
                qs  = input.simulation.qmcsystem
                oldwfn = qs.wavefunction
                newwfn = res.qmcsystem.wavefunction
                if hasattr(oldwfn.determinantset,'multideterminant'):
                    del newwfn.determinantset.slaterdeterminant
                    newwfn.determinantset.multideterminant = oldwfn.determinantset.multideterminant
                    newwfn.determinantset.sposets = oldwfn.determinantset.sposets
                dset = newwfn.determinantset

                if 'jastrows' in newwfn:
                    del newwfn.jastrows
                #end if
                if 'jastrows' in oldwfn:
                    newwfn.jastrows = oldwfn.jastrows
                #end if
                if input.cusp_correction():
                    dset.cuspcorrection = True
                #end if
                if 'orbfile' in result:
                    orb_h5file = result.orbfile
                    if not os.path.exists(orb_h5file) and 'href' in dset:
                        orb_h5file = os.path.join(sim.locdir,dset.href)
                    #end if
                    if not os.path.exists(orb_h5file):
                        self.error('orbital h5 file from convert4qmc does not exist\nlocation checked: {}'.format(orb_h5file))
                    #end if
                    orb_path = os.path.relpath(orb_h5file,self.locdir)
                    dset.href = orb_path
                    detlist = dset.get('detlist')
                    if detlist is not None and 'href' in detlist:
                        detlist.href = orb_path
                    #end if
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
                        js = list(wf.jastrows.values())
                    else:
                        js = []
                    #end if
                    jd = dict()
                    for j in js:
                        jtype = j.type.lower().replace('-','_').replace(' ','_')
                        key = jtype
                        # take care of multiple jastrows of the same type
                        if key in jd:  # use name to distinguish
                            key += j.name
                            if key in jd:  # if still duplicate then error out
                                msg = 'duplicate jastrow in '+self.__class__.__name__
                                self.error(msg)
                            #end if
                        #end if
                        jd[key] = j
                    #end for
                    return jd
                #end def process_jastrow
                if wavefunction is None:
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
            #end try

        elif result_name=='wavefunction':
            if isinstance(sim,Qmcpack):
                opt = QmcpackInput(result.opt_file)
                qs = input.get('qmcsystem')
                qs.wavefunction = opt.qmcsystem.wavefunction.copy()
            elif isinstance(sim,PyscfToAfqmc):
                if not self.input.is_afqmc_input():
                    self.error('incorporating wavefunction from {} is only supported for AFQMC calculations'.format(sim.__class__.__name__))
                #end if
                h5_file =  os.path.relpath(result.h5_file,self.locdir)
                wfn = self.input.simulation.wavefunction
                ham = self.input.simulation.hamiltonian
                wfn.filename = h5_file
                wfn.filetype = 'hdf5'
                if 'filename' not in ham or ham.filename=='MISSING.h5':
                    ham.filename = h5_file
                    ham.filetype = 'hdf5'
                #end if
                if 'xml' in result:
                    xml = QmcpackInput(result.xml)
                    info_new = xml.simulation.afqmcinfo.copy()
                    info = self.input.simulation.afqmcinfo
                    info.set_optional(**info_new)
                    # override particular inputs set by default
                    if 'generation_info' in input._metadata:
                        g = input._metadata.generation_info
                        if 'walker_type' not in g:
                            walker_type = xml.get('walker_type')
                            walkerset = input.get('walkerset')
                            if walker_type is not None and walkerset is not None:
                                walkerset.walker_type = walker_type
                            #end if
                        #end if
                    #end if
                #end if
            else:
                self.error('incorporating wavefunction from '+sim.__class__.__name__+' has not been implemented')
            #end if
        elif result_name=='gc_occupation':
            from pwscf import Pwscf
            from qmcpack_converters import gcta_occupation
            if not isinstance(sim,Pw2qmcpack):
                msg = 'grand-canonical occupation require Pw2qmcpack'
                self.error(msg)
            #endif
            # step 1: extract Fermi energy for each spin from nscf
            nscf = None
            npwdep = 0
            for dep in sim.dependencies:
                if isinstance(dep.sim,Pwscf):
                    nscf = dep.sim
                    npwdep += 1
            if npwdep != 1:
                msg = 'need exactly 1 scf/nscf calculation for Fermi energy'
                msg += '\n found %d' % npwdep
                self.error(msg)
            #end if
            na = nscf.load_analyzer_image()
            Ef_list = na.fermi_energies
            # step 2: analyze ESH5 file for states below Fermi energy
            pa = sim.load_analyzer_image()
            if 'wfh5' not in pa:
              pa.analyze(Ef_list=Ef_list)
              sim.save_analyzer_image(pa)
            #end if
            # step 3: count the number of up/dn electrons at each supertwist
            s1 = self.system.structure
            ntwist = len(s1.kpoints)
            nelecs_at_twist = gcta_occupation(pa.wfh5, ntwist)
            self.nelecs_at_twist = nelecs_at_twist
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
                self.warn('run finished successfully, but output files do not exist')
                self.log(outfiles)
                self.log(os.listdir(self.locdir))
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
                self.twist_average(list(range(len(self.system.structure.kpoints))))
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
        if not self.has_generic_input():
            calctypes = self.input.get_output_info('calctypes')
            opt_run = calctypes!=None and 'opt' in calctypes
            if opt_run:
                opt_file = analyzer.results.optimization.optimal_file
                if opt_file is None:
                    self.failed = True
                #end if
            #end if
            exc_run = 'excitation' in self
            if exc_run:
                exc_failure = False

                edata = self.read_einspline_dat()
                exc_input = self.excitation

                exc_spin,exc_type,exc_spins,exc_types,exc1,exc2 = check_excitation_type(exc_input)

                elns = self.input.get_electron_particle_set()
                
                if exc_type==exc_types.band: 
                    # Band Index 'tw1 band1 tw2 band2'. Eg., '0 45 3 46'
                    # Check that tw1,band1 is no longer in occupied set
                    tw1,bnd1 = exc2.split()[0:2]
                    tw2,bnd2 = exc2.split()[2:4]
                    if exc1 in ('up','down'):
                        spin_channel = exc1
                        dsc = edata[spin_channel]
                        for idx,(tw,bnd) in enumerate(zip(dsc.TwistIndex,dsc.BandIndex)):
                            if tw == int(tw1) and bnd == int(bnd1):
                                # This orbital should no longer be in the set of occupied orbitals
                                if idx<elns.groups[spin_channel[0]].size:
                                    msg  = 'WARNING: You requested \'{}\' excitation of type \'{}\',\n'
                                    msg += '         however, the first orbital \'{} {}\' is still occupied (see einspline file).\n'
                                    msg += '         Please check your input.'
                                    msg = msg.format(spin_channel,exc_input[1],tw1,bnd1)
                                    exc_failure = True
                                #end if
                            elif tw == int(tw2) and bnd == int(bnd2):
                                # This orbital should be in the set of occupied orbitals
                                if idx>=elns.groups[spin_channel[0]].size:
                                    msg  = 'WARNING: You requested \'{}\' excitation of type \'{}\',\n'
                                    msg += '         however, the second orbital \'{} {}\' is not occupied (see einspline file).\n'
                                    msg += '         Please check your input.'
                                    msg = msg.format(spin_channel,exc_input[1],tw2,bnd2)
                                    exc_failure = True
                                #end if
                            #end if
                        #end for
                    else:
                        self.warn('No check for \'{}\' excitation of type \'{}\' was done. When this path is possible, then a check should be written.'.format(exc_input[0],exc_input[1]))
                    #end if
                elif exc_type in (exc_types.energy,exc_types.lowest):
                    # Lowest or Energy Index '-orbindex1 +orbindex2'. Eg., '-4 +5'
                    if exc_type==exc_types.lowest:
                        if exc_spin==exc_spins.down:
                            orb1 = elns.groups.d.size
                        else:
                            orb1 = elns.groups.u.size
                        #end if
                        orb2 = orb1+1 
                    else:
                        orb1 = int(exc_input[1].split()[0][1:])
                        orb2 = int(exc_input[1].split()[1][1:])
                    #end if
                    if exc1 in ('up','down'):

                        spin_channel = exc1
                        nelec = elns.groups[spin_channel[0]].size
                        eigs_spin = edata[spin_channel].Energy

                        # Construct the correct set of occupied orbitals by hand based on
                        # orb1 and orb2 values that were input by the user
                        excited = eigs_spin
                        order = eigs_spin.argsort()
                        ground = excited[order]
                        # einspline orbital ordering for excited state
                        excited = excited[:nelec]
                        # hand-crafted orbital order for excited state

                        # ground can be list or ndarray, but we'll convert it to list
                        # so we can concatenate with list syntax
                        ground = np.asarray(ground).tolist()
                        # After concatenating, convert back to ndarray
                        hc_excited = np.array(ground[:orb1-1]+[ground[orb2-1]]+ground[orb1:nelec])
                            
                        etol = 1e-6
                        if np.abs(hc_excited-excited).max() > etol:
                            msg  = 'WARNING: You requested \'{}\' excitation of type \'{}\',\n'
                            msg += '         however, the second orbital \'{}\' is not occupied (see einspline file).\n'
                            msg += '         Please check your input.'
                            msg = msg.format(spin_channel,exc_input[1],orb1)
                            exc_failure = True
                        #end if

                    elif exc1 in ('singlet','triplet'):
                        wf = self.input.get('wavefunction')
                        occ = wf.determinantset.multideterminant.detlist.csf.occ
                        if occ[int(orb1)-1]!='1':
                            msg  = 'WARNING: You requested \'{}\' excitation of type \'{}\',\n'
                            msg += '         however, this is inconsistent with the occupations in detlist \'{}\'.\n'
                            msg += '         Please check your input.'
                            msg = msg.format(spin_channel,exc_input[1],occ)
                            exc_failure = True
                        #end if
                        if occ[int(orb2)-1]!='1':
                            msg  = 'WARNING: You requested \'{}\' excitation of type \'{}\',\n'
                            msg += '         however, this is inconsistent with the occupations in detlist \'{}\'.\n'
                            msg += '         Please check your input.'
                            msg = msg.format(spin_channel,exc_input[1],occ)
                            exc_failure = True
                        #end if
                    #end if

                else:
                    # The format is: 'gamma vb z cb'
                    if exc1 in ('singlet','triplet'):
                        self.warn('No check for \'{}\' excitation of type \'{}\' was done. When this path is possible, then a check should be written.'.format(exc_input[0],exc_input[1]))
                    else:

                        # assume excitation of form 'gamma vb k cb' or 'gamma vb-1 k cb+1'
                        excitation = exc2.upper().split(' ')
                        k_1, band_1, k_2, band_2 = excitation
                        tilematrix = self.system.structure.tilematrix()
                        
                        wf = self.input.get('wavefunction')
                        if exc_spin==exc_spins.up:
                            sdet =  wf.determinantset.get('updet')
                        else:
                            sdet =  wf.determinantset.get('downdet')
                        #end if
                        from numpy import linalg,where,isclose
                        vb = int(sdet.size / abs(linalg.det(tilematrix))) -1  # Separate for each spin channel
                        cb = vb+1
                        # Convert band_1, band_2 to band indexes
                        bands = [band_1, band_2]
                        for bnum, b in enumerate(bands):
                            b = b.lower()
                            if 'cb' in b:
                                if '-' in b:
                                    b = b.split('-')
                                    bands[bnum] = cb - int(b[1])
                                elif '+' in b:
                                    b = b.split('+')
                                    bands[bnum] = cb + int(b[1])
                                else:
                                    bands[bnum] = cb
                                #end if
                            elif 'vb' in b:
                                if '-' in b:
                                    b = b.split('-')
                                    bands[bnum] = vb - int(b[1])
                                elif '+' in b:
                                    b = b.split('+')
                                    bands[bnum] = vb + int(b[1])
                                else:
                                    bands[bnum] = vb
                                #end if
                            else:
                                QmcpackInput.class_error('{0} in excitation has the wrong formatting'.format(b))
                            #end if
                        #end for
                        band_1, band_2 = bands
                        
                        # Convert k_1 k_2 to wavevector indexes
                        structure = self.system.structure.get_smallest().copy()
                        structure.change_units('A')

                        from structure import get_kpath
                        kpath       = get_kpath(structure=structure)
                        kpath_label = array(kpath['explicit_kpoints_labels'])
                        kpath_rel   = kpath['explicit_kpoints_rel']
                        
                        k1_in = k_1
                        k2_in = k_2
                        if k_1 in kpath_label and k_2 in kpath_label:   
                            k_1 = kpath_rel[where(kpath_label == k_1)][0]
                            k_2 = kpath_rel[where(kpath_label == k_2)][0]

                            kpts = structure.kpoints_unit()
                            found_k1 = False
                            found_k2 = False
                            for knum, k in enumerate(kpts):
                                if isclose(k_1, k).all():
                                    k_1 = knum
                                    found_k1 = True
                                #end if
                                if isclose(k_2, k).all():
                                    k_2 = knum
                                    found_k2 = True
                                #end if
                            #end for
                            if not found_k1 or not found_k2:
                                QmcpackInput.class_error('Requested special kpoint is not in the tiled cell\nRequested "{}", present={}\nRequested "{}", present={}\nAvailable kpoints: {}'.format(k1_in,found_k1,k2_in,found_k2,sorted(set(kpath_label))))
                            #end if
                        else:
                            QmcpackInput.class_error('Excitation wavevectors are not found in the kpath\nlabels requested: {} {}\nlabels present: {}'.format(k_1,k_2,sorted(set(kpath_label))))
                        #end if

                        tw1,bnd1 = (k_1,band_1)
                        tw2,bnd2 = (k_2,band_2)
                        spin_channel = exc1
                        dsc = edata[spin_channel]
                        for idx,(tw,bnd) in enumerate(zip(dsc.TwistIndex,dsc.BandIndex)):
                            if tw == int(tw1) and bnd == int(bnd1):
                                # This orbital should no longer be in the set of occupied orbitals
                                if idx<elns.groups[spin_channel[0]].size:
                                    msg  = 'WARNING: You requested \'{}\' excitation of type \'{}\',\n'
                                    msg += '         however, the first orbital \'{} {}\' is still occupied (see einspline file).\n'
                                    msg += '         Please check your input.'
                                    msg = msg.format(spin_channel,exc_input[1],tw1,bnd1)
                                    exc_failure = True
                                #end if
                            elif tw == int(tw2) and bnd == int(bnd2):
                                # This orbital should be in the set of occupied orbitals
                                if idx>=elns.groups[spin_channel[0]].size:
                                    msg  = 'WARNING: You requested \'{}\' excitation of type \'{}\',\n'
                                    msg += '         however, the second orbital \'{} {}\' is not occupied (see einspline file).\n'
                                    msg += '         Please check your input.'
                                    msg = msg.format(spin_channel,exc_input[1],tw2,bnd2)
                                    exc_failure = True
                                #end if
                            #end if
                        #end for

                #end if

                if exc_failure:
                    self.failed = True
                    self.warn(msg)
                    filename = self.identifier+'_errors.txt'
                    open(os.path.join(self.locdir,filename),'w').write(msg)
                #end if

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
                # write twist info files
                s = self.system.structure
                kweights        = s.kweights.copy()
                kpoints         = s.kpoints.copy()
                kpoints_qmcpack = s.kpoints_qmcpack()
                for file in input.filenames:
                    if file.startswith(self.identifier+'.g'):
                        tokens = file.split('.')
                        twist_index = int(tokens[1].replace('g',''))
                        twist_filename = '{}.{}.twist_info.dat'.format(tokens[0],tokens[1])
                        kw  = kweights[twist_index]
                        kp  = kpoints[twist_index]
                        kpq = kpoints_qmcpack[twist_index]
                        contents = ' {: 16.6f}  {: 16.12f} {: 16.12f} {: 16.12f}  {: 16.12f} {: 16.12f} {: 16.12f}\n'.format(kw,*kp,*kpq)
                        fobj = open(os.path.join(self.locdir,twist_filename),'w')
                        fobj.write(contents)
                        fobj.close()
                    #end if
                #end for
                grand_canonical_twist_average = 'nelecs_at_twist' in self
                if grand_canonical_twist_average:
                    for itwist, qi in enumerate(input.inputs):
                        elecs = self.nelecs_at_twist[itwist]
                        # step 1: resize particlesets
                        nup = elecs[0]
                        ndn = elecs[1]
                        qi.get('u').set(size=nup)
                        qi.get('d').set(size=ndn)
                        # step 2: resize determinants
                        dset = qi.get('determinantset')
                        sdet = dset.slaterdeterminant  # hard-code single det
                        spo_size_map = {}
                        for det in sdet.determinants:
                            nelec = None  # determine from group
                            group = det.get('group')
                            if group == 'u':
                                nelec = nup
                            elif group == 'd':
                                nelec = ndn
                            else:
                                msg = 'need to count number of "%s"' % group
                                self.error(msg)
                            #end if
                            spo_name = det.get('sposet')
                            spo_size_map[spo_name] = nelec
                            det.set(size=nelec)
                        #end for
                        # step 3: resize orbital sets
                        sb = qi.get('sposet_builder')
                        bb = sb.bspline  # hard-code for Bspline orbs
                        assert itwist == bb.twistnum
                        sposets = bb.sposets
                        for spo in sposets:
                            if spo.name in spo_size_map:
                                spo.set(size=spo_size_map[spo.name])
                            #end if
                        #end for
                    #end for
                #end if
            #end if
        #end if
    #end def write_prep

    def read_einspline_dat(self):
        edata = obj()
        import glob
        for einpath in glob.glob(self.locdir+'/einsplin*'):
            ftokens = einpath.split('.')
            fspin = int(ftokens[-5][5])
            if fspin==0:
                spinlab = 'up'
            else:
                spinlab = 'down'
            #end if
            edata[spinlab] = obj()
            with open(einpath) as f:
                data = array(f.read().split()[1:])
                data.shape = len(data)//12,12
                data = data.T
                for darr in data:
                    if darr[0][0]=='K' or darr[0][0]=='E':
                        edata[spinlab][darr[0]] = array(list(map(float,darr[1:])))
                    else:
                        edata[spinlab][darr[0]] = array(list(map(int,darr[1:])))
                    #end if
                #end for
            #end with
        #end for
        return edata
    #end def read_einspline_dat
#end class Qmcpack



def generate_qmcpack(**kwargs):
    sim_args,inp_args = Qmcpack.separate_inputs(kwargs)

    exc = None
    if 'excitation' in inp_args:
        exc = deepcopy(inp_args.excitation)
    #end if

    spp = None
    if 'spin_polarized' in inp_args:
        spp = deepcopy(inp_args.spin_polarized)
    #end if

    if 'input' not in sim_args:
        run_path = None
        if 'path' in sim_args:
            run_path = os.path.join(nexus_core.local_directory,nexus_core.runs,sim_args.path)
        #end if
        inp_args.run_path = run_path
        sim_args.input = generate_qmcpack_input(**inp_args)
    #end if
    qmcpack = Qmcpack(**sim_args)

    if exc is not None:
        qmcpack.excitation = exc
    #end if

    if spp is not None:
        qmcpack.spin_polarized = spp
    #end if

    return qmcpack
#end def generate_qmcpack


def generate_cusp_correction(**kwargs):
    kwargs['input_type']   = 'basic'
    kwargs['bconds']       = 'nnn'
    kwargs['jastrows']     = []
    kwargs['corrections']  = []
    kwargs['calculations'] = []

    sim_args,inp_args = Simulation.separate_inputs(kwargs)

    input = generate_qmcpack_input(**inp_args)

    wf = input.get('wavefunction')
    if not 'determinantset' in wf:
        Qmcpack.class_error('wavefunction does not have determinantset, cannot create cusp correction','generate_cusp_correction')
    #end if
    wf.determinantset.cuspcorrection = True

    sim_args.input = input
    qmcpack = Qmcpack(**sim_args)

    return qmcpack
#end def generate_cusp_correction
