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
from copy import deepcopy
import numpy as np
from .simulation import Simulation, NullSimulationAnalyzer
from .qmcpack_input import (
    QmcpackInput,
    TracedQmcpackInput,
    loop,
    linear,
    cslinear,
    vmc,
    dmc,
    collection,
    determinantset,
    hamiltonian,
    init,
    pairpot,
    bspline_builder,
    generate_qmcpack_input,
    generate_jastrows,
    generate_jastrow,
    generate_jastrow1,
    generate_jastrow2,
    generate_jastrow3,
    generate_opt,
    generate_opts,
    check_excitation_type,
)
from .qmcpack_analyzer import QmcpackAnalyzer
from .qmcpack_converters import Pw2qmcpack, Convert4qmc, Convertpw4qmc, PyscfToAfqmc
from .pyscf_sim import Pyscf
from .developer import DevBase, obj, error, unavailable
from .nexus_base import nexus_core
from .hdfreader import read_hdf
from .unit_converter import convert
from .pwscf import Pwscf
from .xmlreader import XMLreader
try:
    import h5py
except:
    h5py = unavailable('h5py')
#end try


class GCTA(DevBase):
    '''
    This class holds the functionality and data to carry out grand canonical twist averaging in Nexus.
    Throughout the class, the handling of k-points uses unit (crystal) coordinates, which ranges in [0, 1).
    Note that QMCPACK interally uses the range (-0.5, 0.5) for k-points.
    '''
    def __init__(self, input, system, flavor):
        self.flavor = flavor
        self.input = input
        self.system = system
    #end def __init__

    def check_implementation(self, gcta_possible, dependency):
        gcta_flavors = {'safl', 'afl', 'nscf', 'scf'}
        if self.flavor.lower() not in gcta_flavors:
            self.error('GCTA type {} is not recognized. Valid options are {}.'.format(self.flavor, gcta_flavors))
        #end if
        if not gcta_possible:
            self.error('gcta keyword is not yet supported for this workflow. Please contact the developers.')
        #end if
        try:
            symm_kgrid = self.system.generation_info.symm_kgrid
        except:
            symm_kgrid = False
        if (self.flavor.lower() in ['safl', 'afl']) and (symm_kgrid == True):
            self.error('''
                safl and afl are not supported with symm_kgrid = True.
                It is possible to implement the afl and safl algorithms with k-point symmetries
                but it requires significant changes to the current simple implementation
                that strictly uses the Fermi level to set the occupations.
                Please contact the developers if this feature is pressing.
                    ''')
        #end if
        spinor_run = self.input.get('spinor')
        if (self.flavor.lower() == 'safl') and (spinor_run is True):
            self.error('safl is not supported with spinors. Use afl instead.')
        #end if
        if (self.flavor.lower() != 'afl') and (not isinstance(dependency,Pw2qmcpack)):
            self.error('{} flavor of GCTA is only supported with pwscf at the moment.'.format(self.flavor))
        #end if
        twistnum_input = self.input.get('twistnum')
        supercell_nkpoints = len(self.system.structure.kpoints)
        if (twistnum_input is not None) or (supercell_nkpoints == 1):
            self.error('''
                It appears that a single-twist QMC run was attempted using gcta keyword.
                Currently, this is not supported. Please contact the developers if this is needed.''')
        #end if
    #end def check_implementation

    @staticmethod
    def int_kpoint_weight(float_value, atol=1e-8):
        '''
        This function checks if the float k-point weight/norm is close to its integer value. If so, returns the integer value.
        '''
        int_value = round(float_value)
        assert abs(float_value - int_value) < atol, '''
            The k-point weight or norm ({}) is not close to an integer!
            There might be a problem with how the weights were stored.
            Please check the SCF conversion step.
            '''.format(float_value)
        return int_value
    #end def check_kpoint_weight

    def read_eshdf_data(self, filename):
        '''
        Read the ESHDF eigenvalues, k-point info and store the data in the GCTA instance as an attribute
        '''
        def h5_scalar(i):
            value = np.array(i)
            if value.ndim == 0:
                return value.item()
            else:
                return value[0]
        #end def h5_scalar
        h        = read_hdf(filename,view=True)
        nkpoints = h5_scalar(h.electrons.number_of_kpoints)
        if hasattr(h.electrons, 'number_of_spins'):
            nspins   = h5_scalar(h.electrons.number_of_spins) # pwscf collinear
        else:
            nspins   = 1 # convertpw4qmc non-collinear
        data     = obj()
        kweights = []
        for ikpoint in range(nkpoints):
            kp = h.electrons['kpoint_'+str(ikpoint)]
            kw = h5_scalar(kp.weight)
            kweights.append(kw)
            for ispin in range(nspins):
                path = 'electrons/kpoint_{0}/spin_{1}'.format(ikpoint,ispin)
                spin = h.get_path(path)
                eigs = convert(np.array(spin.eigenvalues),'Ha','eV')
                nstates = h5_scalar(spin.number_of_states)
                data[ikpoint,ispin] = obj(
                    eig    = np.array(eigs),
                    kpoint = np.array(kp.reduced_k), # unit (crystal) coordinates for kpoints. The range is [0, 1).
                    kweight = kw,
                    )
            #end for
        #end for
        total_kweight = sum(kweights)
        total_kweight = self.int_kpoint_weight(total_kweight)
        norm_factor = 1.0 / min(kweights) # Multiplicative factor to get integer weights
        res = obj(
            orbfile     = filename,
            nkpoints    = nkpoints,
            total_kw    = total_kweight,
            norm_factor = norm_factor,
            nspins      = nspins,
            nstates     = nstates,
            data        = data,
            )
        self.eig_data = res
    #end def read_eshdf_data

    def unfolded_nelecs(self):
        '''
        Returns the number of electrons in the primitive cell
        '''
        if self.system.folded_system is None:
            n_up = self.system.particles.up_electron.count
            n_dn = self.system.particles.down_electron.count
        else:
            n_up = self.system.folded_system.particles.up_electron.count
            n_dn = self.system.folded_system.particles.down_electron.count
        #end if
        nelecs = n_up + n_dn
        return nelecs
    #end def unfolded_nelecs

    def unfolded_nkpoints(self):
        '''
        Returns the number of unsymmetrized k-points when a supercell is unfolded back to the primitive cell
        '''
        kgrid = np.array(self.system.generation_info.kgrid)
        nkgrid = np.prod(kgrid)
        if self.system.folded_system is None:
            ntile = 1
        else:
            tmatrix = np.array(self.system.structure.tmatrix)
            ntile = np.linalg.det(tmatrix)
        #end if
        nkpoints = round(nkgrid * ntile)
        return nkpoints
    #end def unfolded_nkpoints

    def prim_kpoints(self):
        '''
        Returns the k-points used to build the supercell in unit coordinates
        '''
        if self.system.folded_system is None:
            qmc_kpoints = self.system.structure.kpoints_unit()
        else:
            qmc_kpoints = self.system.folded_system.structure.kpoints_unit()
        #end if
        return qmc_kpoints
    #end def unfolded_nkpoints

    def check_kmesh_size(self):
        '''
        Make sure that NSCF k-points and QMC twists are commensurate for GCTA
        '''
        n_qmc_kpoints = len(self.prim_kpoints())
        n_scf_kpoints = self.eig_data.nkpoints
        assert (n_scf_kpoints == n_qmc_kpoints), '''
            The number of k-points in (N)SCF ({}) and QMC ({}) are not commensurate!
            This is not supported. Please rerun the (N)SCF and conversion steps such
            that the unfolded system contains the same number of k-points in both cases.
            '''.format(n_scf_kpoints, n_qmc_kpoints)
    #end def check_kmesh_size

    def check_kpoint_consistency(self, tol=1e-8):
        '''
        The kpoints expected by the GCTA object and what is found in the self.eig_data.data should be consistent.
        The kpoints in self.eig_data.data are expected to be in unit coordinates. (Conversion: dot(kpoints, inv(kaxes))).
        This function checks if there is 1-to-1 mapping between the GCTA object and the converted data.
        '''
        gcta_kpoints = self.prim_kpoints()
        nkpoints = self.eig_data.nkpoints
        eig_kpoints = []
        for ikpoint in range(nkpoints):
            eig_kpoints.append(self.eig_data.data[ikpoint, 0].kpoint) # 0: only checking the consistency in one spin channel
        #end for
        eig_kpoints = np.array(eig_kpoints)
        # Check if each row of gcta_kpoints exists in eig_kpoints
        for gcta_row in gcta_kpoints:
            if not np.any(np.all(np.isclose(eig_kpoints, gcta_row, atol=tol), axis=1)):
                self.error('''The GCTA k-point {} was not found in the converted data. This is not supposed to happen.
                            Please make sure that the k-points were written in unit coordinates.'''.format(gcta_row))
            #end if
        #end for
    #end def check_kpoint_consistency

    def gcta_converter_kmapping(self, tol=1e-8):
        '''
        The k-points defined by the GCTA object and the k-points written by a converter may have different ordering.
        We need to figure out the mapping between these two so that the k-points fold into correct twists.
        '''
        gcta2conv = {}  # The dictionary that holds the gcta -> converter k-mapping
        gcta_kpoints = self.prim_kpoints()
        nkpoints = self.eig_data.nkpoints
        eig_kpoints = []
        for ikpoint in range(nkpoints):
            eig_kpoints.append(self.eig_data.data[ikpoint, 0].kpoint) # 0: only need one spin channel
        #end for
        eig_kpoints = np.array(eig_kpoints)
        # Check if each row of gcta_kpoints exists in eig_kpoints
        for i, gcta_row in enumerate(gcta_kpoints):
            for k, eig_row in enumerate(eig_kpoints):
                if np.all(np.isclose(gcta_row, eig_row, atol=tol), axis=0):
                    gcta2conv[i] = k
                #end if
            #end for
        #end for
        self.gcta2conv = gcta2conv
    #end def gcta_converter_kmapping

    @staticmethod
    def traceback_dependency(dependency, cls, levels = 1):
        '''
        This function provides limited functionality to go back in dependency by a certain level 
        '''
        if dependency is None:
            error('This function requires a valid dependency. None was given.')
        #end if
        if levels < 1:
            error('Traceback level should be at least one. {} was given.'.format(levels))
        #end if
        current_dep = dependency
        for level in range(levels):
            len_dep = 0
            for dep in current_dep.dependencies:
                if isinstance(dep.sim, cls):
                    found_dep = dep.sim
                    len_dep += 1
                #end if
            #end for
            current_dep = found_dep
            if len_dep != 1:
                error('This function can only traceback using single dependecies! Found {}'.format(len_dep))
            #end if
        #end for
        return current_dep.locdir
    #end def

    @staticmethod
    def pwscf_tot_magnet(filepath):
        file = '{}/pwscf_output/pwscf.xml'.format(filepath)
        xml = XMLreader(file, warn=False).obj
        calculation = xml['qes:espresso']['input']['control_variables']['calculation']['text']
        assert (calculation == 'scf'), 'The total magnetization should be obtained from an SCF run'
        noncolinear = xml['qes:espresso']['input']['spin']['noncolin']['text']
        assert (noncolinear == 'false'), 'Noncollinear calculations are not supported by this function'
        spin_polarized = xml['qes:espresso']['input']['spin']['lsda']['text']
        if spin_polarized == 'true':
            scf_magnet = float(xml['qes:espresso']['output']['magnetization']['total']['text'])
        elif spin_polarized == 'false': # total magnetization is not written for nspin = 1
            scf_magnet = 0.0
        else:
            scf_magnet = None
        #end if
        return scf_magnet
    #end if

    @staticmethod
    def pwscf_fermi(filepath, scf_type):
        file = '{}/pwscf_output/pwscf.xml'.format(filepath)
        xml = XMLreader(file, warn=False).obj
        calculation = xml['qes:espresso']['input']['control_variables']['calculation']['text']
        assert (calculation == scf_type), 'The Fermi level should be obtained from an {} run.'.format(scf_type)
        tot_magnetization = False
        if 'tot_magnetization' in xml['qes:espresso']['input']['bands'].keys():
            tot_magnetization = True
        #end if
        if tot_magnetization == True:
            up_fermi = float(xml['qes:espresso']['output']['band_structure']['two_fermi_energies']['text'].split()[0])
            dn_fermi = float(xml['qes:espresso']['output']['band_structure']['two_fermi_energies']['text'].split()[1])
            fermi_level = np.array([up_fermi, dn_fermi])
        else:
            fermi_level = float(xml['qes:espresso']['output']['band_structure']['fermi_energy']['text'])
        #end if
        fermi_level = convert(fermi_level,'Ha','eV')
        return fermi_level
    #end if

    def adapted_fermi_level(self):
        combined_eigens = []
        data = self.eig_data.data
        norm_factor = self.eig_data.norm_factor # normalization factor to get integer k-weights
        nkpoints = self.eig_data.nkpoints
        nspins = self.eig_data.nspins
        for ispin in range(nspins):
            for ikpoint in range(nkpoints):
                kweight = data[ikpoint,ispin].kweight
                ksym_range = kweight * norm_factor
                ksym_range = self.int_kpoint_weight(ksym_range)
                for ksym in range(ksym_range):
                    combined_eigens.extend(data[ikpoint,ispin].eig)
                #end for
            #end for
        #end for
        spinor_run = self.input.get('spinor')
        if (spinor_run is not True) and (nspins == 1):
            combined_eigens.extend(combined_eigens)
        #end if
        combined_eigens = sorted(combined_eigens)
        nelecs_prim = self.unfolded_nelecs()
        nosym_kpoints = self.unfolded_nkpoints()
        lamda_index = nelecs_prim * nosym_kpoints # The index in the eigenvalue list that produces charge neutral system
        fermi_level = float(combined_eigens[lamda_index-1] + combined_eigens[lamda_index]) / 2
        return fermi_level
    #end def adapted_fermi_level

    def spin_adapted_fermi_level(self, scf_magnet):
        if scf_magnet is None:
            self.error('The reference magnetization in safl can not be None. Please check that the SCF is appropriate.')
        #end if
        combined_eigens = {}
        data = self.eig_data.data
        norm_factor = self.eig_data.norm_factor # normalization factor to get integer k-weights
        nkpoints = self.eig_data.nkpoints
        nspins = self.eig_data.nspins
        for ispin in range(nspins):
            if ispin not in combined_eigens:
                combined_eigens[ispin] = []
            #end if
            for ikpoint in range(nkpoints):
                kweight = data[ikpoint,ispin].kweight
                ksym_range = kweight * norm_factor
                ksym_range = self.int_kpoint_weight(ksym_range)
                for ksym in range(ksym_range):
                    combined_eigens[ispin].extend(data[ikpoint,ispin].eig)
                #end for
            #end for
            combined_eigens[ispin] = sorted(combined_eigens[ispin])
        #end for
        if nspins == 1:
            combined_eigens[1] = combined_eigens[0]
        #end if
        nelecs_prim = self.unfolded_nelecs()
        nosym_kpoints = self.unfolded_nkpoints()
        up_index = round((nelecs_prim + scf_magnet) * nosym_kpoints / 2)
        dn_index = (nelecs_prim * nosym_kpoints) - up_index
        up_fermi = float(combined_eigens[0][up_index-1] + combined_eigens[0][up_index]) / 2
        dn_fermi = float(combined_eigens[1][dn_index-1] + combined_eigens[1][dn_index]) / 2
        fermi_level = np.array([up_fermi, dn_fermi])
        return fermi_level
    #end def adapted_fermi_level

    def set_gcta_occupations(self, fermi_level):
        if fermi_level is None:
            self.error('The Fermi level can not be None. This indicates a bug in {}'.format(self.flavor))
        #end if
        ntwists = len(self.system.structure.kpoints)
        nspins = self.eig_data.nspins
        nstates = self.eig_data.nstates
        gcta2conv = self.gcta2conv
        fermi_levels = fermi_level
        if isinstance(fermi_levels, float):
            fermi_levels = [fermi_levels, fermi_levels]
        #end if
        # kmap is mapping between twists and k-points (internal to gcta, not to be confused by gcta2conv mapping)
        kmap = self.system.structure.kmap()
        if kmap is None:
            kmap = self.system.structure.unique_kpoints()
        #end if
        nelecs_at_twist = []
        for itwist in range(ntwists):
            # calculate nelec for each spin
            nelec_up_dn = []
            for ispin in range(nspins):
                nelec_spin = 0
                for ikpoint in kmap[itwist]:
                    for istate in range(nstates):
                        eig = self.eig_data.data[gcta2conv[ikpoint], ispin].eig[istate]
                        if eig < fermi_levels[ispin]:
                            nelec_spin += 1
                        #end if
                    #end for
                #end for
                nelec_up_dn.append(nelec_spin)
                spinor_run = self.input.get('spinor')
                if (spinor_run is not True) and (nspins == 1):
                    nelec_up_dn.append(nelec_spin)
                #end if
            #end for
            nelecs_at_twist.append(nelec_up_dn)
        #end for
        self.nelecs_at_twist = nelecs_at_twist
    #end set_gcta_occupation

    def sum_charge_twists(self):
        '''
        Returns the net charge of a system with multiple twists (not averaged)
        '''
        n_up = self.system.particles.up_electron.count
        n_dn = self.system.particles.down_electron.count
        n_total = n_up + n_dn
        nelecs_at_twist = self.nelecs_at_twist
        kweights = np.array(self.system.structure.kweights)
        assert (len(kweights) == len(nelecs_at_twist))
        q_sum_twists = 0
        for itwist, nelec_up_dn in enumerate(nelecs_at_twist):
            nelec_twist = sum(nelec_up_dn)
            q_twist = n_total - nelec_twist
            q_sum_twists += q_twist * round(kweights[itwist])
        #end for
        return q_sum_twists
    #end def sum_charge_twists

    def sum_spin_twists(self):
        '''
        Returns the net spin of a system with multiple twists (not averaged)
        '''
        nelecs_at_twist = self.nelecs_at_twist
        kweights = np.array(self.system.structure.kweights)
        assert (len(kweights) == len(nelecs_at_twist))
        spin_sum_twists = 0
        for itwist, nelec_up_dn in enumerate(nelecs_at_twist):
            spin_twist = nelec_up_dn[0] - nelec_up_dn[1]
            spin_sum_twists += spin_twist * round(kweights[itwist])
        #end for
        return spin_sum_twists
    #end def sum_spin_twists

    def check_charge_neutrality(self):
        '''
        Check the net charge of the twist averaged system
        '''
        q_sum_twists = self.sum_charge_twists()
        if (self.flavor.lower() in ['safl', 'afl']) and (q_sum_twists != 0):
            self.error('''
                The sum of charges over all twists is {} electrons!
                This is not supposed to happen for afl or safl!
                Check that the spinor keyword is correctly used in generate_qmcpack.
                Otherwise, there might be a bug in the implementation of gcta.
                '''.format(q_sum_twists))
        #end if
    #end def check_charge_neutrality

    def check_magnetization_accuracy(self, scf_magnet):
        '''
        Check that the net magnetization is close to the reference SCF value
        '''
        if self.flavor.lower() == 'safl':
            nosym_kpoints = self.unfolded_nkpoints()
            spin_sum_twists = self.sum_spin_twists()
            qmc_magnet = spin_sum_twists / nosym_kpoints
            feasible_accuracy = (1.0 / nosym_kpoints) + 1e-8
            error_magnet = abs(qmc_magnet - scf_magnet)
            if error_magnet > feasible_accuracy:
                self.error('''
                    The twist-averaged QMC magnetization ({:.16f}) is not close to the SCF reference value ({:.16f})!
                    This is not supposed to happen for safl. Likely, there is a bug in the implementation of safl.
                    '''.format(qmc_magnet, scf_magnet))
            #end if
        #end if
    #end def check_magnetization_accuracy

    def write_gcta_report(self, locdir, fermi_level, scf_magnet = None):
        spinor_run = self.input.get('spinor')
        nosym_kpoints = self.unfolded_nkpoints()
        q_sum_twists = self.sum_charge_twists()
        qmc_charge = q_sum_twists / nosym_kpoints
        if spinor_run is not True:
            spin_sum_twists = self.sum_spin_twists()
            qmc_magnet = spin_sum_twists / nosym_kpoints
        #end if
        n_up = self.system.particles.up_electron.count
        n_dn = self.system.particles.down_electron.count
        n_total = n_up + n_dn
        nelecs_at_twist = self.nelecs_at_twist
        fermi_level = np.array(fermi_level)
        filepath = '{}/gcta_report.txt'.format(locdir)
        with open(filepath, 'w') as gcta_file:
            # Writing data to a file
            gcta_file.write('SUMMARY FOR GCTA OCCUPATIONS:\n')
            gcta_file.write('==================================================\n')
            gcta_file.write('GCTA Flavor:                    {}\n'.format(self.flavor))
            if fermi_level.size == 1:
                gcta_file.write('Fermi Level [eV]              {:20.16f}\n'.format(fermi_level))
            elif fermi_level.size == 2:
                gcta_file.write('Fermi Level Up [eV]           {:20.16f}\n'.format(fermi_level[0]))
                gcta_file.write('Fermi Level Dn [eV]           {:20.16f}\n'.format(fermi_level[1]))
            else:
                self.error('The number of provided Fermi levels ({}) does not make sense'.format(fermi_level.size))
            #end if
            gcta_file.write('Net Charge:                     {}\n'.format(q_sum_twists))
            gcta_file.write('Net Charge / Prim Cell:       {:20.16f}\n'.format(qmc_charge))
            if spinor_run is not True:
                gcta_file.write('Net Magnetization / Prim Cell:{:20.16f}\n'.format(qmc_magnet))
            #end if
            if scf_magnet is not None:
                gcta_file.write('SCF Magnetization (Reference):{:20.16f}\n'.format(scf_magnet))
            #end if
            gcta_file.write('\n\n')
            if spinor_run is not True:
                gcta_file.write(' TWISTNUM  NELEC_UP  NELEC_DN   CHARGE     SPIN   \n')
            else:
                gcta_file.write(' TWISTNUM    NELEC    CHARGE                      \n')
            #end if
            gcta_file.write('==================================================\n')
            for itwist, nelec_up_dn in enumerate(nelecs_at_twist):
                nelec_twist = sum(nelec_up_dn)
                q_twist = n_total - nelec_twist
                gcta_file.write('{:^10}'.format(itwist))
                gcta_file.write('{:^10}'.format(nelec_up_dn[0]))
                if spinor_run is not True:
                    gcta_file.write('{:^10}'.format(nelec_up_dn[1]))
                #end if
                gcta_file.write('{:^11}'.format(q_twist))
                if spinor_run is not True:
                    spin_twist = nelec_up_dn[0] - nelec_up_dn[1]
                    gcta_file.write('{:^9}'.format(spin_twist))
                #end if
                gcta_file.write('\n')
            #end for
        #end with
        self.log('    See the GCTA occupation report at:  {}'.format(filepath))
    #end def write_gcta_report
#end class GCTA



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
            user_twist_given  = isinstance(twh,htypes) and twh.twist is not None
            user_twist_given |= isinstance(tnh,htypes) and tnh.twistnum is not None
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
            if 'results' not in analyzer or 'optimization' not in analyzer.results:
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
            gcta_possible = False
            if isinstance(sim,Pw2qmcpack) or isinstance(sim,Convertpw4qmc):

                gcta_possible = True
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
                    if system.structure.folded_structure is not None:
                        orb_elem.tilematrix = np.array(system.structure.tmatrix)
                    #end if
                #end if
                defs = obj(
                    #twistnum   = 0,
                    meshfactor = 1.0
                    )
                for var,val in defs.items():
                    if var not in orb_elem:
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
                #end if
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

            elif isinstance(sim,Pyscf):
                sinp = sim.input
                skpoints = None
                if sinp.tiled_kpoints is not None:
                    skpoints = sinp.tiled_kpoints
                elif sinp.kpoints is not None:
                    skpoints = sinp.kpoints
                #end if
                if skpoints is None:
                    self.error('cannot incorporate orbitals from pyscf\nno k-points are present')
                #end if
                nkpoints = len(self.system.structure.kpoints)
                if len(skpoints)!=nkpoints:
                    self.error('cannot incorporate orbitals from pyscf\nwrong number k-points are present\nexpected: {}\npresent: {}'.format(nkpoints,len(skpoints)))
                #end if
                twist_updates = []
                for n,(h5file,kp) in enumerate(zip(result.orb_files,result.kpoints)):
                    filepath = os.path.join(result.location,h5file)
                    tu = obj(
                        twistnum = -1,
                        twist    = tuple(kp),
                        href     = os.path.relpath(filepath,self.locdir),
                        )
                    twist_updates.append(tu)
                #end for
                ds = self.input.get('determinantset')
                ds.twistnum = -1 # set during twist average
                self.twist_average(twist_updates)
            else:
                self.error('incorporating orbitals from '+sim.__class__.__name__+' has not been implemented')
            #end if

            # Activate GCTA occupations if gcta is specified by the user
            gcta_flavor = self.get_optional('gcta', None)

            if (gcta_flavor is not None) and (self.sent_files is not True):
                # Create a GCTA object with deepcopied arguments to avoid interference with Qmcpack instance
                gcta_input = deepcopy(self.input)
                gcta_system = deepcopy(self.system)
                gcta_dependency = deepcopy(sim)
                gcta_locdir = deepcopy(self.locdir)
                gcta_obj = GCTA(gcta_input, gcta_system, gcta_flavor)

                gcta_obj.check_implementation(gcta_possible, gcta_dependency)

                gcta_obj.log('    Reading the eigenvalue and k-point data for GCTA. This might take a while.')
                if isinstance(gcta_dependency,Pw2qmcpack) or isinstance(gcta_dependency,Convertpw4qmc):
                    gcta_obj.read_eshdf_data(h5file)
                else:
                    gcta_obj.error('Reading the eigenvalues for this workflow ({}) is not yet implemented.'.format(gcta_dependency.__class__.__name__))
                #end if

                gcta_obj.check_kmesh_size()
                gcta_obj.check_kpoint_consistency()
                gcta_obj.gcta_converter_kmapping()

                # === Determine the Fermi level ===
                fermi_level = None
                scf_magnet = None

                if gcta_flavor.lower() == 'safl':
                    # We need to get the SCF total magnetization for safl case
                    if isinstance(gcta_dependency,Pw2qmcpack):
                        filepath = gcta_obj.traceback_dependency(gcta_dependency, Pwscf, levels = 2)
                        scf_magnet = gcta_obj.pwscf_tot_magnet(filepath)
                    else:
                        gcta_obj.error('Reading the total magnetization for this workflow ({}) is not yet implemented.'.format(gcta_dependency.__class__.__name__))
                    #end if
                    fermi_level = gcta_obj.spin_adapted_fermi_level(scf_magnet)

                elif gcta_flavor.lower() == 'afl':
                    fermi_level = gcta_obj.adapted_fermi_level()

                elif gcta_flavor.lower() == 'nscf':
                    if isinstance(gcta_dependency,Pw2qmcpack):
                        filepath = gcta_obj.traceback_dependency(gcta_dependency, Pwscf, levels = 1)
                        fermi_level = gcta_obj.pwscf_fermi(filepath, 'nscf')
                    else:
                        gcta_obj.error('Reading the Fermi level for this workflow ({}) is not yet implemented.'.format(gcta_dependency.__class__.__name__))
                    #end if

                elif gcta_flavor.lower() == 'scf':
                    if isinstance(gcta_dependency,Pw2qmcpack):
                        filepath = gcta_obj.traceback_dependency(gcta_dependency, Pwscf, levels = 2)
                        fermi_level = gcta_obj.pwscf_fermi(filepath, 'scf')
                    else:
                        gcta_obj.error('Reading the Fermi level for this workflow ({}) is not yet implemented.'.format(gcta_dependency.__class__.__name__))
                    #end if

                else:
                    gcta_obj.error('GCTA type {} is not recognized.'.format(gcta_flavor))
                # === Finished determining the Fermi level ===

                # Set the twist occupations based on the user-requested Fermi level
                gcta_obj.set_gcta_occupations(fermi_level)

                # Final checks and report
                gcta_obj.check_charge_neutrality()
                gcta_obj.check_magnetization_accuracy(scf_magnet)
                gcta_obj.write_gcta_report(gcta_locdir, fermi_level, scf_magnet)

                # The final GCTA occupations are deepcopied to the Qmcpack instance
                self.nelecs_at_twist = deepcopy(gcta_obj.nelecs_at_twist)
            #end if (GCTA preprocessing is done)

        elif result_name=='jastrow':
            if isinstance(sim,Qmcpack):
                opt_file = result.opt_file
                opt = QmcpackInput(opt_file)
                wavefunction = input.get('wavefunction')
                optwf = opt.qmcsystem.wavefunction
                # handle spinor case
                spinor = input.get('spinor')
                if spinor is not None and spinor:
                    # remove u-d term from optmized jastrow
                    # also set correct cusp condition
                    J2 = optwf.get('J2')
                    if J2 is not None:
                        corr = J2.get('correlation')
                        if 'ud' in corr:
                            del corr.ud
                        if 'uu' in corr:
                            corr.uu.cusp = -0.5
                        #end if
                    #end if
                    J3 = optwf.get('J3')
                    if J3 is not None:
                        corr = J3.get('correlation')
                        if hasattr(corr, 'coefficients'):
                            # For single-species systems, the data structure changes.
                            # In this case, the only J3 term should be 'uu'.
                            # Otherwise, the user might be trying to do something strange.
                            assert 'uu' in corr.coefficients.id, 'Only uu J3 terms are allowed in SOC calculations.'
                        else:
                            j3_ids = []
                            for j3_term in corr:
                                j3_id = j3_term.coefficients.id
                                j3_ids.append(j3_id)
                            #end for
                            for j3_id in j3_ids:
                                if 'ud' in j3_id:
                                    delattr(corr, j3_id)
                                #end if
                            #end for
                    #end if
                #end if
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
                        if jtype not in jopt:
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
                wfn = opt.qmcsystem.wavefunction.copy()
                ovp = 'override_variational_parameters' # name is too long
                if ovp in wfn:
                    href = os.path.join(sim.locdir,wfn[ovp].href)
                    href = os.path.relpath(href,self.locdir)
                    wfn[ovp].href = href
                qs.wavefunction = wfn
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
            from .qmcpack_converters import gcta_occupation
            if not isinstance(sim,Pw2qmcpack):
                msg = 'grand-canonical occupation requires Pw2qmcpack'
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
            # step 2: analyze ESHDF file for states below Fermi energy
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
        elif result_name=='determinantset':
            # This should be removed someday far in the future, 
            # when QMCPACK can actually read its HDF5 files all by itself.
            if isinstance(sim,Convert4qmc):
                wf = input.get('wavefunction')
                if isinstance(wf,collection):
                    wf = wf.get_single('psi0')
                #end if
                spoc_key = None
                if 'sposet_builder' in wf and 'spline' in wf.sposet_builder.type.lower():
                    spoc_key = 'sposet_builder'
                elif 'sposet_builders' in wf and ('bspline' in wf.sposet_builders or 'einspline' in wf.sposet_builders):
                    spoc_key = 'sposet_builders'                    
                elif 'sposet_collection' in wf and 'spline' in wf.sposet_builder.type.lower():
                    spoc_key = 'sposet_builder'
                elif 'sposet_collections' in wf and ('bspline' in wf.sposet_collections or 'einspline' in wf.sposet_collections):
                    spoc_key = 'sposet_collections'                    
                #end if
                if spoc_key is not None:
                    del wf[spoc_key]
                #end if
                c4q_inp = QmcpackInput(result.location)
                ds = c4q_inp.get('determinantset')
                # only one twist is supported by qmcpack right now, so don't need to know it
                if 'twist' in ds:
                    del ds.twist
                #end if
                wf.determinantset = ds
            else:
                self.error('incorporating determinantset from '+sim.__class__.__name__+' has not been implemented')
            #end if                
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
                if sd is not None:
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
            opt_run = calctypes is not None and 'opt' in calctypes
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

                        from .structure import get_kpath
                        kpath       = get_kpath(structure=structure)
                        kpath_label = np.array(kpath['explicit_kpoints_labels'])
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
                        qi.get('u').set(size=nup)
                        if len(elecs) == 2:
                            ndn = elecs[1]
                            qi.get('d').set(size=ndn)
                        #end if
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
                data = np.array(f.read().split()[1:])
                data.shape = len(data)//12,12
                data = data.T
                for darr in data:
                    if darr[0][0]=='K' or darr[0][0]=='E':
                        edata[spinlab][darr[0]] = np.array(list(map(float,darr[1:])))
                    else:
                        edata[spinlab][darr[0]] = np.array(list(map(int,darr[1:])))
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

    gcta = None
    if 'gcta' in inp_args:
        gcta = deepcopy(inp_args.gcta)
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

    if gcta is not None:
        qmcpack.gcta = gcta
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
    if 'determinantset' not in wf:
        Qmcpack.class_error('wavefunction does not have determinantset, cannot create cusp correction','generate_cusp_correction')
    #end if
    wf.determinantset.cuspcorrection = True

    sim_args.input = input
    qmcpack = Qmcpack(**sim_args)

    return qmcpack
#end def generate_cusp_correction
