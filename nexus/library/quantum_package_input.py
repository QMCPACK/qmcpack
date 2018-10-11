##################################################################
##  (c) Copyright 2018-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  quantum_package_input.py                                          #
#    Supports I/O for Quantum Package input data.                    #
#                                                                    #
#  Content summary:                                                  #
#    QuantumPackageInput                                             #
#      SimulationInput class for Quantum Package.                    #
#                                                                    #
#    generate_quantum_package_input                                  #
#      User-facing function to create arbitrary input.               #
#====================================================================#

import os
from generic import obj
from developer import DevBase,log,error
from structure import Structure
from physical_system import PhysicalSystem
from simulation import SimulationInput
from execute import execute


bool_values = dict(T=True,F=False)
bool_values_inv = {True:'T',False:'F'}

def read_qp_value_type(value_filepath):
    f = open(value_filepath,'r')
    svalue = f.read().strip()
    f.close()
    if svalue in bool_values:
        return 'bool'
    else:
        try:
            v = int(svalue)
            return 'int'
        except:
            try:
                v = float(svalue)
                return 'float'
            except:
                return 'str'
            #end try
        #end try
    #end try
#end def read_qp_value_type


def read_qp_value(value_filepath):
    f = open(value_filepath,'r')
    svalue = f.read().strip()
    f.close()
    if svalue in bool_values:
        v = bool_values[svalue]
    else:
        try:
            v = int(svalue)
        except:
            try:
                v = float(svalue)
            except:
                v = svalue
            #end try
        #end try
    #end try
    return v
#end def read_qp_value


def write_qp_value(value_filepath,value):
    if isinstance(value,bool):
        svalue = bool_values_inv[value]
    elif isinstance(value,int):
        svalue = str(value)
    elif isinstance(value,float):
        '{0: 24.15e}'.format(value)
    elif isinstance(value,str):
        svalue = value
    else:
        QuantumPackageInput.class_error('invalid type encountered on write\nattempted to write variable: {0}\nwith type: {1}\nvalid type options: bool,int,float,str'.format(value_filepath,value.__class__.__name__))
    #end if
    f = open(value_filepath,'w')
    f.write(svalue+'\n')
    f.close()
#end def write_qp_value



# quantum package path input specification
input_specification = obj({
    'ao_basis/ao_basis' : str,
    'ao_basis/ao_cartesian' : bool,
    'ao_basis/ao_md5' : str,
    'ao_basis/ao_num' : int,
    'bitmasks/bit_kind' : int,
    'bitmasks/n_int' : int,
    'bitmasks/n_mask_cas' : int,
    'bitmasks/n_mask_gen' : int,
    'davidson/davidson_sze_max' : int,
    'davidson/disk_based_davidson' : bool,
    'davidson/distributed_davidson' : bool,
    'davidson/n_states_diag' : int,
    'davidson/state_following' : bool,
    'davidson/threshold_davidson' : float,
    'determinants/bit_kind' : int,
    'determinants/expected_s2' : float,
    'determinants/mo_label' : str,
    'determinants/n_det' : int,
    'determinants/n_det_max' : int,
    'determinants/n_det_max_jacobi' : int,
    'determinants/n_det_max_property' : int,
    'determinants/n_det_max_stored' : int,
    'determinants/n_int' : int,
    'determinants/n_states' : int,
    'determinants/only_single_double_dm' : bool,
    'determinants/read_wf' : bool,
    'determinants/s2_eig' : bool,
    'determinants/store_full_h_mat' : bool,
    'determinants/target_energy' : float,
    'determinants/threshold_generators' : float,
    'determinants/threshold_selectors' : float,
    'electrons/elec_alpha_num' : int,
    'electrons/elec_beta_num' : int,
    'ezfio/creation' : str,
    'ezfio/library' : str,
    'ezfio/user' : str,
    'full_ci_zmq/energy' : float,
    'full_ci_zmq/energy_pt2' : float,
    'full_ci_zmq/iterative_save' : int,
    'full_ci_zmq/n_iter' : int,
    'hartree_fock/energy' : float,
    'hartree_fock/level_shift' : float,
    'hartree_fock/max_dim_diis' : int,
    'hartree_fock/mo_guess_type' : str,
    'hartree_fock/n_it_scf_max' : int,
    'hartree_fock/no_oa_or_av_opt' : bool,
    'hartree_fock/scf_algorithm' : str,
    'hartree_fock/thresh_scf' : float,
    'hartree_fock/threshold_diis' : float,
    'integrals_bielec/direct' : bool,
    'integrals_bielec/disk_access_ao_integrals' : str,
    'integrals_bielec/disk_access_mo_integrals' : str,
    'integrals_bielec/no_ivvv_integrals' : bool,
    'integrals_bielec/no_vvv_integrals' : bool,
    'integrals_bielec/no_vvvv_integrals' : bool,
    'integrals_bielec/threshold_ao' : float,
    'integrals_bielec/threshold_mo' : float,
    'integrals_monoelec/disk_access_ao_one_integrals' : str,
    'integrals_monoelec/disk_access_mo_one_integrals' : str,
    'mo_basis/ao_md5' : str,
    'mo_basis/mo_label' : str,
    'mo_basis/mo_tot_num' : int,
    'mrpt_utils/do_third_order_1h1p' : bool,
    'nuclei/disk_access_nuclear_repulsion' : str,
    'nuclei/nucl_num' : int,
    'perturbation/correlation_energy_ratio_max' : float,
    'perturbation/do_pt2' : bool,
    'perturbation/pt2_absolute_error' : float,
    'perturbation/pt2_max' : float,
    'perturbation/pt2_relative_error' : float,
    'perturbation/threshold_generators_pt2' : float,
    'perturbation/threshold_selectors_pt2' : float,
    'properties/threshld_two_bod_dm' : float,
    'properties/z_one_point' : float,
    'pseudo/do_pseudo' : bool,
    'pseudo/pseudo_grid_rmax' : float,
    'pseudo/pseudo_grid_size' : int,
    'pseudo/pseudo_klocmax' : int,
    'pseudo/pseudo_kmax' : int,
    'pseudo/pseudo_lmax' : int,
    'qmc/ci_threshold' : float,
    'work/empty' : bool,
    'work/qp_run_address' : str,
    })


# create mapping from variable name to section (directory) name
known_sections = set()
known_variables = set()
variable_section = obj()
section_variables = obj()
for vpath in input_specification.keys():
    secname,varname = vpath.split('/')
    known_sections.add(secname)
    known_variables.add(varname)
    if varname not in variable_section:
        variable_section[varname] = secname
    else:
        vsec = variable_section[varname]
        if isinstance(vsec,str):
            vsec = [vsec,secname]
        else:
            vsec.append(secname)
        #end if
    #end if
    if secname not in section_variables:
        section_variables[secname] = set()
    #end if
    section_variables[secname].add(varname)
#end for
for varname in variable_section.keys():
    vsec = variable_section[varname]
    if isinstance(vsec,list):
        variable_section[varname] = tuple(vsec)
    #end if
#end for


# function to extract and print an updated input_specification based on a list of ezfio directories
def extract_input_specification(*ezfio_paths):
    if len(ezfio_paths)==1 and isinstance(ezfio_paths[0],(list,tuple)):
        ezfio_paths = ezfio_paths[0]
    #end if
    log('\nextracting Quantum Package input specification from ezfio directories')
    typedict = {bool:'bool',int:'int',float:'float',str:'str'}
    new_input_spec = obj()
    for vpath,vtype in input_specification.iteritems():
        new_input_spec[vpath] = typedict[vtype]
    #end for
    for epath in ezfio_paths:
        epath = epath.rstrip('/')
        if not epath.endswith('.ezfio'):
            error('cannot extract input spec from path\ninput path provided is not an ezfio directory\ninput path provided: {0}'.format(epath),'Quantum Package')
        elif not os.path.exists(epath):
            error('cannot extract input spec from path\ninput path provided does not exist\ninput path provided: {0}'.format(epath),'Quantum Package')
        #end if
        log('  extracting from: {0}'.format(epath))
        for path,dirs,files in os.walk(epath):
            for file in files:
                if 'save' not in path:
                    if not file.startswith('.') and not file.endswith('.gz'):
                        filepath = os.path.join(path,file)
                        vtype = read_qp_value_type(filepath)
                        vpath = filepath.replace(epath,'').strip('/')
                        if vpath not in new_input_spec:
                            new_input_spec[vpath] = vtype
                        #end if
                    #end if
                #end if
            #end for
        #end for
    #end for
    log('  extraction complete')

    old_vpaths = set(input_specification.keys())
    new_vpaths = set(new_input_spec.keys())
    if new_vpaths==old_vpaths:
        log('\ninput specification in quantum_package_input.py needs no changes\n')
    else:
        log('\nplease replace input_specification in quantum_package_input.py with the following:\n')
        log('input_specification = obj({')
        s = ''
        for vpath in sorted(new_input_spec.keys()):
            vtype = new_input_spec[vpath]
            s += "    '{0}' : {1},\n".format(vpath,vtype)
        #end for
        s += '    })\n'
        log(s)
    #end if
#end def extract_input_specification



class Section(DevBase):
    None
#end class Section



class QuantumPackageInput(SimulationInput):

    added_keys = '''
        structure
        run_control
        '''.split()

    run_types = set('''
        Gen_Ezfio_from_integral.sh
        H_CORE_guess
        Huckel_guess
        SCF
        SCF_old
        analyze_wf
        check_orthonormality
        cis
        create_ezfio.py
        davidson_slave
        densify_coefmatrix
        diagonalize_restart_and_save_all_states
        diagonalize_restart_and_save_lowest_state
        diagonalize_restart_and_save_one_state
        diagonalize_restart_and_save_two_states
        dump_nto
        dump_one_body_mos
        fci_zmq
        fci_zmq_nos
        four_idx_transform
        guess_doublet
        guess_lowest_state
        guess_singlet
        guess_triplet
        localize_mos
        mo_permutation
        overwrite_with_cas
        print_H_matrix_restart
        print_aos
        print_bitmask
        print_energy
        print_hcc
        print_holes_particles
        print_integrals_ao
        print_integrals_mo
        print_mo_in_space
        print_mulliken
        print_spin_density
        print_wf
        provide_deltarho
        pt2_slave
        pt2_stoch
        pyscf.main
        qmc_create_wf
        qmc_e_curve
        qp_ao_ints
        qp_convert_qmcpack_to_ezfio.py
        read_ao_eri_chunk
        read_integrals_achocol
        read_integrals_achocol2
        read_integrals_ao
        read_integrals_mo
        read_integrals_mo_chocol
        save_HF_determinant
        save_for_qmcchem
        save_for_qmcpack
        save_natorb
        save_only_singles
        save_ortho_mos
        selection_davidson_slave
        selection_slave
        super_ci
        swap_mos
        target_pt2_qmc
        target_pt2_ratio_zmq
        target_pt2_zmq
        test_integrals
        test_two_body_dm
        truncate_wf
        truncate_wf_spin
        '''.split())

    slave_allowed = set('''
        SCF
        cis
        fci_zmq
        '''.split())

    integral_write_allowed = set('''
        SCF
        cis
        fci_zmq
        '''.split())


    def __init__(self,filepath=None):
        self.structure   = None
        self.run_control = obj()
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__


    def present(self,name):
        if name not in known_variables:
            self.error('attempted to check presence of unknown variable "{0}"valid options are: {1}'.format(name,sorted(known_variables)))
        #end if
        secname = variable_section[name]
        return secname in self and name in self[secname]
    #end def present


    def set(self,**kwargs):
        for name,value in kwargs.iteritems():
            if name not in known_variables:
                self.error('cannot set variable\nattempted to set unknown variable "{0}"\nwith value: {1}\nvalid options are: {2}'.format(name,value,sorted(known_variables)))
            #end if
            secname = variable_section[name]
            if secname not in self:
                self[secname] = Section()
            #end if
            self[secname][name] = value
        #end for
    #end def set


    def get(self,name):
        if name not in known_variables:
            self.error('cannot get variable\nattempted to get unknown variable "{0}"\nvalid options are: {1}'.format(name,sorted(known_variables)))
        #end if
        value = None
        secname = variable_section[name]
        if secname in self and name in self[secname]:
            value = self[secname][name]
        #end if
        return value
    #end def get


    def delete(self,name):
        if name not in known_variables:
            self.error('cannot get variable\nattempted to get unknown variable "{0}"\nvalid options are: {1}'.format(name,sorted(known_variables)))
        #end if
        value = None
        secname = variable_section[name]
        if secname in self and name in self[secname]:
            value = self[secname].delete(name)
        #end if
        return value
    #end def delete


    def extract_added_keys(self):
        extra = obj()
        extra.move_from(self,QuantumPackageInput.added_keys)
        return extra
    #end def extract_added_keys


    def restore_added_keys(self,extra):
        extra.move_to(self,QuantumPackageInput.added_keys)
    #end def restore_added_keys


    def read(self,filepath):
        epath = filepath.rstrip('/')
        if not os.path.exists(epath):
            self.error('cannot read input\nprovided ezfio directory does not exist\ndirectory provided:  {0}'.format(epath))
        elif not os.path.isdir(epath):
            self.error('cannot read input\nprovided ezfio path is not a directory\npath provided:  {0}'.format(epath))
        elif not epath.endswith('.ezfio'):
            self.error('cannot read input\nprovided path does not end in an ezfio directory\ndirectory must end with .ezfio\npath provided:  {0}'.format(epath))
        #end if
        for path,dirs,files in os.walk(epath):
            for file in files:
                if 'save' not in path:
                    if not file.startswith('.') and not file.endswith('.gz'):
                        filepath = os.path.join(path,file)
                        v = read_qp_value(filepath)
                        vpath = filepath.replace(epath,'').strip('/')
                        secname,varname = vpath.split('/')
                        if secname not in self:
                            self[secname] = Section()
                        #end if
                        self[secname][varname] = v
                    #end if
                #end if
            #end for
        #end for
    #end def read


    def write(self,filepath=None):
        if filepath is None:
            return str(self)
        #end if

        # get the ezfio path
        epath = filepath.rstrip('/')

        # check that write can occur
        if not epath.endswith('.ezfio'):
            self.error('cannot write input\nprovided path does not end in an ezfio directory\ndirectory must end with .ezfio\npath provided:  {0}'.format(epath))
        #end if
        path,edir = os.path.split(epath)
        if path=='':
            path = './'
        #end if
        if not os.path.exists(path):
            self.error('cannot write input\nattempted to write ezfio directory "{0}" at non-existent destination path\ndestination path: {1}'.format(edir,path))
        #end if

        # if there is no ezfio directory, initialize one
        if not os.path.exists(epath):
            if self.structure is None:
                self.error('cannot write input\nstructure is missing\ninput path provided: {0}'.format(epath))
            elif not isinstance(self.structure,Structure):
                self.error('cannot write input\nstructure must be of type: Structure\ntype provided: {0}\ninput path provided: {1}'.format(self.structure.__class__.__name__,epath))
            #end if
            cwd = os.getcwd()
            os.chdir(path)
            prefix = edir.rsplit('.',1)[0]
            struct_file = prefix+'.xyz'
            self.structure.write_xyz(struct_file)
            command = 'qp_create_ezfio_from_xyz'
            if self.path_exists('ao_basis/ao_basis'):
                command += ' -b '+self.ao_basis.ao_basis
            #end if
            command += ' '+struct_file
            execute(command)
            if not os.path.exists(edir):
                self.error('cannot write input\nezfio creation command failed: {0}\nexecuted at path: {1}\ndirectory {2} not created\nplease source your quantum_package.rc file before running the current script'.format(command,path,edir))
            #end if
            execute('qp_edit -c '+edir)
            os.chdir(cwd)
        #end if

        # write inputs into the ezfio directory/file tree
        extra = self.extract_added_keys()
        for secname,sec in self.iteritems():
            secpath = os.path.join(epath,secname)
            if not os.path.exists(secpath):
                self.error('cannot write input\ninput section path does not exist\nsection path: {0}\nplease ensure that all variables were created previously for this ezfio directory\n(to create all variables, run "qp_edit -c {1}")'.format(secpath,edir))
            #end if
            for varname,val in sec.iteritems():
                vpath = os.path.join(secpath,varname)
                write_qp_value(vpath,val)
            #end for
        #end for
        self.restore_added_keys(extra)

        # take other steps that modify the input, as requested
        if 'run_control' in self:
            rc = self.run_control
            if 'frozen_core' in rc and rc.frozen_core:
                cwd = os.getcwd()
                os.chdir(path)
                execute('qp_set_frozen_core.py '+edir)
                os.chdir(cwd)
            #end if
        #end if

        return ''
    #end def write


    def read_text(self,text,filepath=None):
        self.not_implemented()
    #end def read_text


    def write_text(self,filepath=None):
        self.not_implemented()
    #end def write_text


    def incorporate_system(self,system):
        self.not_implemented()
    #end def incorporate_system


    def check_valid(self,sections=True,variables=True,types=True,exit=True):
        msg = ''

        extra = self.extract_added_keys()

        valid_types = {float:(int,float)}
        if sections:
            for secname,sec in self.iteritems():
                if secname not in known_sections:
                    msg = 'input is invalid\nunknown section encountered\nunknown section provided: {0}\nvalid options are: {1}'.format(secname,sorted(known_sections))
                elif not isinstance(sec,Section):
                    msg = 'input is invalid\ninvalid section type encountered\nsection must be of type: Section\nsection name: {0}\nsection type: {1}\nsection contents: {2}'.format(secname,sec.__class__.__name__,sec)
                #end if
                if len(msg)>0:
                    break
                #end if
                if variables:
                    for varname,var in sec.iteritems():
                        if varname not in known_variables:
                            msg = 'input is invalid\nunknown variable encountered in section "{0}"\nunknown variable: {1}\nvalid options are: {2}'.format(secname,varname,sorted(section_variables[secname]))
                        elif types:
                            vpath = secname+'/'+varname
                            if vpath not in input_specification:
                                msg = 'variable is known but variable path not found in input_specification\nvariable name: {0}\nvariable path: {1}\nthis is a developer error\nplease contact the developers'.format(varname,vpath)
                            else:
                                vtype = input_specification[vpath]
                                if vtype in valid_types:
                                    vtype = valid_types[vtype]
                                #end if
                                if not isinstance(var,vtype):
                                    type_map = {bool:'bool',int:'int',float:'float',str:'str'}
                                    msg = 'input is invalid\nvariable "{0}" in section "{1}" must be of type {2}\ntype provided: {3}\nvalue provided: {4}'.format(varname,secname,type_map[vtype],var.__class__.__name__,var)
                                #end if
                            #end if
                        #end if
                        if len(msg)>0:
                            break
                        #end if
                    #end for
                    if len(msg)>0:
                        break
                    #end if
                #end if
            #end for
        #end if
        is_valid = len(msg)==0
        if not is_valid and exit:
            self.error(msg)
        #end if

        self.restore_added_keys(extra)

        return is_valid
    #end check_valid


    def is_valid(self):
        return self.check_valid(exit=False)
    #end def is_valid
                
#end class QuantumPackageInput



run_inputs = set('''
    prefix
    run_type
    frozen_core
    cis_loop
    sleep
    slave
    postprocess
    '''.split())
gen_inputs = set('''
    system
    defaults
    save_integrals
    validate
    '''.split())
added_inputs = run_inputs | gen_inputs
added_types = obj(
    # run inputs
    prefix         = str,
    run_type       = str,
    frozen_core    = bool,
    cis_loop       = (bool,int),
    sleep          = (int,float),
    slave          = str,
    postprocess    = (tuple,list),
    # gen inputs
    system         = PhysicalSystem,
    defaults       = str,
    save_integrals = bool,
    validate       = bool,
    )
added_required = set('''
    system
    prefix
    sleep
    '''.split())
qp_defaults_version = 'v1'
shared_defaults = obj(
    # run inputs
    postprocess = [],
    # gen inputs
    validate    = True,
    )
qp_defaults = obj(
    none = obj(
        # gen inputs
        save_integrals = False,
        **shared_defaults
        ),
    v1 = obj(
        # run inputs
        sleep          = 30,
        # gen inputs
        save_integrals = True,
        # qp inputs
        n_det_max      = 5000,
        **shared_defaults
        ),
    )
save_ints_defaults = obj(
    disk_access_ao_one_integrals = 'Write',
    disk_access_mo_one_integrals = 'Write',
    disk_access_ao_integrals     = 'Write',
    disk_access_mo_integrals     = 'Write',
    )

def generate_quantum_package_input(**kwargs):

    # make empty input
    qpi = QuantumPackageInput()

    # rewrap keywords and apply defaults
    kw  = obj(**kwargs)
    kw.set_optional(defaults=qp_defaults_version)
    if kw.defaults not in qp_defaults:
        QuantumPackageInput.class_error('cannot generate input\nrequested invalid default set\ndefault set requested: {0}\nvalid options are: {1}'.format(kw.defaults,sorted(qp_defaults.keys())))
    #end if
    kw.set_optional(**qp_defaults[kw.defaults])

    # check for required variables
    req_missing = kw.check_required(added_required,exit=False)
    if len(req_missing)>0:
        QuantumPackageInput.class_error('cannot generate input\nrequired variables are missing\nmissing variables: {0}\nplease supply values for these variables via generate_quantum_package'.format(sorted(req_missing)))
    #end if

    # check types of added variables
    name,vtype = kw.check_types_optional(added_types,exit=False)
    if name is not None:
        QuantumPackageInput.class_error('cannot generate input\nvariable "{0}" has the wrong type\ntype required: {1}\ntype provided: {2}'.format(name,vtype.__name__,kw[name].__class__.__name__))
    #end if

    # separate run inputs from input file variables
    run_kw = kw.extract_optional(run_inputs)
    if run_kw.run_type not in QuantumPackageInput.run_types:
        valid = ''
        for rt in sorted(QuantumPackageInput.run_types):
            valid += '  '+rt+'\n'
        #end for
        QuantumPackageInput.class_error('cannot generate input\ninvalid run_type requested\nrun_type provided: {0}\nvalid options are:\n{1}'.format(run_kw.run_type,valid))
    #end if
    qpi.run_control.set(**run_kw)

    # separate generation inputs from input file variables
    gen_kw = kw.extract_optional(gen_inputs)

    if gen_kw.save_integrals and run_kw.run_type in QuantumPackageInput.integral_write_allowed:
        kw.set_optional(**save_ints_defaults)
    #end if

    # partition inputs into sections and variables
    sections = obj()
    variables = obj()
    for name,value in kw.iteritems():
        is_sec = name in known_sections
        is_var = name in known_variables
        if is_sec and is_var:
            if isinstance(value,(obj,dict)):
                sections[name] = value
            else:
                variables[name] = value
            #end if
        elif is_sec:
            sections[name] = value
        elif is_var:
            variables[name] = value
        else:
            QuantumPackageInput.class_error('cannot generate input\nencountered name that is not known as a section or variable\nunrecognized name provided: {0}\nvalid sections: {1}\nvalid variables: {2}'.format(name,sorted(known_sections),sorted(known_variables)))
        #end if
    #end for

    # assign sections
    for secname,sec in sections.iteritems():
        if isinstance(sec,(obj,dict)):
            sec = Section(sec) # defer checking to check_valid
        #end if
        qpi[secname] = sec
    #end for

    # assign variables to sections
    for varname,var in variables.iteritems():
        if varname not in variable_section:
            QuantumPackageInput.class_error('cannot generate input\nsection cannot be fond for variable provided\nunrecognized variable: {0}'.format(varname))
        #end if
        secname = variable_section[varname]
        if isinstance(secname,tuple):
            QuantumPackageInput.class_error('cannot generate input\nsection cannot be uniquely determined from variable name\nvariable name provided: {0}\npossible sections: {1}\nplease provide this variable directly within on of the input sections listed and try again'.format(varname,secname))
        #end if
        if secname not in qpi:
            qpi[secname] = Section()
        #end if
        sec = qpi[secname][varname] = var
    #end for

    # incorporate atomic and electronic structure
    system = gen_kw.system
    qpi.structure = system.structure.copy()
    if 'electrons' not in qpi:
        qpi.electrons = Section()
    #end if
    nup,ndn = system.particles.electron_counts()
    qpi.electrons.elec_alpha_num = nup
    qpi.electrons.elec_beta_num  = ndn

    # validate the input
    if gen_kw.validate:
        qpi.check_valid()
    #end if

    return qpi
#end def generate_quantum_package_input
