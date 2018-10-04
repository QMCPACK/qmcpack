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
from simulation import SimulationInput


bool_values = dict(T=True,F=False)

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
    def __init__(self,filepath=None):
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__


    def read(self,filepath):
        epath = filepath.rstrip('/')
        if not os.path.exists(epath):
            self.error('provided ezfio directory does not exist\ndirectory provided:  {0}'.format(epath))
        elif not os.path.isdir(epath):
            self.error('provided ezfio path is not a directory\npath provided:  {0}'.format(epath))
        elif not epath.endswith('.ezfio'):
            self.error('provided path does not end in an ezfio directory\ndirectory must end with .ezfio\npath provided:  {0}'.format(epath))
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
                                    msg = 'input is invalid\nvariable "{0}" in section "{1}" must be of type {2}\ntype provided: {3}\nvalue provided: {4}'.format(varname,secname,vtype,var.__class__.__name__,var)
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
        return is_valid
    #end check_valid


    def is_valid(self):
        return self.check_valid(exit=False)
    #end def is_valid
                
#end class QuantumPackageInput



def generate_quantum_package_input(**kwargs):

    # make empty input
    qpi = QuantumPackageInput()

    # partition inputs into sections and variables
    sections = obj()
    variables = obj()
    for name,value in kwargs.iteritems():
        is_sec = name in known_sections
        is_var = name in known_variables
        if is_sec and is_var:
            if isinstance(value,section_types):
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
            sec = Section(sec) # defer checking to validate
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

    # validate the input
    qpi.validate()

    return qpi
#end def generate_quantum_package_input
