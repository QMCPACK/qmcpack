
import testing
from testing import value_eq,object_eq,text_eq,check_object_eq


def format_value(v):
    import numpy as np
    s = ''
    if isinstance(v,np.ndarray):
        pad = 12*' '
        s = 'np.array([\n'
        if len(v.shape)==1:
            s += pad
            for vv in v:
                s += format_value(vv)+','
            #end for
            s = s[:-1]
        else:
            for vv in v:
                s += pad + format_value(list(vv))+',\n'
            #end for
            s = s[:-2]
        #end if
        s += '])'
    elif isinstance(v,(str,np.string_)):
        s = "'"+v+"'"
    else:
        s = str(v)
    #end if
    return s
#end def format_value


def make_serial_reference(qi):
    s = qi.serial()
    ref = '    ref = {\n'
    for k in sorted(s.keys()):
        v = s[k]
        ref +="        '{}' : {},\n".format(k,format_value(v))
    #end for
    ref += '        }\n'
    return ref
#end def make_serial_reference


serial_references = dict()


def generate_serial_references():
    import numpy as np
    from generic import obj
    
    # references for read
    serial_references['h2o.ezfio read'] = {
        'ao_basis/ao_basis' : 'cc-pvtz',
        'ao_basis/ao_cartesian' : False,
        'ao_basis/ao_md5' : 'b0e5878a56051339b81909a4b36ceeef',
        'ao_basis/ao_num' : 65,
        'ao_one_e_ints/io_ao_integrals_e_n' : 'None',
        'ao_one_e_ints/io_ao_integrals_kinetic' : 'None',
        'ao_one_e_ints/io_ao_integrals_overlap' : 'None',
        'ao_one_e_ints/io_ao_integrals_pseudo' : 'None',
        'ao_one_e_ints/io_ao_one_e_integrals' : 'None',
        'ao_two_e_erf_ints/io_ao_two_e_integrals_erf' : 'None',
        'ao_two_e_erf_ints/mu_erf' : 0.5,
        'ao_two_e_ints/direct' : False,
        'ao_two_e_ints/io_ao_two_e_integrals' : 'None',
        'ao_two_e_ints/threshold_ao' : 1e-15,
        'becke_numerical_grid/grid_type_sgn' : 2,
        'davidson/davidson_sze_max' : 15,
        'davidson/disk_based_davidson' : True,
        'davidson/distributed_davidson' : True,
        'davidson/n_det_max_full' : 1000,
        'davidson/n_states_diag' : 4,
        'davidson/only_expected_s2' : True,
        'davidson/state_following' : False,
        'davidson/threshold_davidson' : 1e-10,
        'density_for_dft/damping_for_rs_dft' : 0.5,
        'density_for_dft/density_for_dft' : 'WFT',
        'density_for_dft/no_core_density' : 'full_density',
        'determinants/n_det' : 1,
        'determinants/n_det_max' : 5000,
        'determinants/n_det_print_wf' : 10000,
        'determinants/n_states' : 1,
        'determinants/read_wf' : False,
        'determinants/s2_eig' : True,
        'determinants/target_energy' : 0.0,
        'determinants/threshold_generators' : 0.99,
        'determinants/used_weight' : 1,
        'dft_keywords/correlation_functional' : 'short_range_LDA',
        'dft_keywords/exchange_functional' : 'short_range_LDA',
        'dft_keywords/hf_exchange' : 0.0,
        'dressing/dress_relative_error' : 0.001,
        'dressing/n_it_max_dressed_ci' : 10,
        'dressing/thresh_dressed_ci' : 1e-05,
        'electrons/elec_alpha_num' : 5,
        'electrons/elec_beta_num' : 5,
        'ezfio/creation' : 'Mon Nov  4 08:18:23 EST 2019',
        'ezfio/library' : '/home/user/apps/quantum_package/qp2-2.0.0-beta/external/ezfio',
        'ezfio/user' : 'user',
        'ijkl_ints_in_r3/disk_access_ao_ijkl_r3' : 'None',
        'ijkl_ints_in_r3/disk_access_mo_ijkl_r3' : 'None',
        'mo_one_e_ints/io_mo_integrals_e_n' : 'None',
        'mo_one_e_ints/io_mo_integrals_kinetic' : 'None',
        'mo_one_e_ints/io_mo_integrals_pseudo' : 'None',
        'mo_one_e_ints/io_mo_one_e_integrals' : 'None',
        'mo_two_e_erf_ints/io_mo_two_e_integrals_erf' : 'None',
        'mo_two_e_ints/io_mo_two_e_integrals' : 'None',
        'mo_two_e_ints/no_ivvv_integrals' : False,
        'mo_two_e_ints/no_vvv_integrals' : False,
        'mo_two_e_ints/no_vvvv_integrals' : False,
        'mo_two_e_ints/threshold_mo' : 1e-15,
        'nuclei/disk_access_nuclear_repulsion' : 'None',
        'nuclei/nucl_num' : 3,
        'perturbation/correlation_energy_ratio_max' : 1.0,
        'perturbation/do_pt2' : True,
        'perturbation/h0_type' : 'EN',
        'perturbation/pt2_max' : 0.0001,
        'perturbation/pt2_relative_error' : 0.002,
        'perturbation/variance_max' : 0.0,
        'pseudo/do_pseudo' : False,
        'pseudo/pseudo_grid_rmax' : 10.0,
        'pseudo/pseudo_grid_size' : 1000,
        'pseudo/pseudo_klocmax' : 0,
        'pseudo/pseudo_kmax' : 0,
        'pseudo/pseudo_lmax' : 0,
        'qmcpack/ci_threshold' : 1e-08,
        'rsdft_ecmd/ecmd_functional' : 'short_range_LDA',
        'run_control' : obj(),
        'scf_utils/frozen_orb_scf' : False,
        'scf_utils/level_shift' : 0.0,
        'scf_utils/max_dim_diis' : 15,
        'scf_utils/mo_guess_type' : 'Huckel',
        'scf_utils/n_it_scf_max' : 500,
        'scf_utils/scf_algorithm' : 'DIIS',
        'scf_utils/thresh_scf' : 1e-10,
        'scf_utils/threshold_diis' : 0.0,
        'structure' : None,
        'two_body_dm/ci_threshold' : 1e-05,
        'two_body_dm/mat_mul_svd_vectors' : True,
        'two_body_dm/ontop_approx' : False,
        'two_body_dm/thr_ontop_approx' : 0.001,
        'work/empty' : False,
        }


    # references for generate

    serial_references['h2o.ezfio gen'] = {
        'ao_basis/ao_basis' : 'cc-pvtz',
        'determinants/n_det_max' : 5000,
        'electrons/elec_alpha_num' : 5,
        'electrons/elec_beta_num' : 5,
        'run_control/four_idx_transform' : False,
        'run_control/postprocess' : [],
        'run_control/prefix' : 'h2o',
        'run_control/run_type' : 'scf',
        'run_control/save_for_qmcpack' : False,
        'run_control/save_natorb' : False,
        'run_control/sleep' : 30,
        'structure/axes' : np.array([]),
        'structure/background_charge' : 0,
        'structure/bconds' : np.array([]),
        'structure/center' : np.array([0.0,0.0,0.0]),
        'structure/dim' : 3,
        'structure/elem' : np.array(['O','H','H']),
        'structure/folded_structure' : None,
        'structure/frozen' : None,
        'structure/kaxes' : np.array([]),
        'structure/kpoints' : np.array([]),
        'structure/kweights' : np.array([]),
        'structure/mag' : None,
        'structure/pos' : np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.75716, 0.58626],
            [0.0, 0.75716, -0.58626]]),
        'structure/scale' : 1.0,
        'structure/units' : 'A',
        }

    serial_references['o2.ezfio gen'] = {
        'determinants/n_det_max' : 5000,
        'electrons/elec_alpha_num' : 8,
        'electrons/elec_beta_num' : 8,
        'run_control/four_idx_transform' : True,
        'run_control/postprocess' : [],
        'run_control/prefix' : 'o2',
        'run_control/run_type' : 'fci',
        'run_control/save_for_qmcpack' : False,
        'run_control/save_natorb' : True,
        'run_control/sleep' : 30,
        'structure/axes' : np.array([]),
        'structure/background_charge' : 0,
        'structure/bconds' : np.array([]),
        'structure/center' : np.array([0.0,0.0,0.0]),
        'structure/dim' : 3,
        'structure/elem' : np.array(['O','O']),
        'structure/folded_structure' : None,
        'structure/frozen' : None,
        'structure/kaxes' : np.array([]),
        'structure/kpoints' : np.array([]),
        'structure/kweights' : np.array([]),
        'structure/mag' : None,
        'structure/pos' : np.array([
            [0.0, 0.0, 0.0],
            [1.2074, 0.0, 0.0]]),
        'structure/scale' : 1.0,
        'structure/units' : 'A',
        }

#end def generate_serial_references


def get_serial_references():
    if len(serial_references)==0:
        generate_serial_references()
    #end if
    return serial_references
#end def get_serial_references


def check_vs_serial_reference(qi,name):
    from generic import obj
    sr = obj(get_serial_references()[name])
    sg = qi.serial()
    assert(object_eq(sg,sr))
#end def check_vs_serial_reference



h2o_xyz = '''3

O  0.000000  0.000000  0.000000 
H  0.000000  0.757160  0.586260
H  0.000000  0.757160 -0.586260
'''

o2_xyz = '''2

 O    0.00000000   0.00000000   0.00000000
 O    1.20740000   0.00000000   0.00000000
'''


def test_import():
    from quantum_package_input import QuantumPackageInput
    from quantum_package_input import generate_quantum_package_input
#end def test_import



def test_empty_init():
    from quantum_package_input import QuantumPackageInput
    from quantum_package_input import generate_quantum_package_input

    qi = QuantumPackageInput()

#end def test_empty_init



def test_read():
    import os
    from quantum_package_input import QuantumPackageInput

    tpath = testing.setup_unit_test_output_directory('quantum_package_input','test_generate',file_sets=['h2o.ezfio'])

    ezfio = os.path.join(tpath,'h2o.ezfio')

    qi = QuantumPackageInput(ezfio)

    check_vs_serial_reference(qi,'h2o.ezfio read')
#end def test_read



def test_generate():
    import os
    from generic import obj
    from physical_system import generate_physical_system
    from quantum_package_input import generate_quantum_package_input

    tpath = testing.setup_unit_test_output_directory('quantum_package_input','test_generate')

    # water molecule rhf
    xyz_path = os.path.join(tpath,'H2O.xyz')
    open(xyz_path,'w').write(h2o_xyz)

    system = generate_physical_system(
        structure = xyz_path,
        )

    qi = generate_quantum_package_input(
        system   = system,
        prefix   = 'h2o',
        run_type = 'scf',
        ao_basis = 'cc-pvtz',
        )

    check_vs_serial_reference(qi,'h2o.ezfio gen')

    assert(qi.is_valid())


    # O2 molecule selci
    xyz_path = os.path.join(tpath,'O2.xyz')
    open(xyz_path,'w').write(o2_xyz)

    system = generate_physical_system(
        structure = xyz_path,
        )

    qi = generate_quantum_package_input(
        system             = system,
        prefix             = 'o2',
        run_type           = 'fci',
        n_det_max          = 5000,
        save_natorb        = True,
        four_idx_transform = True,
        )

    check_vs_serial_reference(qi,'o2.ezfio gen')

    assert(qi.is_valid())
#end def test_generate

    
