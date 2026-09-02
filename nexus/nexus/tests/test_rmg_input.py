import pytest
import numpy as np
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.RMG_INPUT)


from importlib.util import find_spec
from . import isolate_nexus_core, register_pseudo_files, TEST_DIR
from ..developer import obj
from ..physical_system import generate_physical_system
from ..rmg_input import RmgInput,generate_rmg_input,input_spec,rmg_modes
from ..testing import value_eq,check_object_eq,dict_serialize
from ..unit_converter import convert

TEST_FILES = {
    "AlN32_input":                                                 TEST_DIR / "test_rmg_input_files/AlN32_input",
    "atomO_polarized_input":                                       TEST_DIR / "test_rmg_input_files/atomO_polarized_input",
    "BlackPhosphorus_Delocalize_both_pp_and_proj_davidson_input":  TEST_DIR / "test_rmg_input_files/BlackPhosphorus_Delocalize_both_pp_and_proj_davidson_input",
    "BlackPhosphorus_Delocalize_both_pp_and_proj_multigrid_input": TEST_DIR / "test_rmg_input_files/BlackPhosphorus_Delocalize_both_pp_and_proj_multigrid_input",
    "BlackPhosphorus_input":                                       TEST_DIR / "test_rmg_input_files/BlackPhosphorus_input",
    "BlackPhosphorus_input_1nongammapoint":                        TEST_DIR / "test_rmg_input_files/BlackPhosphorus_input_1nongammapoint",
    "BlackPhosphorus_input_band":                                  TEST_DIR / "test_rmg_input_files/BlackPhosphorus_input_band",
    "BlackPhosphorus_Localize_both_pp_and_proj_davidson_input":    TEST_DIR / "test_rmg_input_files/BlackPhosphorus_Localize_both_pp_and_proj_davidson_input",
    "BlackPhosphorus_Localize_both_pp_and_proj_multigrid_input":   TEST_DIR / "test_rmg_input_files/BlackPhosphorus_Localize_both_pp_and_proj_multigrid_input",
    "BlackPhosphorus_Localize_pp_only_davidson_input":             TEST_DIR / "test_rmg_input_files/BlackPhosphorus_Localize_pp_only_davidson_input",
    "BlackPhosphorus_Localize_pp_only_multigrid_input":            TEST_DIR / "test_rmg_input_files/BlackPhosphorus_Localize_pp_only_multigrid_input",
    "BlackPhosphorus_Localize_proj_only_davidson_input":           TEST_DIR / "test_rmg_input_files/BlackPhosphorus_Localize_proj_only_davidson_input",
    "BlackPhosphorus_Localize_proj_only_multigrid_input":          TEST_DIR / "test_rmg_input_files/BlackPhosphorus_Localize_proj_only_multigrid_input",
    "C60_input":                                                   TEST_DIR / "test_rmg_input_files/C60_input",
    "CO_qmcpack_semilocal_xml_input":                              TEST_DIR / "test_rmg_input_files/CO_qmcpack_semilocal_xml_input",
    "Diamond16_input":                                             TEST_DIR / "test_rmg_input_files/Diamond16_input",
    "Diamond2_input":                                              TEST_DIR / "test_rmg_input_files/Diamond2_input",
    "Fe_2atom_input":                                              TEST_DIR / "test_rmg_input_files/Fe_2atom_input",
    "graphite_stress_input":                                       TEST_DIR / "test_rmg_input_files/graphite_stress_input",
    "H2O_tddft_electric_field_input":                              TEST_DIR / "test_rmg_input_files/H2O_tddft_electric_field_input",
    "H2O_tddft_point_charge_input":                                TEST_DIR / "test_rmg_input_files/H2O_tddft_point_charge_input",
    "Mg_2atom_input":                                              TEST_DIR / "test_rmg_input_files/Mg_2atom_input",
    "nanotube_80_input":                                           TEST_DIR / "test_rmg_input_files/nanotube_80_input",
    "nanotube_80_input_band":                                      TEST_DIR / "test_rmg_input_files/nanotube_80_input_band",
    "nanotube_80_input_band1":                                     TEST_DIR / "test_rmg_input_files/nanotube_80_input_band1",
    "NiO512_input":                                                TEST_DIR / "test_rmg_input_files/NiO512_input",
    "NiO8_input":                                                  TEST_DIR / "test_rmg_input_files/NiO8_input",
    "PbTiO3_berry_phase_nscf_input":                               TEST_DIR / "test_rmg_input_files/PbTiO3_berry_phase_nscf_input",
    "Pt_bulk_spinorbit_input":                                     TEST_DIR / "test_rmg_input_files/Pt_bulk_spinorbit_input",
    "Pt_bulk_spinorbit_input_band":                                TEST_DIR / "test_rmg_input_files/Pt_bulk_spinorbit_input_band",
    "Si_8atoms_EXX_kpoints_input":                                 TEST_DIR / "test_rmg_input_files/Si_8atoms_EXX_kpoints_input",
    "Si_wannier90_nscf_input":                                     TEST_DIR / "test_rmg_input_files/Si_wannier90_nscf_input",
    "U_bulk_spinorbit_RMG_input":                                  TEST_DIR / "test_rmg_input_files/U_bulk_spinorbit_RMG_input",
    }

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file not found! {file}"


def make_serial_reference(ri):
    s = dict_serialize(ri,dict_type=obj)
    ref = '    ref = {\n'
    for k in sorted(s.keys()):
        v = s[k]
        if isinstance(v,str):
            v = "'"+v+"'"
        #end if
        if not isinstance(v,np.ndarray) or len(v)!=v.size:
            ref +="        '{}' : {},\n".format(k,v)
        else:
            a = 'np.array({})'.format(v)
            a = a.replace('     ','    ,')
            a = a.replace('    ','   ,')
            a = a.replace('   ','  ,')
            a = a.replace('  ',' ,')
            a = a.replace(' ',',')
            a = a.replace(',,,,,','    ,')
            a = a.replace(',,,,','   ,')
            a = a.replace(',,,','  ,')
            a = a.replace(',,',' ,')
            ref +="        '{}' : {},\n".format(k,a)
        #end if
    #end for
    ref += '        }\n'
    return ref
#end def make_serial_reference


serial_references = dict()


def generate_serial_references():
    serial_references['BlackPhosphorus_input'] = {
        'a_length' : 3.3136,
        'atomic_coordinate_type' : 'Absolute',
        'atoms/atoms' : np.array(['P','P','P','P']),
        'atoms/format' : 'movable',
        'atoms/movable' : np.array([True ,True ,True ,True]),
        'atoms/positions' : np.array([[0.    ,  0.     , 0.     ],
                                      [1.6568,  0.     , 1.48364],
                                      [1.6568,  2.13054, 2.18815],
                                      [0.    ,  2.13054, 3.67179]]),
        'b_length' : 10.478,
        'bravais_lattice_type' : 'Orthorhombic Primitive',
        'c_length' : 4.3763,
        'calculation_mode' : 'Quench Electrons',
        'charge_density_mixing' : 0.5,
        'charge_mixing_type' : 'Linear',
        'charge_pulay_order' : 5,
        'charge_pulay_scale' : 0.5,
        'crds_units' : 'Angstrom',
        'description' : 'Black Phosphorus',
        'initial_diagonalization' : True,
        'kohn_sham_mg_levels' : 2,
        'kohn_sham_solver' : 'multigrid',
        'kpoint_distribution' : 1,
        'kpoint_is_shift' : np.array([0,0,0]),
        'kpoint_mesh' : np.array([2,1,2]),
        'length_units' : 'Angstrom',
        'localize_projectors' : False,
        'max_scf_steps' : 20,
        'occupation_electron_temperature_eV' : 0.025,
        'occupation_number_mixing' : 1.0,
        'occupations_type' : 'Fermi Dirac',
        'output_wave_function_file' : 'Waves/wave.out',
        'potential_acceleration_constant_step' : 1.0,
        'potential_grid_refinement' : 2,
        'rms_convergence_criterion' : 1e-07,
        'start_mode' : 'LCAO Start',
        'subdiag_driver' : 'lapack',
        'system_charge' : 0.0,
        'unoccupied_states_per_kpoint' : 10,
        }


    serial_references['BlackPhosphorus_input_band'] = {
        'a_length' : 3.3136,
        'atomic_coordinate_type' : 'Absolute',
        'atoms/atoms' : np.array(['P','P','P','P']),
        'atoms/format' : 'movable',
        'atoms/movable' : np.array([True ,True ,True ,True]),
        'atoms/positions' : np.array([[0.    ,  0.     , 0.     ], 
                                      [1.6568,  0.     , 1.48364],
                                      [1.6568,  2.13054, 2.18815],
                                      [0.    ,  2.13054, 3.67179]]),
        'b_length' : 10.478,
        'bravais_lattice_type' : 'Orthorhombic Primitive',
        'c_length' : 4.3763,
        'calculation_mode' : 'Band Structure Only',
        'charge_density_mixing' : 0.2,
        'charge_mixing_type' : 'Pulay',
        'charge_pulay_order' : 5,
        'charge_pulay_scale' : 0.2,
        'crds_units' : 'Angstrom',
        'description' : 'Black Phosphorus',
        'initial_diagonalization' : True,
        'kohn_sham_mg_levels' : 2,
        'kohn_sham_time_step' : 0.66,
        'kpoint_distribution' : 1,
        'kpoint_is_shift' : np.array([0,0,0]),
        'kpoint_mesh' : np.array([-4 ,2 ,4]),
        'kpoints_bandstructure/counts' : np.array([0,10,10]),
        'kpoints_bandstructure/kpoints' : np.array([[0.5, 0.,  0. ],
                                                    [0. , 0.,  0. ],
                                                    [0. , 0.,  0.5]]),
        'kpoints_bandstructure/labels' : np.array(['X','\\xG','Z']),
        'length_units' : 'Angstrom',
        'localize_projectors' : False,
        'max_scf_steps' : 10,
        'occupation_electron_temperature_eV' : 0.025,
        'occupation_number_mixing' : 1.0,
        'occupations_type' : 'Fermi Dirac',
        'potential_grid_refinement' : 2,
        'processor_grid' : np.array([1,2,1]),
        'projector_mixing' : 0.4,
        'rms_convergence_criterion' : 1e-07,
        'start_mode' : 'LCAO Start',
        'subdiag_driver' : 'lapack',
        'system_charge' : 0.0,
        'unoccupied_states_per_kpoint' : 10,
        'wavefunction_grid' : np.array([20,64,28]),
        }


    serial_references['BlackPhosphorus_Localize_both_pp_and_proj_multigrid_input'] = {
        'a_length' : 3.3136,
        'atomic_coordinate_type' : 'Absolute',
        'atoms/atoms' : np.array(['P','P','P','P']),
        'atoms/format' : 'movable',
        'atoms/movable' : np.array([True ,True ,True ,True]),
        'atoms/positions' : np.array([[0.    ,  0.     , 0.     ],
                                      [1.6568,  0.     , 1.48364],
                                      [1.6568,  2.13054, 2.18815],
                                      [0.    ,  2.13054, 3.67179]]),
        'b_length' : 10.478,
        'bravais_lattice_type' : 'Orthorhombic Primitive',
        'c_length' : 4.3763,
        'calculation_mode' : 'Quench Electrons',
        'charge_density_mixing' : 0.2,
        'charge_mixing_type' : 'Pulay',
        'charge_pulay_order' : 5,
        'charge_pulay_scale' : 0.5,
        'crds_units' : 'Angstrom',
        'description' : 'Black Phosphorus',
        'initial_diagonalization' : True,
        'kohn_sham_mg_levels' : 2,
        'kohn_sham_solver' : 'multigrid',
        'kohn_sham_time_step' : 0.66,
        'kpoint_is_shift' : np.array([0,0,0]),
        'kpoint_mesh' : np.array([2,1,2]),
        'length_units' : 'Angstrom',
        'localize_localpp' : True,
        'localize_projectors' : True,
        'max_scf_steps' : 20,
        'occupation_electron_temperature_eV' : 0.025,
        'occupation_number_mixing' : 1.0,
        'occupations_type' : 'Fermi Dirac',
        'output_wave_function_file' : 'Waves/wave.out',
        'potential_acceleration_constant_step' : 0.0,
        'potential_grid_refinement' : 2,
        'projector_mixing' : 0.4,
        'rms_convergence_criterion' : 1e-07,
        'start_mode' : 'LCAO Start',
        'subdiag_driver' : 'lapack',
        'system_charge' : 0.0,
        'unoccupied_states_per_kpoint' : 10,
        }


    serial_references['C60_input'] = {
        'a_length' : 18.0,
        'alpha' : 0.0,
        'atomic_coordinate_type' : 'Absolute',
        'atoms/atoms' : np.array(['C','C','C','C','C','C','C','C','C','C','C',
                                  'C','C','C','C','C','C','C','C','C','C','C',
                                  'C','C','C','C','C','C','C','C','C','C','C',
                                  'C','C','C','C','C','C','C','C','C','C','C',
                                  'C','C','C','C','C','C','C','C','C','C','C',
                                  'C','C','C','C','C']),
        'atoms/format' : 'movable',
        'atoms/movable' : np.array([True,True,True,True,True,True,True,True,
                                    True,True,True,True,True,True,True,True,
                                    True,True,True,True,True,True,True,True,
                                    True,True,True,True,True,True,True,True,
                                    True,True,True,True,True,True,True,True,
                                    True,True,True,True,True,True,True,True,
                                    True,True,True,True,True,True,True,True,
                                    True,True,True,True]),
        'atoms/positions' : np.array([[15.3293, 12.96  , 19.207 ], 
                                      [17.5363, 12.96  , 17.7857],
                                      [13.7313, 15.1732, 19.2251],
                                      [14.4083, 17.2893, 17.8251],
                                      [11.1292, 14.3262, 19.2527],
                                      [ 9.3055, 15.6378, 17.8881],
                                      [11.1292, 11.5939, 19.2527],
                                      [ 9.3055, 10.2823, 17.8881],
                                      [13.7313, 10.7469, 19.2251],
                                      [14.4083,  8.6308, 17.8251],
                                      [18.2386, 15.1636, 16.3387],
                                      [16.7171, 17.292 , 16.3627],
                                      [18.2386, 10.7565, 16.3387],
                                      [16.717 ,  8.6281, 16.3627],
                                      [19.3641, 14.32  , 13.9962],
                                      [18.9522, 15.6335, 11.7652],
                                      [19.3641, 11.6   , 13.9962],
                                      [18.9522, 10.2865, 11.7652],
                                      [12.5141, 18.6398, 16.4036],
                                      [10.0138, 17.8435, 16.439 ],
                                      [16.2562, 18.6408, 14.0355],
                                      [17.3583, 17.8445, 11.7869],
                                      [13.6568, 19.4834, 14.0698],
                                      [12.2629, 19.4836, 11.85  ],
                                      [ 7.4102, 14.2728, 16.4751],
                                      [ 7.4102, 11.6472, 16.4751],
                                      [ 8.5614, 17.8446, 14.1329],
                                      [ 9.6634, 18.6411, 11.8843],
                                      [ 6.9677, 15.6335, 14.1547],
                                      [ 6.5558, 14.32  , 11.9237],
                                      [10.0138,  8.0765, 16.439 ],
                                      [12.5141,  7.2802, 16.4036],
                                      [ 6.9677, 10.2865, 14.1547],
                                      [ 6.5558, 11.6   , 11.9237],
                                      [ 8.5614,  8.0754, 14.1329],
                                      [ 9.6635,  7.279 , 11.8844],
                                      [13.6568,  6.4367, 14.0698],
                                      [12.2629,  6.4364, 11.85  ],
                                      [16.2562,  7.2793, 14.0355],
                                      [17.3583,  8.0755, 11.7869],
                                      [15.906 , 17.8435,  9.4808],
                                      [13.4057, 18.6398,  9.5163],
                                      [18.5097, 14.2728,  9.4447],
                                      [18.5097, 11.6472,  9.4447],
                                      [16.6144, 15.6378,  8.0317],
                                      [14.7907, 14.3262,  6.6672],
                                      [ 9.2027, 17.2921,  9.5572],
                                      [ 7.6813, 15.1636,  9.5812],
                                      [11.5115, 17.2892,  8.0948],
                                      [12.1885, 15.1732,  6.6948],
                                      [ 7.6813, 10.7564,  9.5812],
                                      [ 9.2027,  8.6279,  9.5572],
                                      [ 8.3836, 12.96  ,  8.1343],
                                      [10.5906, 12.96  ,  6.7129],
                                      [13.4057,  7.2802,  9.5163],
                                      [15.906 ,  8.0765,  9.4808],
                                      [11.5115,  8.6308,  8.0948],
                                      [12.1885, 10.7468,  6.6948],
                                      [16.6144, 10.2823,  8.0317],
                                      [14.7907, 11.5939,  6.6672]]),
        'b_length' : 18.0,
        'beta' : 0.0,
        'bravais_lattice_type' : 'Cubic Primitive',
        'c_length' : 18.0,
        'calculation_mode' : 'Quench Electrons',
        'charge_density_mixing' : 0.7,
        'charge_mixing_type' : 'Broyden',
        'description' : 'C60 test example using Davidson diagonalization',
        'energy_convergence_criterion' : 1e-09,
        'gamma' : 0.0,
        'kohn_sham_mucycles' : 3,
        'kohn_sham_solver' : 'davidson',
        'max_scf_steps' : 20,
        'potential_acceleration_constant_step' : 1.0,
        'potential_grid_refinement' : 2,
        'preconditioner_threshold' : 0.0001,
        'start_mode' : 'LCAO Start',
        'subdiag_driver' : 'scalapack',
        'test_energy' : -343.83860554,
        'unoccupied_states_per_kpoint' : 18,
        'wavefunction_grid' : np.array([48,48,48]),
        'write_data_period' : 50,
        }


    serial_references['Diamond16_input'] = {
        'a_length' : 13.44,
        'alpha' : 0.0,
        'atomic_coordinate_type' : 'Cell Relative',
        'atoms/atoms' : np.array(['C','C','C','C','C','C','C','C','C','C','C','C','C','C','C','C']),
        'atoms/format' : 'moment',
        'atoms/moments' : np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]),
        'atoms/positions' : np.array([[0.   , 0.   , 0.   ],
                                      [0.125, 0.125, 0.125],
                                      [0.   , 0.   , 0.5  ],
                                      [0.125, 0.125, 0.625],
                                      [0.   , 0.5  , 0.   ],
                                      [0.125, 0.625, 0.125],
                                      [0.   , 0.5  , 0.5  ],
                                      [0.125, 0.625, 0.625],
                                      [0.5  , 0.   , 0.   ],
                                      [0.625, 0.125, 0.125],
                                      [0.5  , 0.   , 0.5  ],
                                      [0.625, 0.125, 0.625],
                                      [0.5  , 0.5  , 0.   ],
                                      [0.625, 0.625, 0.125],
                                      [0.5  , 0.5  , 0.5  ],
                                      [0.625, 0.625, 0.625]]),
        'b_length' : 13.44,
        'beta' : 0.0,
        'bravais_lattice_type' : 'Cubic Face Centered',
        'c_length' : 13.44,
        'calculation_mode' : 'Quench Electrons',
        'charge_density_mixing' : 0.7,
        'charge_mixing_type' : 'Broyden',
        'crds_units' : 'Bohr',
        'description' : 'Diamond 16 atom test cell',
        'energy_convergence_criterion' : 1e-12,
        'force_grad_order' : 0,
        'gamma' : 0.0,
        'ionic_time_step' : 10.0,
        'kohn_sham_solver' : 'multigrid',
        'kpoint_mesh' : np.array([2,2,2]),
        'localize_localpp' : False,
        'localize_projectors' : False,
        'max_scf_steps' : 40,
        'occupations_type' : 'Fixed',
        'potential_acceleration_constant_step' : 1.0,
        'potential_grid_refinement' : 2,
        'renormalize_forces' : False,
        'rms_convergence_criterion' : 1e-08,
        'start_mode' : 'LCAO Start',
        'subdiag_driver' : 'lapack',
        'test_energy' : -91.81065082,
        'unoccupied_states_per_kpoint' : 16,
        'wavefunction_grid' : np.array([32,32,32]),
        }


    serial_references['graphite_stress_input'] = {
        'a_length' : 4.64117,
        'alpha' : 0.0,
        'atomic_coordinate_type' : 'Cell Relative',
        'atoms/atoms' : np.array(['C','C','C','C']),
        'atoms/format' : 'basic',
        'atoms/positions' : np.array([[0.16666667, 0.16666667, 0.25      ],
                                      [0.5       , 0.83333333, 0.25      ],
                                      [0.16666667, 0.16666667, 0.75      ],
                                      [0.83333333, 0.5       , 0.75      ]]),
        'b_length' : 4.64117,
        'beta' : 0.0,
        'bravais_lattice_type' : 'Hexagonal Primitive',
        'c_length' : 12.653685887999998,
        'calculation_mode' : 'Relax Structure',
        'cell_movable' : np.array([1,1,0,1,1,0,0,0,1]),
        'cell_relax' : True,
        'charge_density_mixing' : 0.1,
        'charge_mixing_type' : 'Pulay',
        'description' : 'graphite',
        'energy_convergence_criterion' : 1e-09,
        'gamma' : 0.0,
        'kohn_sham_mucycles' : 3,
        'kohn_sham_solver' : 'davidson',
        'kpoint_distribution' : 1,
        'kpoint_is_shift' : np.array([1,1,1]),
        'kpoint_mesh' : np.array([4,4,2]),
        'localize_localpp' : False,
        'localize_projectors' : False,
        'max_md_steps' : 10,
        'max_scf_steps' : 20,
        'occupations_type' : 'Fixed',
        'potential_acceleration_constant_step' : 1.0,
        'potential_grid_refinement' : 2,
        'pseudo_dir' : '/home/luw/Pseudopotentials_collect',
        'start_mode' : 'LCAO Start',
        'stress' : True,
        'subdiag_driver' : 'lapack',
        'test_energy' : -22.98280508,
        'unoccupied_states_per_kpoint' : 16,
        'verbose' : True,
        'wavefunction_grid' : np.array([24,24,64]),
        'write_data_period' : 50,
        }


    serial_references['NiO8_input'] = {
        'Hubbard_U/Ni' : 6.5,
        'a_length' : 7.8811,
        'alpha' : 0.0,
        'atomic_coordinate_type' : 'Cell Relative',
        'atoms/atoms' : np.array(['O','O','O','O','Ni','Ni','Ni','Ni']),
        'atoms/format' : 'movable_moment',
        'atoms/moments' : np.array([0.   ,0.   ,0.   ,0.   ,0.25 ,0.25,-0.25,-0.25]),
        'atoms/movable' : np.array([True ,True ,True ,True ,True ,True ,True ,True]),
        'atoms/positions' : np.array([[0.   , 0. ,   0.5  ],
                                      [0.   , 0.5,   0.   ],
                                      [0.5  , 0. ,   0.   ],
                                      [0.5  , 0.5,   0.5  ],
                                      [0.005, 0. ,   0.   ],
                                      [0.5  , 0.5,   0.   ],
                                      [0.   , 0.5,   0.5  ],
                                      [0.5  , 0. ,   0.5  ]]),
        'b_length' : 7.8811,
        'beta' : 0.0,
        'bravais_lattice_type' : 'Cubic Primitive',
        'c_length' : 7.8811,
        'calculation_mode' : 'Quench Electrons',
        'charge_density_mixing' : 0.25,
        'charge_mixing_type' : 'Broyden',
        'davidson_max_steps' : 15,
        'davidson_multiplier' : 4,
        'description' : 'NiO 8 atom cell in anti-ferromagnetic configuration solved using Davidson diagonalization',
        'energy_convergence_criterion' : 1e-10,
        'force_grad_order' : 0,
        'gamma' : 0.0,
        'kohn_sham_mg_levels' : 3,
        'kohn_sham_solver' : 'davidson',
        'kpoint_is_shift' : np.array([0,0,0]),
        'kpoint_mesh' : np.array([2,2,2]),
        'ldaU_mode' : 'Simple',
        'localize_localpp' : False,
        'localize_projectors' : False,
        'max_scf_steps' : 100,
        'occupations_type' : 'MethfesselPaxton',
        'potential_acceleration_constant_step' : 1.0,
        'potential_grid_refinement' : 2,
        'pseudopotential/pseudos' : np.array(['Ni_oncv.UPF','O_oncv.UPF']),
        'pseudopotential/species' : np.array(['Ni','O']),
        'rms_convergence_criterion' : 1e-08,
        'start_mode' : 'LCAO Start',
        'states_count_and_occupation_spin_down' : '48 1.0 8 0.0',
        'states_count_and_occupation_spin_up' : '48 1.0 8 0.0',
        'subdiag_driver' : 'lapack',
        'test_energy' : -677.71977422,
        'wavefunction_grid' : np.array([36,36,36]),
        'write_data_period' : 5,
        }


    serial_references['Pt_bulk_spinorbit_input'] = {
        'a_length' : 7.42,
        'alpha' : 0.0,
        'atomic_coordinate_type' : 'Absolute',
        'atoms/atoms' : np.array(['Pt']),
        'atoms/format' : 'full_spin',
        'atoms/movable' : np.array([[ True,  True,  True]]),
        'atoms/positions' : np.array([[0., 0., 0.]]),
        'atoms/spin_phi' : np.array([0.]),
        'atoms/spin_ratio' : np.array([0.]),
        'atoms/spin_theta' : np.array([90.]),
        'b_length' : 7.42,
        'beta' : 0.0,
        'bravais_lattice_type' : 'Cubic Face Centered',
        'c_length' : 7.42,
        'calculation_mode' : 'Quench Electrons',
        'charge_density_mixing' : 0.5,
        'charge_mixing_type' : 'Broyden',
        'compressed_infile' : False,
        'compressed_outfile' : False,
        'description' : 'atom_Pt_pp',
        'energy_convergence_criterion' : 1e-09,
        'gamma' : 0.0,
        'kohn_sham_mucycles' : 3,
        'kohn_sham_solver' : 'davidson',
        'kpoint_is_shift' : np.array([1,1,1]),
        'kpoint_mesh' : np.array([4,4,4]),
        'localize_localpp' : False,
        'localize_projectors' : False,
        'max_scf_steps' : 100,
        'noncollinear' : True,
        'occupation_electron_temperature_eV' : 0.2,
        'occupations_type' : 'MethfesselPaxton',
        'potential_acceleration_constant_step' : 1.0,
        'potential_grid_refinement' : 2,
        'pseudo_dir' : './',
        'pseudopotential/pseudos' : np.array(['Pt.rel-pbe-n-rrkjus.UPF']),
        'pseudopotential/species' : np.array(['Pt']),
        'spinorbit' : True,
        'start_mode' : 'LCAO Start',
        'states_count_and_occupation' : '10 1.0 10 0.0',
        'subdiag_driver' : 'lapack',
        'test_energy' : -45.09801811,
        'wavefunction_grid' : np.array([32,32,32]),
        'write_data_period' : 10,
        'write_qmcpack_restart' : False,
        }


    serial_references['Pt_bulk_spinorbit_input_band'] = {
        'a_length' : 7.42,
        'alpha' : 0.0,
        'atomic_coordinate_type' : 'Absolute',
        'atoms/atoms' : np.array(['Pt']),
        'atoms/format' : 'full_spin',
        'atoms/movable' : np.array([[ True,  True,  True]]),
        'atoms/positions' : np.array([[0., 0., 0.]]),
        'atoms/spin_phi' : np.array([0.]),
        'atoms/spin_ratio' : np.array([0.]),
        'atoms/spin_theta' : np.array([90.]),
        'b_length' : 7.42,
        'beta' : 0.0,
        'bravais_lattice_type' : 'Cubic Face Centered',
        'c_length' : 7.42,
        'calculation_mode' : 'Band Structure Only',
        'charge_density_mixing' : 0.5,
        'charge_mixing_type' : 'Broyden',
        'compressed_infile' : False,
        'compressed_outfile' : False,
        'description' : 'atom_Pt_pp',
        'energy_convergence_criterion' : 1e-09,
        'gamma' : 0.0,
        'kohn_sham_mucycles' : 3,
        'kohn_sham_solver' : 'davidson',
        'kpoint_distribution' : 1,
        'kpoint_is_shift' : np.array([1,1,1]),
        'kpoint_mesh' : np.array([-4 ,4 ,4]),
        'kpoints_bandstructure/counts' : np.array([1,10]),
        'kpoints_bandstructure/kpoints' : np.array([[0. , 0. , 0. ],
                                                    [0.5, 0.5, 0. ]]),
        'kpoints_bandstructure/labels' : np.array(['G','X']),
        'localize_localpp' : False,
        'localize_projectors' : False,
        'max_scf_steps' : 10,
        'noncollinear' : True,
        'occupation_electron_temperature_eV' : 0.2,
        'occupations_type' : 'Fermi Dirac',
        'potential_acceleration_constant_step' : 1.0,
        'potential_grid_refinement' : 2,
        'pseudo_dir' : './',
        'pseudopotential/pseudos' : np.array(['Pt.rel-pbe-n-rrkjus.UPF']),
        'pseudopotential/species' : np.array(['Pt']),
        'spinorbit' : True,
        'start_mode' : 'LCAO Start',
        'states_count_and_occupation' : '10 1.0 10 0.0',
        'subdiag_driver' : 'lapack',
        'wavefunction_grid' : np.array([32,32,32]),
        'write_data_period' : 10,
        'write_qmcpack_restart' : False,
        }


    serial_references['Si_8atoms_EXX_kpoints_input'] = {
        'a_length' : 10.2,
        'alpha' : 0.0,
        'atomic_coordinate_type' : 'Cell Relative',
        'atoms/atoms' : np.array(['Si','Si','Si','Si','Si','Si','Si','Si']),
        'atoms/format' : 'full_spin',
        'atoms/movable' : np.array([[ True , True,  True],
                                    [ True , True,  True],
                                    [ True , True,  True],
                                    [ True , True,  True],
                                    [ True , True,  True],
                                    [ True , True,  True],
                                    [ True , True,  True],
                                    [ True , True,  True]]),
        'atoms/positions' : np.array([[0.  , 0.  , 0.  ],
                                      [0.5 , 0.5 , 0.  ],
                                      [0.  , 0.5 , 0.5 ],
                                      [0.5 , 0.  , 0.5 ],
                                      [0.25, 0.25, 0.25],
                                      [0.75, 0.75, 0.25],
                                      [0.25, 0.75, 0.75],
                                      [0.75, 0.25, 0.75]]),
        'atoms/spin_phi' : np.array([0.,0.,0.,0.,0.,0.,0.,0.]),
        'atoms/spin_ratio' : np.array([0.,0.,0.,0.,0.,0.,0.,0.]),
        'atoms/spin_theta' : np.array([90.,90.,90.,90.,90.,90.,90.,90.]),
        'b_length' : 10.2,
        'beta' : 0.0,
        'bravais_lattice_type' : 'Cubic Primitive',
        'c_length' : 10.2,
        'calculation_mode' : 'Quench Electrons',
        'charge_density_mixing' : 0.5,
        'charge_mixing_type' : 'Broyden',
        'compressed_infile' : False,
        'compressed_outfile' : False,
        'description' : 'Si bulk',
        'dos_broading' : 0.5,
        'dos_method' : 'Gaussian',
        'energy_convergence_criterion' : 1e-09,
        'exchange_correlation_type' : 'gaupbe',
        'exx_mode' : 'Local fft',
        'exxdiv_treatment' : 'none',
        'gamma' : 0.0,
        'input_wave_function_file' : 'Waves/wave.out',
        'kohn_sham_mucycles' : 3,
        'kohn_sham_solver' : 'davidson',
        'kpoint_is_shift' : np.array([0,0,0]),
        'kpoint_mesh' : np.array([4,4,4]),
        'localize_localpp' : False,
        'localize_projectors' : False,
        'max_scf_steps' : 100,
        'occupations_type' : 'Fixed',
        'output_wave_function_file' : 'Waves/wave.out',
        'potential_acceleration_constant_step' : 1.0,
        'potential_grid_refinement' : 2,
        'pseudopotential/pseudos' : np.array(['Si_ONCV_PBE_sr.upf']),
        'pseudopotential/species' : np.array(['Si']),
        'start_mode' : 'LCAO Start',
        'states_count_and_occupation' : '16 2.0 8 0.0',
        'subdiag_driver' : 'lapack',
        'test_energy' : -33.43040801,
        'wavefunction_grid' : np.array([24,24,24]),
        'write_data_period' : 10,
        'write_pseudopotential_plots' : True,
        'write_qmcpack_restart' : False,
        'x_gamma_extrapolation' : False,
        }


    serial_references['CO_qmcpack_semilocal_xml_input'] = {
        'atomic_coordinate_type' : 'Absolute',
        'atoms/atoms' : np.array(['C','O']),
        'atoms/format' : 'movable',
        'atoms/movable' : np.array([True,True]),
        'atoms/positions' : np.array([[ 7.9316601,9.,9.],
                                      [10.0683399,9.,9.]]),
        'calculation_mode' : 'Quench Electrons',
        'charge_mixing_type' : 'Broyden',
        'davidson_multiplier' : 2,
        'description' : 'CO dimer length test',
        'energy_convergence_criterion' : 1e-14,
        'force_grad_order' : 0,
        'ionic_time_step' : 20.0,
        'kohn_sham_solver' : 'davidson',
        'lattice_vector' : np.array([18.,0.,0.,0.,18.,0.,0.,0.,18.]),
        'max_md_steps' : 100,
        'max_scf_steps' : 30,
        'occupations_type' : 'Fixed',
        'pseudopotential/pseudos' : np.array(['C.ccECP.xml','O.ccECP.xml']),
        'pseudopotential/species' : np.array(['C','O']),
        'qmc_nband' : 6,
        'relax_max_force' : 1e-05,
        'relax_method' : 'LBFGS',
        'renormalize_forces' : False,
        'rms_convergence_criterion' : 1e-10,
        'semilocal_projectors' : 10,
        'start_mode' : 'LCAO Start',
        'test_energy' : -21.70721649,
        'unoccupied_states_per_kpoint' : 4,
        'use_bessel_projectors' : True,
        'wavefunction_grid' : np.array([128,128,128]),
        'write_qmcpack_restart' : True,
        }


    h2o_tddft_common = {
        'a_length' : 24.0,
        'atomic_coordinate_type' : 'Absolute',
        'atoms/atoms' : np.array(['O','H','H']),
        'atoms/format' : 'movable',
        'atoms/movable' : np.array([True,True,True]),
        'atoms/positions' : np.array([[9., 9.,      9.65808 ],
                                      [9.,10.4362, 10.762926],
                                      [9., 7.5638, 10.762926]]),
        'b_length' : 24.0,
        'bravais_lattice_type' : 'Orthorhombic Primitive',
        'c_length' : 24.0,
        'calculation_mode' : 'TDDFT',
        'crds_units' : 'Bohr',
        'dipole_correction' : np.array([1,1,1]),
        'energy_convergence_criterion' : 1e-14,
        'input_tddft_file' : 'Wave/wave_tddft',
        'input_wave_function_file' : 'Wave/wave',
        'length_units' : 'Bohr',
        'max_scf_steps' : 100,
        'occupation_electron_temperature_eV' : 0.1,
        'occupation_number_mixing' : 0.3,
        'occupations_type' : 'Fermi Dirac',
        'output_tddft_file' : 'Wave/wave_tddft',
        'output_wave_function_file' : 'Wave/wave',
        'poisson_solver' : 'pfft',
        'pseudopotential/pseudos' : np.array(['../H.UPF','../O.UPF']),
        'pseudopotential/species' : np.array(['H','O']),
        'restart_tddft' : False,
        'rms_convergence_criterion' : 1e-10,
        'start_mode' : 'LCAO Start',
        'tddft_steps' : 10000,
        'unoccupied_states_per_kpoint' : 50,
        'verbose' : True,
        'wavefunction_grid' : np.array([64,64,64]),
        }

    serial_references['H2O_tddft_electric_field_input'] = dict(
        h2o_tddft_common,
        description = 'H2O TDDFT electric-field perturbation',
        electric_field_tddft = np.array([0.,0.,0.00272]),
        tddft_mode = 'electric field',
        )

    serial_references['H2O_tddft_point_charge_input'] = dict(
        h2o_tddft_common,
        description = 'H2O TDDFT point-charge perturbation',
        tddft_mode = 'point charge',
        tddft_qgau = 1.0,
        tddft_qpos = np.array([12.,12.,18.]),
        )


    serial_references['PbTiO3_berry_phase_nscf_input'] = {
        'BerryPhase' : True,
        'BerryPhaseDirection' : 2,
        'a_length' : 7.3699,
        'atomic_coordinate_type' : 'Cell Relative',
        'atoms/atoms' : np.array(['Pb','Ti','O','O','O']),
        'atoms/format' : 'basic',
        'atoms/positions' : np.array([[0. ,0. ,0.01],
                                      [0.5,0.5,0.5 ],
                                      [0. ,0.5,0.5 ],
                                      [0.5,0.5,0.  ],
                                      [0.5,0. ,0.5 ]]),
        'b_length' : 7.3699,
        'bravais_lattice_type' : 'Cubic Primitive',
        'c_length' : 7.3699,
        'calculation_mode' : 'NSCF',
        'charge_density_mixing' : 0.2,
        'charge_mixing_type' : 'Broyden',
        'description' : 'PbTiO3',
        'dos_method' : 'Gaussian',
        'energy_convergence_criterion' : 1e-07,
        'internal_pseudo_type' : 'sg15',
        'kohn_sham_mucycles' : 3,
        'kohn_sham_solver' : 'davidson',
        'kpoint_distribution' : 1,
        'kpoint_is_shift' : np.array([1,1,0]),
        'kpoint_mesh' : np.array([4,4,6]),
        'localize_localpp' : False,
        'localize_projectors' : False,
        'max_scf_steps' : 40,
        'output_wave_function_file' : '/dev/null',
        'potential_grid_refinement' : 2,
        'preconditioner_threshold' : 0.0001,
        'start_mode' : 'LCAO Start',
        'subdiag_driver' : 'scalapack',
        'wavefunction_grid' : np.array([16,16,16]),
        'write_data_period' : 50,
        }


    si_kvalues = np.arange(4,dtype=float)/4
    si_kpoints = np.array([
        (kx,ky,kz)
        for kx in si_kvalues
        for ky in si_kvalues
        for kz in si_kvalues
        ])
    serial_references['Si_wannier90_nscf_input'] = {
        'a_length' : 10.2,
        'atomic_coordinate_type' : 'Cell Relative',
        'atoms/atoms' : np.array(['Si','Si']),
        'atoms/format' : 'full_spin',
        'atoms/movable' : np.array([[True,True,True],[True,True,True]]),
        'atoms/positions' : np.array([[0.75,0.75,0.75],[0.,0.,0.]]),
        'atoms/spin_phi' : np.array([0.,0.]),
        'atoms/spin_ratio' : np.array([0.,0.]),
        'atoms/spin_theta' : np.array([0.,0.]),
        'b_length' : 10.2,
        'bravais_lattice_type' : 'Cubic Face Centered',
        'c_length' : 10.2,
        'calculation_mode' : 'NSCF',
        'crds_units' : 'Bohr',
        'description' : 'Silicon.txt',
        'exchange_correlation_type' : 'AUTO_XC',
        'internal_pseudo_type' : 'sg15',
        'kohn_sham_solver' : 'davidson',
        'kpoint_distribution' : 1,
        'kpoint_is_shift' : np.array([0,0,0]),
        'kpoint_mesh' : np.array([4,4,4]),
        'kpoints/kpoints' : si_kpoints,
        'kpoints/weights' : np.full(64,0.015625),
        'lattice_units' : 'Bohr',
        'num_wanniers' : 12,
        'potential_grid_refinement' : 2,
        'start_mode' : 'LCAO Start',
        'subdiag_driver' : 'auto',
        'time_reversal' : False,
        'unoccupied_states_per_kpoint' : 8,
        'use_symmetry' : 0,
        'wannier90' : True,
        'wannier90_scdm' : -1,
        'wavefunction_grid' : np.array([32,32,32]),
        'write_pseudopotential_plots' : False,
        }

#end def generate_serial_references


def get_serial_references():
    if len(serial_references)==0:
        generate_serial_references()
    #end if
    return serial_references
#end def get_serial_references


def check_vs_serial_reference(gi,name):
    sr = obj(get_serial_references()[name])
    sg = dict_serialize(gi,dict_type=obj)
    assert(check_object_eq(sg,sr))
#end def check_vs_serial_reference


def test_empty_init():
    ri = RmgInput()
#end test_empty_init


def test_input_spec():
    documented_sections = input_spec.section_order[:15]
    documented_count = sum(len(input_spec.section_contents[s]) for s in documented_sections)
    assert(documented_count==271)
    assert(len(input_spec.keywords)==284)

    added_keywords = set([
        'AFM',
        'BerryPhase',
        'BerryPhaseCycle',
        'BerryPhaseDirection',
        'STM_bias',
        'STM_height',
        'adaptive_cmix',
        'adaptive_convergence',
        'afd_cfac',
        'all_electron_parm',
        'davidson_1stage_ortho',
        'davidson_2stage_ortho',
        'davidson_premg',
        'drho_precond',
        'drho_precond_q0',
        'drho_precond_type',
        'electric_field',
        'electric_field_tddft',
        'epsg_guard',
        'freeze_ldaU_steps',
        'gpu_managed_memory',
        'internal_pseudo_type',
        'kpoint_units',
        'kpoints',
        'kpoints_bandstructure',
        'lambda_max',
        'lambda_min',
        'ldau_mixing',
        'ldau_mixing_type',
        'ldau_pulay_order',
        'ldau_pulay_refresh',
        'ldau_pulay_scale',
        'ldos_end_grid',
        'ldos_start_grid',
        'prolong_order',
        'qmc_nband',
        'resta_beta',
        'semilocal_projectors',
        'sts_end_grid',
        'sts_start_grid',
        'subdiag_groups',
        'tddft_frequency',
        'tddft_gpu',
        'tddft_noscf',
        'tddft_start_state',
        'test_bond_length',
        'test_bond_length_tolerance',
        'test_steps',
        'test_steps_tolerance',
        'tetra_method',
        'use_bessel_projectors',
        'use_block_diag',
        'use_cmix',
        'use_energy_correction',
        'use_gpu_fd',
        'use_rmm_diis',
        ])
    assert(added_keywords <= set(input_spec.keywords.keys()))

    assert(input_spec.keywords.AFM.key_type=='boolean')
    assert(input_spec.keywords.adaptive_cmix.max_value==10.0)
    assert(input_spec.keywords.electric_field.key_type=='double array')
    assert(input_spec.keywords.tddft_mode.allowed==set(
        ['electric field','point charge','vector potential']))
    assert('Example: 2 means system is missing two electrons' in
           input_spec.keywords.system_charge.description)
    assert(rmg_modes.full_mode('nscf')=='NSCF')
    assert(rmg_modes.full_mode('stm')=='STM')

#end def test_input_spec


def test_hubbard_u_records():
    text = '''
        Hubbard_U = "
        Ni 6.5 3d 0.1 0.2 0.3
        Mn 4.0 3d
        "
        '''
    ri = RmgInput()
    ri.read_text(text)
    assert(ri.Hubbard_U.Ni.U==6.5)
    assert(ri.Hubbard_U.Ni.orbital=='3d')
    assert(value_eq(ri.Hubbard_U.Ni.J,np.array([0.1,0.2,0.3])))
    assert(ri.Hubbard_U.Mn.U==4.0)
    assert(ri.Hubbard_U.Mn.orbital=='3d')

    ri_roundtrip = RmgInput()
    ri_roundtrip.read_text(ri.write_text())
    assert(check_object_eq(ri_roundtrip,ri))

    ri_generated = RmgInput()
    ri_generated.assign(Hubbard_U=obj(Ni=(6.5,'3d',0.1,0.2,0.3)))
    assert(ri_generated.Hubbard_U.Ni.U==6.5)
    assert(ri_generated.Hubbard_U.Ni.orbital=='3d')
    assert(value_eq(ri_generated.Hubbard_U.Ni.J,np.array([0.1,0.2,0.3])))

#end def test_hubbard_u_records


def test_kpoints_bandstructure_records():
    text = '''
        kpoints_bandstructure = "
        0.00 0.00 0.00  0 G
        0.50 0.00 0.00 20 M
        0.25 0.25 0.00 20 K
        0.00 0.00 0.00 20 G
        "
        '''
    ri = RmgInput()
    ri.read_text(text)
    assert(ri.kpoints_bandstructure.kpoints.shape==(4,3))
    assert(value_eq(ri.kpoints_bandstructure.counts,np.array([0,20,20,20])))
    assert(value_eq(ri.kpoints_bandstructure.labels,np.array(['G','M','K','G'])))

    ri_roundtrip = RmgInput()
    ri_roundtrip.read_text(ri.write_text())
    assert(check_object_eq(ri_roundtrip,ri))

#end def test_kpoints_bandstructure_records


def test_atoms_spin_ratio_records():
    text = '''
        atoms = "
        H 0.0 0.0 0.0 1 1 1  0.5
        H 1.0 0.0 0.0 1 1 1 -0.5
        "
        '''
    ri = RmgInput()
    ri.read_text(text)
    assert(ri.atoms.format=='spin_ratio')
    assert(value_eq(ri.atoms.spin_ratio,np.array([0.5,-0.5])))

    ri_roundtrip = RmgInput()
    ri_roundtrip.read_text(ri.write_text())
    assert(check_object_eq(ri_roundtrip,ri))

#end def test_atoms_spin_ratio_records


def test_read():
    infiles_read = {}
    for infile in TEST_FILES:
        ri_read = RmgInput(TEST_FILES[infile])
        assert(ri_read.is_valid())
        infiles_read[infile] = ri_read
    #end for

    input_files_check = (
        "BlackPhosphorus_input",
        "BlackPhosphorus_input_band",
        "BlackPhosphorus_Localize_both_pp_and_proj_multigrid_input",
        "C60_input",
        "CO_qmcpack_semilocal_xml_input",
        "Diamond16_input",
        "graphite_stress_input",
        "H2O_tddft_electric_field_input",
        "H2O_tddft_point_charge_input",
        "NiO8_input",
        "PbTiO3_berry_phase_nscf_input",
        "Pt_bulk_spinorbit_input",
        "Pt_bulk_spinorbit_input_band",
        "Si_8atoms_EXX_kpoints_input",
        "Si_wannier90_nscf_input",
        )

    # print out the reference text (used in generate_serial_references)
    #for infile in input_files_check:
    #    ri_read = infiles_read[infile]
    #    print()
    #    print(infile)
    #    print(80*'=')
    #    print(make_serial_reference(ri_read))
    ##end if

    for infile in input_files_check:
        ri_read = infiles_read[infile]
        check_vs_serial_reference(ri_read,infile)
    #end for

#end def test_read


def test_run_mode_input_coverage():
    mode_inputs = {
        'BlackPhosphorus_input'          : 'Quench Electrons',
        'BlackPhosphorus_input_band'     : 'Band Structure Only',
        'graphite_stress_input'           : 'Relax Structure',
        'Si_wannier90_nscf_input'         : 'NSCF',
        'H2O_tddft_electric_field_input'  : 'TDDFT',
        }
    for name,mode in mode_inputs.items():
        assert(RmgInput(TEST_FILES[name]).calculation_mode==mode)
    #end for

    wannier_nscf = RmgInput(TEST_FILES['Si_wannier90_nscf_input'])
    assert(wannier_nscf.wannier90)
    assert(len(wannier_nscf.kpoints.kpoints)==64)

    berry_nscf = RmgInput(TEST_FILES['PbTiO3_berry_phase_nscf_input'])
    assert(berry_nscf.calculation_mode=='NSCF')
    assert(berry_nscf.BerryPhase)

    electric_tddft = RmgInput(TEST_FILES['H2O_tddft_electric_field_input'])
    point_charge_tddft = RmgInput(TEST_FILES['H2O_tddft_point_charge_input'])
    assert(electric_tddft.tddft_mode=='electric field')
    assert(point_charge_tddft.tddft_mode=='point charge')

    qmcpack_orbitals = RmgInput(TEST_FILES['CO_qmcpack_semilocal_xml_input'])
    assert(qmcpack_orbitals.write_qmcpack_restart)
    assert(qmcpack_orbitals.qmc_nband==6)
    assert(qmcpack_orbitals.use_bessel_projectors)
    assert(qmcpack_orbitals.semilocal_projectors==10)
    assert(all(p.endswith('.xml') for p in qmcpack_orbitals.pseudopotential.pseudos))

#end def test_run_mode_input_coverage



def test_write(tmp_path):
    new_input_files = {
        'CO_qmcpack_semilocal_xml_input',
        'H2O_tddft_electric_field_input',
        'H2O_tddft_point_charge_input',
        'PbTiO3_berry_phase_nscf_input',
        'Si_wannier90_nscf_input',
        }
    input_files_write = tuple(TEST_FILES.keys())
    assert(new_input_files <= set(input_files_write))

    for infile in input_files_write:
        write_file = tmp_path / infile

        ri_read = RmgInput(TEST_FILES[infile])

        ri_read.write(write_file)

        ri_write = RmgInput(write_file)

        assert(check_object_eq(ri_write,ri_read))
    #end for

#end def test_write



@isolate_nexus_core
def test_generate():
    register_pseudo_files([
        'Ni_oncv.UPF','O_oncv.UPF','Pt.rel-pbe-n-rrkjus.UPF'
        ])

    # recreate 'BlackPhosphorus_input'
    infile = 'BlackPhosphorus_input'

    shared_inputs = obj(
        description               = 'Black Phosphorus',
        calculation_mode          = 'Quench Electrons',
        charge_density_mixing     = 0.5,
        charge_mixing_type        = 'Linear',
        charge_pulay_order        = 5,
        charge_pulay_scale        = 0.5,
        initial_diagonalization   = True,
        kohn_sham_mg_levels       = 2,
        kohn_sham_solver          = 'multigrid',
        localize_projectors       = False,
        max_scf_steps             = 20,
        occupation_electron_temperature_eV = 0.025,
        occupation_number_mixing  = 1.0,
        occupations_type          = 'Fermi Dirac',
        output_wave_function_file = 'Waves/wave.out',
        potential_acceleration_constant_step = 1.0,
        potential_grid_refinement = 2,
        rms_convergence_criterion = 1e-07,
        start_mode                = 'LCAO Start',
        subdiag_driver            = 'lapack',
        system_charge             = 0.0,
        unoccupied_states_per_kpoint = 10,
        bravais_lattice_type      = 'Orthorhombic Primitive',
        length_units              = 'Angstrom',
        a_length                  = 3.3136,
        b_length                  = 10.478,
        c_length                  = 4.3763,
        kpoint_distribution       = 1,
        kpoint_is_shift           = (0,0,0),
        kpoint_mesh               = (2,1,2),
        atomic_coordinate_type    = 'Absolute',
        crds_units                = 'Angstrom',
        )

    ri = generate_rmg_input(
        atoms                     = '''
         P   0.0            0.0            0.0            1
         P   1.6568         0.0            1.48364        1
         P   1.6568         2.13054        2.18815        1
         P   0.0            2.13054        3.67179        1
        ''',
        **shared_inputs
        )
    check_vs_serial_reference(ri,infile)

    ri = generate_rmg_input(
        atoms                     = obj(
            format    = 'movable',
            atoms     = ['P','P','P','P'],
            positions = [[0.    ,  0.     , 0.     ],
                         [1.6568,  0.     , 1.48364],
                         [1.6568,  2.13054, 2.18815],
                         [0.    ,  2.13054, 3.67179]],
            movable   = [True ,True ,True ,True],
            ),
        **shared_inputs
        )
    check_vs_serial_reference(ri,infile)


    # recreate 'NiO8_input'
    infile = 'NiO8_input'

    shared_inputs = obj(
        description               = 'NiO 8 atom cell in anti-ferromagnetic configuration solved using Davidson diagonalization',
        calculation_mode          = 'Quench Electrons',
        alpha                     = 0.0,
        beta                      = 0.0,
        charge_density_mixing     = 0.25,
        charge_mixing_type        = 'Broyden',
        davidson_max_steps        = 15,
        davidson_multiplier       = 4,
        energy_convergence_criterion = 1e-10,
        force_grad_order          = 0,
        gamma                     = 0.0,
        kohn_sham_mg_levels       = 3,
        kohn_sham_solver          = 'davidson',
        ldaU_mode                 = 'Simple',
        localize_localpp          = False,
        localize_projectors       = False,
        max_scf_steps             = 100,
        occupations_type          = 'MethfesselPaxton',
        potential_acceleration_constant_step = 1.0,
        potential_grid_refinement = 2,
        rms_convergence_criterion = 1e-08,
        start_mode                = 'LCAO Start',
        subdiag_driver            = 'lapack',
        test_energy               = -677.71977422,
        wavefunction_grid         = (36,36,36),
        write_data_period         = 5,
        bravais_lattice_type      = 'Cubic Primitive',
        a_length                  = 7.8811,
        b_length                  = 7.8811,
        c_length                  = 7.8811,
        kpoint_is_shift           = (0,0,0),
        kpoint_mesh               = (2,2,2),
        atomic_coordinate_type    = 'Cell Relative',
        )

    ri = generate_rmg_input(
        states_count_and_occupation_spin_down = '48 1.0 8 0.0',
        states_count_and_occupation_spin_up   = '48 1.0 8 0.0',
        Hubbard_U                 = 'Ni 6.5',
        pseudopotential           = '''
            Ni Ni_oncv.UPF
            O  O_oncv.UPF
            ''',
        atoms                     = '''
            O     0.000000     0.000000     0.500000    1    0.0
            O     0.000000     0.500000     0.000000    1    0.0
            O     0.500000     0.000000     0.000000    1    0.0
            O     0.500000     0.500000     0.500000    1    0.0
            Ni    0.005000     0.000000     0.000000    1    0.25
            Ni    0.500000     0.500000     0.000000    1    0.25
            Ni    0.000000     0.500000     0.500000    1   -0.25
            Ni    0.500000     0.000000     0.500000    1   -0.25
            ''',
        **shared_inputs
        )
    check_vs_serial_reference(ri,infile)

    ri = generate_rmg_input(
        states_count_and_occupation_spin_down = '48 1.0 8 0.0',
        states_count_and_occupation_spin_up = '48 1.0 8 0.0',
        Hubbard_U                 = obj(Ni=6.5),
        pseudopotential           = obj(
            species = ['Ni','O'],
            pseudos = ['Ni_oncv.UPF','O_oncv.UPF'],
            ),
        atoms                     = obj(
            format    = 'movable_moment',
            atoms     = ['O','O','O','O','Ni','Ni','Ni','Ni'],
            positions = [[0.   , 0. ,   0.5  ],
                         [0.   , 0.5,   0.   ],
                         [0.5  , 0. ,   0.   ],
                         [0.5  , 0.5,   0.5  ],
                         [0.005, 0. ,   0.   ],
                         [0.5  , 0.5,   0.   ],
                         [0.   , 0.5,   0.5  ],
                         [0.5  , 0. ,   0.5  ]],
            movable = [True,True,True,True,True,True,True,True],
            moments = [0.,0.,0.,0.,0.25,0.25,-0.25,-0.25],
            ),
        **shared_inputs
        )
    check_vs_serial_reference(ri,infile)

    if (
        find_spec("spglib") is not None
        and find_spec("seekpath") is not None
        and find_spec("scipy") is not None
        ):
        nio8 = generate_physical_system(
            units     = 'B',
            axes      = 7.8811*np.identity(3),
            elem      = ['O','O','O','O','Ni','Ni','Ni','Ni'],
            mag       = [0,0,0,0,.25,.25,-.25,-.25],
            posu      = [[0.   , 0. ,   0.5  ],
                         [0.   , 0.5,   0.   ],
                         [0.5  , 0. ,   0.   ],
                         [0.5  , 0.5,   0.5  ],
                         [0.005, 0. ,   0.   ],
                         [0.5  , 0.5,   0.   ],
                         [0.   , 0.5,   0.5  ],
                         [0.5  , 0. ,   0.5  ]],
            Ni        = 18,
            O         = 6,
            )
        nio8.structure.freeze(negate=True)
        s_trans,rmg_inputs,R,tmatrix,bv = nio8.structure.rmg_transform(all_results=True)
        nio8.structure = s_trans
        assert(value_eq(R,np.eye(3,dtype=float)))
        assert(tmatrix is None)
        keys = 'bravais_lattice_type','a_length','b_length','c_length','wavefunction_grid'
        for k in keys:
            del shared_inputs[k]
        d = dict(**rmg_inputs)
        d.update(**shared_inputs)
        ri = generate_rmg_input(
            Hubbard_U       = obj(Ni=6.5),
            virtual_frac    = 1./6,
            wf_grid_spacing = 0.22,
            pseudos         = ['Ni_oncv.UPF','O_oncv.UPF'],
            system          = nio8,
            **d
            )
        assert(value_eq(ri.length_units,'Bohr'))
        del ri.length_units
        check_vs_serial_reference(ri,infile)
    #end if


    # recreate 'Pt_bulk_spinorbit_input_band'
    infile = 'Pt_bulk_spinorbit_input_band'

    shared_inputs = obj(
        description               = 'atom_Pt_pp',
        calculation_mode          = 'Band Structure Only',
        alpha                     = 0.0,
        beta                      = 0.0,
        charge_density_mixing     = 0.5,
        charge_mixing_type        = 'Broyden',
        compressed_infile         = False,
        compressed_outfile        = False,
        energy_convergence_criterion = 1e-09,
        gamma                     = 0.0,
        kohn_sham_mucycles        = 3,
        kohn_sham_solver          = 'davidson',
        localize_localpp          = False,
        localize_projectors       = False,
        max_scf_steps             = 10,
        noncollinear              = True,
        occupation_electron_temperature_eV = 0.2,
        occupations_type          = 'Fermi Dirac',
        potential_acceleration_constant_step = 1.0,
        potential_grid_refinement = 2,
        pseudo_dir                = './',
        spinorbit                 = True,
        start_mode                = 'LCAO Start',
        subdiag_driver            = 'lapack',
        wavefunction_grid         = (32,32,32),
        write_data_period         = 10,
        write_qmcpack_restart     = False,
        bravais_lattice_type      = 'Cubic Face Centered',
        a_length                  = 7.42,
        b_length                  = 7.42,
        c_length                  = 7.42,
        kpoint_distribution       = 1,
        kpoint_is_shift           = (1,1,1),
        kpoint_mesh               = (-4, 4, 4),
        atomic_coordinate_type    = 'Absolute',
        )

    ri = generate_rmg_input(
        states_count_and_occupation = '10 1.0 10 0.0',
        kpoints_bandstructure     = '''
            0.0   0.0   0.0   1   G
            0.5   0.5   0.0   10  X
            ''',
        pseudopotential           = '''
            Pt    Pt.rel-pbe-n-rrkjus.UPF
            ''',
        atoms                     = '''
            Pt   0.0   0.0   0.0   1 1  1  0.0  90.0 00.0
            ''',
        **shared_inputs
        )
    check_vs_serial_reference(ri,infile)

    ri = generate_rmg_input(
        states_count_and_occupation = '10 1.0 10 0.0',
        kpoints_bandstructure = obj(
            counts  = [1,10],
            kpoints = [[0. , 0. , 0. ],
                       [0.5, 0.5, 0. ]],
            labels  = ['G','X'],
            ),
        pseudopotential       = obj(
            pseudos   = ['Pt.rel-pbe-n-rrkjus.UPF'],
            species   = ['Pt'],
            ),
        atoms                 = obj(
            format     = 'full_spin',
            atoms      = ['Pt'],
            positions  = [[0., 0., 0.]],
            movable    = [[ True,  True,  True]],
            spin_ratio = [0.],
            spin_phi   = [0.],
            spin_theta = [90.],
            ),
        **shared_inputs
        )
    check_vs_serial_reference(ri,infile)



    # test diamond generation

    a = 3.57 # A
    a = convert(a,'A','B')

    shared_inputs = obj(
        # nexus inputs
        input_type = 'generic',
        # control options
        calculation_mode       = 'Quench Electrons',
        compressed_infile      = False,
        compressed_outfile     = False,
        description            = 'diamond',
        energy_convergence_criterion = 1.0e-09,
        max_scf_steps          = 100,
        #start_mode             = 'Restart From File',
        write_data_period      = 10,
        # cell parameter options
        atomic_coordinate_type = 'Cell Relative',
        potential_grid_refinement = 2,
        # pseudopotential related options
        localize_localpp       = False,
        localize_projectors    = False,
        # kohn sham solver options
        kohn_sham_mucycles     = 3,
        kohn_sham_solver       = 'davidson',
        # orbital occupation options
        occupations_type       = 'Fixed',
        # charge density mixing options
        charge_density_mixing  = 0.5,
        charge_mixing_type     = 'Broyden',
        potential_acceleration_constant_step = 1.0,
        # diagonalization options
        subdiag_driver         = 'lapack',
        ## testing options
        #test_energy            = -11.32982439,
        # miscellaneous options
        kpoint_distribution    = 8,
        )


    # manual specification of rmg system (no physical system object)
    mi = generate_rmg_input(
        bravais_lattice_type   = 'Cubic Face Centered',
        a_length               = a,
        b_length               = a,
        c_length               = a,
        alpha                  = 0.0,
        beta                   = 0.0,
        gamma                  = 0.0,
        atoms                  = '''
          C      0.250   0.250   0.250 1 1 1  0.0000   0.00   0.00
          C      0.000   0.000   0.000 1 1 1  0.0000   0.00   0.00
          ''',
        kpoint_mesh            = (4,4,4),
        kpoint_is_shift        = (0,0,0),
        states_count_and_occupation = '4 2.0 4 0.0',
        wavefunction_grid           = (32,32,32),
        **shared_inputs
        )

    # generated from system
    system = generate_physical_system(
        structure = 'diamond',
        cell      = 'prim',
        C         = 4,
        units     = 'A',
        kgrid     = (4,4,4),
        kshift    = (0,0,0),
        )

    gi = generate_rmg_input(
        system          = system,
        virtual_frac    = 1.0,
        wf_grid_spacing = 0.15,
        **shared_inputs
        )

    # generated from system w/ kgrid override
    ki = generate_rmg_input(
        system          = system,
        virtual_frac    = 1.0,
        wf_grid_spacing = 0.15,
        kpoint_mesh     = (4,4,4),
        kpoint_is_shift = (0,0,0),
        **shared_inputs
        )


    # check shared keys match
    mkeys = set(mi.keys())
    gkeys = set(gi.keys())
    kkeys = set(ki.keys())

    skeys = mkeys & gkeys & kkeys

    skeys_ref = set([
        'atomic_coordinate_type', 'atoms', 'calculation_mode', 
        'charge_density_mixing', 'charge_mixing_type', 'compressed_infile', 
        'compressed_outfile', 'description', 'energy_convergence_criterion', 
        'kohn_sham_mucycles', 'kohn_sham_solver', 'kpoint_distribution', 
        'localize_localpp', 'localize_projectors', 'max_scf_steps', 
        'occupations_type', 'potential_acceleration_constant_step', 
        'potential_grid_refinement', 'states_count_and_occupation', 
        'subdiag_driver', 'wavefunction_grid', 'write_data_period'])

    assert(skeys==skeys_ref)

    mr = obj()
    gr = obj()
    kr = obj()
    for ri,rk,rr in [(mi,mkeys,mr),(gi,gkeys,gr),(ki,kkeys,kr)]:
        for k in rk-skeys:
            rr[k] = ri[k]
            del ri[k]
        #end for
    #end for

    assert(check_object_eq(gi,ki))
    del mi.atoms
    del gi.atoms
    assert(check_object_eq(mi,gi))


    # check that residual (differing) keys match reference
    mr_ref = obj(
        a_length        = 6.746322294401746,
        alpha           = 0.0,
        b_length        = 6.746322294401746,
        beta            = 0.0,
        bravais_lattice_type = 'Cubic Face Centered',
        c_length        = 6.746322294401746,
        gamma           = 0.0,
        kpoint_is_shift = np.array([0, 0, 0],dtype=int),
        kpoint_mesh     = np.array([4, 4, 4],dtype=int),
        )
    assert(check_object_eq(mr,mr_ref))

    kr_ref = obj(
        kpoint_is_shift = np.array([0, 0, 0],dtype=int),
        kpoint_mesh     = np.array([4, 4, 4],dtype=int),
        lattice_units   = 'Bohr',
        lattice_vector  = np.array([[3.37316115, 3.37316115, 0.        ],
                                    [0.        , 3.37316115, 3.37316115],
                                    [3.37316115, 0.        , 3.37316115]],
                                   dtype=float),
        )
    assert(check_object_eq(kr,kr_ref))

    gr_ref = obj(
        lattice_units   = 'Bohr',
        lattice_vector  = np.array([[3.37316115, 3.37316115, 0.        ],
                                    [0.        , 3.37316115, 3.37316115],
                                    [3.37316115, 0.        , 3.37316115]],
                                   dtype=float),
        kpoints = obj(
            kpoints         = np.array(
                [[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                 [ 2.50000000e-01,  3.09994390e-17,  3.24386335e-18],
                 [ 5.00000000e-01, -2.33137274e-17, -2.33137274e-17],
                 [ 7.50000000e-01,  5.63983306e-18,  5.63983306e-18],
                 [ 4.14583178e-17,  2.50000000e-01,  1.16568637e-17],
                 [ 2.50000000e-01,  2.50000000e-01,  0.00000000e+00],
                 [ 5.00000000e-01,  2.50000000e-01,  3.24386335e-18],
                 [ 7.50000000e-01,  2.50000000e-01,  3.62891808e-17],
                 [ 2.33137274e-17,  5.00000000e-01,  2.33137274e-17],
                 [ 2.50000000e-01,  5.00000000e-01,  1.16568637e-17],
                 [ 5.00000000e-01,  5.00000000e-01,  0.00000000e+00],
                 [ 7.50000000e-01,  5.00000000e-01,  3.24386335e-18],
                 [-5.63983306e-18,  7.50000000e-01, -5.63983306e-18],
                 [ 2.50000000e-01,  7.50000000e-01,  2.12678489e-17],
                 [ 5.00000000e-01,  7.50000000e-01, -3.24386335e-18],
                 [ 7.50000000e-01,  7.50000000e-01,  0.00000000e+00],
                 [ 3.24386335e-18,  1.16568637e-17,  2.50000000e-01],
                 [ 2.50000000e-01,  3.15405006e-17,  2.50000000e-01],
                 [ 5.00000000e-01,  3.24386335e-18,  2.50000000e-01],
                 [ 7.50000000e-01,  3.62891808e-17,  2.50000000e-01],
                 [ 3.15405006e-17,  2.50000000e-01,  2.50000000e-01],
                 [ 2.50000000e-01,  2.50000000e-01,  2.50000000e-01],
                 [ 5.00000000e-01,  2.50000000e-01,  2.50000000e-01],
                 [ 7.50000000e-01,  2.50000000e-01,  2.50000000e-01],
                 [ 1.16568637e-17,  5.00000000e-01,  2.50000000e-01],
                 [ 2.50000000e-01,  5.00000000e-01,  2.50000000e-01],
                 [ 5.00000000e-01,  5.00000000e-01,  2.50000000e-01],
                 [ 7.50000000e-01,  5.00000000e-01,  2.50000000e-01],
                 [ 2.12678489e-17,  7.50000000e-01,  2.50000000e-01],
                 [ 2.50000000e-01,  7.50000000e-01,  2.50000000e-01],
                 [ 5.00000000e-01,  7.50000000e-01,  2.50000000e-01],
                 [ 7.50000000e-01,  7.50000000e-01,  2.50000000e-01],
                 [-2.33137274e-17,  2.33137274e-17,  5.00000000e-01],
                 [ 2.50000000e-01,  1.16568637e-17,  5.00000000e-01],
                 [ 5.00000000e-01,  0.00000000e+00,  5.00000000e-01],
                 [ 7.50000000e-01,  3.24386335e-18,  5.00000000e-01],
                 [ 3.24386335e-18,  2.50000000e-01,  5.00000000e-01],
                 [ 2.50000000e-01,  2.50000000e-01,  5.00000000e-01],
                 [ 5.00000000e-01,  2.50000000e-01,  5.00000000e-01],
                 [ 7.50000000e-01,  2.50000000e-01,  5.00000000e-01],
                 [ 0.00000000e+00,  5.00000000e-01,  5.00000000e-01],
                 [ 2.50000000e-01,  5.00000000e-01,  5.00000000e-01],
                 [ 5.00000000e-01,  5.00000000e-01,  5.00000000e-01],
                 [ 7.50000000e-01,  5.00000000e-01,  5.00000000e-01],
                 [-3.24386335e-18,  7.50000000e-01,  5.00000000e-01],
                 [ 2.50000000e-01,  7.50000000e-01,  5.00000000e-01],
                 [ 5.00000000e-01,  7.50000000e-01,  5.00000000e-01],
                 [ 7.50000000e-01,  7.50000000e-01,  5.00000000e-01],
                 [ 5.63983306e-18, -5.63983306e-18,  7.50000000e-01],
                 [ 2.50000000e-01,  2.12678489e-17,  7.50000000e-01],
                 [ 5.00000000e-01, -3.24386335e-18,  7.50000000e-01],
                 [ 7.50000000e-01, -3.15405006e-17,  7.50000000e-01],
                 [ 3.62891808e-17,  2.50000000e-01,  7.50000000e-01],
                 [ 2.50000000e-01,  2.50000000e-01,  7.50000000e-01],
                 [ 5.00000000e-01,  2.50000000e-01,  7.50000000e-01],
                 [ 7.50000000e-01,  2.50000000e-01,  7.50000000e-01],
                 [ 3.24386335e-18,  5.00000000e-01,  7.50000000e-01],
                 [ 2.50000000e-01,  5.00000000e-01,  7.50000000e-01],
                 [ 5.00000000e-01,  5.00000000e-01,  7.50000000e-01],
                 [ 7.50000000e-01,  5.00000000e-01,  7.50000000e-01],
                 [-3.15405006e-17,  7.50000000e-01,  7.50000000e-01],
                 [ 2.50000000e-01,  7.50000000e-01,  7.50000000e-01],
                 [ 5.00000000e-01,  7.50000000e-01,  7.50000000e-01],
                 [ 7.50000000e-01,  7.50000000e-01,  7.50000000e-01]],
                dtype=float),
            weights         = np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
                                        1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
                                        1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
                                        1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
                                        1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],
                                       dtype=float),
            ),
        )

    assert(check_object_eq(gr,gr_ref,atol=1e-12))


    # Regenerate the run-mode inputs from their parsed values. This follows the
    # same keyword assignment path as generating them explicitly from scratch,
    # while keeping the larger k-point and atomic tables readable in this test.
    run_mode_inputs = (
        'CO_qmcpack_semilocal_xml_input',
        'H2O_tddft_electric_field_input',
        'H2O_tddft_point_charge_input',
        'PbTiO3_berry_phase_nscf_input',
        'Si_wannier90_nscf_input',
        )
    for infile in run_mode_inputs:
        ri_read = RmgInput(TEST_FILES[infile])
        ri_generated = generate_rmg_input(
            defaults = 'none',
            **dict(ri_read.items())
            )
        check_vs_serial_reference(ri_generated,infile)
    #end for

#end def test_generate
