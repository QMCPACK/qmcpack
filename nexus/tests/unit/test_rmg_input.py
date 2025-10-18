
import testing
from testing import failed
from testing import divert_nexus_log,restore_nexus_log
from testing import value_eq,object_eq,check_object_eq
import versions


associated_files = dict()

input_files = '''
    AlN32_input
    atomO_polarized_input
    BlackPhosphorus_Delocalize_both_pp_and_proj_davidson_input
    BlackPhosphorus_Delocalize_both_pp_and_proj_multigrid_input
    BlackPhosphorus_input
    BlackPhosphorus_input_1nongammapoint
    BlackPhosphorus_input_band
    BlackPhosphorus_Localize_both_pp_and_proj_davidson_input
    BlackPhosphorus_Localize_both_pp_and_proj_multigrid_input
    BlackPhosphorus_Localize_pp_only_davidson_input
    BlackPhosphorus_Localize_pp_only_multigrid_input
    BlackPhosphorus_Localize_proj_only_davidson_input
    BlackPhosphorus_Localize_proj_only_multigrid_input
    C60_input
    Diamond16_input
    Diamond2_input
    Fe_2atom_input
    graphite_stress_input
    Mg_2atom_input
    nanotube_80_input
    nanotube_80_input_band
    nanotube_80_input_band1
    NiO512_input
    NiO8_input
    Pt_bulk_spinorbit_input
    Pt_bulk_spinorbit_input_band
    Si_8atoms_EXX_kpoints_input
    U_bulk_spinorbit_RMG_input
    '''.split()
    



def get_files():
    return testing.collect_unit_test_file_paths('rmg_input',associated_files)
#end def get_files



def make_serial_reference(ri):
    import numpy as np
    s = ri.serial()
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
    import numpy as np

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

#end def generate_serial_references


def get_serial_references():
    if len(serial_references)==0:
        generate_serial_references()
    #end if
    return serial_references
#end def get_serial_references


def check_vs_serial_reference(gi,name):
    from generic import obj
    sr = obj(get_serial_references()[name])
    sg = gi.serial()
    assert(check_object_eq(sg,sr))
#end def check_vs_serial_reference


def test_files():
    filenames = input_files
    files = get_files()
    assert(set(files.keys())==set(filenames))
#end def test_files



def test_import():
    from rmg_input import RmgInput
#end def test_import



def test_empty_init():
    from rmg_input import RmgInput
    ri = RmgInput()
#end test_empty_init



def test_read():
    from rmg_input import RmgInput
    
    files = get_files()

    infiles_read = {}
    for infile in input_files:
        ri_read = RmgInput(files[infile])
        assert(ri_read.is_valid())
        infiles_read[infile] = ri_read
    #end for

    input_files_check = '''
        BlackPhosphorus_input
        BlackPhosphorus_input_band
        BlackPhosphorus_Localize_both_pp_and_proj_multigrid_input
        C60_input
        Diamond16_input
        graphite_stress_input
        NiO8_input
        Pt_bulk_spinorbit_input
        Pt_bulk_spinorbit_input_band
        Si_8atoms_EXX_kpoints_input
        '''.split()

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



def test_write():
    import os
    from rmg_input import RmgInput

    tpath = testing.setup_unit_test_output_directory('rmg_input','test_write')
    
    files = get_files()

    for infile in input_files:
        write_file = os.path.join(tpath,infile)

        ri_read = RmgInput(files[infile])

        ri_read.write(write_file)

        ri_write = RmgInput(write_file)

        assert(check_object_eq(ri_write,ri_read))
    #end for

#end def test_write



def test_generate():
    import numpy as np
    from generic import obj
    from unit_converter import convert
    from physical_system import generate_physical_system
    from rmg_input import generate_rmg_input

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

    if versions.spglib_available and versions.seekpath_available:
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
        shared_inputs.delete('bravais_lattice_type','a_length','b_length','c_length','wavefunction_grid')
        ri = generate_rmg_input(
            Hubbard_U       = obj(Ni=6.5),
            virtual_frac    = 1./6,
            wf_grid_spacing = 0.22,
            pseudos         = ['Ni_oncv.UPF','O_oncv.UPF'],
            system          = nio8,
            **obj(rmg_inputs,shared_inputs)
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

#end def test_generate
