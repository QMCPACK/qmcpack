#! /usr/bin/env python3

import numpy as np

from unit_converter import convert

from nexus import settings,job,obj,run_project
from nexus import generate_physical_system
from nexus import generate_rmg


settings(
    pseudo_dir = '../../qmcpack/pseudopotentials',
    results    = '',
    sleep      = 3,
    machine    = 'ws8',
    )


# lattice constant for diamond
a = 3.57 # A
a = convert(a,'A','B')



# RMG inputs shared for both manual and generated cell specification
shared_inputs = obj(
    # nexus inputs
    identifier             = 'scf',
    job                    = job(cores=8,app='rmg-cpu'),
    pseudos                = ['C.BFD.xml'],
    input_type             = 'generic',
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
    kpoint_is_shift        = (0,0,0),
    kpoint_mesh            = (4,4,4),
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
    #test_energy            = -11.39547539 # reference energy for validation
    # miscellaneous options
    kpoint_distribution    = 8,
    )


# manual specification of rmg system (no physical system object)
scf_man = generate_rmg(
    path       = 'diamond2/scf_man',
    # rmg specific inputs
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
    states_count_and_occupation = '4 2.0 4 0.0',
    wavefunction_grid           = (32,32,32),
    **shared_inputs
    )



# generated specification of rmg system
system = generate_physical_system(
    structure = 'diamond',
    cell      = 'prim',
    C         = 4,
    units     = 'A',
    )

scf_gen = generate_rmg(
    path            = 'diamond2/scf_gen',
    system          = system,
    virtual_frac    = 1.0,
    wf_grid_spacing = 0.15,
    **shared_inputs
    )



run_project()
