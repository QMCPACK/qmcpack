#! /usr/bin/env python3
"""RMG + QMCPACK example: SCF, NSCF on the supercell twists, J1/J2 optimization and twist-averaged DMC."""

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_rmg
from nexus import generate_qmcpack

settings(
    pseudo_dir = '../../qmcpack/pseudopotentials',
    results    = '',
    sleep      = 3,
    machine    = 'ws8',
    )

# 8-atom cubic diamond cell with a 2x2x2 twist grid
system = generate_physical_system(
    units    = 'A',
    axes     = '''1.785   1.785   0.000
                  0.000   1.785   1.785
                  1.785   0.000   1.785''',
    elem_pos = '''
               C  0.0000  0.0000  0.0000
               C  0.8925  0.8925  0.8925
               ''',
    tiling   = [[ 1, -1,  1],
                [ 1,  1, -1],
                [-1,  1,  1]],
    kgrid    = (2,2,2),
    kshift   = (0,0,0),
    C        = 4,
    )

# RMG inputs shared by the SCF and NSCF runs.  RMG runs in the primitive
# cell.  The NSCF run restarts from the per-processor restart files of the
# SCF run, so both runs must use the same wavefunction grid, the same
# number of MPI ranks and the same kpoint_distribution.
rmg_job = job(cores=4,app='rmg-cpu')
shared_inputs = obj(
    job                       = rmg_job,
    system                    = system,
    pseudos                   = ['C.BFD.xml'],
    input_type                = 'generic',
    wf_grid_spacing           = 0.15,
    atomic_coordinate_type    = 'Cell Relative',
    exchange_correlation_type = 'pbe',
    localize_localpp          = False,
    localize_projectors       = False,
    kohn_sham_solver          = 'davidson',
    subdiag_driver            = 'lapack',
    occupations_type          = 'Fixed',
    kpoint_distribution       = 4,
    )

scf = generate_rmg(
    identifier                   = 'scf',
    path                         = 'diamond/scf',
    calculation_mode             = 'Quench Electrons',
    kpoint_mesh                  = (4,4,4),
    kpoint_is_shift              = (0,0,0),
    energy_convergence_criterion = 1e-9,
    max_scf_steps                = 100,
    charge_density_mixing        = 0.5,
    charge_mixing_type           = 'Broyden',
    **shared_inputs
    )

# The NSCF run computes the orbitals at the k-points of the twists (no
# kpoint_mesh, Nexus writes the k-point list) and writes the ESHDF orbital
# file Waves/wave.out.h5 for QMCPACK.
nscf = generate_rmg(
    identifier            = 'nscf',
    path                  = 'diamond/nscf_twist',
    calculation_mode      = 'NSCF',
    write_qmcpack_restart = True,
    dependencies          = (scf,'charge_density'),
    **shared_inputs
    )

opt = generate_qmcpack(
    identifier   = 'opt',
    path         = 'diamond/optJ12',
    job          = job(cores=8,threads=8,app='qmcpack'),
    input_type   = 'basic',
    system       = system,
    pseudos      = ['C.BFD.xml'],
    twistnum     = 0,
    J1           = True,
    J2           = True,
    qmc          = 'opt',
    cycles       = 6,
    samples      = 25600,
    dependencies = (nscf,'orbitals'),
    )

# twist-averaged DMC: one MPI rank per twist (8 twists)
qmc = generate_qmcpack(
    identifier       = 'dmc',
    path             = 'diamond/dmc_twist',
    job              = job(cores=8,threads=1,app='qmcpack'),
    input_type       = 'basic',
    system           = system,
    pseudos          = ['C.BFD.xml'],
    J1               = True,
    J2               = True,
    qmc              = 'dmc',
    walkers_per_rank = 32,
    eq_dmc           = True,
    timestep         = 0.02,
    dependencies     = [(nscf,'orbitals'),
                        (opt,'jastrow')],
    )

run_project()
