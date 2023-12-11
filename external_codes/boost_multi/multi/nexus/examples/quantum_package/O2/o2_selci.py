#! /usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_quantum_package

settings(
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16',
    qprc          = '/your/path/to/quantum_package.rc',
    )

# define run details
qp_job = job(cores=16,threads=16)

# read in structure for oxygen dimer
dimer = generate_physical_system(
    structure = './O2.xyz',
    net_spin  = 2,
    )

# path, job, system details are shared across runs
qp_shared = dict(
    path   = 'O_dimer_selected_CI',
    job    = qp_job,
    system = dimer,
    prefix = 'fci', # single shared ezfio
    )

# run Hartree-Fock
scf = generate_quantum_package(
    identifier            = 'scf',
    run_type              = 'scf',
    ao_basis              = 'aug-cc-pvdz',
    io_ao_two_e_integrals = 'Write',
    four_idx_transform    = True,
    **qp_shared
    )

# initial selected CI run
fci0 = generate_quantum_package(
    identifier         = 'fci0',
    run_type           = 'fci',
    n_det_max          = 5000,
    save_natorb        = True,
    four_idx_transform = True,
    dependencies       = (scf,'other'),
    **qp_shared
    )

# final selected CI based on natural orbitals
fci = generate_quantum_package(
    identifier    = 'fci',
    run_type      = 'fci',
    n_det_max     = 5000,
    dependencies  = (fci0,'other'),
    **qp_shared
    )

run_project()
