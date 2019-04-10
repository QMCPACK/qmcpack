#! /usr/bin/env python

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_quantum_package
from nexus import generate_convert4qmc
from nexus import generate_qmcpack

settings(
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws12',
    qprc          = '/home/j1k/apps/quantum_package/qp2-2.0.0-beta/quantum_package.rc',
    )

# define run details
qp_job  = job(cores=12,threads=12)
c4q_job = job(cores=1)
qmc_job = job(cores=12,threads=12)

# read in structure for oxygen dimer
dimer = generate_physical_system(
    structure = './O2.xyz',
    net_spin  = 2,
    )

# path, job, system details are shared across runs
qp_shared = dict(
    path   = 'O_dimer/selci',
    job    = qp_job,
    system = dimer,
    )

# run Hartree-Fock
scf = generate_quantum_package(
    identifier            = 'scf',
    run_type              = 'scf',
    ao_basis              = 'cc-pvtz',
    io_mo_two_e_integrals = 'Write',
    **qp_shared
    )

fit = generate_quantum_package(
    identifier    = 'fit',
    run_type      = 'four_idx_transform',
    prefix        = 'fci', 
    dependencies  = (scf,'other'),
    **qp_shared
    )

# singles calcs in case HF converged to wrong state
cis = generate_quantum_package(
    identifier            = 'cis',
    run_type              = 'cis',
    cis_loop              = True,
    frozen_core           = True,
    io_mo_two_e_integrals = 'Write',
    prefix                = 'fci',
    dependencies          = (fit,'other'),
    **qp_shared
    )

# initial selected CI run
fci0 = generate_quantum_package(
    identifier    = 'fci0',
    run_type      = 'fci',
    converge_dets = True,
    dependencies  = (cis,'other'),
    **qp_shared
    )

# save natural orbitals from selCI run
sno = generate_quantum_package(
    identifier    = 'sno',
    run_type      = 'save_natorb',
    dependencies  = (fci0,'other'),
    **qp_shared
    )

# final selected CI based on natural orbitals
fci = generate_quantum_package(
    identifier    = 'fci',
    run_type      = 'fci',
    converge_dets = True,
    dependencies  = (sno,'other'),
    **qp_shared
    )

# write wavefunction output for QMCPACK
swf = generate_quantum_package(
    identifier    = 'swf',
    run_type      = 'save_for_qmcpack',
    dependencies  = (fci,'other'),
    **qp_shared
    )

# convert orbitals and multidet to QMCPACK format
c4q = generate_convert4qmc(
    identifier   = 'c4q',
    path         = 'O_dimer/selci',
    job          = c4q_job,
    hdf5         = True,
    dependencies = (swf,'orbitals'),
    )

# run VMC with QMCPACK
qmc = generate_qmcpack(
    identifier   = 'vmc',
    path         = 'O_dimer/vmc',
    job          = qmc_job,
    system       = dimer,
    jastrows     = [],
    qmc          = 'vmc',
    dependencies = (c4q,'orbitals'),
    )

run_project()
