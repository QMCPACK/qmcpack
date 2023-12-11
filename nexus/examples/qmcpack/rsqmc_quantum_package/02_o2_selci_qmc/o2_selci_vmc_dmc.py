#! /usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_quantum_package
from nexus import generate_convert4qmc
from nexus import generate_cusp_correction
from nexus import generate_qmcpack

# note: you must source the QP config file before running this script
#   source /home/ubuntu/apps/qp2/quantum_package.rc

settings(
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16',
    qprc          = '/home/ubuntu/apps/qp2/quantum_package.rc',
    )

# define run details
qp_job  = job(cores=16,threads=16)
c4q_job = job(cores=1)
qmc_job = job(cores=16,threads=16)

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
    identifier       = 'fci',
    run_type         = 'fci',
    n_det_max        = 5000,
    save_for_qmcpack = True,
    dependencies     = (fci0,'other'),
    **qp_shared
    )

# convert orbitals and multidet to QMCPACK format
c4q = generate_convert4qmc(
    identifier   = 'c4q',
    path         = 'O_dimer/selci',
    job          = c4q_job,
    hdf5         = True,
    dependencies = (fci,'orbitals'),
    )

# calculate cusp correction
cc = generate_cusp_correction(
    identifier   = 'cusp',
    path         = 'O_dimer/cuspcorr',
    job          = qmc_job,
    system       = dimer,
    dependencies = (c4q,'orbitals'),
    )

# optimize 2-body Jastrow
optJ2 = generate_qmcpack(
    block           = True,
    identifier      = 'opt',
    path            = 'O_dimer/optJ2',
    job             = qmc_job,
    system          = dimer,
    J2              = True,
    qmc             = 'opt',
    minmethod       = 'oneshiftonly',
    init_cycles     = 3,
    init_minwalkers = 0.1,
    cycles          = 3,
    samples         = 25600,
    dependencies    = [(c4q,'orbitals'),
                       (cc,'cuspcorr')],
    )

# optimize 3-body Jastrow
optJ3 = generate_qmcpack(
    block           = True,
    identifier      = 'opt',
    path            = 'O_dimer/optJ3',
    job             = qmc_job,
    system          = dimer,
    J3              = True,
    qmc             = 'opt',
    minmethod       = 'oneshiftonly',
    init_cycles     = 3,
    init_minwalkers = 0.1,
    cycles          = 3,
    samples         = 51200,
    dependencies    = [(c4q,'orbitals'),
                       (cc,'cuspcorr'),
                       (optJ2,'jastrow')],
    )

# run VMC with QMCPACK
qmc = generate_qmcpack(
    block        = True,
    identifier   = 'vmc',
    path         = 'O_dimer/vmc',
    job          = qmc_job,
    system       = dimer,
    jastrows     = [],
    qmc          = 'vmc',
    blocks       = 800,
    steps        = 100,
    dependencies = [(c4q,'orbitals'),
                    (cc,'cuspcorr'),
                    (optJ3,'jastrow')],
    )

# run DMC with QMCPACK
qmc = generate_qmcpack(
    block        = True,
    identifier   = 'dmc',
    path         = 'O_dimer/dmc',
    job          = qmc_job,
    system       = dimer,
    jastrows     = [],
    qmc          = 'dmc',
    vmc_samples  = 1024, # dmc walkers drawn from vmc
    eq_dmc       = True, # run intermediate dmc w/ larger timestep
    dependencies = [(c4q,'orbitals'),
                    (cc,'cuspcorr'),
                    (optJ3,'jastrow')],
    )

run_project()
