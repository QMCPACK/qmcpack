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
    qprc          = \
'/home/ubuntu/apps/qp2/quantum_package.rc',
    )

scf_job = job(cores=16,threads=16)
c4q_job = job(cores=1)
qmc_job = job(cores=16,threads=16)

system = generate_physical_system(
    structure = 'H2O.xyz',
    )

# perform Hartree-Fock
scf = generate_quantum_package(
    identifier       = 'hf',      # log output goes to hf.out
    path             = 'H2O/hf',  # directory to run in
    job              = scf_job,
    system           = system,
    prefix           = 'h2o',     # create/use h2o.ezfio
    run_type         = 'scf',     # qprun scf h2o.ezfio
    ao_basis         = 'cc-pvtz', # use cc-pvtz basis
    save_for_qmcpack = True,      # write h5 file for qmcpack
    )

# convert orbitals to QMCPACK format
c4q = generate_convert4qmc(
    identifier   = 'c4q',
    path         = 'H2O/hf',
    job          = c4q_job,
    hdf5         = True,          # use hdf5 format
    dependencies = (scf,'orbitals'),
    )

# calculate cusp correction
cc = generate_cusp_correction(
    identifier   = 'cusp',
    path         = 'H2O/cuspcorr',
    job          = qmc_job,
    system       = system,
    dependencies = (c4q,'orbitals'),
    )

# optimize 2-body Jastrow
optJ2 = generate_qmcpack(
    block           = True,
    identifier      = 'opt',
    path            = 'H2O/optJ2',
    job             = qmc_job,
    system          = system,
    J2              = True,
    J2_rcut         = 8.0,
    qmc             = 'opt',          # use opt defaults
    minmethod       = 'oneshiftonly', # adjust for oneshift
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
    path            = 'H2O/optJ3',
    job             = qmc_job,
    system          = system,
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
    path         = 'H2O/vmc',
    job          = qmc_job,
    system       = system,
    jastrows     = [],
    qmc          = 'vmc',    # use vmc defaults
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
    path         = 'H2O/dmc',
    job          = qmc_job,
    system       = system,
    jastrows     = [],
    qmc          = 'dmc',    # use dmc defaults
    eq_dmc       = True,     # add equilibration run
    dependencies = [(c4q,'orbitals'),
                    (cc,'cuspcorr'),
                    (optJ3,'jastrow')],
    )

run_project()
