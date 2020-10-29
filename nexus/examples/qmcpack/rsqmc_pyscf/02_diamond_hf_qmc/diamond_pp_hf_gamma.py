#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
from nexus import generate_convert4qmc
from nexus import generate_qmcpack

settings(
    pseudo_dir = '../../pseudopotentials',
    results    = '',
    sleep      = 3,
    machine    = 'ws16',
    )

system = generate_physical_system(
    units    = 'A',
    axes     = '''1.785   1.785   0.000
                  0.000   1.785   1.785
                  1.785   0.000   1.785''',
    elem_pos = '''
               C  0.0000  0.0000  0.0000
               C  0.8925  0.8925  0.8925
               ''',
    tiling   = (2,1,1),
    kgrid    = (1,1,1),
    kshift   = (0,0,0),
    C        = 4,
    )

scf = generate_pyscf(
    identifier = 'scf',                      # log output goes to scf.out
    path       = 'diamond/scf',              # directory to run in
    job        = job(serial=True,threads=16),# pyscf must run w/o mpi
    template   = './scf_template.py',        # pyscf template file
    system     = system,
    cell       = obj(                        # used to make Cell() inputs
        basis         = 'bfd-vdz',
        ecp           = 'bfd',
        drop_exponent = 0.1,
        verbose       = 5,
        ),
    save_qmc   = True,                # save wfn data for qmcpack
    )

c4q = generate_convert4qmc(
    identifier   = 'c4q',
    path         = 'diamond/scf',
    job          = job(cores=1),
    no_jastrow   = True,
    hdf5         = True,              # use hdf5 format
    dependencies = (scf,'orbitals'),
    )

opt = generate_qmcpack(
    block           = True,
    identifier      = 'opt',
    path            = 'diamond/optJ2',
    job             = job(cores=16,threads=4,app='qmcpack'),
    input_type      = 'basic',
    system          = system,
    pseudos         = ['C.BFD.xml'],
    corrections     = [],
    J2              = True,
    qmc             = 'opt',
    minmethod       = 'oneshiftonly', # adjust for oneshift
    init_cycles     = 3,
    init_minwalkers = 0.1,
    cycles          = 3,
    samples         = 25600,
    dependencies    = (c4q,'orbitals'),
    )

qmc = generate_qmcpack(
    block        = True,
    identifier   = 'vmc',
    path         = 'diamond/vmc',
    job          = job(cores=16,threads=4,app='qmcpack'),
    input_type   = 'basic',
    system       = system,
    pseudos      = ['C.BFD.xml'],
    corrections  = [],
    qmc          = 'vmc',
    dependencies = [(c4q,'orbitals'),
                    (opt,'jastrow')],
    )

qmc = generate_qmcpack(
    block        = True,
    identifier   = 'dmc',
    path         = 'diamond/dmc',
    job          = job(cores=16,threads=4,app='qmcpack'),
    input_type   = 'basic',
    system       = system,
    pseudos      = ['C.BFD.xml'],
    corrections  = [],
    qmc          = 'dmc',
    vmc_samples  = 800,
    eq_dmc       = True,
    dependencies = [(c4q,'orbitals'),
                    (opt,'jastrow')],
    )

run_project()
