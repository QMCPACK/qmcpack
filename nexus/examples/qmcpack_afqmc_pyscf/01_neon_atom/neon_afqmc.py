#! /usr/bin/env python

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
from nexus import generate_pyscf_to_afqmc
from nexus import generate_qmcpack


settings(
    results    = '',
    sleep      = 3,
    machine    = 'ws16',
    )

qmcpack = '/home/j1k/apps/qmcpack/qmcpack/build_afqmc_intel/bin/qmcpack'

system = generate_physical_system(
    units    = 'A',
    elem_pos = 'Ne 0 0 0',
    )

scf = generate_pyscf(
    identifier = 'scf',               # log output goes to scf.out
    path       = 'rhf',               # directory to run in
    job        = job(serial=True),    # run pyscf serially         
    template   = './scf_template.py', # pyscf template file
    system     = system,
    mole       = obj(                 # used to make Mole() inputs
        verbose  = 4,
        basis    = 'aug-cc-pvdz',
        ),
    checkpoint = True,
    )

p2a = generate_pyscf_to_afqmc(
    identifier         = 'p2a',
    path               = 'rhf',
    job                = job(serial=True),
    cholesky_threshold = 1e-5,
    dependencies       = (scf,'wavefunction'),
    )

qmc = generate_qmcpack(
    identifier   = 'qmc',
    path         = 'afqmc',
    job          = job(cores=1,app=qmcpack),
    input_type   = 'basic_afqmc',
    dependencies = (p2a,'wavefunction'),
    )

run_project()
