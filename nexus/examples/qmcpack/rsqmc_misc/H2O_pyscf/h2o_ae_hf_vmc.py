#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
from nexus import generate_convert4qmc
from nexus import generate_qmcpack

settings(
    results    = '',
    sleep      = 3,
    machine    = 'ws16',
    )

system = generate_physical_system(
    structure = 'H2O.xyz',
    )

scf = generate_pyscf(
    identifier = 'scf',               # log output goes to scf.out
    path       = 'h2o_ae_hf',         # directory to run in
    job        = job(serial=True),    # pyscf must run serially         
    template   = './scf_template.py', # pyscf template file
    system     = system,
    mole       = obj(                 # used to make Mole() inputs
        verbose  = 5,
        basis    = 'ccpvtz',
        symmetry = True,
        ),
    save_qmc   = True,                # save wfn data for qmcpack
    )

c4q = generate_convert4qmc(
    identifier   = 'c4q',
    path         = 'h2o_ae_hf',
    job          = job(cores=1),
    no_jastrow   = True,
    dependencies = (scf,'orbitals'),
    )

qmc = generate_qmcpack(
    identifier   = 'vmc',
    path         = 'h2o_ae_hf',
    job          = job(cores=10),
    system       = system,
    input_type   = 'basic',
    jastrows     = [],
    qmc          = 'vmc',
    blocks       = 8000,
    steps        = 100,
    dependencies = (c4q,'orbitals'),
    )

run_project()
