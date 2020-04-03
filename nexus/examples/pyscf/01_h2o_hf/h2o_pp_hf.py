#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf

settings(
    results = '',
    sleep   = 3,
    machine = 'ws16',
    )

system = generate_physical_system(
    structure = 'H2O.xyz',
    O         = 6,                    # Zeff=6 for BFD oxygen
    H         = 1,                    # Zeff=1 for BFD hydrogen
    )

scf = generate_pyscf(
    identifier = 'scf',               # log output goes to scf.out
    path       = 'h2o_pp_hf',         # directory to run in
    job        = job(serial=True),    # pyscf must run serially         
    template   = './scf_template.py', # pyscf template file
    system     = system,
    mole       = obj(                 # used to make Mole() inputs
        basis    = 'bfd-vtz',
        ecp      = 'bfd',
        symmetry = True,
        ),
    )

run_project()
