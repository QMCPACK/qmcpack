#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
from nexus import generate_pyscf_to_afqmc
from nexus import generate_qmcpack

settings(
    results = '',
    sleep   = 3,
    machine = 'ws16',
    )

system = generate_physical_system(
    units    = 'A',
    elem_pos = 'Ne 0 0 0',
    )

scf = generate_pyscf(
    identifier = 'scf',
    path       = 'rhf',
    job        = job(serial=True),
    template   = './scf_template.py',
    system     = system,
    mole       = obj(
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
    job          = job(cores=1,app='qmcpack'),
    input_type   = 'basic_afqmc',
    dependencies = (p2a,'wavefunction'),
    )

run_project()
