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
    elem_pos = 'C 0 0 0',
    net_spin = 2,
    )

scf = generate_pyscf(
    identifier = 'scf',
    path       = 'uhf',
    job        = job(serial=True),
    template   = './scf_template.py',
    system     = system,
    mole       = obj(
        verbose  = 4,
        basis    = 'cc-pvtz',
        ),
    checkpoint = True,
    )

p2a = generate_pyscf_to_afqmc(
    identifier         = 'p2a',
    path               = 'uhf',
    job                = job(serial=True),
    cholesky_threshold = 1e-5,
    ao                 = True,
    verbose            = True,
    dependencies       = (scf,'wavefunction'),
    )

qmc = generate_qmcpack(
    identifier   = 'qmc',
    path         = 'afqmc',
    job          = job(cores=1,app='qmcpack'),
    system       = system,
    input_type   = 'basic_afqmc',
    blocks       = 1000,
    timestep     = 0.01,
    dependencies = (p2a,'wavefunction'),
    )

run_project()
