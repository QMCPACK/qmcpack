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

a = 3.6

system = generate_physical_system(
    units  = 'A',
    axes   = [[   0, a/2, a/2 ],
              [ a/2,   0, a/2 ],
              [ a/2, a/2,   0 ]],
    elem   = ('C','C'),
    pos    = [[   0,   0,   0 ],
              [ a/4, a/4, a/4 ]],
    tiling = (2,2,2),
    kgrid  = (1,1,1),
    kshift = (0,0,0),
    )

scf = generate_pyscf(
    identifier = 'scf',               
    path       = 'rhf',               
    job        = job(serial=True),    
    template   = './scf_template.py', 
    system     = system,
    cell       = obj(
        basis   = 'gth-szv',
        pseudo  = 'gth-pade',
        mesh    = [25,25,25],
        verbose = 5,
        ),
    checkpoint = True,
    )

p2a = generate_pyscf_to_afqmc(
    identifier         = 'p2a',
    path               = 'rhf',
    job                = job(serial=True),
    cholesky_threshold = 1e-5,
    ao                 = True,
    kpoint             = True,
    verbose            = True,
    dependencies       = (scf,'wavefunction'),
    )

qmc = generate_qmcpack(
    identifier   = 'qmc',
    path         = 'afqmc',
    job          = job(cores=1,app='qmcpack'),
    system       = system,
    input_type   = 'basic_afqmc',
    blocks       = 100,
    timestep     = 0.01,
    dependencies = (p2a,'wavefunction'),
    )

run_project()
