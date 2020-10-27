#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
from nexus import generate_pyscf_to_afqmc
from nexus import generate_qmcpack
from qmcpack_input import back_propagation,onerdm

settings(
    results = '',
    sleep   = 3,
    machine = 'ws16',
    )

system = generate_physical_system(
    units    = 'A',
    elem_pos = '''
               C  0.00000000  0.00000000  0.00000000
               H  0.00000000  0.00000000  1.10850000
               H  1.04510382  0.00000000 -0.36950000
               H -0.52255191  0.90508646 -0.36950000
               H -0.52255191 -0.90508646 -0.36950000
               ''',
    )

scf = generate_pyscf(
    identifier = 'scf',               
    path       = 'rhf',               
    job        = job(serial=True),    
    template   = './scf_template.py', 
    system     = system,
    mole       = obj(                 
        basis    = 'sto-3g',
        ),
    checkpoint = True,
    )

p2a = generate_pyscf_to_afqmc(
    identifier         = 'p2a',
    path               = 'rhf',
    job                = job(serial=True),
    cholesky_threshold = 1e-5,
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
    estimators   = [
        back_propagation(
            ortho      = 1,
            naverages  = 4,
            block_size = 2,
            nsteps     = 200,
            onerdm     = onerdm(),
            )
        ],
    dependencies = (p2a,'wavefunction'),
    )

run_project()
