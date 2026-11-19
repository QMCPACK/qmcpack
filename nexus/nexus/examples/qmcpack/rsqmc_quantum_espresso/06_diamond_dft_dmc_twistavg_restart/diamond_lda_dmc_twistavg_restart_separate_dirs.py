#! /usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc,dmc

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
    tiling   = [[ 1, -1,  1],
                [ 1,  1, -1],
                [-1,  1,  1]],
    kgrid    = (2,2,2),
    kshift   = (0,0,0),
    C        = 4,
    )

scf = generate_pwscf(
    identifier  = 'scf',
    path        = 'diamond/scf',
    job         = job(cores=16,app='pw.x'),
    input_type  = 'generic',
    calculation = 'scf',
    input_dft   = 'lda',
    ecutwfc     = 200,
    conv_thr    = 1e-8,
    system      = system,
    pseudos     = ['C.BFD.upf'],
    kgrid       = (4,4,4),
    kshift      = (0,0,0),
    )

nscf = generate_pwscf(
    identifier   = 'nscf',
    path         = 'diamond/nscf_twist',
    job          = job(cores=16,app='pw.x'),
    input_type   = 'generic',
    calculation  = 'nscf',
    input_dft    = 'lda',
    ecutwfc      = 200,
    conv_thr     = 1e-8,
    system       = system,
    pseudos      = ['C.BFD.upf'],
    nosym        = True,
    dependencies = (scf,'charge_density'),
    )

conv = generate_pw2qmcpack(
    identifier   = 'conv',
    path         = 'diamond/nscf_twist',
    job          = job(cores=16,app='pw2qmcpack.x'),
    write_psir   = False,
    dependencies = (nscf,'orbitals'),
    )

qmc1 = generate_qmcpack(
    identifier   = 'dmc',
    path         = 'diamond/dmc_twist_initial',
    job          = job(cores=16,threads=2,app='qmcpack'),
    input_type   = 'basic',
    driver       = 'batched',
    system       = system,
    pseudos      = ['C.BFD.xml'],
    J2           = True,
    calculations = [
        vmc(
            total_walkers = 256,
            blocks        = 20,
            steps         = 10,
            checkpoint    = 0,
            ),
        dmc(
            total_walkers = 256,
            blocks        = 20,
            steps         = 10,
            timestep      = 0.01,
            checkpoint    = 0,
            ),
        ],
    dependencies = (conv,'orbitals'),
    )

qmc2 = generate_qmcpack(
    identifier   = 'dmc_restart',
    path         = 'diamond/dmc_twist_restart',
    job          = job(cores=16,threads=2,app='qmcpack'),
    input_type   = 'basic',
    driver       = 'batched',
    system       = system,
    pseudos      = ['C.BFD.xml'],
    J2           = True,
    calculations = [
        dmc(
            total_walkers = 256,
            blocks        = 20,
            steps         = 10,
            timestep      = 0.01,
            checkpoint    = 0,
            ),
        ],
    dependencies = [(conv,'orbitals'),
                    (qmc1,'restart')],
    )

run_project()
