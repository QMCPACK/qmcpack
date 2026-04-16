#! /usr/bin/env python3

'''
This example is similar to diamond_lda_vmc.py (legacy drivers) 
but for the batched drivers.
'''

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
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
    tiling   = [[ 1, -1,  1],
                [ 1,  1, -1],
                [-1,  1,  1]],
    kgrid    = (1,1,1),
    kshift   = (0,0,0),
    C        = 4,
    )

scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'diamond/scf',
    job          = job(cores=16,app='pw.x'),
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    system       = system,
    pseudos      = ['C.BFD.upf'],
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    )

nscf = generate_pwscf(
    identifier   = 'nscf',
    path         = 'diamond/nscf',
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
    path         = 'diamond/nscf',
    job          = job(cores=16,app='pw2qmcpack.x'),
    write_psir   = False,
    dependencies = (nscf,'orbitals'),
    )

opt = generate_qmcpack(
    identifier   = 'opt',
    path         = 'diamond/optJ2_batched',
    job          = job(cores=16,threads=4,app='qmcpack'),
    input_type   = 'basic',
    system       = system,
    pseudos      = ['C.BFD.xml'],
    J2           = True,
    qmc          = 'opt',
    driver       = 'batched',
    cycles       = 6,
    samples      = 51200,
    dependencies = (conv,'orbitals'),
    )

qmc = generate_qmcpack(
    identifier       = 'vmc',
    path             = 'diamond/vmc_batched',
    job              = job(cores=16,threads=4,app='qmcpack'),
    input_type       = 'basic',
    system           = system,
    pseudos          = ['C.BFD.xml'],
    J2               = True,
    qmc              = 'vmc',
    driver           = 'batched',
    walkers_per_rank = 4,  # specify this way to get 1 walker per thread
    dependencies     = [(conv,'orbitals'),
                        (opt,'jastrow')],
    )

run_project()
