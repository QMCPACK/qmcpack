#! /usr/bin/env python3
from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pwscf, generate_hp

import pdb

settings(
    pseudo_dir = '../pseudopotentials',
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
    C        = 4,
    )


scf1 = generate_pwscf(
    identifier   = 'scf',
    path         = 'diamond/scf',
    job          = job(cores=16,app='/Users/ksu/Software/qe-7.1/build_mpi/bin/pw.x'),
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    system       = system,
    pseudos      = ['C.BFD.upf'],
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    nogamma      = True,
    hubbard      = {'V' : {('C-2p', 'C-2p'): 1e-8, ('C-2p', 'C-2s'): 1e-8}},
)

hp = generate_hp(
    nq1          = 2,
    nq2          = 2,
    nq3          = 2,
    lmin         = 0,
    job          = job(cores=16,app='/Users/ksu/Software/qe-7.1/build_mpi/bin/hp.x'),
    path         = 'diamond/scf',
    dependencies = (scf1, 'other')
)

# scf2 = generate_pwscf(
#     identifier   = 'scf2',
#     path         = 'diamond/scf',
#     job          = job(cores=16,app='/Users/ksu/Software/qe-7.1/build_mpi/bin/pw.x'),
#     input_type   = 'generic',
#     calculation  = 'scf',
#     input_dft    = 'lda', 
#     ecutwfc      = 200,   
#     conv_thr     = 1e-8, 
#     system       = system,
#     pseudos      = ['C.BFD.upf'],
#     kgrid        = (4,4,4),
#     kshift       = (0,0,0),
#     hubbard_proj = 'atomic',
#     hubbard      = {'U' : {'C-2p':3.0},
#                     'J' : {'C-2p':2.0},
#                     'V' : {('C-2p', 'C-2s'): 1.0}},
#     nogamma      = True,
#     )
run_project()
