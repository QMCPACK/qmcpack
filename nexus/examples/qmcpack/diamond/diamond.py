#! /usr/bin/env python

from nexus import settings,Job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc

settings(
    pseudo_dir    = '../pseudopotentials',
    status_only   = 0,
    #generate_only   = False,
    # Complicated setting only so examples can be run in test harness.
    # For real runs, use the plain setting of 'generate_only' above.
    generate_only   = globals().get('override_generate_only_setting',False),
    sleep         = 3,
    machine       = 'ws16'
    )

dia16 = generate_physical_system(
    units  = 'A',
    axes   = [[ 1.785,  1.785,  0.   ],
              [ 0.   ,  1.785,  1.785],
              [ 1.785,  0.   ,  1.785]],
    elem   = ['C','C'],
    pos    = [[ 0.    ,  0.    ,  0.    ],
              [ 0.8925,  0.8925,  0.8925]],
    tiling = (2,2,2),
    kgrid  = (1,1,1),
    kshift = (0,0,0),
    C      = 4
    )
              
scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'diamond/scf',
    job          = Job(cores=16,app='pw.x'),
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    nosym        = True,
    wf_collect   = True,
    system       = dia16,
    pseudos      = ['C.BFD.upf'], 
    )

conv = generate_pw2qmcpack(
    identifier   = 'conv',
    path         = 'diamond/scf',
    job          = Job(cores=1,app='pw2qmcpack.x'),
    write_psir   = False,
    dependencies = (scf,'orbitals')
    )

qmc = generate_qmcpack(
    identifier   = 'vmc',
    path         = 'diamond/vmc',
    job          = Job(cores=16,threads=4,app='qmcapp'),
    input_type   = 'basic',
    system       = dia16,
    pseudos      = ['C.BFD.xml'],
    jastrows     = [],
    calculations = [
        vmc(
            walkers     =   1,
            warmupsteps =  20,
            blocks      = 200,
            steps       =  10,
            substeps    =   2,
            timestep    =  .4
            )
        ],
    dependencies = (conv,'orbitals')
    )

run_project(scf,conv,qmc)
