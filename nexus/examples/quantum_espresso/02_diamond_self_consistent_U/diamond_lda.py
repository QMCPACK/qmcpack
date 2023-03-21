#! /usr/bin/env python3
import sys 
sys.path = ['/Users/ksu/Documents/GitHub/qmcpack/nexus/lib', 
            '/Users/ksu/Documents/GitHub/qmcpack/nexus/examples/quantum_espresso/01_diamond_scf',
            '/Users/ksu/Documents/GitHub/qmcpack/nexus', 
            # '/Users/ksu/Software/qmcpack/latest/nexus/lib', 
            '/Users/ksu/Documents/GitHub/pyqmc/pyqmc', 
            '/Users/ksu/Documents/GitHub/pyscf/pyscf', 
            '/Users/ksu/Documents/GitHub/qmcpack/nexus/examples/quantum_espresso/01_diamond_scf',
            '/usr/local/anaconda3/lib/python39.zip',
            '/usr/local/anaconda3/lib/python3.9',
            '/usr/local/anaconda3/lib/python3.9/lib-dynload',
            '/Users/ksu/.local/lib/python3.9/site-packages',
            '/usr/local/anaconda3/lib/python3.9/site-packages',
            '/usr/local/anaconda3/lib/python3.9/site-packages/aeosa']
from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pwscf, generate_hp

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

sims = []
for step in range(4):
    if step > 0:
        hubbard_result = (sims[-1], 'hubbard_parameters')
        hubbard      = None
    else:
        hubbard_result = []
        hubbard      = {'V' : {('C-2p', 'C-2p'): 1e-8}}

    scf = generate_pwscf(
        identifier   = 'scf',
        path         = 'diamond/scf_step_{}'.format(step),
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
        hubbard      = hubbard,
        dependencies = hubbard_result
    )
    sims.append(scf)

    hp = generate_hp(
        nq1          = 2,
        nq2          = 2,
        nq3          = 2,
        lmin         = 0,
        job          = job(cores=16,app='/Users/ksu/Software/qe-7.1/build_mpi/bin/hp.x'),
        path         = 'diamond/scf_step_{}'.format(step),
        dependencies = (sims[-1], 'other')
    )
    sims.append(hp)

run_project()
