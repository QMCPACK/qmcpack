#! /usr/bin/env python3
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

num_steps = 4
sims = []
for step in range(num_steps):
    if step > 0:
        hubbard_result = (sims[-1], 'hubbard_parameters')
        hubbard      = None
    else:
        hubbard_result = []
        hubbard      = {'V' : {('C-2p', 'C-2p'): 1e-8}}
        # Other use examples
        # hubbard      = {'V' : {('C-2p', 'C-2p'): [{'indices':(1,2), 'value':1e-8},
        #                                           {'indices':(2,1), 'value':1e-8}]}}
        # hubbard      = {'V' : {('C-2p', 'C-2p'): [{'radius' : 5.0, 'value':1e-8}]}} # radius is in units of Bohr
        # hubbard      = {'U':{'C-2p': 1.0},
        #                 'V':{('C-2p', 'C-2p'): 1e-8}}
    #end if

    scf = generate_pwscf(
        identifier   = 'scf',
        path         = 'diamond/scf_step_{}'.format(step),
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
        nogamma      = True,
        hubbard      = hubbard,
        dependencies = hubbard_result
    )
    scf.show_input()
    sims.append(scf)

    hp = generate_hp(
        nq1          = 2,
        nq2          = 2,
        nq3          = 2,
        lmin         = 0,
        job          = job(cores=16,app='hp.x'),
        path         = 'diamond/scf_step_{}'.format(step),
        dependencies = (sims[-1], 'other')
    )
    sims.append(hp)

run_project()
