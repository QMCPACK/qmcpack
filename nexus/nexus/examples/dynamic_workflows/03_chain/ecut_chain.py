#! /usr/bin/env python3

'''
Same diamond / BFD setup and stepping as the ecut part of example 02. 
'''

import sys
from nexus import settings, job, workflow_manager
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus.dynamic_workflows.pwscf import Ecut


settings(
    results    = '',
    pseudo_dir = '../../qmcpack/pseudopotentials',
    machine    = 'ws8',
    dynamic    = True,
    )


system = generate_physical_system(
    structure = 'diamond',
    cell      = 'prim',
    kgrid     = (1, 1, 1),
    C         = 4,
    )

chain = Ecut(start=50, tol=1e-4, max_runs=10)


def make_scf(params):
    ecut = params['ecutwfc']
    return generate_pwscf(
        identifier = 'scf',
        path       = f'ecut_{ecut}',
        job        = job(cores=4),
        system     = system,
        pseudos    = ['C.BFD.upf'],
        input_type = 'generic',
        ecutwfc    = ecut,
        kgrid      = (1, 1, 1),
        dynamic_id = f'ecut{ecut}',
        requires   = 'none',
        )
#end def make_scf

wm = workflow_manager()
decision = chain.drive(make_scf, wm, products='energy')
print()
print('status :', decision.status)
print('ecut   :', decision.params['ecutwfc'])
print('products:', decision.products)
if not decision.completed:
    sys.exit(1)
#end if
