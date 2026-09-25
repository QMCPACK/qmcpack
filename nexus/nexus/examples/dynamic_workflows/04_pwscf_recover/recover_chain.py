#! /usr/bin/env python3

'''
Planewave cutoff chain with recoverable PWscf jobs.
'''

import sys
from nexus import settings, job, workflow_manager
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus.dynamic_workflows.pwscf import Ecut
from nexus.dynamic_workflows.pwscf.error_handler import PwscfErrorHandler

settings(
    results    = '',
    pseudo_dir = '../../qmcpack/pseudopotentials',
    machine    = 'ws8',
    dynamic    = True,
    verbose    = True,
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
        identifier       = 'scf',
        path             = f'ecut_{ecut}',
        job              = job(cores=4),
        system           = system,
        pseudos          = ['C.BFD.upf'],
        input_type       = 'generic',
        ecutwfc          = ecut,
        kgrid            = (1, 1, 1),
        electron_maxstep = 2, # Added to simulate failure
        mixing_beta      = 0.7,
        dynamic_id       = f'ecut{ecut}',
        requires         = 'none',
        error_handling   = True,
        # error_handling = 3,
        # error_handling = PwscfErrorHandler(max_runs=10, mixing_factor=0.3)
        )
#end def make_scf

wm = workflow_manager()
decision = chain.drive(make_scf, wm, products='energy')
print('status :', decision.status)
print('ecut   :', decision.params['ecutwfc'])
print('products:', decision.products)
