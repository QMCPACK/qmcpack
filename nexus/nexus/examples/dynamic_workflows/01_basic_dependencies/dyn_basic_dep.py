#! /usr/bin/env python3

import numpy as np
from nexus import settings,job
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import workflow_manager


settings(
    results    = '',
    pseudo_dir = '../../qmcpack/pseudopotentials',
    machine    = 'ws12',
    dynamic    = True,
    )


system = generate_physical_system(
    structure = 'diamond',
    cell      = 'prim',
    kgrid     = (2,1,1),
    C         = 4,
    )


# small random displacement
s = system.structure
disp = 0.003*np.random.randn(*s.pos.shape)
disp = np.abs(disp)
disp = np.dot(disp,s.axes)
s.pos += disp


relax = generate_pwscf(
    identifier = 'relax',
    path       = 'relax',
    job        = job(cores=8),
    system     = system,
    pseudos    = ['C.BFD.upf'],
    input_type = 'relax',
    input_dft  = 'pbe',
    ecut       = 200,
    conv_thr   = 1e-6,
    kgrid      = (2,1,1),
    kshift     = (1,1,1),
    # combined [path,identifier,dynamic_id] must be unique
    dynamic_id = 'relax1',
    # requires nothing (runs immediately)
    requires   = 'none', 
    )


scf = generate_pwscf(
    identifier = 'scf',
    path       = 'scf',
    job        = job(cores=4),
    system     = system,
    pseudos    = ['C.BFD.upf'],
    input_type = 'generic',
    ecutwfc    = 200,
    kgrid      = (2,2,1),
    nosym      = True,
    dynamic_id = 'scf1',
    # needs structure prior to run
    requires   = 'structure',
    )


nscf = generate_pwscf(
    identifier = 'nscf',
    path       = 'nscf',
    job        = job(cores=4),
    system     = system,
    pseudos    = ['C.BFD.upf'],
    input_type = 'generic',
    calculation = 'nscf',
    ecutwfc    = 200,
    dynamic_id = 'nscf1',
    # needs structure and charge density prior to run
    requires   = ('charge_density','structure'),
    )


# Workflow manager oversees dynamic workflow
wm = workflow_manager()


# Poll to submit/check runs
max_polls = 500
for i in range(max_polls):
    print('\npoll ' + str(i))

    if relax.succ:
        print('\nRelaxation completed successfully,\npassing on structure.')
        scf.structure  = relax.structure
        nscf.structure = relax.structure

    if scf.succ:
        print('\nSCF completed successfully,\npassing on charge density.')
        nscf.charge_density = scf.charge_density

    if relax.fail or scf.fail or nscf.done:
        if nscf.succ:
            print('\nNSCF completed successfully.')
        break

    # Runs execute immediately at poll if ready.
    # Poll should go at loops end.
    wm.poll(1) # Sleep for 1 second between polls.


# Check end status
if i==max_polls-1:
    print('\nMax polls exceeded without runs finishing.')
elif relax.fail:
    print('\nRelaxation failed.')
elif scf.fail:
    print('\nSCF failed.')
elif nscf.fail:
    print('\nSCF failed.')
else:
    print('\nAll runs completed successfully!')
