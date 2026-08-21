#! /usr/bin/env python3

'''
k-point mesh walk using ``Kgrid`` (walk mode).

Same diamond / BFD setup and stepping as the k-grid half of
``02_energy_convergence`` (start 1x1x1, +1 per axis, successive
|ΔE| ≤ 1e-3 Ry).  ``drive`` owns the poll loop; ``Kgrid`` only
decides the next ``kgrid`` and when the target is reached.
``generate_pwscf`` stays in this script.

Set ``ecutwfc`` to the cutoff from ``ecut_walk.py`` (example 02
uses 330 Ry).  50 Ry is used here so this script can run on its
own.  For a spacing-based start instead of ``start=1``::

    walk = Kgrid(spacing=0.5, structure=system, increment=1, tol=1e-3)
'''

import sys
from nexus import settings, job, workflow_manager
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus.dynamic_workflows.pwscf import Kgrid


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

ecutwfc = 50
walk    = Kgrid(start=1, increment=1, tol=1e-3, max_runs=10)


def make_scf(params):
    kgrid = params['kgrid']
    n1, n2, n3 = kgrid
    return generate_pwscf(
        identifier = 'scf',
        path       = f'kgrid_{n1}{n2}{n3}',
        job        = job(cores=4),
        system     = system,
        pseudos    = ['C.BFD.upf'],
        input_type = 'generic',
        ecutwfc    = ecutwfc,
        kgrid      = kgrid,
        dynamic_id = f'kgrid{n1}{n2}{n3}',
        requires   = 'none',
        )
#end def make_scf

wm = workflow_manager()
decision = walk.drive(make_scf, wm, products='energy')
print()
print('status :', decision.status)
print('kgrid  :', decision.params['kgrid'])
print('products:', decision.products)
if not decision.completed:
    sys.exit(1)
#end if
