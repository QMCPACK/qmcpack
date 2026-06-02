#! /usr/bin/env python3

from nexus import settings,job,workflow_manager
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack


settings(
    results    = '',
    pseudo_dir = '../../qmcpack/pseudopotentials',
    machine    = 'ws8',
    dynamic    = True,
    )


def gen_system(tiling=None,kgrid=(1,1,1)):
    system = generate_physical_system(
        structure = 'diamond',
        cell      = 'prim',
        tiling    = tiling,
        kgrid     = kgrid,
        C         = 4,
        )
    return system


def gen_qe(run_type   = 'scf',
           dynamic_id = 'qe',
           path       = None,
           system     = None,
           ecutwfc    = None,
           nkgrid     = 1,
           ):
    assert run_type in ('scf','nscf')
    if run_type=='scf':
        kgrid = 3*[nkgrid]
        extra = dict(kgrid=kgrid,requires='none')
    elif run_type=='nscf':
        extra = dict(nosym=True,requires='charge_density')
    if system is None:
        system = gen_system()
    qe = generate_pwscf(
        identifier = run_type,
        path       = path,
        job        = job(cores=4),
        system     = system,
        pseudos    = ['C.BFD.upf'],
        input_type = 'generic',
        ecutwfc    = ecutwfc,
        dynamic_id = dynamic_id,
        **extra
        )
    return qe


# start the workflow manager
wm = workflow_manager()


# determine planewave energy cutoff
ecut      = 50  
tol       = 1e-4
converged = False
ecuts     = []
path_str  = lambda ecut: '01_ecut_conv/ecut_'+str(ecut)
qe = gen_qe(path=path_str(ecut),ecutwfc=ecut)
while not converged:
    if qe.succ:
        energies.append(qe.products.energy)
        ecuts.append(ecut)
        if len(energies)<2 or abs(energies[-1]-energies[-2])>tol:
            # dynamic portion
            ecut = int(((1.5*ecut)//10)*10)
            qe = gen_qe(path=path_str(ecut),ecutwfc=ecut)
        else:
            print('Converged!!!  Final energy cutoff: '+str(ecut))
            converged = True
    elif qe.fail:
        print('\nQE run failed!!!')
        break
    wm.poll(1)


# determine k-point grid, using the converged energy cutoff
energies  = []
nkgrid    = 1
tol       = 1e-3
converged = False
nkgrids   = []
path_str  = lambda nk: '02_kgrid_conv/kgrid_{0}{0}{0}'.format(nk)
qe = gen_qe(path=path_str(nkgrid),ecutwfc=ecut,nkgrid=nkgrid)
while not converged:
    if qe.succ:
        energies.append(qe.products.energy)
        nkgrids.append(nkgrid)
        if len(energies)<2 or abs(energies[-1]-energies[-2])>tol:
            # dynamic portion
            nkgrid += 1
            qe = gen_qe(path=path_str(nkgrid),ecutwfc=ecut,nkgrid=nkgrid)
        else:
            print('Converged!!!  Final kgrid: {0}x{0}x{0}'.format(nkgrid))
            converged = True
    elif qe.fail:
        print('\nQE run failed!!!')
        break
    wm.poll(1)



