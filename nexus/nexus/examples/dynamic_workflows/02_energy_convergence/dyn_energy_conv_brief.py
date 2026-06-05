#! /usr/bin/env python3

from nexus import settings,job,workflow_manager
from nexus import generate_physical_system
from nexus import generate_pwscf

'''
A simple type of dynamic workflow is to automatically determine 
converged parameter values.  In DFT, two such cases are convergence 
of the total energy with respect to increasing planewave energy 
cutoff and to increasingly large k-poing grids.  

This example first iteratively finds a converged planewave energy 
cutoff for a diamond primitive cell using a BFD potential for carbon 
and terminating when successive iterations produce total energies 
within 1e-4 Ry of each other.  The resulting energy cutoff is then 
fed into another iterative convergence procedure for the k-point grid, 
which stops when a tolerance of 1e-3 Ry has been met.  The converged 
values for the energy cutoff and k-point grid are 330 Ry and 7x7x7, 
respectively.  

From this example, it is clear how the workflow can be extended to 
include subsequent supercell expansion, Jastrow factor optimization, 
and VMC/DMC runs with total energies converged below a specified 
minimum statistical errorbar.
'''

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
    if nkgrid is None:
        nkgrid = 1
        path = '01_ecut_conv/ecut_'+str(ecutwfc)
    else:
        path = '02_kgrid_conv/kgrid_{0}{0}{0}'.format(nkgrid)
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
energies  = []
ecut      = 50  
tol       = 1e-4
converged = False
ecuts     = []
qe = gen_qe(ecutwfc=ecut)
while not converged:
    if qe.succ:
        energies.append(qe.products.energy)
        ecuts.append(ecut)
        if len(energies)<2 or abs(energies[-1]-energies[-2])>tol:
            # dynamic portion
            ecut = int(((1.5*ecut)//10)*10)
            qe = gen_qe(ecutwfc=ecut)
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
qe = gen_qe(ecutwfc=ecut,nkgrid=nkgrid)
while not converged:
    if qe.succ:
        energies.append(qe.products.energy)
        nkgrids.append(nkgrid)
        if len(energies)<2 or abs(energies[-1]-energies[-2])>tol:
            # dynamic portion
            nkgrid += 1
            qe = gen_qe(ecutwfc=ecut,nkgrid=nkgrid)
        else:
            print('Converged!!!  Final kgrid: {0}x{0}x{0}'.format(nkgrid))
            converged = True
    elif qe.fail:
        print('\nQE run failed!!!')
        break
    wm.poll(1)



