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


# convenience function to make QE sims
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


# create the dynamic workflow manager
wm = workflow_manager()



#==============================================#
# Planewave energy cutoff convergence workflow #
#                                              #
#  Iteratively increase the cutoff until the   #
#  energy difference between iterations falls  #
#  below a tolerance threshold.                #
#==============================================#

# convenience function for log printing
def print_progress(ecut,energies):
    print()
    print(50*'=')
    print('  ecuts   :',ecuts)
    print('  energies:',energies)
    if len(energies)>=2:
        dE = abs(energies[-1]-energies[-2])
        print('  dE = {:7.5f} , tol = {:7.5f}'.format(dE,tol))
    print(50*'=')
    print()


print(3*'\n'+70*'*')
print('Determining converged PW energy cutoff')
print(70*'*'+'\n')
energies  = []     # total energy history
ecut      = 50     # initial energy cutoff
tol       = 1e-4   # tolerance for convergence (Ry)
converged = False
ecuts = []
path_str = lambda ecut: '01_ecut_conv/ecut_'+str(ecut)

# Start with a single QE simulation, ecut = 50 Ry
qe = gen_qe(path=path_str(ecut),ecutwfc=ecut)

# Convergence iterations
while not converged:
    print('poll')

    if qe.succ:
        # Current QE run succeeded
        # Capture the total energy and the energy cutoff
        energies.append(qe.products.energy)
        ecuts.append(ecut)
        # Check if the tolerance has been met
        if len(energies)<2 or abs(energies[-1]-energies[-2])>tol:
            # If not, increase the energy cutoff by 50%
            ecut = int(((1.5*ecut)//10)*10)
            # Run subsequent QE with higher cutoff
            qe = gen_qe(path=path_str(ecut),ecutwfc=ecut)
            print_progress(ecut,energies)
        else:
            # When converged, report the final energy cutoff
            print_progress(ecut,energies)
            print('\n'+50*'*')
            print('Converged!!!  Final energy cutoff: '+str(ecut))
            print(50*'*')
            converged = True # to exit polling loop
    elif qe.fail:
        print('\nQE run failed!!!')
        break

    # Run simulations actively upon poll
    wm.poll(1)



#==============================================#
# k-point grid convergence workflow            #
#                                              #
#  Iteratively increase the k-point grid until #
#  the energy difference between iterations    #
#  falls below a tolerance threshold.          #
#==============================================#

# convenience function for log printing
def print_progress(nkgrid,energies):
    print()
    print(50*'=')
    print('  nkgrids :',nkgrids)
    print('  energies:',energies)
    if len(energies)>=2:
        dE = abs(energies[-1]-energies[-2])
        print('  dE = {:6.4f} , tol = {:6.4f}'.format(dE,tol))
    print(50*'=')
    print()

print(3*'\n'+70*'*')
print('Determining converged k-point grid')
print(70*'*'+'\n')
energies  = []
nkgrid    = 1
tol       = 1e-3
converged = False
nkgrids   = []
path_str  = lambda nk: '02_kgrid_conv/kgrid_{0}{0}{0}'.format(nk)
qe = gen_qe(path=path_str(nkgrid),ecutwfc=ecut,nkgrid=nkgrid)
while not converged:
    print('poll')
    if qe.succ:
        energies.append(qe.products.energy)
        nkgrids.append(nkgrid)
        if len(energies)<2 or abs(energies[-1]-energies[-2])>tol:
            nkgrid += 1
            qe = gen_qe(path=path_str(nkgrid),ecutwfc=ecut,nkgrid=nkgrid)
            print_progress(nkgrid,energies)
        else:
            print_progress(nkgrid,energies)
            print('\n'+50*'*')
            print('Converged!!!  Final kgrid: {0}x{0}x{0}'.format(nkgrid))
            print(50*'*')
            converged = True # to exit polling loop
    elif qe.fail:
        print('\nQE run failed!!!')
        break
    wm.poll(1)



