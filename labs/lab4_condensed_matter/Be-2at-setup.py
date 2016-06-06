#! /usr/bin/env python

# import Nexus functions
from nexus import settings,job,run_project,get_machine,obj
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc,loop,linear,dmc

# Nexus settings
settings(
    pseudo_dir    = './pseudopotentials',
    runs          = '',
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 10,
    machine       = 'vesta',
    account       = 'QMCPACK-Training'
    )
 
# allow max of one job at a time (lab only)
vesta = get_machine('vesta')
vesta.queue_size = 1

# locations of pwscf, pw2qmcpack and qmcpack executables
pwscf      = '/soft/applications/qmcpack/Binaries/pw.x'
pw2qmcpack = '/soft/applications/qmcpack/Binaries/pw2qmcpack.x'
qmcpack    = '/soft/applications/qmcpack/Binaries/qmcpack'

# run directory and pseudopotentials
directory = 'bcc-beryllium' # directory to perform runs
dft_pps   = ['Be.ncpp']   # pwscf pseudopotentials
qmc_pps   = ['Be.xml']   # qmcpack pseudopotentials

# job details
dft_job = job(cores=16,hours=2,queue="qmcpack",app=pwscf)
p2q_job = job(cores=1,hours=2,queue="qmcpack",app=pw2qmcpack)
qmc_job = job(nodes=32,hours=2,threads=16,queue="qmcpack",app=qmcpack)

# specify k-point grids
kgrids = [(2,2,2),(3,3,3)]

sims = []
first = True
for kgrid in kgrids:
    ks = '{0}{1}{2}'.format(*kgrid)

    # create conventional cell tiled from primitive one
    bcc_Be = generate_physical_system(
        lattice    = 'cubic',
        cell       = 'primitive',
        centering  = 'I',
        atoms      = 'Be',
        constants  = 3.490,
        units      = 'A',
        net_charge = 0,
        net_spin   = 0,
        Be         = 2,
        tiling     = [[a,b,c],[d,e,f],[g,h,i]],
        kgrid      = kgrid,
        kshift     = (.5,.5,.5)
        )

    # scf run to generate converged charge density
    if first:
        scf = generate_pwscf(
            identifier     = 'scf',
            path           = directory+'/scf',
            job            = dft_job,
            input_type     = 'scf',
            system         = bcc_Be,
            spin_polarized = False,
            pseudos        = dft_pps,
            input_dft      = 'lda', 
            ecut           = 200,
            conv_thr       = 1e-8, 
            mixing_beta    = .7,
            nosym          = False,
            wf_collect     = False,
            kgrid          = (8,8,8),
            kshift         = (0,0,0)
            )
        sims.append(scf)
    #end if

    # nscf run to generate orbitals
    nscf = generate_pwscf(
        identifier     = 'nscf',
        path           = directory+'/nscf_'+ks,
        job            = dft_job,
        input_type     = 'nscf',
        system         = bcc_Be,
        spin_polarized = False,
        pseudos        = dft_pps,
        input_dft      = 'lda', 
        ecut           = 200,
        conv_thr       = 1e-8, 
        mixing_beta    = .7,
        tprnfor        = False,
        tstress        = False,
        nosym          = True,
        wf_collect     = True,
        dependencies   = (scf,'charge_density')
        )
    sims.append(nscf)

    # conversion step to create h5 file with orbitals
    p2q = generate_pw2qmcpack(
        identifier   = 'p2q',
        path         = directory+'/nscf_'+ks,
        job          = p2q_job,
        write_psir   = False,
        dependencies = (nscf,'orbitals'),
        )
    sims.append(p2q)
    first = False
#end for

# write input files and submit jobs
run_project(sims)

