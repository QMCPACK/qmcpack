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
qmcpack    = '/soft/applications/qmcpack/Binaries/qmcpack_comp'

# run directory and pseudopotentials
directory = 'bcc-beryllium' # directory to perform runs
dft_pps   = ['Be.ncpp']   # pwscf pseudopotentials
qmc_pps   = ['Be.xml']   # qmcpack pseudopotentials

# job details
dft_job = job(cores=16,minutes=10,queue="qmcpack",app=pwscf)
p2q_job = job(cores=1,minutes=10,queue="qmcpack",app=pw2qmcpack)
qmc_job = job(nodes=32,minutes=10,threads=16,queue="qmcpack",app=qmcpack)

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

    p2q = generate_pw2qmcpack(
        identifier   = 'p2q',
        path         = directory+'/nscf_'+ks,
        job          = p2q_job,
        write_psir   = False,
        #dependencies = (nscf,'orbitals'),
        )
    sims.append(p2q)

    if first:
        # optimization inputs
        linopt1 = linear(
            energy               = 0.0,
            unreweightedvariance = 1.0,
            reweightedvariance   = 0.0,
            timestep             = 0.4,
            samples              = 8192,
            warmupsteps          = 50,
            blocks               = 64,
            substeps             = 4,
            nonlocalpp           = True,
            usebuffer            = True,
            walkers              = 1,
            minwalkers           = 0.5,
            maxweight            = 1e9, 
            usedrift             = True,
            minmethod            = 'quartic',
            beta                 = 0.0,
            exp0                 = -16,
            bigchange            = 15.0,
            alloweddifference    = 1e-4,
            stepsize             = 0.2,
            stabilizerscale      = 1.0,
            nstabilizers         = 3
            )

        linopt2 = linopt1.copy()
        linopt2.samples = 16384

        # optimization run
        opt = generate_qmcpack(
            identifier   = 'opt',
            path         = directory+'/opt',
            job          = qmc_job,
            input_type   = 'basic',
            system       = bcc_Be,
            twistnum     = 0,
            bconds       = 'ppp',
            pseudos      = qmc_pps,
            jastrows     = [('J1','bspline',8),
                            ('J2','bspline',8)],
            calculations = [loop(max=4,qmc=linopt1),
                            loop(max=4,qmc=linopt2)],
            dependencies = (p2q,'orbitals')
            )
        sims.append(opt)
    #end if

    # DMC run
    ntwists = kgrid[0]*kgrid[1]*kgrid[2]
    qmc = generate_qmcpack(
        identifier   = 'dmc',
        path         = directory+'/dmc_'+ks,
        job          = job(nodes=ntwists,hours=2,threads=16,queue="qmcpack",app=qmcpack,queue="R.qmc"),
        input_type   = 'basic',
        system       = bcc_Be,
        bconds       = 'ppp',
        pseudos      = qmc_pps,
        jastrows     = [],            
        calculations = [
            vmc(
                walkers     =   1,
                warmupsteps =  30,
                blocks      =  20,
                steps       =  10,
                substeps    =   2,
                timestep    =  .4,
                samples     = 2048
                ),
            dmc(
                warmupsteps   =  20, 
                blocks        =  100,
                steps         =  10,
                timestep      = 0.04,
                nonlocalmoves = True
                )
            ],
        dependencies = [(p2q,'orbitals'),(opt,'jastrow')]
        )
    sims.append(qmc)
    
    first = False
#end for

# write input files and submit jobs
run_project(sims)

