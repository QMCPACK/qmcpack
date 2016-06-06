#! /usr/bin/env python

# import Nexus functions
from nexus import settings,job,get_machine,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf,generate_pw2qmcpack
from nexus import generate_qmcpack,loop,linear,vmc,dmc

# Nexus settings
settings(
    pseudo_dir    = './pseudopotentials',
    runs          = '',
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    #machine       = 'ws4',
    machine       = 'vesta',
    account       = 'QMCPACK-Training',
    ) 

# specify job details
if settings.machine.startswith('ws'):    # running on workstation
    dftjob = job(cores=4,app='pw.x')
    p2qjob = job(cores=1,app='pw2qmcpack.x')
    qmcjob = job(cores=4,app='qmcpack')
else:                                    # running on Vesta
    appdir = '/soft/applications/qmcpack/Binaries/'
    dftjob  = job(nodes=32,threads= 1,queue='qmcpack',hours=1,app=appdir+'pw.x')
    p2qjob  = job(cores= 1,threads= 1,queue='qmcpack',hours=1,app=appdir+'pw2qmcpack.x')
    qmcjob  = job(nodes=32,threads=16,queue='qmcpack',hours=1,app=appdir+'qmcpack')

    vesta = get_machine('vesta') # allow one job at a time (lab only)
    vesta.queue_size = 2
#end if


# create DFT, OPT, & DMC sim's for each bond length
sims = []
scales = [1.00]
for scale in scales:
    directory = 'scale_'+str(scale)

    # make stretched/compressed dimer
    dimer = generate_physical_system(
        type       = 'dimer',
        dimer      = ('O','O'),
        separation = 1.2074*scale,
        Lbox       = 10.0,  # use 15.0 or so for production
        units      = 'A',
        net_spin   = 2,
        O          = 6,
        )

    # describe DFT run
    dft = generate_pwscf(
        identifier   = 'dft',
        path         = directory,
        job          = dftjob,
        system       = dimer,
        input_type   = 'scf',
        pseudos      = ['O.BFD.upf'],
        input_dft    = 'lda', 
        ecut         = 300,
        conv_thr     = 1e-7, 
        mixing_beta  = .7,
        nosym        = True,
        wf_collect   = True,
        )
    sims.append(dft)

    # describe orbital conversion run
    p2q = generate_pw2qmcpack(
        identifier   = 'p2q',
        path         = directory,
        job          = p2qjob,
        write_psir   = False,
        dependencies = (dft,'orbitals'),
        )
    sims.append(p2q)

    # describe optimization run
    if scale==1.00:                   # use eqm. Jastrow for all bond lengths
        opt = generate_qmcpack(
            identifier   = 'opt',
            path         = directory,
            job          = qmcjob,
            system       = dimer,
            input_type   = 'basic',
            pseudos      = ['O.BFD.xml'],
            bconds       = 'nnn',
            jastrows     = [('J1','bspline',8,5.0), # 1 & 2 body Jastrows
                            ('J2','bspline',8,10.0)],
            calculations = [
                loop(max=12,
                     qmc=linear(
                        energy               = 0.0,
                        unreweightedvariance = 1.0,
                        reweightedvariance   = 0.0,
                        timestep             = 0.3,
                        samples              = 61440,
                        warmupsteps          = 50,
                        blocks               = 200,
                        substeps             = 1,
                        nonlocalpp           = True,
                        usebuffer            = True,
                        walkers              = 1,
                        minwalkers           = 0.5,
                        maxweight            = 1e9, 
                        usedrift             = False,
                        minmethod            = 'quartic',
                        beta                 = 0.025,
                        exp0                 = -16,
                        bigchange            = 15.0,
                        alloweddifference    = 1e-4,
                        stepsize             = 0.2,
                        stabilizerscale      = 1.0,
                        nstabilizers         = 3,
                        )
                     )
                ],
            dependencies = (p2q,'orbitals'),
            )
        sims.append(opt)
    #end if

    # describe DMC run
    qmc = generate_qmcpack(
        identifier   = 'qmc',
        path         = directory,
        job          = qmcjob,
        system       = dimer,
        input_type   = 'basic',
        pseudos      = ['O.BFD.xml'],
        bconds       = 'nnn',
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
                warmupsteps   = 100, 
                blocks        = 400,
                steps         =  32,
                timestep      = 0.01,
                nonlocalmoves = True,
                )
            ],
        dependencies = [(p2q,'orbitals'),(opt,'jastrow') ],
        )
    sims.append(qmc)
#end for

# execute all simulations
run_project(sims) 
