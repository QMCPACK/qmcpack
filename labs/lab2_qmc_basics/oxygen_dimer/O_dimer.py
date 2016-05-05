#! /usr/bin/env python

# import project suite functions
from project import settings,Job,get_machine,run_project
from project import generate_physical_system
from project import generate_qmcpack,loop,linear,vmc,dmc

# project suite settings
settings(
    pseudo_dir    = './pseudopotentials',
    runs          = '',
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'vesta',
    account       = 'QMC_2014_training'
    ) 

 # allow max of one job at a time (lab only)
vesta = get_machine('vesta')
vesta.queue_size = 1

# specify job parameters
qmcjob = Job(    
    nodes   = 32,
    threads = 16,
    hours   = 1, 
    app     = '/soft/applications/qmcpack/build_XL_real/bin/qmcapp'
    )

# specify optimization parameters
linopt1 = linear(
    energy               = 0.0,
    unreweightedvariance = 1.0,
    reweightedvariance   = 0.0,
    timestep             = 0.4,
    samples              = 5120, # opt w/ 5120 samples
    warmupsteps          = 50,
    blocks               = 200,
    substeps             = 1,
    nonlocalpp           = True,
    usebuffer            = True,
    walkers              = 1,
    minwalkers           = 0.5,
    maxweight            = 1e9, 
    usedrift             = True,
    minmethod            = 'quartic',
    beta                 = 0.025,
    exp0                 = -16,
    bigchange            = 15.0,
    alloweddifference    = 1e-4,
    stepsize             = 0.2,
    stabilizerscale      = 1.0,
    nstabilizers         = 3
    )

linopt2 = linopt1.copy()  
linopt2.samples = 20480 # opt w/ 20000 samples

linopt3 = linopt1.copy()
linopt3.samples = 40960 # opt w/ 40000 samples

opt_calcs = [loop(max=8,qmc=linopt1), # loops over opt's
             loop(max=6,qmc=linopt2),
             loop(max=4,qmc=linopt3)]

# specify DMC parameters
qmc_calcs = [
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
        nonlocalmoves = True
        )
    ]

# create opt & DMC sim's for each bond length
sims = []
scale = 1.00
directory = 'scale_'+str(scale)

# make stretched/compressed dimer
dimer = generate_physical_system(
    type       = 'dimer',
    dimer      = ('O','O'),
    separation = 1.2074*scale,
    Lbox       = 15.0,
    units      = 'A',
    net_spin   = 2,
    O          = 6
    )

# describe optimization run
opt = generate_qmcpack(
    identifier   = 'opt',
    path         = directory,
    system       = dimer,
    job          = qmcjob,
    input_type   = 'basic',
    pseudos      = ['O.BFD.xml'],
    orbitals_h5  = 'O2.pwscf.h5',
    bconds       = 'nnn',
    jastrows     = [('J1','bspline',8,4.5), # 1 & 2 body Jastrows
                    ('J2','pade',0.5,0.5)],
    calculations = opt_calcs
    )
sims.append(opt)

# describe DMC run
qmc = generate_qmcpack(
    identifier   = 'qmc',
    path         = directory,
    system       = dimer,
    job          = qmcjob,
    input_type   = 'basic',
    pseudos      = ['O.BFD.xml'],
    orbitals_h5  = 'O2.pwscf.h5',
    bconds       = 'nnn',
    jastrows     = [],            
    calculations = qmc_calcs,
    dependencies = (opt,'jastrow') # Jastrows result from opt
    )
sims.append(qmc)

# execute all simulations
run_project(sims) 
