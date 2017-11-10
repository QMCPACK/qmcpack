#! /usr/bin/env python

# import nexus functions
from nexus import settings,Job,get_machine,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack
from nexus import loop,linear,vmc,dmc

# nexus settings
settings(
    pseudo_dir    = '../pseudopotentials',
    runs          = '',
    results       = '',
    status_only   = 0,
    #generate_only   = False,
    # Complicated setting only so examples can be run in test harness.
    # For real runs, use the plain setting of 'generate_only' above.
    generate_only   = globals().get('override_generate_only_setting',False),
    sleep         = 3,
    machine       = 'ws16',
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
    type       = 'dimer',      # dimer selected
    dimer      = ('O','O'),    # atoms in dimer
    separation = 1.2074*scale, # dimer bond length
    Lbox       = 15.0,         # box size
    units      = 'A',          # Angstrom units
    net_spin   = 2,            # Nup-Ndown = 2
    O          = 6             # O has 6 val. electrons
    )

# describe scf run
scf = generate_pwscf(
    identifier   = 'scf',
    path         = directory,
    system       = dimer,
    job          = Job(cores=16),
    input_type   = 'scf',
    pseudos      = ['O.BFD.upf'],
    input_dft    = 'lda',
    ecut         = 200,
    conv_thr     = 1e-7,
    mixing_beta  = .7,
    nosym        = True,
    wf_collect   = True
    )
sims.append(scf)

# describe orbital conversion
p2q = generate_pw2qmcpack(
    identifier   = 'p2q',
    path         = directory,
    job          = Job(cores=1),
    write_psir   = False,
    dependencies = (scf,'orbitals')
    )
sims.append(p2q)
    
# describe optimization run
opt = generate_qmcpack(
    identifier   = 'opt',
    path         = directory,
    system       = dimer,
    job          = Job(cores=16,app='qmcapp'),
    input_type   = 'basic',
    pseudos      = ['O.BFD.xml'],
    bconds       = 'nnn',
    jastrows     = [('J1','bspline',8,4.5), # 1 & 2 body Jastrows
                    ('J2','pade',0.5,0.5)],
    calculations = opt_calcs,
    dependencies = (p2q,'orbitals')
    )
sims.append(opt)

# describe DMC run
qmc = generate_qmcpack(
    identifier   = 'qmc',
    path         = directory,
    system       = dimer,
    job          = Job(cores=16,app='qmcapp'),
    input_type   = 'basic',
    pseudos      = ['O.BFD.xml'],
    bconds       = 'nnn',
    jastrows     = [],            
    calculations = qmc_calcs,
    dependencies = [(p2q,'orbitals'),
                    (opt,'jastrow')]
    )
sims.append(qmc)

# execute all simulations
run_project(sims) 
