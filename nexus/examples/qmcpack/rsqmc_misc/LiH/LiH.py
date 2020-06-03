#! /usr/bin/env python3
# LiH crystal with Quantum ESPRESSO orbitals

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc,loop,linear,dmc

# General Settings (Directories For I/O, Machine Type, etc.)
settings(
    pseudo_dir      = 'pseudopotentials',
    runs            = 'runs',
    results         = 'results',
    sleep           = 3,
    generate_only   = 0,
    status_only     = 0,
    machine         = 'ws1',
    )

# Executables (Indicate Path If Needed)
pwscf               = 'pw.x'
pw2qmcpack          = 'pw2qmcpack.x'
qmcpack             = 'qmcpack'

# Pseudopotentials
dft_pps             = ['Li.TN-DF.upf','H.TN-DF.upf']
qmc_pps             = ['Li.pp.data','H.pp.data']

# job Definitions (MPI Tasks, MP Threading, PBS Queue, Time, etc.)
scf_job             = job(app=pwscf,serial=True)
nscf_job            = job(app=pwscf,serial=True)
p2q_job             = job(app=pw2qmcpack,serial=True)
opt_job             = job(threads=4,app=qmcpack,serial=True)
dmc_job             = job(threads=4,app=qmcpack,serial=True)

# System To Be Simulated
rocksalt_LiH = generate_physical_system(
    lattice         = 'cubic',
    cell            = 'primitive',
    centering       = 'F',
    atoms           = ('Li','H'),
    basis           = [[0.0,0.0,0.0],
                       [0.5,0.5,0.5]],
    basis_vectors   = 'conventional',
    constants       = 7.1,
    units           = 'B',
    kgrid           = (1,1,1),
    kshift          = (0,0,0),
    net_charge      = 0,
    net_spin        = 0,
    Li              = 1,
    H               = 1,
    )

# DFT SCF To Generate Converged Density
scf = generate_pwscf(
    identifier      = 'scf',
    path            = '.',
    job             = scf_job,
    input_type      = 'scf',
    system          = rocksalt_LiH,
    pseudos         = dft_pps,
    ecut            = 450,
    ecutrho         = 1800,
    conv_thr        = 1.0e-10,
    mixing_beta     = 0.7,
    kgrid           = (7,7,7),
    kshift          = (1,1,1),
    )

# DFT NSCF To Generate Wave Function At Specified K-points
nscf = generate_pwscf(
    identifier      = 'nscf',
    path            = '.',
    job             = nscf_job,
    input_type      = 'nscf',
    system          = rocksalt_LiH,
    pseudos         = dft_pps,
    ecut            = 450,
    ecutrho         = 1800,
    conv_thr        = 1.0e-10,
    mixing_beta     = 0.7,
    dependencies    = (scf,'charge-density'),
    )

# Convert DFT Wavefunction Into HDF5 File For QMCPACK
p2q = generate_pw2qmcpack(
    identifier      = 'p2q',
    path            = '.',
    job             = p2q_job,
    write_psir      = False,
    dependencies    = (nscf,'orbitals'),
    )

# QMC Optimization Parameters - Coarse Sampling Set
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
    nstabilizers         = 3,
    )

# QMC Optimization Parameters - Finer Sampling Set
linopt2 = linopt1.copy()
linopt2.samples = 16384

# QMC Optimization
opt = generate_qmcpack(
    identifier      = 'opt',
    path            = '.',
    job             = opt_job,
    input_type      = 'basic',
    system          = rocksalt_LiH,
    twistnum        = 0,
    bconds          = 'ppp',
    pseudos         = qmc_pps,
    jastrows        = [('J1','bspline',8),
                       ('J2','bspline',8)],
    calculations    = [loop(max=4,qmc=linopt1),
                       loop(max=4,qmc=linopt2)],
    dependencies    = (p2q,'orbitals'),
    )
pp = opt.input.get('pseudos')
pp.Li.set(
    format  = 'casino',
    l_local = 's',
    nrule   = 2,
    lmax    = 2,
    cutoff  = 2.19,
    )
pp.H.set(
    format  = 'casino',
    l_local = 's',
    nrule   = 2,
    lmax    = 2,
    cutoff  = 0.50,
    )

# QMC VMC/DMC With Optimized Jastrow Parameters
qmc = generate_qmcpack(
    identifier      = 'dmc',
    path            = '.',
    job             = dmc_job,
    input_type      = 'basic',
    system          = rocksalt_LiH,
    pseudos         = qmc_pps,
    bconds          = 'ppp',
    jastrows        = [],
    calculations    = [
        vmc(
            walkers              = 1,
            samplesperthread     = 64,
            stepsbetweensamples  = 1,
            substeps             = 5,
            warmupsteps          = 100,
            blocks               = 1,
            timestep             = 1.0,
            usedrift             = False,
           ),
        dmc(
            minimumtargetwalkers = 128,
            reconfiguration      = 'no',
            warmupsteps          = 100,
            timestep             = 0.005,
            steps                = 10,
            blocks               = 200,
            nonlocalmoves        = True,
           )
        ],
    dependencies   = [(p2q,'orbitals'),(opt,'jastrow')],
    )
pp = qmc.input.get('pseudos')
pp.Li.set(
    format  = 'casino',
    l_local = 's',
    nrule   = 2,
    lmax    = 2,
    cutoff  = 2.37,
    )
pp.H.set(
    format  = 'casino',
    l_local = 's',
    nrule   = 2,
    lmax    = 2,
    cutoff  = 0.50,
    )

run_project()
