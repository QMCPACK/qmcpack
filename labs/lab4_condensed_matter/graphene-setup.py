#! /usr/bin/env python

# import Nexus functions
from nexus import settings,job,run_project,get_machine
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
    sleep         = 3,
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
directory = 'graphene' # directory to perform runs
dft_pps   = ['C.BFD.upf']   # pwscf pseudopotentials
qmc_pps   = ['C.BFD.xml']   # qmcpack pseudopotentials

# job details
dft_job = job(nodes=4,minutes=10,queue="qmcpack",app=pwscf)
p2q_job = job(cores=1,minutes=10,queue="qmcpack",app=pw2qmcpack)
qmc_job = job(nodes=32,minutes=10,threads=16,queue="qmcpack",app=qmcpack)


# create 2 atom sheet of graphene



graphene = generate_physical_system(
    axes      = [[9.30501148, 0.00000000,   0.0000000],
                 [-4.6525058, 8.05837632,   0.0000000],
                 [0.00000000, 0.00000000,  15.0000000]],
    elem      = ['C','C','C','C','C','C','C','C'],
    pos       = [[ 0.00000000, 0.00000000, 7.50000000],
                 [ 2.32625287, 1.34306272, 7.50000000],
                 [ 4.65250574, 0.00000000, 7.50000000],
                 [ 6.97875861, 1.34306272, 7.50000000],
                 [-2.32625290, 4.02918816, 7.50000000],
                 [-0.00000003, 5.37225088, 7.50000000],
                 [ 2.32625284, 4.02918816, 7.50000000],
                 [ 4.65250571, 5.37225088, 7.50000000]],

    units     = 'B',
    kgrid     = (1,1,1),
    kshift    = (0,0,0),
    net_charge= 0,
    net_spin  = 0,
    C         = 4
    )

sims = []

# scf run to generate converged charge density
scf = generate_pwscf(
    identifier   = 'scf',
    path         = directory+'/scf',
    job          = dft_job,
    input_type   = 'scf',
    system       = graphene,
    pseudos      = dft_pps,
    input_dft    = 'lda', 
    ecut         = 200,
    conv_thr     = 1e-8, 
    mixing_beta  = .7,
    nosym        = False,
    wf_collect   = False,
    kgrid        = (4,4,1),
    kshift       = (0,0,0)
    )
sims.append(scf)

# nscf run to generate orbitals
nscf = generate_pwscf(
    identifier   = 'nscf',
    path         = directory+'/nscf',
    job          = dft_job,
    input_type   = 'nscf',
    system       = graphene,
    pseudos      = dft_pps,
    input_dft    = 'lda', 
    ecut         = 200,
    conv_thr     = 1e-8, 
    mixing_beta  = .7,
    tprnfor      = False,
    tstress      = False,
    nosym        = True,
    wf_collect   = True,
    dependencies = (scf,'charge_density')
    )
sims.append(nscf)

# conversion step to create h5 file with orbitals
p2q = generate_pw2qmcpack(
    identifier   = 'p2q',
    path         = directory+'/nscf',
    job          = p2q_job,
    write_psir   = False,
    dependencies = (nscf,'orbitals')
    )
sims.append(p2q)

# optimization inputs
linopt1 = linear(
    energy               = 0.0,
    unreweightedvariance = 1.0,
    reweightedvariance   = 0.0,
    timestep             = 0.4,
    samples              = 16384,
    warmupsteps          = 40,
    blocks               = 32,
    steps                = 16,
    substeps             = 1,
    nonlocalpp           = False,
    usebuffer            = False,
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
linopt2.samples = 32768
linopt2.nonlocalpp = True
linopt2.usebuffer = True
linopt2.unreweightedvariance = 0.2
linopt2.energy = 0.8

# optimization run
opt = generate_qmcpack(
    identifier   = 'opt',
    path         = directory+'/opt',
    job          = qmc_job,
    input_type   = 'basic',
    system       = graphene,
    bconds       = 'ppn',
    pseudos      = qmc_pps,
    jastrows     = [('J1','bspline',8),
                    ('J2','bspline',8,'coeff',[8*[0],8*[0]])],
    calculations = [loop(max=4,qmc=linopt1),
                    loop(max=4,qmc=linopt2)],
    dependencies = (p2q,'orbitals')
    )
sims.append(opt)

run_project(sims)

