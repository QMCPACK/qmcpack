#!/usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system,read_structure
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc,dmc
from nexus import ppset,obj
from qmcpack_input import spindensity

settings(
    pseudo_dir    = './pseudopotentials',
    results       = '',
    sleep         = 3,
    machine       = 'ws16',
    status_only   = 0,
    generate_only = 0,
    )

ppset(
    label   = 'ccECP-soft',
    pwscf   = 'Fe.ccECP-soft.upf'.split(),
    qmcpack = 'Fe.ccECP-soft.xml'.split(),
    )
pseudos = 'ccECP-soft'

# inputs shared by all dft calculations
dft_shared = obj(
    input_type       = 'generic',
    input_dft        = 'lda',
    ecutwfc          = 400,
    bandfac          = 1.5,
    occupations      = 'smearing',
    smearing         = 'fermi-dirac',
    degauss          = 1e-3,
    mixing_beta      = 0.4,
    nspin            = 2,
    hubbard          = {'U':{'Fe-3d': 5.5}},
    hubbard_proj     = 'ortho-atomic',
    start_mag        = obj(Fe=0.1),
    kshift           = (0,0,0),
    pseudos          = pseudos,
    )

# inputs shared by all qmc calculations
qmc_shared = obj(
    driver         = 'batched',
    input_type     = 'basic',
    meshfactor     = 1.00,
    spin_polarized = True,
    pseudos        = pseudos,
    )

# define the physical system
system = generate_physical_system(
    structure  = 'Fe_bcc.xsf',
    kgrid      = (6,6,6),
    kshift     = (0,0,0),
    net_spin   = 6, # The net_spin to be used in Jastrow optimizations is chosen as an integer close to the SCF value of 5.66.
    net_charge = 0,
    symm_kgrid = False,
    tiling     = (1,1,1),
    Fe         = 16,
    )

basepath = './'

scf = generate_pwscf(
    identifier        = 'scf',
    path              = basepath + 'scf',
    job               = job(cores=8, app='pw.x'),
    calculation       = 'scf',
    system            = system,
    kgrid             = (24,24,24),
    **dft_shared
    )

nscf = generate_pwscf(
    identifier        = 'nscf',
    path              = basepath + 'nscf', # A different path is currently required for nscf since gcta-safl requires scf magnetization, which gets overwritten otherwise
    job               = job(cores=16, app='pw.x'),
    calculation       = 'nscf',
    system            = system,
    nosym             = True,
    dependencies      = (scf,'charge_density'),
    **dft_shared
    )

p2q = generate_pw2qmcpack(
    identifier        = 'p2q',
    path              = basepath + 'nscf',
    job               = job(cores=16,app='pw2qmcpack.x'),
    write_psir        = False,
    dependencies      = (nscf,'orbitals'),
    )

optJ12 = generate_qmcpack(
    identifier           = 'optJ12',
    path                 = basepath + 'optJ12',
    job                  = job(cores=16,threads=4,app='qmcpack'),
    system               = system,
    twistnum             = 0,
    qmc                  = 'opt',
    J2                   = True,
    minmethod            = 'oneshift',
    init_cycles          = 5,
    cycles               = 5,
    walkers_per_rank     = 8,
    warmupsteps          = 10,
    samples              = 2e5,
    blocks               = 100,
    minwalkers           = 0.3,
    corrections          = [],
    dependencies         = (p2q,'orbitals'),
    **qmc_shared
    )

J3_rcut = system.structure.rwigner()
optJ123 = generate_qmcpack(
    identifier           = 'optJ123',
    path                 = basepath + 'optJ123',
    job                  = job(cores=16,threads=4,app='qmcpack'),
    system               = system,
    twistnum             = 0,
    qmc                  = 'opt',
    J2                   = True,
    J3                   = True,
    J3_rcut              = J3_rcut,
    minmethod            = 'oneshift',
    init_cycles          = 0,
    cycles               = 10,
    walkers_per_rank     = 8,
    warmupsteps          = 10,
    samples              = 2e5,
    blocks               = 100,
    minwalkers           = 0.5,
    corrections          = [],
    dependencies         = [(p2q,'orbitals'), (optJ12,'jastrow')],
    **qmc_shared
    )

sdens = spindensity(
    grid = (100,100,100),
    )

qmc = generate_qmcpack(
    #skip_submit          = True,
    identifier           = 'qmc',
    path                 = basepath + 'qmc',
    job                  = job(cores=12,threads=4,app='qmcpack_complex'),
    system               = system,
    gcta                 = 'safl', # This is the only keyword needed to activate the GCTA occupations
    estimators           = [sdens],
    calculations         = [
    vmc(
        walkers_per_rank = 8,
        warmupsteps      = 100,
        blocks           = 20,
        steps            = 50,
        timestep         = 0.3,
        ),
    dmc( # equilibration
        walkers_per_rank = 8,
        timestep         = 0.02,
        warmupsteps      = 20,
        blocks           = 20,
        steps            = 20,
        nonlocalmoves    = 'v3',
        ),
    dmc( # main dmc
        walkers_per_rank = 8,
        timestep         = 0.005,
        warmupsteps      = 100,
        blocks           = 100,
        steps            = 200,
        nonlocalmoves    = 'v3',
        ),
    ],
    dependencies         = [(p2q,'orbitals'), (optJ123,'jastrow')],
    **qmc_shared
    )

run_project()
