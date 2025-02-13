#!/usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc,dmc
from nexus import ppset,obj

settings(
    pseudo_dir    = './pseudopotentials',
    results       = '',
    sleep         = 3,
    machine       = 'ws16',
    status_only   = 0,
    generate_only = 1,
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
    bandfac          = 1.5,    # Adjust as needed for the 1RDM orbital basis
    occupations      = 'smearing',
    smearing         = 'fermi-dirac',
    degauss          = 1e-3,
    mixing_beta      = 0.4,
    nspin            = 2,
    hubbard          = {'U':{'Fe-3d': 2.0}},
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

system = generate_physical_system(
    structure  = 'Fe_bcc.xsf',
    kgrid      = (1,1,1),
    kshift     = (0,0,0),
    net_spin   = 6,
    net_charge = 0,
    symm_kgrid = False,
    tiling     = (1,1,1),
    Fe         = 16,
    )

basepath = './'

scf = generate_pwscf(
    identifier        = 'scf',
    path              = basepath + 'scf',
    job               = job(cores=16, app='pw.x'),
    calculation       = 'scf',
    system            = system,
    kgrid             = (24,24,24),
    **dft_shared
    )

nscf = generate_pwscf(
    identifier        = 'nscf',
    path              = basepath + 'nscf',
    job               = job(cores=16, app='pw.x'),
    calculation       = 'nscf',
    system            = system,
    nosym             = True,
    nogamma           = True,
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
    warmupsteps          = 10,
    init_cycles          = 5,
    cycles               = 5,
    walkers_per_rank     = 256,
    samples              = 1e5,
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
    warmupsteps          = 10,
    init_cycles          = 0,
    cycles               = 10,
    walkers_per_rank     = 256,
    samples              = 1e5,
    blocks               = 100,
    minwalkers           = 0.4,
    corrections          = [],
    dependencies         = [(p2q,'orbitals'), (optJ12,'jastrow')],
    **qmc_shared
    )


#===== Spin density =====
from qmcpack_input import spindensity
sdens = spindensity(
    dr = (0.05, 0.05, 0.05),   # Bohr units
    #grid = (100,100,100),  # Alternative to dr, an NxNxN grid can be specified
    )

#===== Energy density =====
from qmcpack_input import generate_energydensity
edens = generate_energydensity(
    coord = 'cartesian',
    grid  = (100, 100, 100),
)

#===== Momentum distribution =====
from qmcpack_input import momentumdistribution
mom_dist = momentumdistribution(
    samples = 40,
    kmax    = 8.0,
)

#===== One body density matrix =====
from qmcpack_input import onebodydensitymatrices, sposet
nbnd = nscf.input.system.nbnd   # Total number of bands solved in KS-DFT for each spin channel (occupied + virtual)
# For magnetic systems such as iron, the spin-up and spin-down 1RDMs are distinct.
# Currently, the 1RDM estimator can only project onto a single specified basis.
# The estimator below uses the spin-up orbitals as the basis in the 1RDM for both spin channels.
# To get the proper spin-down 1RDM, a seperate QMC run with spin-down basis (spindataset = 1) is needed at the moment.
dm_est = onebodydensitymatrices(
    type          = 'OneBodyDensityMatrices',
    energy_matrix = False,
    check_overlap = False,
    integrator    = 'density',
    evaluator     = 'matrix',
    samples       = 10,
    basis         = sposet(type='bspline', size=nbnd, spindataset=0), # Uses the spin up channel (0) KS-DFT orbitals as the DM basis
)

# In production runs containing estimators, large numbers (100+) of blocks are preferred for robust statistical analysis.
# More steps per block than specified here are additionally preferred for computational efficiency due to the I/O cost of estimators.
# Finally, note that for large numbers of blocks, the stat.h5 file's disk footprint can get considerably large.
# This is especially true for dense grids, say spindensity with (300x300x300).
qmc = generate_qmcpack(
    #skip_submit          = True,
    identifier           = 'qmc',
    path                 = basepath + 'qmc',
    job                  = job(cores=16,threads=4,app='qmcpack'),
    system               = system,
    twistnum             = 0,
    estimators           = [sdens, edens, mom_dist, dm_est],  # Requested estimators. These will be run in all QMC series.
    calculations         = [
    vmc(
        walkers_per_rank = 256,
        warmupsteps      = 100,
        blocks           = 200, # Pure VMC estimators will be needed to extrapolate the mixed DMC estimators
        steps            = 50,
        timestep         = 0.3,
        ),
    dmc( # equilibration
        walkers_per_rank = 256,
        timestep         = 0.02,
        warmupsteps      = 20,
        blocks           = 20,
        steps            = 20,
        nonlocalmoves    = 'v3',
        ),
    dmc( # main dmc
        walkers_per_rank = 256,
        timestep         = 0.005,
        warmupsteps      = 100,
        blocks           = 200,
        steps            = 200,
        nonlocalmoves    = 'v3',
        ),
    ],
    dependencies         = [(p2q,'orbitals'), (optJ123,'jastrow')],
    **qmc_shared
    )

run_project()
