#! /usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc
from structure import *

'''
This nexus example shows a variety of ways that excitations can be specified.
'''

settings(
    pseudo_dir    = '../../pseudopotentials',
    runs          = './runs_excitation_alternatives'
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16'
    )

#Input structure
dia = generate_physical_system(
    units  = 'A',
    axes   = [[ 1.785,  1.785,  0.   ],
              [ 0.   ,  1.785,  1.785],
              [ 1.785,  0.   ,  1.785]],
    elem   = ['C','C'],
    pos    = [[ 0.    ,  0.    ,  0.    ],
              [ 0.8925,  0.8925,  0.8925]],
    C      = 4
    )

kg = dia.structure.kgrid_from_kspacing(0.3) # Get SCF kmesh from k-spacing

scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'diamond/scf',
    job          = job(nodes=1,app='pw.x',hours=1),
    input_type   = 'generic',
    calculation  = 'scf',
    nspin        = 2,
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    nosym        = False,
    wf_collect   = False,
    system       = dia,
    tot_magnetization = 0,
    kgrid        = kg,
    kshift       = (0,0,0),
    pseudos      = ['C.BFD.upf'], 
    )

nscf = generate_pwscf(
    identifier   = 'nscf',
    path         = 'diamond/nscf',
    job          = job(nodes=1,app='pw.x',hours=1),
    input_type   = 'generic',
    calculation  = 'nscf',
    input_dft    = 'lda',
    ecutwfc      = 200,
    nspin        = 2,
    conv_thr     = 1e-8,
    nosym        = True,
    wf_collect   = True,
    system       = dia,
    nbnd         = 8,      #a sensible nbnd value can be given
    verbosity    = 'high', #verbosity must be set to high
    pseudos      = ['C.BFD.upf'],
    dependencies = (scf,'charge_density'),
    )

conv = generate_pw2qmcpack(
    identifier   = 'conv',
    path         = 'diamond/nscf',
    job          = job(cores=1,app='pw2qmcpack.x', hours = 1),
    write_psir   = False,
    dependencies = (nscf,'orbitals'),
    )

opt = generate_qmcpack(
    identifier      = 'opt',
    path            = 'diamond/opt',
    job             = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type      = 'basic',
    system          = dia,
    pseudos         = ['C.BFD.xml'],
    twistnum        = 0,
    J2              = True,         # Add a 2-body B-spline Jastrow
    spin_polarized  = True,
    qmc             = 'opt',        # Do a wavefunction optimization
    minmethod       = 'oneshift',   # Optimization algorithm (assumes energy minimization)
    init_cycles     = 4,            # First 4 iterations allow large parameter changes
    cycles          = 10,           # 8 subsequent iterations with smaller parameter changes
    warmupsteps     = 8,            # First 8 steps are not recorded
    blocks          = 100,          # Number of blocks to write in the .scalar.dat file
    timestep        = 0.4,          # MC step size (nothing to do with time for VMC)
    init_minwalkers = 0.01,         # Smaller values -> bigger parameter change
    minwalkers      = 0.5,          #
    samples         = 5000,         # VMC samples per iteration
    use_nonlocalpp_deriv = False,
    dependencies    = (conv,'orbitals'),
    )

################################################################################
############ Ground State at Gamma #############################################
################################################################################
qmc_ground = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_ground',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    twistnum        = 0,
    system         = dia,
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

################################################################################
############ Single Determinant Excitations ####################################
################################################################################

# In each of the following 4 examples, an optical excitation is performed in the up-channel
# corresponding to the homo-lumo gap at the gamma k-point. All 4 examples lead to the same 
# excitation, but show the various ways that the excitation can be specfified

# up channel, gamma vb gamma cb
qmc_optical = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical_up_g-vb-g-cb',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    twistnum        = 0,
    excitation     = ['up', 'gamma vb gamma cb'],
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

# up channel, band index 
qmc_optical = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical_up_band-index',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    twistnum        = 0,
    excitation     = ['up', '0 3 0 4'],
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

# up channel, energy index
qmc_optical = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical_up_energy-index',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    twistnum        = 0,
    excitation     = ['up', '-4 +5'],
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

# up channel, lowest index
qmc_optical = generate_qmcpack(
    skip_submit    = 0,
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical_up_lowest',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    twistnum        = 0,
    excitation     = ['up', 'lowest'],
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

################################################################################
############ Triplet Excitations ###############################################
################################################################################

# In each of the following 2 examples, an optical excitation is performed for a triplet state
# corresponding to the homo-lumo gap at the gamma k-point. Both examples lead to the same 
# excitation, but show the various ways that the excitation can be specfified

# triplet, energy index
qmc_optical = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical_triplet_energy-index',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    twistnum        = 0,
    excitation     = ['triplet', '-4 +5'],
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

# triplet, lowest
qmc_optical = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical_triplet_lowest',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    twistnum        = 0,
    excitation     = ['triplet', 'lowest'],
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

################################################################################
############ Singlet Excitations ###############################################
################################################################################

# In each of the following 2 examples, an optical excitation is performed for a singlet state
# corresponding to the homo-lumo gap at the gamma k-point. Both examples lead to the same 
# excitation, but show the various ways that the excitation can be specfified

# singlet, energy index
qmc_optical = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical_singlet_energy-index',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    twistnum        = 0,
    excitation     = ['singlet', '-4 +5'],
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

# singlet, lowest
qmc_optical = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical_singlet_lowest',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    twistnum        = 0,
    excitation     = ['singlet', 'lowest'],
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 2400,
            steps       =  25,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

run_project()
