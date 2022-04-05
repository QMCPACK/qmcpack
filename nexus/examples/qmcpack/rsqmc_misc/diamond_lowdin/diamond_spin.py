#! /usr/bin/env python

from nexus import settings
from nexus import job
from nexus import run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_projwfc
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack
from nexus import vmc

from structure import *

from qmcpack_input import dm1b
from qmcpack_input import sposet

settings(
    pseudo_dir    = '../../pseudopotentials',
    runs          = 'runs_spin',
    results       = '',
    status_only   = 0,
    generate_only = 0,
    skip_submit   = 0,
    sleep         = 3,
    machine       = 'ws4'
    )

dia16 = generate_physical_system(
    units  = 'A',
    axes   = [[ 1.785,  1.785,  0.   ],
              [ 0.   ,  1.785,  1.785],
              [ 1.785,  0.   ,  1.785]],
    elem   = ['C','C'],
    pos    = [[ 0.    ,  0.    ,  0.    ],
              [ 0.8925,  0.8925,  0.8925]],
    tiling = (1,1,1),
    C      = 4
    )
              
# k-mesh used for density
scf_kg = dia16.structure.kgrid_from_kspacing(0.5) # Get SCF kmesh from k-spacing

# twist-mesh used for qmc
dia16.structure.add_symmetrized_kmesh(kgrid=(2,2,2),kshift=(0,0,0))


number_of_ks_orbs = 11

scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'scf',
    job          = job(cores=1,app='pw.x',hours=1),
    input_type   = 'generic',
    calculation  = 'scf',
    nspin        = 2,
    tot_magnetization = 0,
    nbnd         = number_of_ks_orbs,
    input_dft    = 'lda',
    ecutwfc      = 200,
    conv_thr     = 1e-8,
    nosym        = False,
    wf_collect   = False,
    system       = dia16,
    kgrid        = scf_kg,
    kshift       = (0,0,0),
    pseudos      = ['C.BFD.upf'],
    )

nscf = generate_pwscf(
    identifier   = 'nscf',
    path         = 'nscf',
    job          = job(cores=1,app='pw.x',hours=1),
    input_type   = 'generic',
    calculation  = 'nscf',
    input_dft    = 'lda',
    ecutwfc      = 200,
    nspin        = 2,
    tot_magnetization = 0,
    conv_thr     = 1e-8,
    nosym        = True,
    wf_collect   = True,
    system       = dia16,
    nbnd         = number_of_ks_orbs,
    verbosity    = 'high', #verbosity must be set to high
    pseudos      = ['C.BFD.upf'],
    dependencies = (scf,'charge_density'),
    )

# To obtain the overlaps between the Bloch states and atomic orbitals,
# projwfc.x needs to be run. The overlaps will be stored in:
# pwscf_output/pwscf.save/atomic_proj.xml
# WARNING: Always check the the <OVERLAPS> element is written to atomic_proj.xml
#          Sometimes QE will not write <OVERLAPS> if running on >1 core.
pwf = generate_projwfc(
    identifier      = 'pwf',
    path            = 'nscf',
    job             = job(cores=1,app='projwfc.x',hours=1),
    lwrite_overlaps = True,
    lsym            = False,
    dependencies    = (nscf,'other')
    )

# Generate orbital h5 file
conv = generate_pw2qmcpack(
    identifier   = 'conv',
    path         = 'nscf',
    job          = job(cores=1,app='pw2qmcpack.x',hours=1),
    write_psir   = False,
    dependencies = (nscf,'orbitals'),
    )

# Define 1RDM Parameters
dm_estimator = dm1b(
        energy_matrix = False,
        integrator    = 'uniform_grid',
        points        = 6,
        scale         = 1.0,
        basis         = sposet(type='bspline',size=number_of_ks_orbs,spindataset=0),
        evaluator     = 'matrix',
        center        = (0,0,0),
        check_overlap = False,
        )

down_dm_estimator = dm1b(
        energy_matrix = False,
        integrator    = 'uniform_grid',
        points        = 6,
        scale         = 1.0,
        basis         = sposet(type='bspline',size=number_of_ks_orbs,spindataset=1),
        evaluator     = 'matrix',
        center        = (0,0,0),
        check_overlap = False,
        )


qmc = generate_qmcpack(
    identifier   = 'vmc_1rdm_noJ',
    path         = 'vmc_1rdm_noJ',
    job          = job(cores=3,app='qmcpack_complex',hours=1),
    input_type   = 'basic',
    system       = dia16,
    pseudos      = ['C.BFD.xml'],
    estimators   = [dm_estimator],
    jastrows     = [],
    calculations = [
        vmc(
            walkers     =   1,
            warmupsteps =  20,
            blocks      = 200,
            steps       =  10,
            substeps    =   2,
            timestep    =  .4
            )
        ],
    dependencies = (conv,'orbitals'),
    )

qmc = generate_qmcpack(
    identifier   = 'vmc_1rdm_down_noJ',
    path         = 'vmc_1rdm_down_noJ',
    job          = job(cores=3,app='qmcpack_complex',hours=1),
    input_type   = 'basic',
    system       = dia16,
    pseudos      = ['C.BFD.xml'],
    estimators   = [down_dm_estimator],
    jastrows     = [],
    calculations = [
        vmc(
            walkers     =   1,
            warmupsteps =  20,
            blocks      = 200,
            steps       =  10,
            substeps    =   2,
            timestep    =  .4
            )
        ],
    dependencies = (conv,'orbitals'),
    )
run_project()

