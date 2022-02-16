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

from qmcpack_input import spindensity

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

# Generate orbital h5 file
conv = generate_pw2qmcpack(
    identifier   = 'conv',
    path         = 'nscf',
    job          = job(cores=1,app='pw2qmcpack.x',hours=1),
    write_psir   = False,
    dependencies = (nscf,'orbitals'),
    )


qmc = generate_qmcpack(
    identifier   = 'vmc',
    path         = 'vmc',
    job          = job(cores=3,app='qmcpack_complex',hours=1),
    input_type   = 'basic',
    system       = dia16,
    pseudos      = ['C.BFD.xml'],
    estimators   = [spindensity(dr=(0.5,0.5,0.5))],
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

