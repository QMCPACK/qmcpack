#! /usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc
from structure import *

settings(
    pseudo_dir    = '../../pseudopotentials',
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

#Standardized Primitive cell -- run rest of the calculations on this cell
dia2_structure   = get_primitive_cell(structure=dia.structure)['structure']
#get_band_tiling requires "SeeK-path" python3 libraries

dia2 = generate_physical_system(
    structure    = dia2_structure,
    kgrid  = (1,1,1), 
    kshift = (0,0,0), # Assumes we study transitions from Gamma. For non-gamma tilings, use kshift appropriately
    tiling = [3,1,3],
    C            = 4,
    )

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
    nosym        = True,
    wf_collect   = True,
    system       = dia2,
    tot_magnetization = 0,
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
    system       = dia2,
    nbnd         = 8,      #a sensible nbnd value can be given 
    verbosity    = 'high', #verbosity must be set to high
    pseudos      = ['C.BFD.upf'], 
    dependencies = (scf, 'charge_density'),
)

conv = generate_pw2qmcpack(
    identifier   = 'conv',
    path         = 'diamond/nscf',
    job          = job(cores=1,app='pw2qmcpack.x', hours = 1),
    write_psir   = False,
    dependencies = (nscf,'orbitals'),
    )

qmc = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia2,
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            walkers     =  16,
            warmupsteps =  20,
            blocks      = 1000,
            steps       =  10,
            substeps    =   2,
            timestep    =  .4
            )
        ],
    dependencies = (conv,'orbitals'),
    )

qmc_optical = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'diamond/vmc_optical',
    job            = job(cores=16,threads=16,app='qmcpack', hours = 1),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia2,
    excitation     = ['up', '0 3 4 4'], #
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            walkers     =  16,
            warmupsteps =  20,
            blocks      = 1000,
            steps       =  10,
            substeps    =   2,
            timestep    =  .4
            )
        ],
    dependencies = (conv,'orbitals'),
    )

run_project()
