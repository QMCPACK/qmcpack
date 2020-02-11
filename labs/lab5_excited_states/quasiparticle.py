#! /usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,loop,linear,vmc


settings(
    pseudo_dir    = './pseudopotentials',
    results       = '',
    runs          = 'diamond_runs',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16'
    )

#Input structure
dia_inputs = dict(
    units     = 'A',
    axes      = [[ 1.785,  1.785,  0.   ],
                 [ 0.   ,  1.785,  1.785],
                 [ 1.785,  0.   ,  1.785]],
    elem      = ['C','C'],
    pos       = [[ 0.    ,  0.    ,  0.    ],
                 [ 0.8925,  0.8925,  0.8925]],
    use_prim  = True,    # Use SeeK-path library to identify prim cell
    tiling    = [3,1,3], # Tile the cell
    kgrid     = (1,1,1), 
    kshift    = (0,0,0), # Assumes we study transitions from Gamma. For non-gamma tilings, use kshift appropriately
    C         = 4,
    )

dia = generate_physical_system(
    **dia_inputs
    )

scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'scf',
    job          = job(cores=16,app='pw.x'),
    input_type   = 'generic',
    calculation  = 'scf',
    nspin        = 2,
    input_dft    = 'lda', 
    ecutwfc      = 200,   
    conv_thr     = 1e-8, 
    nosym        = True, # Important in this case.
                         # With nosym=False a 1eV larger gap is found.
                         # Issue does not appear in larger, e.g. 3x3x3 cells.
    wf_collect   = False,
    system       = dia,
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    tot_magnetization = 0,
    pseudos      = ['C.BFD.upf'], 
    )

nscf = generate_pwscf(
    identifier   = 'nscf',
    path         = 'nscf',
    job          = job(cores=16,app='pw.x'),
    input_type   = 'generic',
    calculation  = 'nscf',
    input_dft    = 'lda', 
    ecutwfc      = 200,
    nspin        = 2,   
    conv_thr     = 1e-8,
    nosym        = True,
    nosym_evc    = True,
    noinv        = True,
    wf_collect   = True,
    system       = dia,
    nbnd         = 8,      #a sensible nbnd value can be given 
    verbosity    = 'high', #verbosity must be set to high
    pseudos      = ['C.BFD.upf'], 
    dependencies = (scf,'charge_density'),
    )


conv = generate_pw2qmcpack(
    identifier   = 'conv',
    path         = 'nscf',
    job          = job(cores=1,app='pw2qmcpack.x'),
    write_psir   = False,
    dependencies = (nscf,'orbitals'),
    )

opt = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'opt',
    path           = 'opt',
    job            = job(cores=16,threads=16,app='qmcpack_complex'),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    pseudos        = ['C.BFD.xml'],
    jastrows       = [('J1','bspline',8),('J2','bspline',8)],
    calculations   = [
        loop(max=4,
             qmc = linear(
                 minmethod   = 'OneShiftOnly',
                 minwalkers  = 0.5,
                 samples     = 5000,
                 walkers     = 1,
                 warmupsteps =  20,
                 blocks      =  25,
                 substeps    =   5,
                 timestep    =  .4,
                 ),
             ),
        ],
    dependencies = (conv,'orbitals'),
    )

qmc_ground = generate_qmcpack(
    det_format     = 'old',
    identifier     = 'vmc',
    path           = 'vmc_ground',
    job            = job(cores=16,threads=16,app='qmcpack_complex'),
    input_type     = 'basic',
    spin_polarized = True,
    system         = dia,
    pseudos        = ['C.BFD.xml'],
    jastrows       = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 800,
            steps       =   5,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )

qmc_minus = generate_qmcpack(
    det_format   = 'old',
    identifier   = 'vmc',
    path         = 'vmc_-e',
    job          = job(cores=16,threads=16,app='qmcpack_complex'),
    input_type   = 'basic',
    spin_polarized = True,
    system       = dia,
    pseudos      = ['C.BFD.xml'],
    jastrows     = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 800,
            steps       =   5,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )
u = qmc_minus.input.get('u')
u.size-=1
updet = qmc_minus.input.get('updet')
updet.size-=1

qmc_plus = generate_qmcpack(
    det_format   = 'old',
    identifier   = 'vmc',
    path         = 'vmc_+e',
    job          = job(cores=16,threads=16,app='qmcpack_complex'),
    input_type   = 'basic',
    spin_polarized = True,
    system       = dia,
    pseudos      = ['C.BFD.xml'],
    jastrows     = [],
    calculations   = [
        vmc(
            warmupsteps =  20,
            blocks      = 800,
            steps       =   5,
            substeps    =   2,
            timestep    =  .4,
            )
        ],
    dependencies = [(conv,'orbitals'),
                    (opt,'jastrow')],
    )
u = qmc_plus.input.get('u')
u.size+=1
updet = qmc_plus.input.get('updet')
updet.size+=1


run_project()
