#! /usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack,vmc


settings(
    pseudo_dir    = './pseudopotentials',
    results       = '',
    runs          = 'diamond_runs',  
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16'
    )

dia = generate_physical_system(
    units     = 'A',
    axes      = [[ 1.785,  1.785,  0.   ],
                 [ 0.   ,  1.785,  1.785],
                 [ 1.785,  0.   ,  1.785]],
    elem      = ['C','C'],
    pos       = [[ 0.    ,  0.    ,  0.    ],
                 [ 0.8925,  0.8925,  0.8925]],
    use_prim  = True, # use SeeK-path library to identify prim cell
    add_kpath = True, # use SeeK-path library to get prim k-path
    C         = 4
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
    nosym        = True,
    wf_collect   = True,
    system       = dia,
    kgrid        = (4,4,4),
    kshift       = (0,0,0),
    tot_magnetization = 0,
    pseudos      = ['C.BFD.upf'], 
    )

band = generate_pwscf(
    identifier   = 'nscf',
    path         = 'band',
    job          = job(cores=16,app='pw.x'),
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
    dependencies = (scf, 'charge_density'),
    )

run_project()

p = band.load_analyzer_image()
p.plot_bandstructure()
print "VBM:\n{0}".format(p.bands.vbm) 
print "CBM:\n{0}".format(p.bands.cbm)

#print 'CBM, Delta'
#print p.bands.up[51]
#print 'Delta prime'
#print p.bands.up[46]
