#! /usr/bin/env python3

from nexus import *

settings(
    pseudo_dir    = '../../pseudopotentials',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16'
    )

relax_job = job(cores=16,app='pw.x')
scf_job   = job(cores=16,app='pw.x')

dia16 = generate_physical_system(
    structure = './d16vac.POSCAR',
    C         = 4
    )

relax = generate_pwscf(      
    identifier   = 'relax', 
    path         = 'diamond_vacancy/relax',
    job          = relax_job,
    input_type   = 'generic',
    calculation  = 'relax',
    ion_dynamics = 'bfgs',
    input_dft    = 'lda',        
    ecutwfc      = 35,
    conv_thr     = 1e-6, 
    system       = dia16,            
    pseudos      = ['C.BFD.upf'], 
    kgrid        = (2,2,2),                
    kshift       = (0,0,0),              
    )
              
scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'diamond_vacancy/scf',
    job          = scf_job,
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda', 
    ecutwfc      = 75,   
    conv_thr     = 1e-7, 
    system       = dia16,
    pseudos      = ['C.BFD.upf'], 
    kgrid        = (2,2,2),                
    kshift       = (0,0,0),              
    dependencies = (relax,'structure')
    )

run_project()
