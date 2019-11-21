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
dia2_kpath       = get_kpath(structure=dia2_structure)
#get_band_tiling and get_kpath require "SeeK-path" python3 libraries

dia2 = generate_physical_system(
    structure    = dia2_structure,
    kgrid  = (4,4,4),
    kshift = (0,0,0),
    C            = 4,
    )

scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'diamond/scf',
    job          = job(nodes=1,app='pw.x', hours = 1),
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
#K-path of the standardized primitive cell
dia2_structure.clear_kpoints()
dia2_kpts = generate_physical_system(
    structure    = dia2_structure,
    kpoints      = dia2_kpath['explicit_kpoints_abs_inv_B'],
    C            = 4,
    )

band = generate_pwscf(
    identifier   = 'nscf',
    path         = 'diamond/band',
    job          = job(nodes=1,app='pw.x', hours = 1),
    input_type   = 'generic',
    calculation  = 'nscf',
    input_dft    = 'lda', 
    ecutwfc      = 200,
    nspin        = 2,   
    conv_thr     = 1e-8,
    nosym        = True,
    wf_collect   = True,
    system       = dia2_kpts,
    nbnd         = 8,      #a sensible nbnd value can be given 
    verbosity    = 'high', #verbosity must be set to high
    pseudos      = ['C.BFD.upf'], 
    dependencies = (scf, 'charge_density'),
    )

run_project()

if band.finished:
    from pwscf_analyzer import PwscfAnalyzer
    p = PwscfAnalyzer(band)
    p.analyze()
    p.plot_bandstructure()
    print("VBM: {0}".format(p.bands.vbm))
    print("CBM: {0}".format(p.bands.cbm))
#end if
