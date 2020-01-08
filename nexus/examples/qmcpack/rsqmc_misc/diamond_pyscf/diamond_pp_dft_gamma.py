#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
from nexus import generate_convert4qmc
from nexus import generate_qmcpack

settings(
    pseudo_dir = '../../pseudopotentials',
    results    = '',
    sleep      = 3,
    machine    = 'ws16',
    )

system = generate_physical_system(
    units    = 'A',
    axes     = '''1.785   1.785   0.000
                  0.000   1.785   1.785
                  1.785   0.000   1.785''',
    elem_pos = '''
               C  0.0000  0.0000  0.0000
               C  0.8925  0.8925  0.8925
               ''',
    kgrid    = (1,1,1),
    kshift   = (0,0,0),
    C        = 4,
    )

scf = generate_pyscf(
    identifier = 'scf',                      # log output goes to scf.out
    path       = 'diamond_pp_dft_gamma',     # directory to run in
    job        = job(serial=True,threads=16),# pyscf must run w/o mpi
    template   = './dft_template.py',        # pyscf template file
    system     = system,
    cell       = obj(                        # used to make Cell() inputs
        basis         = 'bfd-vdz',
        ecp           = 'bfd',
        drop_exponent = 0.1,
        verbose       = 5,
        ),
    save_qmc   = True,                # save wfn data for qmcpack
    )

c4q = generate_convert4qmc(
    identifier   = 'c4q',
    path         = 'diamond_pp_dft_gamma',
    job          = job(cores=1),
    no_jastrow   = True,
    dependencies = (scf,'orbitals'),
    )

qmc = generate_qmcpack(
    identifier   = 'vmc',
    path         = 'diamond_pp_dft_gamma',
    job          = job(cores=16),
    system       = system,
    pseudos      = ['C.BFD.xml'],
    input_type   = 'basic',
    jastrows     = [],
    corrections  = [],
    qmc          = 'vmc',
    dependencies = (c4q,'orbitals'),
    )

run_project()
