#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf

settings(
    results = '',
    sleep   = 3,
    machine = 'ws16',
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

    # Use Hartree-Fock methods
    path       = 'diamond_pp_hf_gamma',      # directory to run in
    template   = './hf_template.py',         # pyscf hf template file

    # Alternative using DFT.
    # Uncomment the lines below and comment out the lines above to run.
    # path       = 'diamond_pp_dft_gamma',      # directory to run in
    # template   = './dft_template.py',         # pyscf dft template file

    job        = job(serial=True,threads=16),# pyscf must run w/o mpi
    system     = system,
    cell       = obj(                        # used to make Cell() inputs
        basis         = 'bfd-vdz',
        ecp           = 'bfd',
        drop_exponent = 0.1,
        verbose       = 5,
        ),
    )

run_project()
