#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_gamess

settings(
    pseudo_dir    = './pseudopotentials',
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16',
    ericfmt       = '/your/path/to/ericfmt.dat'
    )

gms_job = job(cores=16,app='gamess.x')

h2o = generate_physical_system(
    # full atomic/electronic structure
    elem        = ['O','H','H'], 
    pos         = [[0.000000, 0.000000, 0.000000],
                   [0.000000,-0.757160, 0.586260],
                   [0.000000, 0.757160, 0.586260]],
    units       = 'A', # Angstroms
    net_spin    = 0,   # multiplicity=1, nup-ndown=0
    O           = 6,   # Zeff=6 for BFD ECP
    H           = 1,   # Zeff=1 for BFD ECP
    # C2v symmetry structure
    folded_elem = ['O','H'],     
    folded_pos  = [[0.000000, 0.000000, 0.000000],
                   [0.000000, 0.757160, 0.586260]],
    )

rhf = generate_gamess(
    identifier = 'rhf',
    path       = 'pp_cisd',
    job        = gms_job,
    system     = h2o,
    pseudos    = ['O.BFD_V5Z.gms','H.BFD_V5Z_ANO.gms'],
    scftyp     = 'rohf',
    runtyp     = 'energy',
    exetyp     = 'run',
    ispher     = 1,
    maxit      = 200,
    memory     = 150000000,
    dirscf     = True,
    guess      = 'huckel',
    symmetry   = 'Cnv 2',
    prtmo      = True,
    )

cisd = generate_gamess(
    identifier = 'cisd',
    path       = 'pp_cisd',
    job        = gms_job,
    system     = h2o,
    pseudos    = ['O.BFD_V5Z.gms','H.BFD_V5Z_ANO.gms'],
    scftyp     = 'none',
    cityp      = 'guga',
    runtyp     = 'energy',
    exetyp     = 'run',
    ispher     = 1,
    maxit      = 200,
    memory     = 150000000,
    dirscf     = True,
    symmetry   = 'Cnv 2',
    cidrt = obj(
        group  = 'C2v',
        nfzc   = 0,
        ndoc   = 4,
        nalp   = 0,
        nval   = 60,
        nprt   = 2,
        istsym = 1,
        iexcit = 2,
        mxnint = 500000,
        ),
    gugdia = obj(
        prttol = 0.001,
        cvgtol = 1.0e-5,
        itermx = 1000,
        ),
    dependencies = (rhf,'orbitals')
    )

run_project()
