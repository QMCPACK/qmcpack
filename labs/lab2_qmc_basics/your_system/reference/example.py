#! /usr/bin/env python

# import project suite functions
from project import settings,Job,run_project,get_machine
from project import generate_physical_system
from project import generate_pwscf
from project import generate_pw2qmcpack
from project import generate_qmcpack,vmc,loop,linear,dmc

# project suite settings
settings(
    pseudo_dir    = './pseudopotentials',
    runs          = '',
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'vesta',
    account       = 'QMC_2014_training'
    )
 
# allow max of one job at a time (lab only)
vesta = get_machine('vesta')
vesta.queue_size = 1

# locations of pwscf, pw2qmcpack and qmcpack executables
pwscf      = '/soft/applications/qmcpack/DFT_Binaries/pw.x'
pw2qmcpack = '/soft/applications/qmcpack/DFT_Binaries/pw2qmcpack.x'
qmcpack    = '/soft/applications/qmcpack/build_XL_real/bin/qmcapp'


# details of your physical system (diamond conventional cell below)
my_project_name = 'diamond_vmc'   # directory to perform runs
my_dft_pps      = ['C.BFD.upf']   # pwscf pseudopotentials
my_qmc_pps      = ['C.BFD.xml']   # qmcpack pseudopotentials

#  generate your system
#    units      :  'A'/'B' for Angstrom/Bohr
#    axes       :  simulation cell axes in cartesian coordinates (a1,a2,a3)
#    elem       :  list of atoms in the system
#    pos        :  corresponding atomic positions in cartesian coordinates
#    kgrid      :  Monkhorst-Pack grid
#    kshift     :  Monkhorst-Pack shift (between 0 and 0.5)
#    net_charge :  system charge in units of e
#    net_spin   :  # of up spins - # of down spins
#    C = 4      :  (pseudo) carbon has 4 valence electrons
my_system = generate_physical_system(
    units      = 'A',
    axes       = [[ 3.57000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [ 0.00000000e+00, 3.57000000e+00, 0.00000000e+00],
                  [ 0.00000000e+00, 0.00000000e+00, 3.57000000e+00]],
    elem       = ['C','C','C','C','C','C','C','C'],
    pos        = [[ 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [ 8.92500000e-01, 8.92500000e-01, 8.92500000e-01],
                  [ 0.00000000e+00, 1.78500000e+00, 1.78500000e+00],
                  [ 8.92500000e-01, 2.67750000e+00, 2.67750000e+00],
                  [ 1.78500000e+00, 0.00000000e+00, 1.78500000e+00],
                  [ 2.67750000e+00, 8.92500000e-01, 2.67750000e+00],
                  [ 1.78500000e+00, 1.78500000e+00, 0.00000000e+00],
                  [ 2.67750000e+00, 2.67750000e+00, 8.92500000e-01]],
    kgrid      = (1,1,1),
    kshift     = (0,0,0),
    net_charge = 0,
    net_spin   = 0,
    C          = 4       # one line like this for each atomic species
    )

my_bconds       = 'ppp'  #  ppp/nnn for periodic/open BC's in QMC
                         #  if nnn, center atoms about (a1+a2+a3)/2


sims = []

# scf run to generate orbitals
scf = generate_pwscf(
    identifier   = 'scf',
    path         = my_project_name,
    job          = Job(nodes=32,hours=2,app=pwscf),
    input_type   = 'scf',
    system       = my_system,
    pseudos      = my_dft_pps,
    input_dft    = 'lda', 
    ecut         = 200,   # PW energy cutoff in Ry
    conv_thr     = 1e-8, 
    mixing_beta  = .7,
    nosym        = True,
    wf_collect   = True
    )
sims.append(scf)

# conversion step to create h5 file with orbitals
p2q = generate_pw2qmcpack(
    identifier   = 'p2q',
    path         = my_project_name,
    job          = Job(cores=1,hours=2,app=pw2qmcpack),
    write_psir   = False,
    dependencies = (scf,'orbitals')
    )
sims.append(p2q)

# vmc run with qmcpack
qmc = generate_qmcpack(
    identifier   = 'vmc',
    path         = my_project_name,
    job          = Job(nodes=32,hours=2,threads=16,app=qmcpack),
    input_type   = 'basic',
    system       = my_system,
    bconds       = my_bconds,
    pseudos      = my_qmc_pps,
    jastrows     = [('J1','bspline',8,4.0),
                    ('J2','bspline',8    )], 
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
    dependencies = (p2q,'orbitals')
    )
sims.append(qmc)

# write input files and submit jobs
run_project(sims)




## additional runs for optimization and DMC
##
##   to use, comment out "run_project" above and uncomment below
#linopt1 = linear(
#    energy               = 0.0,
#    unreweightedvariance = 1.0,
#    reweightedvariance   = 0.0,
#    timestep             = 0.4,
#    samples              = 5000,
#    warmupsteps          = 50,
#    blocks               = 200,
#    substeps             = 1,
#    nonlocalpp           = True,
#    usebuffer            = True,
#    walkers              = 1,
#    minwalkers           = 0.5,
#    maxweight            = 1e9, 
#    usedrift             = True,
#    minmethod            = 'quartic',
#    beta                 = 0.025,
#    exp0                 = -16,
#    bigchange            = 15.0,
#    alloweddifference    = 1e-4,
#    stepsize             = 0.2,
#    stabilizerscale      = 1.0,
#    nstabilizers         = 3
#    )
#
#linopt2 = linopt1.copy()
#linopt2.samples = 20000
#
#opt = generate_qmcpack(
#    identifier   = 'opt',
#    path         = my_project_name,
#    job          = Job(nodes=32,hours=2,threads=16,app=qmcpack),
#    input_type   = 'basic',
#    system       = my_system,
#    bconds       = my_bconds,
#    pseudos      = my_qmc_pps,
#    jastrows     = [('J1','bspline',8,5.0),
#                    ('J2','bspline',8    )],
#    calculations = [loop(max=8,qmc=linopt1),
#                    loop(max=4,qmc=linopt2)],
#    dependencies = (p2q,'orbitals')
#    )
#sims.append(opt)
#
#qmc = generate_qmcpack(
#    identifier   = 'qmc',
#    path         = my_project_name,
#    job          = Job(nodes=32,hours=2,threads=16,app=qmcpack),
#    input_type   = 'basic',
#    system       = my_system,
#    bconds       = my_bconds,
#    pseudos      = my_qmc_pps,
#    jastrows     = [],            
#    calculations = [
#        vmc(
#            walkers     =   1,
#            warmupsteps =  30,
#            blocks      =  20,
#            steps       =  10,
#            substeps    =   2,
#            timestep    =  .4,
#            samples     = 2048
#            ),
#        dmc(
#            warmupsteps   =  20, 
#            blocks        = 200,
#            steps         =  10,
#            timestep      = 0.01,
#            nonlocalmoves = True
#            )
#        ],
#    dependencies = [(p2q,'orbitals'),(opt,'jastrow')]
#    )
#sims.append(qmc)
#
#run_project(sims)

