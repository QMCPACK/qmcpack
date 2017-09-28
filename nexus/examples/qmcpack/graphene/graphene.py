#! /usr/bin/env python

from nexus import settings,Job,run_project
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack
from nexus import loop,linear,vmc,dmc


# general settings for nexus
settings(
    pseudo_dir    = '../pseudopotentials',# directory with all pseudopotentials
    sleep         = 3,                    # check on runs every 'sleep' seconds
    #generate_only   = False,
    # Complicated setting only so examples can be run in test harness.
    # For real runs, use the plain setting of 'generate_only' above.
    generate_only   = globals().get('override_generate_only_setting',False),
    status_only   = 0,                    # only show status of runs
    machine       = 'ws16',               # local machine is 16 core workstation
    )



# generate the graphene physical system
graphene = generate_physical_system(
    lattice   = 'hexagonal',      # hexagonal cell shape
    cell      = 'primitive',      # primitive cell
    centering = 'P',              # primitive basis centering
    constants = (2.462,10.0),     # a,c constants
    units     = 'A',              # in Angstrom
    atoms     = ('C','C'),        # C primitive atoms
    basis     = [[ 0  , 0  , 0],  # basis vectors
                 [2./3,1./3, 0]],
    tiling    = (2,2,1),          # tiling of primitive cell
    kgrid     = (1,1,1),          # Monkhorst-Pack grid
    kshift    = (.5,.5,.5),       # and shift
    C         = 4                 # C has 4 valence electrons
    ) 


# list of simulations in workflow
sims = []

# scf run produces charge density
scf = generate_pwscf(
    # nexus inputs
    identifier   = 'scf',           # identifier/file prefix
    path         = 'graphene/scf',  # directory for scf run
    job          = Job(cores=16),   # run on 16 cores
    pseudos      = ['C.BFD.upf'],   # pwscf PP file
    system       = graphene,        # run graphene
    # input format selector
    input_type   = 'scf',           # scf, nscf, relax, or generic
    # pwscf input parameters
    input_dft    = 'lda',           # dft functional
    ecut         =  150 ,           # planewave energy cutoff (Ry)
    conv_thr     =  1e-6,           # scf convergence threshold (Ry)
    mixing_beta  =    .7,           # charge mixing factor
    kgrid        = (8,8,8),         # MP grid of primitive cell
    kshift       = (1,1,1),         #  to converge charge density
    wf_collect   = False,           # don't collect orbitals
    use_folded   = True             # use primitive rep of graphene
    )
sims.append(scf)  

# nscf run to produce orbitals for jastrow optimization
nscf_opt = generate_pwscf(
    # nexus inputs
    identifier   = 'nscf',          # identifier/file prefix      
    path         = 'graphene/nscf_opt', # directory for nscf run       
    job          = Job(cores=16),   # run on 16 cores             
    pseudos      = ['C.BFD.upf'],   # pwscf PP file               
    system       = graphene,        # run graphene                
    # input format selector                                      
    input_type   = 'nscf',          # scf, nscf, relax, or generic
    # pwscf input parameters
    input_dft    = 'lda',           # dft functional
    ecut         =  150 ,           # planewave energy cutoff (Ry)
    conv_thr     =  1e-6,           # scf convergence threshold (Ry)
    mixing_beta  =    .7,           # charge mixing factor
    nosym        = True,            # don't symmetrize k-points
    use_folded   = True,            # use primitive rep of graphene
    wf_collect   = True,            # write out orbitals
    kgrid        = (1,1,1),         # single k-point for opt
    kshift       = (0,0,0),         # gamma point
    # workflow dependencies
    dependencies = (scf,'charge_density')
    )
sims.append(nscf_opt)

# orbital conversion job for jastrow optimization
p2q_opt = generate_pw2qmcpack(
    # nexus inputs
    identifier   = 'p2q',
    path         = 'graphene/nscf_opt',
    job          = Job(cores=1),
    # pw2qmcpack input parameters
    write_psir   = False,
    # workflow dependencies
    dependencies = (nscf_opt,'orbitals')
    )
sims.append(p2q_opt)

# Jastrow optimization
opt = generate_qmcpack(
    # nexus inputs
    identifier   = 'opt',           # identifier/file prefix
    path         = 'graphene/opt',  # directory for opt run
    job          = Job(cores=16,app='qmcapp'),
    pseudos      = ['C.BFD.xml'],   # qmcpack PP file
    system       = graphene,        # run graphene
    # input format selector   
    input_type   = 'basic',
    # qmcpack input parameters
    corrections  = [], 
    jastrows     = [('J1','bspline',8),   # 1 body bspline jastrow
                    ('J2','bspline',8)],  # 2 body bspline jastrow
    calculations = [
        loop(max = 6,                        # No. of loop iterations
             qmc = linear(                   # linearized optimization method
                energy               =  0.0, # cost function
                unreweightedvariance =  1.0, #   is all unreweighted variance
                reweightedvariance   =  0.0, #   no energy or r.w. var. 
                timestep             =  0.5, # vmc timestep (1/Ha)
                warmupsteps          =  100, # MC steps before data collected 
                samples              = 16000,# samples used for cost function 
                stepsbetweensamples  =   10, # steps between uncorr. samples
                blocks               =   10, # ignore this  
                minwalkers           =   0.1,#  and this
                bigchange            =  15.0,#  and this
                alloweddifference    =  1e-4 #  and this, for now
                )
             )        
        ],
    # workflow dependencies
    dependencies = (p2q_opt,'orbitals')        
    )
sims.append(opt)


# nscf run to produce orbitals for final dmc
nscf = generate_pwscf(
    # nexus inputs
    identifier   = 'nscf',          # identifier/file prefix      
    path         = 'graphene/nscf', # directory for nscf run       
    job          = Job(cores=16),   # run on 16 cores             
    pseudos      = ['C.BFD.upf'],   # pwscf PP file               
    system       = graphene,        # run graphene                
    # input format selector                                      
    input_type   = 'nscf',          # scf, nscf, relax, or generic
    # pwscf input parameters
    input_dft    = 'lda',           # dft functional
    ecut         =  150 ,           # planewave energy cutoff (Ry)
    conv_thr     =  1e-6,           # scf convergence threshold (Ry)
    mixing_beta  =    .7,           # charge mixing factor
    nosym        = True,            # don't symmetrize k-points
    use_folded   = True,            # use primitive rep of graphene
    wf_collect   = True,            # write out orbitals
    # workflow dependencies
    dependencies = (scf,'charge_density')
    )
sims.append(nscf)

# orbital conversion job for final dmc
p2q = generate_pw2qmcpack(
    # nexus inputs
    identifier   = 'p2q',
    path         = 'graphene/nscf',
    job          = Job(cores=1),
    # pw2qmcpack input parameters
    write_psir   = False,
    # workflow dependencies
    dependencies = (nscf,'orbitals')
    )
sims.append(p2q)
    
# final dmc run
qmc = generate_qmcpack( 
    # nexus inputs
    identifier   = 'qmc',           # identifier/file prefix       
    path         = 'graphene/qmc',  # directory for dmc run       
    job          = Job(cores=16,app='qmcapp'),
    pseudos      = ['C.BFD.xml'],   # qmcpack PP file
    system       = graphene,        # run graphene
    # input format selector                                      
    input_type   = 'basic',
    # qmcpack input parameters
    corrections  = [],              # no finite size corrections
    jastrows     = [],              # overwritten from opt
    calculations = [                # qmcpack input parameters for qmc
        vmc(                        # vmc parameters 
            timestep      = 0.5,    # vmc timestep (1/Ha)
            warmupsteps   = 100,    # No. of MC steps before data is collected
            blocks        = 200,    # No. of data blocks recorded in scalar.dat
            steps         =  10,    # No. of steps per block
            substeps      =   3,    # MC steps taken w/o computing E_local
            samplesperthread = 40   # No. of dmc walkers per thread
            ),                      
        dmc(                        # dmc parameters
            timestep      = 0.01,   # dmc timestep (1/Ha)
            warmupsteps   =  50,    # No. of MC steps before data is collected
            blocks        = 400,    # No. of data blocks recorded in scalar.dat
            steps         =   5,    # No. of steps per block
            nonlocalmoves = True    # use Casula's T-moves
            ),                      #  (retains variational principle for NLPP's)
        ],
    # workflow dependencies
    dependencies = [(p2q,'orbitals'),
                    (opt,'jastrow')]
    )



# nexus monitors all runs
run_project(sims)



# print out the total energy
performed_runs = not settings.generate_only and not settings.status_only
if performed_runs:
    # get the qmcpack analyzer object
    # it contains all of the statistically analyzed data from the run
    qa = qmc.load_analyzer_image()
    # get the local energy from dmc.dat
    le = qa.dmc[1].dmc.LocalEnergy  # dmc series 1, dmc.dat, local energy
    #  print the total energy for the 8 atom system
    print 'The DMC ground state energy for graphene is:'
    print '    {0} +/- {1} Ha'.format(le.mean,le.error)
#end if



