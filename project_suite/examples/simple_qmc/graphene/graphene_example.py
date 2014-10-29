#! /usr/bin/env python

from nexus import settings,ProjectManager,Job
from nexus import generate_physical_system
from nexus import loop,linear,vmc,dmc
from qmcpack_calculations import standard_qmc


#general settings for nexus
settings(
    pseudo_dir    = './pseudopotentials',# directory with all pseudopotentials
    sleep         = 3,                   # check on runs every 'sleep' seconds
    generate_only = 0,                   # only make input files
    status_only   = 0,                   # only show status of runs
    machine       = 'node16',            # local machine is 16 core workstation
    )



#generate the graphene physical system
graphene = generate_physical_system(
    structure = 'graphite_aa',  # graphite keyword
    cell      = 'hex',          # hexagonal cell shape
    tiling    = (2,2,1),        # tiling of primitive cell
    constants = (2.462,10.0),   # a,c constants
    units     = 'A',            # in Angstrom
    kgrid     = (1,1,1),        # Monkhorst-Pack grid
    kshift    = (.5,.5,.5),     # and shift
    C         = 4               # C has 4 valence electrons
    ) 


#generate the simulations for the qmc workflow
qsims = standard_qmc(
    # subdirectory of runs 
    directory       = 'graphene_test',
    # description of the physical system
    system          = graphene,
    pseudos         = ['C.BFD.upf',  # pwscf PP file
                       'C.BFD.xml'], # qmcpack PP file
    # job parameters
    scfjob          = Job(cores=16), # cores to run scf 
    nscfjob         = Job(cores=16), # cores to run non-scf 
    optjob          = Job(cores=16), # cores for optimization 
    qmcjob          = Job(cores=16), # cores for qmc
    # dft parameters (pwscf)
    functional      = 'lda',         # dft functional
    ecut            =  150 ,         # planewave energy cutoff (Ry)
    conv_thr        =  1e-6,         # scf convergence threshold (Ry)
    mixing_beta     =    .7,         # charge mixing factor
    scf_kgrid       = (8,8,8),       # MP grid of primitive cell
    scf_kshift      = (1,1,1),       #  to converge charge density
    # qmc wavefunction parameters (qmcpack)
    meshfactor      = 1.0,           # bspline grid spacing, larger is finer
    jastrows        = [
        dict(type     = 'J1',        # 1-body
             function = 'bspline',   # bspline jastrow
             size     = 8),          # with 8 knots
        dict(type     = 'J2',        # 2-body
             function = 'bspline',   # bspline jastrow
             size     = 8)           # with 8 knots
        ],
    # opt parameters (qmcpack)
    perform_opt     = True,     # produce optimal jastrows
    block_opt       = False,    # if true, ignores opt and qmc
    skip_submit_opt = False,    # if true, writes input files, does not run opt
    opt_kpoint      = 'L',      # supercell k-point for the optimization
    opt_calcs       = [         # qmcpack input parameters for opt
        loop(max = 4,                        # No. of loop iterations
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
             ),
        loop(max = 4,
             qmc = linear(                   # same as above, except
                energy               =  0.5, # cost function
                unreweightedvariance =  0.0, #   is 50/50 energy and r.w. var.
                reweightedvariance   =  0.5, 
                timestep             =  0.5,  
                warmupsteps          =  100, 
                samples              = 64000,# and there are more samples 
                stepsbetweensamples  =   10, 
                blocks               =   10,   
                minwalkers           =   0.1, 
                bigchange            =  15.0,
                alloweddifference    =  1.0e-4
                )
             )
        ],
    # qmc parameters (qmcpack)
    block_qmc       = False,    # if true, ignores qmc
    skip_submit_qmc = False,    # if true, writes input file, does not run qmc
    qmc_calcs       = [         # qmcpack input parameters for qmc
        vmc(                      # vmc parameters 
            timestep      = 0.5,  # vmc timestep (1/Ha)
            warmupsteps   = 100,  # No. of MC steps before data is collected
            blocks        = 200,  # No. of data blocks recorded in scalar.dat
            steps         =  10,  # No. of steps per block
            substeps      =   3,  # MC steps taken w/o computing E_local
            samplesperthread = 40 # No. of dmc walkers per thread
            ),                    
        dmc(                      # dmc parameters
            timestep      = 0.01, # dmc timestep (1/Ha)
            warmupsteps   =  50,  # No. of MC steps before data is collected
            blocks        = 400,  # No. of data blocks recorded in scalar.dat
            steps         =   5,  # No. of steps per block
            nonlocalmoves = True  # use Casula's T-moves
            ),                    #  (retains variational principle for NLPP's)
        ],
    # return a list or object containing simulations
    return_list = False
    )


#the project manager monitors all runs
pm = ProjectManager()  

# give it the simulation objects
pm.add_simulations(qsims.list()) 

# run all the simulations    
pm.run_project()  



# print out the total energy
performed_runs = not settings.generate_only and not settings.status_only
if performed_runs:
    # get the qmcpack analyzer object
    # it contains all of the statistically analyzed data from the run
    qa = qsims.qmc.load_analyzer_image()
    # get the local energy from dmc.dat
    le = qa.dmc[1].dmc.LocalEnergy  # dmc series 1, dmc.dat, local energy
    #  print the total energy for the 8 atom system
    print 'The DMC ground state energy for graphene is:'
    print '    {0} +/- {1} Ha'.format(le.mean,le.error)
#end if



