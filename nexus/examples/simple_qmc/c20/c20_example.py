#! /usr/bin/env python

from nexus import settings,ProjectManager,Job
from nexus import Structure,PhysicalSystem
from nexus import loop,linear,vmc,dmc
from qmcpack_calculations import basic_qmc


#general settings for nexus
settings(
    pseudo_dir    = './pseudopotentials',# directory with all pseudopotentials
    sleep         = 3,                   # check on runs every 'sleep' seconds
    generate_only = 0,                   # only make input files
    status_only   = 0,                   # only show status of runs
    machine       = 'node16',            # local machine is 16 core workstation
    )



#generate the C20 physical system
# specify the xyz file
structure_file = 'c20.cage.xyz'
# make an empty structure object
structure = Structure()
# read in the xyz file
structure.read_xyz(structure_file)
# place a bounding box around the structure
structure.bounding_box(
    box   = 'cubic',         # cube shaped cell
    scale = 1.5              # 50% extra space
    )
# make it a gamma point cell
structure.add_kmesh(
    kgrid      = (1,1,1),    # Monkhorst-Pack grid
    kshift     = (0,0,0)     # and shift
    )
# add electronic information
c20 = PhysicalSystem(
    structure = structure,   # C20 structure
    net_charge = 0,          # net charge in units of e
    net_spin   = 0,          # net spin in units of e-spin
    C          = 4           # C has 4 valence electrons
    ) 


#generate the simulations for the qmc workflow
qsims = basic_qmc(
    # subdirectory of runs 
    directory       = 'c20_test',
    # description of the physical system
    system          = c20,
    pseudos         = ['C.BFD.upf',  # pwscf PP file
                       'C.BFD.xml'], # qmcpack PP file
    # job parameters
    scfjob          = Job(cores=16), # cores to run scf 
    optjob          = Job(cores    = 16,        # cores for optimization 
                          app_name = 'qmcapp'), # use real-valued qmcpack
    qmcjob          = Job(cores    = 16,        # cores for qmc
                          app_name = 'qmcapp'), # use real-valued qmcpack
    # dft parameters (pwscf)
    functional      = 'lda',         # dft functional
    ecut            =  150 ,         # planewave energy cutoff (Ry)
    conv_thr        =  1e-6,         # scf convergence threshold (Ry)
    mixing_beta     =    .7,         # charge mixing factor
    # qmc wavefunction parameters (qmcpack)
    bconds          = 'nnn',         # open boundary conditions
    meshfactor      = 1.0,           # bspline grid spacing, larger is finer
    jastrows        = [
        dict(type     = 'J1',        # 1-body
             function = 'bspline',   # bspline jastrow
             size     = 8,           # with 8 knots
             rcut     = 6.0),        # and a radial cutoff of 6 bohr
        dict(type     = 'J2',        # 2-body
             function = 'bspline',   # bspline jastrow
             size     = 8,           # with 8 knots
             rcut     = 8.0),        # and a radial cutoff of 8 bohr
        ],
    # opt parameters (qmcpack)
    perform_opt     = True,     # produce optimal jastrows
    block_opt       = False,    # if true, ignores opt and qmc
    skip_submit_opt = False,    # if true, writes input files, does not run opt
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
    #  print the total energy for the 20 atom system
    print 'The DMC ground state energy for C20 is:'
    print '    {0} +/- {1} Ha'.format(le.mean,le.error)
#end if





