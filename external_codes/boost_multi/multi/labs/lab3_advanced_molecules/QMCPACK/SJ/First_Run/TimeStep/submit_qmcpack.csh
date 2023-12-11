#!/bin/csh
#COBALT -A QMC_2014_training 
#COBALT -t 60
#COBALT -n 32 

# location of codes
set qmc = "/soft/applications/qmcpack/build_XL_real/bin/qmcapp"
set convert = "/soft/applications/qmcpack/build_XL_real/bin/convert4qmc"

### MODIFY FROM HERE ###

###########################
#####     QMCPACK     #####
###########################

# general 
set nodes = 32 
set threads = 16
set tasks = 1

set qmcpack_input = "vmc_dmc.xml" 
set qmcpack_output = "optm.output" 

# uncomment for single-determinant
runjob --np ${nodes} -p ${tasks} --block $COBALT_PARTNAME --verbose=INFO --envs OMP_NUM_THREADS=$threads : $qmc ${qmcpack_input} > ${qmcpack_output} 


