#!/bin/csh
#COBALT -A QMC_2014_training 
#COBALT -O output_$jobid
#COBALT -t 10
#COBALT -n 1 

# location of codes
set qmc = "/soft/applications/qmcpack/build_XL_real/bin/qmcapp"
set convert = "/soft/applications/qmcpack/build_XL_real/bin/convert4qmc"

### MODIFY FROM HERE ###

###########################
#####      GAMESS     #####
###########################

# general 
set nodes = 1
set threads = 1
set tasks = 1  

# Single Det
set output_orbitals = "SCF.GAMESS.output"

# Multi-Det
set output_ci = "CI.GAMESS.output"
set thres = 0.0075
set norb = 40

# uncomment for single-determinant
runjob --np ${nodes} -p ${tasks} --block $COBALT_PARTNAME --verbose=INFO --envs OMP_NUM_THREADS=$threads : $convert -gamessAscii ${output_orbitals} -add3BodyJ

# uncomment for multi-determinant 
#runjob --np ${nodes} -p ${tasks} --block $COBALT_PARTNAME --verbose=INFO --envs OMP_NUM_THREADS=$threads : $convert -gamessAscii ${output_orbitals} -ci ${output_ci} -threshold ${thres} -readInitialGuess ${norb} -add3BodyJ

mv sample.Gaussian-G2.xml H2O.wfs.xml
mv sample.Gaussian-G2.ptcl.xml H2O.ptcl.xml
