#!/bin/bash
#COBALT -q debug-cache-quad
#COBALT -A PSFMAT_2 
#COBALT -n 8  
#COBALT -t 60
#COBALT -O CI_TESTs 
#COBALT --attrs mcdram=cache:numa=quad

##COBALT -q debug-cache-quad
export ATP_ENABLED=1

file_prefix=qmc_ref-SD-C4_AE_Mol_QP
#file_prefix=qmc_ref-MD-C4_AE_Mol_Excited_QP
#file_prefix=qmc_ref-MD-C4_AE_Mol_Ground_QPO

exe=/home/abenali/THETA/qmcpack/build_real/bin/qmcpack

NCORES=64
HT=4
NTHREADS=$((NCORES * HT))

aprun -n 8 -N 1 -cc depth -d $NTHREADS -j $HT $exe ${file_prefix}.in.xml > ${file_prefix}.output 

