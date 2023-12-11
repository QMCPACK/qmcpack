#!/bin/bash

# TEST FOR QP2
source /soft/applications/quantum_package/quantum_package.rc

#Create the input file for Quantum Package 2 (QP2) using a modified  Basis Set for CC-pv5z not including the H orbitals for Li atom
qp_create_ezfio -b 'cc-pv5z' LiH.xyz  

###Run Hartree Fock 
mpirun -n 1 qp_run scf LiH.ezfio &> LiH.qp.out 


qp_run save_for_qmcpack LiH.ezfio &> LiH.dump 

convert4qmc -orbitals QP2QMCACK.h5 -nojastrow 

