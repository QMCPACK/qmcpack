#!/bin/bash

# TEST FOR QP2
source ${QP2}/qp2/quantum_package.rc

#Create the input file for Quantum Package 2 (QP2) using a modified  Basis Set for CC-pv5z not including the H orbitals for Li atom
qp_create_ezfio -b cc-pvdz C2.xyz  

###Run Hartree Fock 
mpirun -n 1 qp_run scf C2.ezfio 

echo "2" > C2.ezfio/determinants/n_states
echo "200" > C2.ezfio/determinants/n_det_max

mpirun -n 2 qp_run fci C2.ezfio &> C2.fci

qp_run save_for_qmcpack C2.ezfio 
mv QP2QMCACK.h5 C2.h5
convert4qmc -orbitals C2.h5 -multidet C2.h5 -addCusp -nojastrow 

