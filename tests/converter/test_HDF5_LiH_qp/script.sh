#!/bin/bash

#On cooley, one needs to source quantum_package
source /soft/applications/quantum_package/quantum_package.rc

#Create the input file for Quantum Package(QP) using a modified  Basis Set for CC-pv5z not including the H orbitals for Li atoms as they are not handled in the AOS implementation
qp_create_ezfio_from_xyz LiH.xyz -b 'cc-pv5z' 

###Run Hartree Fock in parallel with 10 nodes (1 Master + 9 slaves)
#Master needs one node only While the slave nodes will attach to this job
mpirun -n 1 qp_run SCF LiH.ezfio &> LiH.qp.out &
sleep 10
mpirun -n 9 qp_run -slave qp_ao_ints LiH.ezfio &> LiH-slave.out 
wait


qp_run save_for_qmcpack LiH.ezfio &> LiH.dump 

convert4qmc -QP C2-Dump-1e-3.out -nojastrow -hdf5 

