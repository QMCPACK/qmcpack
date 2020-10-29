#!/bin/bash

#One needs to source quantum_package
source ${PATH_TO_QP2}/qp2/quantum_package.rc

#Create the input file for Quantum Package 2 (QP2) using the BFD ECP and the corresponding Basis Set for CC-pvTz 
qp_create_ezfio C2.xyz -b cc-pvtz_ecp_bfd -p bfd

#Forcing Quantum Package to create every variable in the Input. (Avoid using VIM to edit input)
qp_edit -c C2.ezfio

# Save the AO orbitals (save time for restart)
sed -i s/"None"/"Write"/ C2.ezfio/ao_one_e_ints/io_ao_one_e_integrals       

# Save the MO orbitals (save time for restart)
sed -i s/"None"/"Write"/ C2.ezfio/mo_two_e_ints/io_mo_two_e_integrals 

# Set the max determinants to 6M instead of 1M
echo "             6000000" > C2.ezfio/determinants/n_det_max

###Run Hartree Fock
mpirun -n 1 qp_run scf C2.ezfio &> C2.HF.out 

#run the 4 index transformation (MO transformation) Serial run
qp_run four_idx_transform C2.ezfio 


###Run CIPSI in parallel with 10 nodes (1 Master + 9 slaves)
#Master needs one node only While the slave nodes will attach to this job
mpirun -n 1 qp_run fci C2.ezfio &> C2.FCI.out &
sleep 10
mpirun -n 9 qp_run -s fci C2.ezfio &> C2.FCI-slave.out 
wait

#Truncate the determinants to keep only weights of 10^-3 and higher. 
echo "1e-3" >  C2.ezfio/qmcpack/ci_threshold 
mpirun -n 1 qp_run truncate_wf_spin C2.ezfio &> C2.Trunc-1e-3.out


qp_run save_for_qmcpack C2.ezfio > C2-Dump-1e-3.out 
mv QP2QMCPACK.h5 C2.h5
convert4qmc -orbitals C2.h5  -multidets C2.h5  -add3BodyJ -prefix C2 

